#!/usr/bin/env python
'''
A wrapper for RepeatMasker - breaks sequences into chunks of 50Mbp (with 10kb overlap)
and runs RepeatMasker on each chunk, running in parallel but each RepeatMasker job is run
as a single thread. - see https://github.com/Dfam-consortium/RepeatMasker/issues/274
implemented just using standard python libraries.
'''
import argparse
import shutil
import subprocess
import multiprocessing
import os
import sys
import tempfile

# Canonical RepeatMasker .out header (2 columns lines + 1 blank). Only used as
# a fallback when NO chunk produced a .out at all (a genome with zero library
# hits — pathological; rDNA normally always matches). Exact spacing is
# immaterial because downstream parsers skip the first two lines and there are
# no data rows in that case.
RM_OUT_HEADER_LINES = [
    "   SW   perc perc perc  query      position in query    matching  repeat      position in repeat\n",
    "score   div. del. ins.  sequence   begin end   (left)   repeat    class/family begin end (left)  ID\n",
    "\n",
]


def read_fasta_sequence_size(fasta_file):
    """Read size of sequence into dictionary"""
    fasta_dict = {}
    with open(fasta_file, 'r') as f:
        for line in f:
            if line[0] == '>':
                header = line.strip().split()[0][1:]  # remove part of name after space
                fasta_dict[header] = 0
            else:
                fasta_dict[header] += len(line.strip())
    return fasta_dict

def read_single_fasta_to_dictionary(fh):
    """
    Read fasta file into dictionary
    :param fh:
    :return:
    fasta_dict
    """
    fasta_dict = {}
    for line in fh:
        if line[0] == '>':
            header = line.strip().split()[0][1:]  # remove part of name after space
            fasta_dict[header] = []
        else:
            fasta_dict[header] += [line.strip()]
    fasta_dict = {k: ''.join(v) for k, v in fasta_dict.items()}
    return fasta_dict


def iter_fasta_records(fasta_file):
    """Yield ``(header, sequence)`` one record at a time.

    Only the current record is held in memory (peak = the largest single
    sequence), not the whole file — the streaming replacement for
    ``read_single_fasta_to_dictionary``, which materialised the entire
    (genome-sized) FASTA in RAM (~90 GB on a 90 Gbp assembly). The header is the
    token up to the first whitespace, matching ``read_fasta_sequence_size`` /
    ``read_single_fasta_to_dictionary`` so the same ``matching_table`` keys line
    up. The file handle is closed on exit (the old callers leaked it).
    """
    header = None
    chunks = []
    with open(fasta_file, 'r') as fh:
        for line in fh:
            if line and line[0] == '>':
                if header is not None:
                    yield header, ''.join(chunks)
                header = line.strip().split()[0][1:]
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, ''.join(chunks)


def split_fasta_to_chunks(fasta_file, chunk_size=100000000, overlap=100000, temp_dir=None):
    """
    Split fasta file to chunks, sequences longe than chuck size are split to overlaping
    peaces. If sequences are shorter, chunck with multiple sequences are created.
    :param fasta_file:

    :param fasta_file:
    :param chunk_size:
    :param overlap:
    :return:
    fasta_file_split
    matching_table
    """
    min_chunk_size = chunk_size * 2
    print("analyzing fasta file")
    fasta_dict = read_fasta_sequence_size(fasta_file)
    print("fasta file loaded")
    print("number of sequences: ", len(fasta_dict))
    # calculates ranges for splitting of fasta files and store them in list
    matching_table = []
    fasta_file_split = tempfile.NamedTemporaryFile(delete=False, dir=temp_dir).name
    print("calculating ranges for splitting of fasta files")
    for header, size in fasta_dict.items():
        if size > min_chunk_size:
            number_of_chunks = int(size / chunk_size)
            adjusted_chunk_size = int(size / number_of_chunks)
            for i in range(number_of_chunks):
                start = i * adjusted_chunk_size
                end = ((i + 1) *
                       adjusted_chunk_size
                       + overlap) if i + 1 < number_of_chunks else size
                new_header = header + '_' + str(i)
                matching_table.append([header, i, start, end, new_header])
        else:
            new_header = header + '_0'
            matching_table.append([header, 0, 0, size, new_header])
    print("splitting fasta file")
    # Index matching_table rows by original header once. This was previously a
    # per-header `[x for x in matching_table if x[0]==header]` scan inside the
    # loop -> O(headers x chunks), pathological on highly fragmented assemblies.
    rows_by_header = {}
    for row in matching_table:
        rows_by_header.setdefault(row[0], []).append(row)
    # Stream the FASTA one record at a time and write its chunk slices — the old
    # read_single_fasta_to_dictionary held the entire genome in RAM (~90 GB on a
    # 90 Gbp assembly). Records stream in file order, so the output is identical.
    with open(fasta_file_split, 'w') as fh_out:
        for header, sequence in iter_fasta_records(fasta_file):
            for header2, i, start, end, new_header in rows_by_header.get(header, []):
                fh_out.write('>' + new_header + '\n')
                fh_out.write(sequence[start:end] + '\n')
    return fasta_file_split, matching_table

def make_temp_files(number_of_files):
    """
    Make named temporary files, file will not be deleted upon exit!
    :param number_of_files:
    :return:
    filepaths
    """
    temp_files = []
    for i in range(number_of_files):
        temp_files.append(tempfile.NamedTemporaryFile(delete=False).name)
        os.remove(temp_files[-1])
    return temp_files


def get_new_header_and_coordinates(header, start, end, matching_table):
    """
    Get new header and coordinates for sequence
    :param header:
    :param start:
    :param end:
    :param matching_table:
    :return:
    new_header
    new_start
    new_end
    """
    matching_table_part = [x for x in matching_table if x[0] == header]
    new_coords = []
    for chunk in matching_table_part:
        if chunk[2] <= start < chunk[3]:
            new_header = chunk[4]
            new_start = start - chunk[2]
            new_end = end - chunk[2]
            new_sequence_length = chunk[3] - chunk[2]
            new_coords.append([new_header, new_start, new_end, new_sequence_length])
    return new_coords


def get_original_header_and_coordinates(new_header, new_start, new_end, matching_table):
    """
    Get original header and coordinates for sequence
    :param new_header:
    :param new_start:
    :param new_end:
    :param matching_table:
    :return:
    original_header
    original_start
    original_end
    """
    matching_table_part = [x for x in matching_table if x[4] == new_header]
    real_chunk_size = matching_table_part[0][3] - matching_table_part[0][2]
    ori_header = matching_table_part[0][0]
    start = matching_table_part[0][2]
    ori_start = new_start + start
    ori_end = new_end + start
    return ori_header, ori_start, ori_end, real_chunk_size



def recalculate_rm_out_coordinates(rm_out_file, matching_table_dict,
                                   rm_out_file_recalculated):

    """
    Recalculate RepeatMasker output coordinates
    :param rm_out_file:
    :param matching_table:
    :return:
    rm_out_file_recalculated

    out format is :
   SW   perc perc perc  query                                  position in query              matching                                              repeat                                                   position in repeat
score   div. del. ins.  sequence                               begin    end          (left)   repeat                                                class/family                                         begin   end    (left)       ID

14724    1.8  0.0  0.0  ENA|CALLYH020000001|CALLYH020000001.1       285     2092   (316552) + Contig5_17                                            rDNA_45S/18S                                               1   1808     (0)       1
 1059   10.7  0.0  1.3  ENA|CALLYH020000001|CALLYH020000001.1      2326     2486   (316158) + Contig2_19                                            rDNA_45S/5.8S                                              2    160     (4)       2
24611    6.4  0.1  0.1  ENA|CALLYH020000001|CALLYH020000001.1      2709     6097   (312547) + 24_25S_PS_CL2Contig6_rc/3985-7377                     rDNA_45S/25S                                               4   3393     (0)       3
  241   15.8  0.0  5.3  ENA|CALLYH020000001|CALLYH020000001.1      6223     6282   (312362) + 24_25S_PS_CL2Contig6_rc/3985-7377                     rDNA_45S/25S                                            3337   3393     (0)       4
   33    5.2  0.0  0.0  ENA|CALLYH020000001|CALLYH020000001.1      7758     7797   (310847) + (CGGGTG)n                                             Simple_repeat                                              1     40     (0)       5
14688    1.8  0.0  0.1  ENA|CALLYH020000001|CALLYH020000001.1     10564    12373   (306271) + Contig5_17                                            rDNA_45S/18S                                               1   1808     (0)       6
 1059   10.7  0.0  1.3  ENA|CALLYH020000001|CALLYH020000001.1     12607    12767   (305877) + Contig2_19                                            rDNA_45S/5.8S                                              2    160     (4)       7


query sequence and begin and end positions will be recalculated



    """

    with open(rm_out_file, 'r') as fh_in:
        with open(rm_out_file_recalculated, 'w') as fh_out:
            # first thre lines are header, third line must be empty and fourth line
            # contain first data line
            # skip firt two lines
            fh_in.readline()
            fh_in.readline()
            # check if third line is empty
            if fh_in.readline().strip() != '':
                raise ValueError('Third line in RepeatMasker output file is not empty')
            for line in fh_in:
                line = line.strip()
                if line == '':
                    continue
                fields = line.split()
                query = fields[4]
                begin = int(fields[5])
                end = int(fields[6])

                ori_header = matching_table_dict[query][0]
                ori_start = matching_table_dict[query][2] + begin
                ori_end = matching_table_dict[query][2] + end
                fields[4] = ori_header
                fields[5] = str(ori_start)
                fields[6] = str(ori_end)
                fh_out.write(' '.join(fields) + '\n')

def repeatmasker(fasta_file, library, workdir='.', sensitivity='default'):
    """
    Run RepeatMasker on fasta file
    :param library:
    :param fasta_file:
    :param threads:
    :param workdir:
    :param out:
    :return:
    """
    # run RepeatMasker

    sensitivity_options = {'quick': '-q', 'default': '', 'rush': '-qq'}

    command = ['RepeatMasker', '-pa', '1', '-dir', workdir,
               '-lib', library, '-xsmall', '-e ncbi', '-no_is', fasta_file,
               sensitivity_options[sensitivity]]
    print("----------")
    print(' '.join(command))
    os.system(' '.join(command))

    # get RepeatMasker output file
    rm_out_file = fasta_file + '.out'
    # get RepeatMasker gff3 file
    return rm_out_file

def split_fasta_to_files(fasta_file, workdir, chunk_size=5000000):
    """Pack chunk records into multi-FASTA files of <= chunk_size bp each.

    Previously this wrote ONE file per chunk record, so every small scaffold
    became its own RepeatMasker invocation — and each RM invocation
    re-preprocesses the entire library. On a fragmented assembly (thousands of
    kilobase scaffolds) that is thousands of redundant library preprocessings
    and a near-serial run. Here, consecutive chunk records are binned together
    up to chunk_size, collapsing those tiny scaffolds into a handful of
    multi-FASTA files.

    Chunk boundaries are UNCHANGED (the big-sequence pieces from
    split_fasta_to_chunks, ~chunk_size each, still land in a bin of their own);
    only the file grouping changes. RepeatMasker masks each record in a
    multi-FASTA independently, so the per-record .out rows — and therefore the
    final annotation — are identical; only the number of RM processes (and the
    redundant per-process library preprocessing) drops.
    """
    fasta_files = []
    cur_fh = None
    cur_size = 0
    # Stream one record at a time (bin-packing needs no random access) — the old
    # read_single_fasta_to_dictionary held the entire genome-sized split file in
    # RAM (~90 GB). Records stream in file order, so the bins are identical.
    for header, sequence in iter_fasta_records(fasta_file):
        seq_len = len(sequence)
        # Close the current bin when adding this record would push it over
        # chunk_size. A record >= chunk_size simply lands in a bin of its own.
        if cur_fh is not None and cur_size > 0 and cur_size + seq_len > chunk_size:
            cur_fh.close()
            cur_fh = None
        if cur_fh is None:
            path = tempfile.NamedTemporaryFile(delete=False, dir=workdir).name
            fasta_files.append(path)
            cur_fh = open(path, 'w')
            cur_size = 0
        cur_fh.write('>' + header + '\n')
        cur_fh.write(sequence + '\n')
        cur_size += seq_len
    if cur_fh is not None:
        cur_fh.close()
    return fasta_files

def setup_rmblast_culling_shim(culling_limit, workdir):
    """Point RMBLAST_DIR at a culling shim so every per-chunk RepeatMasker worker
    passes ``-culling_limit N`` to rmblastn (no RepeatMasker/container patch).

    The shim is built by the shared ``rmblast_culling_shim`` helper (also used by
    the tidecluster_reannotate rule). Returns the shim dir, or None when culling
    is off / rmblastn can't be located — in which case the run proceeds unchanged
    (the optimisation is never a correctness requirement). The RMBLAST_DIR env var
    must be set before the multiprocessing.Pool below so the forked workers
    inherit it.
    """
    # rmblast_culling_shim.py lives beside this script (scripts/); make it
    # importable in both the repo and container contexts.
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from rmblast_culling_shim import make_culling_shim
    shim_dir = make_culling_shim(culling_limit, workdir)
    if shim_dir is not None:
        os.environ['RMBLAST_DIR'] = shim_dir
        print(F"RepeatMasker culling enabled: -culling_limit {int(culling_limit)} "
              F"(RMBLAST_DIR shim: {shim_dir})")
    return shim_dir

def main():
    # parse arguments
    parser = argparse.ArgumentParser(description='RepeatMasker wrapper')
    parser.add_argument('-f','--fasta', help='Input fasta file', required=True)
    parser.add_argument('-d', '--workdir', help='Working directory', default='.')

    parser.add_argument('-p', '--threads', help='Number of threads', default=1, type=int)
    parser.add_argument('-o', '--out', help='output file (required)', required=True)
    parser.add_argument('-l', '--library', help='RepeatMasker library', required=True)
    parser.add_argument('-s', '--sensitivity', help='Sensitivity, possible values are '
                                                    'quick, default and rush',
                        required=False,default='default')
    parser.add_argument('--culling-limit', dest='culling_limit', type=int,
                        default=0,
                        help='rmblastn -culling_limit value (0 = off, the '
                             'default). N>0 caps redundant HSPs per locus via an '
                             'RMBLAST_DIR shim; 2 is the validated recommendation '
                             '(~3x faster, <1%% masked-bp cost on a redundant '
                             'de-novo library).')
    args = parser.parse_args()

    # create working directory if not exists
    if not os.path.exists(args.workdir):
        os.makedirs(args.workdir)

    tmp_dir = F"{args.workdir}/rm_tmp"
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # Optionally inject -culling_limit into every per-chunk rmblastn call. Must
    # be set up before the multiprocessing.Pool below so the RMBLAST_DIR env var
    # is inherited by the forked RepeatMasker workers.
    shim_dir = setup_rmblast_culling_shim(args.culling_limit, args.workdir)

    rm_chunk_size = 5000000
    fasta_parts, matching_table = split_fasta_to_chunks(args.fasta,
                                                        chunk_size=rm_chunk_size,
                                                        overlap=2000,
                                                        temp_dir=tmp_dir)

    matching_table_dict = {x[4]: x for x in matching_table}

    # Pack chunk records into <= rm_chunk_size multi-FASTA files so small
    # scaffolds share one RepeatMasker invocation (see split_fasta_to_files).
    fasta_parts_files = split_fasta_to_files(fasta_parts, tmp_dir, rm_chunk_size)

    print("Number of fasta parts: ", len(fasta_parts_files))

    # Dispatch order = largest part first (LPT / longest-processing-time-first
    # scheduling). Each RepeatMasker worker is single-threaded (culling pins
    # rmblastn to -num_threads 1 for determinism; see rmblast_culling_shim.py),
    # so with a few big chromosome-sized parts and many small scaffolds, file
    # order would leave a big part running alone on one core at the tail while
    # the other cores idle. Starting the big parts first lets the small ones
    # backfill the tail, minimising makespan. chunksize=1 so the pool hands out
    # one part at a time in this order (the default batches the iterable, which
    # would defeat LPT). Sort key is (-size, path) — fully deterministic.
    # NB: this only reorders *dispatch*; the recalculation and concatenation
    # loops below still walk fasta_parts_files in the original genome order, so
    # the merged .out is byte-identical regardless of scheduling.
    scheduled_parts = sorted(fasta_parts_files,
                             key=lambda f: (-os.path.getsize(f), f))
    rm_arguments = list(zip(scheduled_parts,
                            [args.library] * len(scheduled_parts),
                            [tmp_dir] * len(scheduled_parts),
                            [args.sensitivity] * len(scheduled_parts)))

    with multiprocessing.Pool(args.threads) as pool:
        pool.starmap(repeatmasker, rm_arguments, chunksize=1)

    # recalculate coordinates in RepeatMasker output files. RepeatMasker emits
    # no .out for a chunk with zero hits; treat a missing .out as contributing
    # no rows (the original code crashed with FileNotFoundError here on
    # fragmented inputs).
    for fasta_part in fasta_parts_files:
        rm_out = fasta_part + '.out'
        out_recalculated = fasta_part + '.out.recalculated'
        if not os.path.exists(rm_out):
            open(out_recalculated, 'w').close()
            continue
        recalculate_rm_out_coordinates(rm_out, matching_table_dict,
                                       out_recalculated)
    # concatenate RepeatMasker output files
    with open(args.out, 'w') as fh_out:
        # RM header (first three lines) from the first chunk that produced a
        # .out; fall back to the canonical header only if no chunk had hits.
        header_lines = None
        for fasta_part in fasta_parts_files:
            rm_out = fasta_part + '.out'
            if os.path.exists(rm_out):
                with open(rm_out, 'r') as fh_in:
                    header_lines = [fh_in.readline() for _ in range(3)]
                break
        fh_out.writelines(header_lines if header_lines is not None
                          else RM_OUT_HEADER_LINES)

        for fasta_part in fasta_parts_files:
            with open(fasta_part + '.out.recalculated', 'r') as fh_in:
                for line in fh_in:
                    fh_out.write(line)
    # remove tmp_dir, use python function to remove directory
    shutil.rmtree(tmp_dir)
    # remove the culling shim (if any)
    if shim_dir is not None and os.path.isdir(shim_dir):
        shutil.rmtree(shim_dir, ignore_errors=True)





if __name__ == '__main__':
    main()