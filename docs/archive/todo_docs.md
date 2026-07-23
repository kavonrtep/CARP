# Bugfix:
- In Density plot in html report when satellited density from TideCluster are shown, labels on side show only TRC_1(?bp)  - only question mark without actuall value


# Feature/maybe fix
- LTR_RT_TR  structural category - this structural category should be included as separate track in density plots
- this should be also shown in Full classification table in Complete TEs count as value in parathesis
- so it wll shou cound and how many are in LTR-RT-TR structure - thios must be explain in legend what this category actually mean
- potential fix - verify that the overlaping LTR does not count twice when reporting composition
- For child which are part of LTR_RT_TR, it should have attibute which will make it clear that this feature is part of LTR_RT_TR 
- Examble dataset/analysis output is `/mnt/ceph/454_data/TideCluster/tmp/run-000116/mnt/ceph/454_data/TideCluster/tmp/run-000116` - this finished analysis is suitable test case for improved reporting as it contain number LTR_RT_TR elements

# Improvements 
## File cleanup audit:
- Cleaning up unneccessary file (but not to breake anything downsterm - check manifest)
- list of dir/files which could be removed after analysis is completed:
  - `DANTE/DANTE_filtered.gff3.tmp.gff3`
  - `DANTE_LTR/library/mmseqs2/mmseqs_all_seqs.fasta`
  - `DANTE_LTR/library/mmseqs2/partitioned_s900_w1000.fasta`
  - `DANTE_TIR/DANTE_TIR.RData`
  - in DANTE_TIR_FALLBACK - `TPase_5prime_alignment.tsv` and `TPase_3prime_alignment.tsv` files
  - in RepeatMasker dir - `genome_cleaned.fasta`
  - possible to delete? - `DANTE_LTR/library/TE*.fasta`   - verify that it is safe
  - `TideCluster/default/TideCluster_kite/kitehor.periodogram1` - verify that it is safe
  - `TRC_X.fasta_YY.kmers` file in _tarean directories - verify it is safe
- Cleaning should be dafaul but it should be possible to run it with 'keep all'/debug mode when files are kept
- Any other significan(large, or large number of small) files which can be safelly removed? We can create minimalistic/maximalistic cleanup variant. 
- Deletion should be probably done on level of run script - not snakemake? and only if run is without errors

## Evaluate result from two version - include so discrepancy in LINE
Same genome, distinct LINE quantity, total repeat content is the same
- /mnt/ceph/454_data/TideCluster/tmp/run-000082/mnt/ceph/454_data/TideCluster/tmp/run-000082
- /mnt/ceph/454_data/TideCluster/tmp/run-000097/mnt/ceph/454_data/TideCluster/tmp/run-000097
Identify source of this difference if there is not underlying bug/problem
