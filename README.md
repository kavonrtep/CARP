# Assembly Repeat Annotation Pipeline


## Requirements 
Singularity is required to use the container. Singularity can be installed using conda environment. 

```bash
conda create -n singularity3 -c conda-forge "singularity>=3.6"
conda activate singularity3
```

## Quick Start
Singularity image (.sif file) can be downloaded from https://github.com/kavonrtep/assembly_repeat_annotation_pipeline/releases

Format of config.yaml file is as follows:

```yaml
genome_fasta: data/CEN6_ver_220406_part.fasta
output_dir: output
custom_library: data/pisum_custom_library.fasta
tandem_repeat_library: data/FabTR_all_sequences_210901.db.RM_format.fasta
```

Lines specifying `custom_library` and `tandem_repeat_library` are optional. File `custom_library` is 
used by RepeatMasker fo similarity based annotation. Sequences must be as FASTA with 
sequence IDs in format `>repeatname#class/subclass`
Classification categories are: 

```text
Class_I/LINE
Class_I/LTR/Ty1_copia
Class_I/LTR/Ty1_copia/Ale
Class_I/LTR/Ty1_copia/Angela
Class_I/LTR/Ty1_copia/Bianca
Class_I/LTR/Ty1_copia/Ikeros
Class_I/LTR/Ty1_copia/Ivana
Class_I/LTR/Ty1_copia/SIRE
Class_I/LTR/Ty1_copia/TAR
Class_I/LTR/Ty1_copia/Tork
Class_I/LTR/Ty3_gypsy
Class_I/LTR/Ty3_gypsy/chromovirus
Class_I/LTR/Ty3_gypsy/chromovirus/CRM
Class_I/LTR/Ty3_gypsy/chromovirus/Reina
Class_I/LTR/Ty3_gypsy/chromovirus/Tekay
Class_I/LTR/Ty3_gypsy/non-chromovirus/OTA
Class_I/LTR/Ty3_gypsy/non-chromovirus/OTA/Athila
Class_I/LTR/Ty3_gypsy/non-chromovirus/OTA/Tat
Class_I/LTR/Ty3_gypsy/non-chromovirus/OTA/Tat/Ogre
Class_I/LTR/Ty3_gypsy/non-chromovirus/OTA/Tat/Retand
Class_II/Subclass_1/TIR/EnSpm_CACTA
Class_II/Subclass_1/TIR/hAT
Class_II/Subclass_1/TIR/MITE
Class_II/Subclass_1/TIR/MITE/Stowaway
Class_II/Subclass_1/TIR/MuDR_Mutator
Class_II/Subclass_1/TIR/PIF_Harbinger
Class_II/Subclass_1/TIR/Tc1_Mariner
Class_II/Subclass_2/Helitron
rDNA_45S/18S
rDNA_45S/25S
rDNA_45S/5.8S
rDNA_5S/5S
```
It is possible to include additional categories to repeat library. Note that `Subclass_1` repeats are used to clean `Class_I/LTR` library which is built by DANTE_LTR thus it is critical to use correct codes for `Class_II/Subclass_1` mobile elements in repeat library.

File `custom_library` is used by TideCluster to annotate discovered tandem repeats based
on the similarity. Format is the same as above repeat database. E.g. 
`>sequence_id/Satellite/PisTR-B`

To run annotation pipeline, execute following command:

```bash
singularity run -B /path/to/ -B $PWD assembly_repeat_annotation_pipeline.sif -c config.yaml -t 20
```
Parameter `-t` specifies the number of threads to use. Singularity parameter `-B` is used to bind the input and output directories to the container. Without this parameter, the container will not be able to access the input and output files. File `config.yaml` must be also in directory which is accessible to the container. In the example above this is the current directory `$PWD`. 


## Running pipeline on metacentrum
Use [./scripts/annotate_repeats_metacentrum.sh](./scripts/annotate_repeats_metacentrum.sh) script to run the pipeline on metacentrum. Adjust paths to the input files , output directory and singularity image. 

 

 



## Output structure
TODO

## Build the container

To build the container, run the following command:

```bash
SINGULARITY=`which singularity`
sudo $SINGULARITY build assembly_repeat_annotation_pipeline_0.4.0.sif Singularity
```

