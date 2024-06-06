# assembly_repeat_annotation_pipeline




## Build the container

To build the container, run the following command:

```bash
SINGULARITY=`which singularity`
sudo $SINGULARITY build assembly_repeat_annotation_pipeline.sif Singularity
```



 RepeatMasker -pa {threads}  $genome_absolute_path  -lib $library_absolute_path -dir $PWD -s -xsmall -e ncbi -no_is