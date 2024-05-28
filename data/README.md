# custom library example:

concatenate from multiple files, remove string "All/repeat/mobile_element/" with ""

```bash
SRC_DIR="/mnt/raid/454_data/Pisum_pangenome/assemblies/JI2822__2023-06-22/analysis/RM_LIBS/"
cat ${SRC_DIR}/PST_Zaba_Stowaway_Psmar_from_old_papers \
 ${SRC_DIR}/RE_RM_LINE_library_short_names.fasta \
 ${SRC_DIR}/rDNA_45S_Pisum.RM_format.fasta \
 ${SRC_DIR}/RM_lib_custom_Class_II_and_pararetrovirus.fasta | \
    sed 's/All\/repeat\/mobile_element\///' > pisum_custom_library.fasta
 grep ">" pisum_custom_library.fasta | cut -d "#" -f 2 | sort | uniq -c 
```

## Composition of library:
  165 Class_II/Subclass_1/TIR/EnSpm_CACTA
      5 Class_II/Subclass_1/TIR/hAT
      2 Class_II/Subclass_1/TIR/MITE
     17 Class_II/Subclass_1/TIR/MITE/Stowaway
    459 Class_II/Subclass_1/TIR/MuDR_Mutator
      1 Class_II/Subclass_1/TIR/Tc1_Mariner
     16 Class_II/Subclass_2/Helitron
     63 Class_I/LINE
      6 Class_I/pararetrovirus
      1 rDNA_45S/18S
      1 rDNA_45S/25S
      1 rDNA_45S/5_8S
      1 rDNA_45S/IGS
      1 rDNA_45S/ITS1
      1 rDNA_45S/ITS2

Library should include repeats but without Tandem Repeats and LTR-RTs. Tandem repeat library is annotated with TideCluster and LTR-RTs library is annotated with DANTE_LTR. 


# tandem repeat library example:

```bash
cp /mnt/ceph/454_data/vicieae/databases/satellites/db/FabTR_all_sequences_210901.db.RM_format FabTR_all_sequences_210901.db.RM_format.fasta
```