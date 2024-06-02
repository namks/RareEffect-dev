# RareEffect

RareEffect is a novel method to estimate the variant-level effect size (beta) and region-level (gene-level) heritability for rare variants.
RareEffect is integrated into SAIGE and is available in SAIGE version 1.3.5 or higher.

## Procedure

1. First, run SAIGE step 1. The output file from SAIGE step 1 (`.rda` file) will be used as the input for RareEffect.
2. The inputs for RareEffect are as follows:

  * Genotype file (plink bed/bim/fam)
  * Group file (same format as the group file for SAIGE-GENE+)
  * SAIGE step 1 output rda file

3. Execute the main function of RareEffect. An example command is provided below.

```
Rscript run_RareEffect.R \
    --rdaFile=$PATH_TO_RDA_FILE \
    --chrom=$CHR \
    --geneName=$GENE_NAME \
    --groupFile=$PATH_TO_GROUP_FILE \
    --traitType=binary \
    --bedFile=$PATH_TO_BED_FILE \
    --bimFile=$PATH_TO_BIM_FILE \
    --famFile=$PATH_TO_FAM_FILE \
    --macThreshold=10 \
    --collapseLoF=FALSE \
    --collapsemis=FALSE \
    --collapsesyn=FALSE \
    --outputPrefix=$PATH_TO_OUTPUT_FILE
```

4. RareEffect generates two output files: variant-level effect size and region-level heritability.
5. Using the variant-level effect size from RareEffect, you can calculate the polygenic score in another cohort.

```
Rscript calculate_RareEffect_PRS.R \
    --
```
