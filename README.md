# RareEffect

RareEffect is a novel method to estimate the variant-level effect size (beta) and region-level (gene-level) heritability for rare variants.
RareEffect is integrated into SAIGE and is available in SAIGE version released after June 13, 2024 or later.

## System Requirements

RareEffect (integrated in SAIGE) is available on Linux. (We tested RareEffect integrated in SAIGE released on June 13, 2024 (commit hash: 7d7831b) on Ubuntu 20.04 LTS.)
All software dependencies (including version numbers) are specified in `RareEffect.yml`.

## Installation

We recommend the installation using the conda environment. The installation process took about 10 minutes.

1. Download and install Anaconda (or miniconda).
2. Create a conda environment using the following command:

```
conda env create -f RareEffect.yml
```

This command will make a conda environment named `SAIGE_RE`.

3. Activate `SAIGE_RE` environment and set system variables.

```
conda activate SAIGE_RE
FLAGPATH=`which python | sed 's|/bin/python$||'`
export LDFLAGS="-L${FLAGPATH}/lib"
export CPPFLAGS="-I${FLAGPATH}/include"
```

4. Clone SAIGE source code from the GitHub repository using the following command:

```
src_branch=main
repo_src_url=https://github.com/saigegit/SAIGE
git clone --depth 1 -b $src_branch $repo_src_url
```

5. Install required packages using the following command:

```
Rscript ./SAIGE/extdata/install_packages.R
```

6. Install SAIGE package using the following command. This command should be executed from the directory containing the SAIGE code.

```
R CMD INSTALL .
```

## Demo (example data)

1. The inputs for RareEffect are as follows:

  * Genotype file (plink bed/bim/fam)
  * Group file (same format as the group file for SAIGE-GENE+)
  * SAIGE step 1 output (`.rda` file)

You can use the example step 1 output in SAIGE package located in the following path. (In the real data analysis, you need to run SAIGE step 1 using your data.)

```
extdata/output/example_quantitative.rda
```
 
2. Execute the main function of RareEffect. An example command is provided below. This step will be finished in 1 second for the example data.

```
Rscript run_RareEffect.R \
    --rdaFile output/example_quantitative.rda \
    --chrom 1 \
    --geneName GENE1 \
    --groupFile input/group_new_snpid.txt \
    --traitType quantitative \
    --bedFile input/genotype_100markers.bed \
    --bimFile input/genotype_100markers.bim \
    --famFile input/genotype_100markers.fam \
    --macThreshold 10 \
    --collapseLoF FALSE \
    --outputPrefix output/RareEffect_example_output
```

  * rdaFile: SAIGE step 1 output file
  * chrom: chromosome number (eg. 11)
  * geneName: name of gene to analyze (eg. APOC3)
  * groupFile: group file containing variant list and functional annotations.
  * traitType: type of trait to be analyzed, `binary` or `quantitative`
  * bedFile, bimFile, famFile: genotype file
  * macThreshold: minor allele count threshold for ultra-rare variant collapsing
  * collapseLoF: if true, RareEffect collapses all LoF variants into one super-variant like Burden test (regardless of their MAC)
  * outputPrefix: path to output

3. You will get two output files: variant-level effect size and region-level heritability.

  * `RareEffect_example_output_effect.txt`
    
```
variant effect PEV
rs31 3.37458398173721e-07 9.99984127282007e-07
rs33 -2.01027399753112e-07 9.99983127313752e-07
lof_UR -5.5315652973989e-07 9.99955129043553e-07
rs17 -1.3901101330302e-07 9.99988116933531e-07
mis_UR 7.00908607301673e-07 9.99989116910765e-07
```

  * `RareEffect_example_output_h2.txt`

```
LoF mis syn all
7.80000038242669e-08 2.30000021572006e-08 0 1.01000005981467e-07
```

4. Using the variant-level effect size from RareEffect, you can calculate the polygenic scores. This step will be finished in 1 second for the example data.

```
Rscript calculate_RareEffect_PRS.R \
    --effectFile output/RareEffect_example_output_effect.txt \
    --bedFile input/genotype_100markers.bed \
    --bimFile input/genotype_100markers.bim \
    --famFile input/genotype_100markers.fam \
    --groupFile input/group_new_snpid.txt \
    --geneName GENE1 \
    --outputPrefix output/RareEffect_example_PRS_output
```

  * effectFile: variant-level effect size output from RareEffect main function
  * variantListFile (optional): list of common and rare, but non-ultra-rare variants. These variants will not be used for ultra-rare variant scoring.

  * Individual RareEffect PRS
```
IID	PRS
1a1	0
1a2	0
1a3	0
1a4	0
1a5	0
1a6	0
1a7	0
1a8	0
1a9	0
1a10	0
1a11	0
1a12	0
1a13	0
1a14	0
1a15	0
1a16	-1.3901101330302e-07
1a17	0
1a18	0
1a19	-1.3901101330302e-07
1a20	0
```

## How to run the software on your own data

1. RareEffect main function

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
    --outputPrefix=$PATH_TO_OUTPUT_FILE
```

2. PRS calculation

```
Rscript calculate_RareEffect_PRS.R \
    --effectFile $PATH_TO_EFFECT_SIZE_FILE \
    --bedFile $PATH_TO_BED_FILE \
    --bimFile $PATH_TO_BIM_FILE \
    --famFile $PATH_TO_FAM_FILE \
    --groupFile $PATH_TO_GROUP_FILE \
    --geneName $GENE_NAME \
    --variantListFile $PATH_TO_VARIANT_LIST_FILE \
    --outputPrefix $PATH_TO_OUTPUT_FILE
```

## Example output from the real data

We estimated variant-level effect size and gene-level heritability for genome-wide significant genes (or top 10 genes) using 393,247 White British individuals in UKB WES data.
The example output files can be downloaded at:
  * Variant-level effect size: https://storage.googleapis.com/leelabsg/RareEffect/RareEffect_effect_size.zip
  * Gene-level signed heritability: https://storage.googleapis.com/leelabsg/RareEffect/RareEffect_h2.zip
