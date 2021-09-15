# Integrating "Targeted Analysis of sequencing Reads for GenoTyping" (TARGT) Genotype Calls of Classical HLA Genes with the 1000 Genomes Project

This repo contains code associated with Pierini et al 202X. If you only use the code in this repo please cite TBD and if you use the updated TGP VCF file please cite  TBD and [Pierini, F., Nutsua, M., Böhme, L. _et al._ Targeted analysis of polymorphic loci from low-coverage shotgun sequence data allows accurate genotyping of HLA genes in historical human populations. _Sci Rep_ __10__, 7339 (2020)](https://doi.org/10.1038/s41598-020-64312-w).





## Generating the TGP VCF File

Before integrating the genotype calls from the TARGT pipeline we first need to download the phase 3 realse of the TGP chromosome 6 VCF file and it's corresponding index, which can be done with the following code:

```bash
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
```

and:

```bash
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi
```

For our purposes we are only interested in the protein coding regions of chromosome 6, using [BCFtools](https://samtools.github.io/bcftools/bcftools.html) we can filter the TGP VCF file to only contain our regions of interest with the following code (you can find the `.bed` file in the `resources` directory):

```bash
bcftools view -R chr6_ProtCod_GRCh37hg19_nonoverlapping_May2021.bed -Oz -o tgp_chr6_ProtCod_GRCh37hg19_nonoverlapping_unfiltered.vcf.gz ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
```

We then can index `tgp_chr6_ProtCod_GRCh37hg19_nonoverlapping_unfiltered.vcf.gz` using [Tabix](http://www.htslib.org/doc/tabix.html) with the following code:

```bash
tabix -p vcf tgp_chr6_ProtCod_GRCh37hg19_nonoverlapping_unfiltered.vcf.gz
```

Lastly, to make sure my computer doesn't yell at me during the analyses I will subset `tgp_chr6_ProtCod_GRCh37hg19_nonoverlapping_unfiltered.vcf.gz` from the start position of the HLA-A gene to the end position of HLA-DQB1 (the positions were defined by [Ensembl](https://useast.ensembl.org/index.html)) with the following code:

```bash
bcftools view -r 6:29909037-32636160 -Oz -o tgp_hla_ProtCod_GRCh37hg19_nonoverlapping_unfiltered.vcf.gz tgp_chr6_ProtCod_GRCh37hg19_nonoverlapping_unfiltered.vcf.gz
```





## TARGT Table Conversion and QC

I have preformatted the TARGT tables (which can be found in the `targt_tables` directory) to include the necessary header line information as specified by the [VCF documentation](https://samtools.github.io/hts-specs/VCFv4.3.pdf) here I make a couple of assumptions that should not effect down stream analysis:

* Since we called genotypes using the TARGT pipeline we do not have variant IDs associated with our new calls, so I annotate the `ID` field as missing using a `"."`

* Since the quality of each call was assed in the TARGT pipeline I annotate the filter statuts in the `FILTER` field as `PASS` for every site

In converting the TARGT genotype calls to the VCF format I make two additional assumptions that **_COULD_** effect downstream analysis:

* Since the TGP data is imputed I assume that missing genotype calls are homozygous for the reference allele
* At position `32552155` all samples are homozygous for the reference allele except for individual `HG03009` who has a deletion on one of their chromosomes, because I was hesitant to manually change the refernce allele from `G` to the proper format of `TG` and since every other sample is homozygous for the `G` allele I also assumed that individual `HG03009` is homozygous for the reference allele `G`

I believe I chose a conservative approach, but we can easily pivot as necessary. Now that we have the assumptions out of the way let's run:

```bash
python3 targt_tables_qc.py
```

which will output three files that we will now go over.

### `targt_hla_gt_calls_no_header.vcf`

Is a headerless (but otherwise perfectly formatted) VCF file consiting of **ONLY** the sites called by the TARGT pipeline.

### `targt_sample_validation.txt`

Is the output from making sure that all samples in the `targt_hla_gt_calls_no_header.vcf` are indentical and in the same order as the TGP data. If nothing is wrong with the samples or their order (and there is not) then the file will contain the following message:

```basic
All good Dave!
```

### `targt_ref_allele_validation.txt`

Is the output of checking that something did not go wrong when the TARGT pipeline was determining the reference allele, since both data sets were aligned to the `hg19` reference assembly then theoretically the intersect of sites between the TARGT and TGP data should have the same refernce allele. This file has four columns in the following format:

```basic
position	TGP_ref_allele	TARGT_ref_allele	agreement
```

where a value of `1` in the agreemnet column would represent that the reference alleles are the same between data sets and a value of `0` indicates that the refernce allele is different between the two data sets. There are six positions where the reference alleles are in disagrement, all due to a deletion event. I will illustrate this point by looking at position `31239150` using the following code:

```bash
zcat tgp_hla_ProtCod_GRCh37hg19_nonoverlapping_unfiltered.vcf.gz | awk '$2 == "31239150" { print }'
```

which will output the following entry (for brevity I will just show you the first nine columns of the VCF file):

```basic
6	31239150	rs373379283	CA	C	100	PASS	AC=179;AF=0.0357428;AN=5008;NS=2504;DP=7582;EAS_AF=0.001;AMR_AF=0.0274;AFR_AF=0.0817;EUR_AF=0.0447;SAS_AF=0.0061;AA=?|A|-|unsure;VT=INDEL;EX_TARGET	GT
```

The reference allele for this position is `CA` and the alternative allele is `C_` where `_` denotes the deletion of the `A` allele. In all six cases the TARGT pipeline called the deletion event the reference allele. We should dicuss how we want to adress this. The easiest route is a manual fix since it only effects six positions.





## Integrating the TARGT Classical HLA Calls into the TGP VCF File

To actually generate the new VCF file simply run:

```bash
python3 update_hla_vcf.py
```

which will create two files. The first file `unannotated_replaced_hla_exons.vcf` is simply an intemediary VCF file that can be deleted and the second file `annotated_replaced_hla_exons.vcf` is proably the reason why you are looking at this repo. Lastly, I will gzip and index `annotated_replaced_hla_exons.vcf` using Tabix:

```bash
bgzip annotated_replaced_hla_exons.vcf
tabix -p vcf annotated_replaced_hla_exons.vcf.gz
```

### `annotated_replaced_hla_exons.vcf.gz`

I added the following `INFO` flags to the header `annotated_replaced_hla_exons.vcf.gz`:

```basic
##INFO=<ID=TARGT,Number=0,Type=Flag,Description="indicates that the genotype call was produced by the TARGT pipeline"
##INFO=<ID=REPLACED_GT,Number=0,Type=Flag,Description="indicates that the TGP genotype call was replaced by TARGT pipeline genotype call"
##INFO=<ID=NEW_GT,Number=0,Type=Flag,Description="indicates a new genotype call introduced by TARGT pipeline genotype call"
```

I added these info flags such that we can subset `annotated_replaced_hla_exons.vcf.gz` in three unique ways using BCFtools.

Scenario 1: Generate a VCF file of all TARGT calls.

```bash
bcftools view -i 'INFO/TARGT=1' -Oz -o all_targt_calls.vcf.gz annotated_replaced_hla_exons.vcf.gz
```

Scenario 2: Generate a VCF file of only the new calls introudced by the TARGT pipeline.

```bash
bcftools view -i 'INFO/NEW_GT=1' -Oz -o only_new_targt_calls.vcf.gz annotated_replaced_hla_exons.vcf.gz
```

Scenario 3: Generate a VCF file that contains only sites that were replaced with TARGT calls.

```bash
bcftools view -i 'INFO/REPLACED_GT=1' -Oz -o only_replaced_targt_calls.vcf.gz annotated_replaced_hla_exons.vcf.gz
```





## TARGT Results

Using the filtering scheme described in scenario 1 we can use one of my favorite python packages [`scikit-allel`](scikit-allel) to analyze the TARGT calls we just introduced into our updated VCF file.

```python
### Dependicies ###
import allel
import numpy as np

### Preprocessing ###
# Load the TARGT VCF file.
targt_hla_vcf = 'all_targt_calls.vcf.gz'

# Load the data into scikit-allel.
targt_hla_a = allel.read_vcf(
        targt_hla_vcf,
        fields='calldata/GT',
        region='6:29909037-29913661',
        )
targt_hla_b = allel.read_vcf(
        targt_hla_vcf,
        fields='calldata/GT',
        region='6:31321649-31324965',
        )
targt_hla_c = allel.read_vcf(
        targt_hla_vcf,
        fields='calldata/GT',
        region='6:31236526-31239907',
        )
targt_hla_drb1 = allel.read_vcf(
        targt_hla_vcf,
        fields='calldata/GT',
        region='6:32546546-32557625',
        )
targt_hla_dqb1 = allel.read_vcf(
        targt_hla_vcf,
        fields='calldata/GT',
        region='6:32627244-32636160',
        )

# Convert the genotype information into a genotype matrix.
targt_hla_a_gt = allel.GenotypeArray(targt_hla_a['calldata/GT'])
targt_hla_b_gt = allel.GenotypeArray(targt_hla_b['calldata/GT'])
targt_hla_c_gt = allel.GenotypeArray(targt_hla_c['calldata/GT'])
targt_hla_drb1_gt = allel.GenotypeArray(targt_hla_drb1['calldata/GT'])
targt_hla_dqb1_gt = allel.GenotypeArray(targt_hla_dqb1['calldata/GT'])

### Analysis ###
# Count the total number of introduced sites per locus.
targt_hla_a_tot = targt_hla_a_gt.shape[0]
targt_hla_b_tot = targt_hla_b_gt.shape[0]
targt_hla_c_tot = targt_hla_c_gt.shape[0]
targt_hla_drb1_tot = targt_hla_drb1_gt.shape[0]
targt_hla_dqb1_tot = targt_hla_dqb1_gt.shape[0]
# Count the number of invariant sites per locus.
targt_hla_a_invar = (np.count_nonzero(targt_hla_a_gt.count_alleles(), axis=1) == 1).sum()
targt_hla_a_invar = (np.count_nonzero(targt_hla_a_gt.count_alleles(), axis=1) == 1).sum()
targt_hla_b_invar = (np.count_nonzero(targt_hla_b_gt.count_alleles(), axis=1) == 1).sum()
targt_hla_c_invar = (np.count_nonzero(targt_hla_c_gt.count_alleles(), axis=1) == 1).sum()
targt_hla_drb1_invar = (np.count_nonzero(targt_hla_drb1_gt.count_alleles(), axis=1) == 1).sum()
targt_hla_dqb1_invar = (np.count_nonzero(targt_hla_dqb1_gt.count_alleles(), axis=1) == 1).sum()
# Count the number of bi-allelic sites per locus.
targt_hla_a_bi = (np.count_nonzero(targt_hla_a_gt.count_alleles(), axis=1) == 2).sum()
targt_hla_a_bi = (np.count_nonzero(targt_hla_a_gt.count_alleles(), axis=1) == 2).sum()
targt_hla_b_bi = (np.count_nonzero(targt_hla_b_gt.count_alleles(), axis=1) == 2).sum()
targt_hla_c_bi = (np.count_nonzero(targt_hla_c_gt.count_alleles(), axis=1) == 2).sum()
targt_hla_drb1_bi = (np.count_nonzero(targt_hla_drb1_gt.count_alleles(), axis=1) == 2).sum()
targt_hla_dqb1_bi = (np.count_nonzero(targt_hla_dqb1_gt.count_alleles(), axis=1) == 2).sum()
# Count the number of multi-allelic sites per locus.
targt_hla_a_multi = (np.count_nonzero(targt_hla_a_gt.count_alleles(), axis=1) > 2).sum()
targt_hla_b_multi = (np.count_nonzero(targt_hla_b_gt.count_alleles(), axis=1) > 2).sum()
targt_hla_c_multi = (np.count_nonzero(targt_hla_c_gt.count_alleles(), axis=1) > 2).sum()
targt_hla_drb1_multi = (np.count_nonzero(targt_hla_drb1_gt.count_alleles(), axis=1) > 2).sum()
targt_hla_dqb1_multi = (np.count_nonzero(targt_hla_dqb1_gt.count_alleles(), axis=1) > 2).sum()
```

To make things a little bit easier I compiled the results into the following table:

|  Locus   | Total Sites | Invariant Sites | Bi-allelic Sites | Multi-allelic Sites |
| :------: | :---------: | :-------------: | :--------------: | :-----------------: |
|  HLA-A   |     787     |       366       |       342        |         79          |
|  HLA-B   |     791     |       689       |        83        |         19          |
|  HLA-C   |     796     |       736       |        49        |         11          |
| HLA-DRB1 |     270     |       178       |        67        |         25          |
| HLA-DQB1 |     270     |       210       |        51        |          9          |















