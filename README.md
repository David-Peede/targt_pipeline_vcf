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

For our purposes we are only interested in the exons for the classical HLA class 1 and class 2 genes, using [BCFtools](https://samtools.github.io/bcftools/bcftools.html) we can filter the TGP VCF file to only contain our regions of interest with the following code (you can find the `.bed` file in the `resources` directory):

```bash
bcftools view -R hla_a_c_b_drb1_dqb1_exons.bed -Oz -o tgp_hla_a_c_b_drb1_dqb1_exons_unfiltered.vcf.gz ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
```

Lastly, we can index `tgp_hla_a_c_b_drb1_dqb1_exons_unfiltered.vcf.gz` using [Tabix](http://www.htslib.org/doc/tabix.html) with the following code:

```bash
tabix -p vcf tgp_hla_a_c_b_drb1_dqb1_exons_unfiltered.vcf.gz
```





## TARGT Table Conversion and QC

I have preformatted the TARGT tables (which can be found in the `targt_tables` directory) to include the necessary header line information as specified by the [VCF documentation](https://samtools.github.io/hts-specs/VCFv4.3.pdf) here I make a couple of assumptions that should not effect down stream analysis:

* Since we called genotypes using the TARGT pipeline we do not have variant IDs associated with our new calls, so I annotate the `ID` field as missing using a `"."`

* Since the quality of each call was assed in the TARGT pipeline I annotate the filter statuts in the `FILTER` field as `PASS` for every site

In converting the TARGT genotype calls to the VCF format you will notice that there are some missing genotype calls. I always denote missing genotype information with a `"."` and I do not do any subsequent filtering of missing genotype calls. However the `targt_tables_qc.py` will produce a file that has site level annotations of which individual(s) are missing genotype information.

Now that we have the assumptions out of the way let's run:

```bash
python3 targt_tables_qc.py
```

which will output four files that we will now go over.

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
position		TGP_ref_allele		TARGT_ref_allele		agreement
```

where a value of `1` in the agreemnet column would represent that the reference allele strings are **IDENTICAL** between data sets and a value of `0` indicates that the refernce allele strings are **NOT IDENTICAL** between the two data sets. There are two positions (`31324493` and `31324496`) where the reference alleles are not encoded the same way, but the actual reference alleles are the same! This discordance is all due to how VCF files encode deletion events—see section 5 in the [VCF documentation](https://samtools.github.io/hts-specs/VCFv4.3.pdf) for a more thorough explanation of how VCF files encode different types of variants. I will illustrate this point by looking at our two positions using the following code:

```bash
zcat tgp_hla_a_c_b_drb1_dqb1_exons_unfiltered.vcf.gz | awk '$2 == "31324493" { print }'
zcat tgp_hla_a_c_b_drb1_dqb1_exons_unfiltered.vcf.gz | awk '$2 == "31324496" { print }'
```

which will output the following entry (for brevity I will just show you the first nine columns of the VCF file):

```basic
6	31324493	rs576010607	CA	C	100	PASS	AC=330;AF=0.0658946;AN=5008;NS=2504;DP=6820;EAS_AF=0.0913;AMR_AF=0.0418;AFR_AF=0.0348;EUR_AF=0.0805;SAS_AF=0.0838;AA=|||unknown(NO_COVERAGE);VT=INDEL;EX_TARGET	GT

6	31324496	rs540530530	GT	G	100	PASS	AC=330;AF=0.0658946;AN=5008;NS=2504;DP=6940;EAS_AF=0.0913;AMR_AF=0.0418;AFR_AF=0.0348;EUR_AF=0.0805;SAS_AF=0.0838;AA=|||unknown(NO_COVERAGE);VT=INDEL;EX_TARGET	GT
```

The reference allele for position `31324493` is `CA` and the alternative allele is `C_` where `_` denotes the deletion of the `A` allele. Similarly the reference allele for position `31324496` is `GT` and the alternative allele is `G_` where `_` denotes the deletion of the `T` allele. From section 5 in the [VCF documentation](https://samtools.github.io/hts-specs/VCFv4.3.pdf) it is then clear to see that the actual reference allele for these two positions are `C` and `G` respectively thus all positions are concordant, which is exactly what we would expect!

### `targt_missing_data_qc.txt`

Is the output of checking missing genotype calls per haplotype for every site produced by the TARGT pipeline. This file has four columns in the following format:

``` basic
POS		num_missing_calls		haps_with_missing_calls		locus
```

where `num_missing_calls` refers to the total number of chromosomes with missing genotype calls, `haps_with_missing_calls`refers to the comma seperated list of the specific haplotype(s) that the missing genotype call is on, and `locus` refers to the HLA gene that site is falls on. Since each loci varies in the number of missing genotype calls I do not filter the VCF file by missing information and allow the researcher to make that call—haha "call", get it? That pun was 100% unintentional.





## Integrating the TARGT Classical HLA Calls into the TGP VCF File

To actually generate the new VCF file simply run:

```bash
python3 update_hla_vcf.py
```

which will create two files. The first file `unannotated_replaced_hla_exons.vcf` is simply an intemediary VCF file that can be deleted and the second file `annotated_replaced_hla_exons.vcf` is proably the reason why you are looking at this repo. Lastly, I will bgzip and index `annotated_replaced_hla_exons.vcf` using Tabix:

```bash
bgzip annotated_replaced_hla_exons.vcf
tabix -p vcf annotated_replaced_hla_exons.vcf.gz
```

### `annotated_replaced_hla_exons.vcf.gz`

I added the following `INFO` flags to the header `annotated_replaced_hla_exons.vcf.gz` (which can be found in the `resources` directory along with its indexed file):

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
import vcf_functions.py as vf

# Load VCF files.
original_tgp_hla_vcf = 'tgp_hla_a_c_b_drb1_dqb1_exons_unfiltered.vcf.gz'
all_targt_hla_vcf = 'all_targt_calls.vcf.gz'
new_targt_hla_vcf = 'only_new_targt_calls.vcf.gz'
replaced_targt_hla_vcf = 'only_replaced_targt_calls.vcf.gz'

# Analyze the variant type data per locus.
vf.hla_exon_report(original_tgp_hla_vcf)
vf.hla_exon_report(all_targt_hla_vcf)
vf.hla_exon_report(new_targt_hla_vcf)
vf.hla_exon_report(replaced_targt_hla_vcf)

# Create site dictionary for the original TGP calls and the calls replaced by the TARGT pipeline.
original_tgp_hla_dicc = vf.vcf_site_info(original_tgp_hla_vcf)
replaced_targt_hla_dicc = vf.vcf_site_info(replaced_targt_hla_vcf)

# Analyze the the variant identities before and after replacement.
vf.compare_replaced_calls(original_tgp_hla_dicc, replaced_targt_hla_dicc)
```

So first we will look at the original calls for our HLA exons from `tgp_hla_a_c_b_drb1_dqb1_exons_unfiltered.vcf.gz` (__NOTE__: since the data was imputed the TGP VCF file only contains variable sites):

| Locus             |   Total Sites |   Invariant Sites |   Bi-Allelic Sites |   Multi-Allelic Sites |
|:-------------------:|:---------------:|:-------------------:|:--------------------:|:-----------------------:|
| HLA-A (Exon 1)    |            45 |                 0 |                 40 |                     5 |
| HLA-A (Exon 2)    |            37 |                 0 |                 28 |                     9 |
| HLA-B (Exon 1)    |            36 |                 0 |                 27 |                     9 |
| HLA-B (Exon 2)    |            57 |                 0 |                 48 |                     9 |
| HLA-C (Exon 1)    |            31 |                 0 |                 25 |                     6 |
| HLA-C (Exon 2)    |            25 |                 0 |                 23 |                     2 |
| HLA-DRB1 (Exon 1) |            24 |                 0 |                 24 |                     0 |
| HLA-DQB1 (Exon 1) |            33 |                 0 |                 33 |                     0 |

Next let's look at how the calls changed once we included the TARGT calls from `all_targt_calls.vcf.gz`:

| Locus             |   Total Sites |   Invariant Sites |   Bi-Allelic Sites |   Multi-Allelic Sites |
|:-------------------:|:---------------:|:-------------------:|:--------------------:|:-----------------------:|
| HLA-A (Exon 1)    |           270 |                59 |                172 |                    39 |
| HLA-A (Exon 2)    |           276 |                69 |                167 |                    40 |
| HLA-B (Exon 1)    |           276 |               229 |                 35 |                    12 |
| HLA-B (Exon 2)    |           270 |               216 |                 47 |                     7 |
| HLA-C (Exon 1)    |           276 |               240 |                 26 |                    10 |
| HLA-C (Exon 2)    |           270 |               246 |                 23 |                     1 |
| HLA-DRB1 (Exon 1) |           270 |               177 |                 68 |                    25 |
| HLA-DQB1 (Exon 1) |           270 |               210 |                 51 |                     9 |

So for the HLA-A locus the TARGT calls introduced quite a few informative polymorphic sites and we introduced a modest amount of new informative polymorphic sites for the HLA-DRB1 and HLA-DQB1 loci! However, for the HLA-B and HLA-C loci the majority of the introduced TARGT calls are invariant. Now let's look at the relative contributions of new TARGT calls vs replaced TARGT calls. First let's look at the breakdown of new calls introduced by the TARGT pipeline from `only_new_targt_calls.vcf.gz`:

| Locus             |   Total Sites |   Invariant Sites |   Bi-Allelic Sites |   Multi-Allelic Sites |
|:-------------------:|:---------------:|:-------------------:|:--------------------:|:-----------------------:|
| HLA-A (Exon 1)    |           225 |                59 |                152 |                    14 |
| HLA-A (Exon 2)    |           239 |                68 |                153 |                    18 |
| HLA-B (Exon 1)    |           240 |               226 |                 13 |                     1 |
| HLA-B (Exon 2)    |           213 |               212 |                  1 |                     0 |
| HLA-C (Exon 1)    |           245 |               233 |                 12 |                     0 |
| HLA-C (Exon 2)    |           245 |               243 |                  2 |                     0 |
| HLA-DRB1 (Exon 1) |           246 |               177 |                 48 |                    21 |
| HLA-DQB1 (Exon 1) |           237 |               209 |                 19 |                     9 |

So the majority of new calls introduced by the TARGT pipeline appear to be invariant, however, for the HLA-A and HLA-DRB1 loci it looks like the majority of informative polymorphic sites are new calls! Let's compare this to the calls from the TARGT pipeline that replaced the original HLA calls from `only_replaced_targt_calls.vcf.gz`:

| Locus             |   Total Sites |   Invariant Sites |   Bi-Allelic Sites |   Multi-Allelic Sites |
|:-------------------:|:---------------:|:-------------------:|:--------------------:|:-----------------------:|
| HLA-A (Exon 1)    |            45 |                 0 |                 20 |                    25 |
| HLA-A (Exon 2)    |            37 |                 1 |                 14 |                    22 |
| HLA-B (Exon 1)    |            36 |                 3 |                 22 |                    11 |
| HLA-B (Exon 2)    |            57 |                 4 |                 46 |                     7 |
| HLA-C (Exon 1)    |            31 |                 7 |                 14 |                    10 |
| HLA-C (Exon 2)    |            25 |                 3 |                 21 |                     1 |
| HLA-DRB1 (Exon 1) |            24 |                 0 |                 20 |                     4 |
| HLA-DQB1 (Exon 1) |            33 |                 1 |                 32 |                     0 |

Interestingly, every single original position was replaced with TARGT calls! However, if we look at only the replaced calls we are actually losing informative information since for every loci we are not only introducing bi-allelic site but also multi-allelic and invariant sites. Lastly, lets look at the identity of the calls we replaced, which is the output of `vf.compare_replaced_calls(original_tgp_hla_dicc, replaced_targt_hla_dicc)`:

| Locus             |   Identical Variant |   Bi-allelic -> Invariant |   Bi-allelic -> Bi-allelic |   Bi-allelic -> Multi-allelic |   Multi-allelic -> Invariant |   Multi-allelic -> Bi-allelic |   Multi-allelic -> Multi-allelic |
|:-------------------:|:---------------------:|:---------------------------:|:----------------------------:|:-------------------------------:|:------------------------------:|:-------------------------------:|:----------------------------------:|
| HLA-A (Exon 1)    |                  17 |                         0 |                          4 |                            20 |                            0 |                             0 |                                4 |
| HLA-A (Exon 2)    |                  14 |                         1 |                          2 |                            13 |                            0 |                             0 |                                7 |
| HLA-B (Exon 1)    |                  21 |                         3 |                          0 |                             3 |                            1 |                             0 |                                8 |
| HLA-B (Exon 2)    |                  45 |                         4 |                          0 |                             0 |                            2 |                             0 |                                6 |
| HLA-C (Exon 1)    |                  16 |                         7 |                          0 |                             4 |                            0 |                             0 |                                4 |
| HLA-C (Exon 2)    |                  20 |                         3 |                          0 |                             0 |                            1 |                             0 |                                1 |
| HLA-DRB1 (Exon 1) |                  20 |                         0 |                          0 |                             4 |                            0 |                             0 |                                0 |
| HLA-DQB1 (Exon 1) |                  32 |                         1 |                          0 |                             0 |                            0 |                             0 |                                0 |













