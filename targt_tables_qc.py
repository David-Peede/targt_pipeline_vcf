### Dependicies ###
import gzip
import numpy as np
import pandas as pd
import re
import vcf_functions as vf


### Load Required Data ###
# TARGT tables.
hla_a_table = '1Kfull_HLASeqA.csv'
hla_c_table = '1Kfull_HLASeqC.csv'
hla_b_table = '1Kfull_HLASeqB.csv'
drb1_table = '1Kfull_HLASeqDRB1.csv'
dqb1_table = '1Kfull_HLASeqDQB1.csv'

# TGP unfiltered protein coding HLA VCF file.
tgp_hla_vcf = 'tgp_hla_ProtCod_GRCh37hg19_nonoverlapping_unfiltered.vcf.gz'


### TARGT Table Conversion & QC ###
# First we convert the individual TARGT tables into one and output the resulting
# headerless VCF file to our working directory.
_ = vf.combine_vcf_tables(
        hla_a_table,
        hla_c_table,
        hla_b_table,
        drb1_table,
        dqb1_table,
        )

# We now load the the headerless VCF file from our working directory as a text file.
targt_hla_vcf_table = 'targt_hla_gt_calls_no_header.vcf'

# Before we can QC we need to generate the sample arrays and site dictionaries for each
# data set.
tgp_header, tgp_samps, tgp_dicc = vf.tgp_site_info(tgp_hla_vcf)
targt_samps, targt_dicc = vf.targt_site_info(targt_hla_vcf_table)

# Make sure that the sample IDs are identical and in the correct order.
vf.validate_sample_indicies(tgp_samps, targt_samps)

# Make sure that the both data sets are using the same reference allele at each site.
vf.validate_ref_alleles(tgp_dicc, targt_dicc)

##########################################################################################
# Note: You will need to inspect the QC files prior to running update_hla_vcf.py
##########################################################################################