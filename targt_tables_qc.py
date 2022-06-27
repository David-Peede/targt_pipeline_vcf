### Dependicies ###
import gzip
import numpy as np
import pandas as pd
import re
import vcf_functions as vf
import sys


### Load Required Data ###
# TARGT tables.
hla_a_table = f'{sys.argv[1]}full_HLASeqA.csv'
hla_c_table = f'{sys.argv[1]}full_HLASeqC.csv'
hla_b_table = f'{sys.argv[1]}full_HLASeqB.csv'
drb1_table  = f'{sys.argv[1]}full_HLASeqDRB1.csv'
dqb1_table  = f'{sys.argv[1]}full_HLASeqDQB1.csv'

# Original unfiltered protein coding HLA VCF file.
org_hla_vcf = f'{sys.argv[1]}_hla_a_c_b_drb1_dqb1_exons_unfiltered.vcf.gz'


### TARGT Table Conversion & QC ###
# First we convert the individual TARGT tables into one and output the resulting
# headerless VCF file to our working directory.
_ = vf.combine_vcf_tables(
        hla_a_table,
        hla_c_table,
        hla_b_table,
        drb1_table,
        dqb1_table,
        sys.argv[1],
        )

# We now load the the headerless VCF file from our working directory as a text file.
targt_hla_vcf_table = f'targt_{sys.argv[1]}_hla_gt_calls_no_header.vcf'

# Before we can QC we need to generate the sample arrays and site dictionaries for each
# data set.
org_header, org_samps, org_dicc = vf.org_site_info(org_hla_vcf)
targt_samps, targt_dicc = vf.targt_site_info(targt_hla_vcf_table)

# Make sure that the sample IDs are identical and in the correct order.
vf.validate_sample_indicies(org_samps, targt_samps, sys.argv[1])

# Make sure that the both data sets are using the same reference allele at each site.
vf.validate_ref_alleles(org_dicc, targt_dicc, sys.argv[1])

# Lastly we will track how many and which haplotypes are missing per site.
_ = vf.combine_missing_info_qcs(
        hla_a_table,
        hla_c_table,
        hla_b_table,
        drb1_table,
        dqb1_table,
        sys.argv[1],
        )

##########################################################################################
# Note: You will need to inspect the QC files prior to running update_hla_vcf.py
##########################################################################################
