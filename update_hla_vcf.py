### Dependicies ###
import gzip
import numpy as np
import pandas as pd
import re
import vcf_functions as vf
import sys


### Load Required Data ###
# Original unfiltered protein coding HLA VCF file.
org_hla_vcf = sys.argv[2]

# QC'ed headerless VCF file from the TARGT pipeline.
targt_hla_vcf_table = f'targt_{sys.argv[1]}_hla_gt_calls_no_header.vcf'


### Generate the New HLA VCF File ###
# First we need to grab the header from the original VCF file and create site
# dictionaries for each data set.
org_header, _, org_dicc = vf.org_site_info(org_hla_vcf)
_, targt_dicc = vf.targt_site_info(targt_hla_vcf_table)

# Next we need to create the updated VCF file and save it to our working directory.
vf.create_new_vcf(org_header, org_dicc, targt_dicc, sys.argv[1])

# We now load the the unannotated VCF file from our working directory as a text file.
unannotated_vcf = f'unannotated_{sys.argv[1]}_replaced_hla_exons.vcf'

# Lastly we annotate the the unannotated VCF file and save it to our working directory.
vf.annotate_new_vcf(unannotated_vcf, org_dicc, targt_dicc, sys.argv[1])

##########################################################################################
# Note: You will need to gzip the resulting annotated VCF file with bgzip if you would like
#       to index it using tabix.
##########################################################################################
