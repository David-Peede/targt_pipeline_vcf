### Dependicies ###
import gzip
import numpy as np
import pandas as pd
import re
import vcf_functions as vf


### Load Required Data ###
# TGP unfiltered protein coding HLA VCF file.
tgp_hla_vcf = 'tgp_hla_ProtCod_GRCh37hg19_nonoverlapping_unfiltered.vcf.gz'

# QC'ed headerless VCF file from the TARGT pipeline.
targt_hla_vcf_table = 'targt_hla_gt_calls_no_header.vcf'


### Generate the New HLA VCF File ###
# First we need to grab the header from the TGP VCF file and create site dictionaries for
# each data set.
tgp_header, _, tgp_dicc = vf.tgp_site_info(tgp_hla_vcf)
_, targt_dicc = vf.targt_site_info(targt_hla_vcf_table)

# Next we need to create the updated VCF file and save it to our working directory.
vf.create_new_vcf(tgp_header, tgp_dicc, targt_dicc)

# We now load the the unannotated VCF file from our working directory as a text file.
unannotated_vcf = 'unannotated_replaced_hla_exons.vcf'

# Lastly we annotate the the unannotated VCF file and save it to our working directory.
vf.annotate_new_vcf(unannotated_vcf, tgp_dicc, targt_dicc)

##########################################################################################
# Note: You will need to gzip the resulting annotated VCF file with bgzip if you would like
#       to index it using tabix.
##########################################################################################