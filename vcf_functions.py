### Dependicies ###
import gzip
import numpy as np
import pandas as pd
import re


### Define Functions ###
def table_2_vcf(targt_table):
    """
    ###########################################################################
    INPUT: One csv produced by the TARGT pipeline formatted with the
           appropiate columns specified by the VCF format but with unformatted
           indivual genotype calls.
    ---------------------------------------------------------------------------
    OUTPUT: A pandas data frame in VCF format with indivual genotype calls
            properly formatted.
    ###########################################################################
    """
    # Intialize an empty list for later.
    call = []
    # Read in the csv file as a pandas dataframe.
    table = pd.read_csv(targt_table)
    # Convert the reference allele to 0 in all samples.
    for col in table.columns[9:]:
        table.loc[table[col] == table['REF'], col] = '0'
    # Convert all missing calls to '.'
    for col in table.columns[9:]:
        table.loc[table[col] == 'X', col] = '.'
    # Creat a list of every REF and ALT allele per site.
    for indv in range(table.index.size):
        call.append(table.iloc[indv, 9:].unique().tolist())
    # Remove 0s from polymorphic sites.
    for lists in call:
        if len(lists) == 1:
            lists
        else:
            lists[:] = (value for value in lists if value != '0')
    # Remove .s from polymorphic sites.
    for genos in call:
        if len(genos) == 1:
            genos
        else:
            genos[:] = (value for value in genos if value != '.')
    # Merge alternate alleles in multiallelic sites.
    for calls in call:
        if len(calls) == 1:
            calls
        else:
            multi = ','.join(calls)
            calls.append(multi)
    # Remove all extra elements from every list.
    for values in call:
        if len(values) == 1:
            values
        else:
            del values[:-1]
    # Flatten the alternate allele list.
    alt = [val[0] for val in call]
    # Convert the 0's to .'s for all invariant sites.
    alt_col = ['.' if site == '0' else site for site in alt]
    # Replace the nan values with the alternative alleles.
    table.ALT = alt_col
    # Convert alternate alleles to 1s, 2s, and 3s.
    alternates = table['ALT'].str.split(',', expand=True).rename(
            columns={0:'ALT1', 1:'ALT2', 2:'ALT3'})
    vcf_table = pd.concat([table, alternates], axis=1)
    # Ensure that missing calls aren't encoded as an alternate allele.
    vcf_table.loc[vcf_table['ALT1'] == '.', 'ALT1'] = None
    for col in vcf_table.columns[9:-3]:
        vcf_table.loc[vcf_table[col] == vcf_table['ALT1'], col] = '1'
        vcf_table.loc[vcf_table[col] == vcf_table['ALT2'], col] = '2'
        vcf_table.loc[vcf_table[col] == vcf_table['ALT3'], col] = '3'
    vcf_table.drop(['ALT1', 'ALT2', 'ALT3'], axis=1, inplace=True)
    # Merge individual genotypes into the proper format.
    for col in vcf_table.columns[9:]:
        if col[-2] != '.':
            vcf_table[col] = vcf_table[col].str.cat(
                    vcf_table[str(col)+'.1'], sep='|')
    # Remove unwanted columns.
    for col in vcf_table.columns:
        if col[-2] == '.':
            vcf_table.drop(col, axis=1, inplace=True)
    return vcf_table


def combine_vcf_tables(
        hla_a_targt_table,
        hla_c_targt_table,
        hla_b_targt_table,
        hla_drb1_targt_table,
        hla_dqb1_targt_table,
        ):
    """
    ###########################################################################
    INPUT: One csv per classical HLA locus produced by the TARGT pipeline
           formatted with the appropiate columns specified by the VCF format
           but with unformatted indivual genotype calls.
    ---------------------------------------------------------------------------
    OUTPUT: Saves a headerless classical HLA VCF file to your working directory.
    ###########################################################################
    """
    # Convert each TARGT table to a vcf-esque pandas data frame.
    hla_a_vcf_df = table_2_vcf(hla_a_targt_table)
    hla_c_vcf_df = table_2_vcf(hla_c_targt_table)
    hla_b_vcf_df = table_2_vcf(hla_b_targt_table)
    drb1_vcf_df = table_2_vcf(hla_drb1_targt_table)
    dqb1_vcf_df = table_2_vcf(hla_dqb1_targt_table)
    # Concatenate the individual vcf-esque pandas data frames into one data frame.
    targt_hla_vcf_df = pd.concat(
            [hla_a_vcf_df, hla_c_vcf_df, hla_b_vcf_df, drb1_vcf_df, dqb1_vcf_df],
            axis=0, sort=False, ignore_index=True,
            )
    # Save the combined TARGT data frame as a headerless vcf file.
    targt_hla_vcf_df.to_csv(
            'targt_hla_gt_calls_no_header.vcf',
            sep='\t',
            index=False,
            )
    return targt_hla_vcf_df


def missing_info_qc(targt_table):
    """
    ###########################################################################
    INPUT: One csv produced by the TARGT pipeline formatted with the
           appropiate columns specified by the VCF format but with unformatted
           indivual genotype calls.
    ---------------------------------------------------------------------------
    OUTPUT: A pandas data frames with the number of missing calls per site and
            the haplotypes that the missing call is from.
    ###########################################################################
    """
    # Intialize an empty list for the missing haplotypes.
    missing_haps = []
    # Read in the csv file as a pandas dataframe.
    qc_df = pd.read_csv(targt_table)
    # Create a dictionary from the pandas dataframe to parse.
    qc_df_dict = dict(list(qc_df.groupby(qc_df.index)))
    # For every row in the dataframe figure out which haplotypes having missing
    # genotype calls and append them to the 'missing_haps' list.
    for idx, df in qc_df_dict.items():
        inds = df.columns[(df == 'X').any()]
        missing_haps.append(','.join(inds.to_list()))
    # Idependently count how many missing genotype calls there are per site.
    qc_df['num_missing_calls'] = qc_df.apply(
            lambda row: sum(row[8:] == 'X'), axis=1,
            )
    # Add the identity of the missing genotype calls to the dataframe.
    qc_df['haps_with_missing_calls'] = missing_haps
    # Get rid of all the surperfulous data.
    qc_df = qc_df[['POS', 'num_missing_calls', 'haps_with_missing_calls']]
    return qc_df


def combine_missing_info_qcs(
        hla_a_targt_table,
        hla_c_targt_table,
        hla_b_targt_table,
        hla_drb1_targt_table,
        hla_dqb1_targt_table,
        ):
    """
    ###########################################################################
    INPUT: One csv per classical HLA locus produced by the TARGT pipeline
           formatted with the appropiate columns specified by the VCF format
           but with unformatted indivual genotype calls.
    ---------------------------------------------------------------------------
    OUTPUT: Saves a missing data qc report to your working directory.
    ###########################################################################
    """
    # Convert each TARGT table to a missing data qc pandas dataframe.
    hla_a_qc_df = missing_info_qc(hla_a_targt_table)
    hla_a_qc_df['locus'] = 'HLA-A'
    hla_c_qc_df = missing_info_qc(hla_c_targt_table)
    hla_c_qc_df['locus'] = 'HLA-C'
    hla_b_qc_df = missing_info_qc(hla_b_targt_table)
    hla_b_qc_df['locus'] = 'HLA-B'
    drb1_qc_df = missing_info_qc(hla_drb1_targt_table)
    drb1_qc_df['locus'] = 'HLA-DRB1'
    dqb1_qc_df = missing_info_qc(hla_dqb1_targt_table)
    dqb1_qc_df['locus'] = 'HLA-DQB1'
    # Concatenate the individual missing data qc pandas dataframes into one
    # data frame.
    targt_hla_qc_df = pd.concat(
            [hla_a_qc_df, hla_c_qc_df, hla_b_qc_df, drb1_qc_df, dqb1_qc_df],
            axis=0, sort=False, ignore_index=True,
            )
    # Save the combined missing data qc dataframe as a txt file.
    targt_hla_qc_df.to_csv(
            'targt_missing_data_qc.txt',
            sep='\t',
            index=False,
            )
    return targt_hla_qc_df


def tgp_site_info(tgp_vcf):
    """
    ###########################################################################
    INPUT: A gzipped VCF file from the TGP.
    ---------------------------------------------------------------------------
    OUTPUT
        header: Each line of the TGP VCF header (list)
        tgp_sample_array: The TGP samples in their proper order (array)
        tgp_dicc: The TGP VCF line information, reference allele, and ancestral
                  call for every site in the TGP VCF file (dictionary)
    ###########################################################################
    """
    # Iterate through every line in the TGF vcf file.
    with gzip.open(tgp_vcf, 'rt') as data:
        # Intialize a list for the header lines.
        header = []
        # Intialize dictionary.
        tgp_dicc = {}
        for line in data:
            # Append all the header lines to the header list.
            if line.startswith('##'):
                header.append(line)
                # Grab the line of the vcf that has the sample IDs.
            elif line.startswith('#'):
                header.append(line)
                spline = line.split()
                # Save the sample array for the validation step.
                tgp_sample_array = np.asarray(spline[9:])
            else:
                spline = line.split()
                # Grab the position.
                pos = spline[1]
                # Grab the refernce allele.
                ref = spline[3]
                # Grab the info field.
                info = spline[7]
                # If the site has information about the ancestral call grab it.
                if 'AA' in info:
                    anc_allele_search = re.search('AA=(\S+)VT', info)
                    anc_allele_group = anc_allele_search.group()
                    anc_allele_groups = anc_allele_group.split(';')
                    anc_allele = ';'+anc_allele_groups[0]
                    # If the site doesn't have informtaion about the ancestral
                    # call then insert a missing AA call.
                else:
                    anc_allele = ';AA=.|||'
                tgp_dicc[int(pos)] = {
                        'site': line,
                        'ref': ref,
                        'aa': anc_allele,
                        }
    return header, tgp_sample_array, tgp_dicc

def targt_site_info(targt_headerless_vcf):
    """
    ###########################################################################
    INPUT: A headerless classical HLA VCF file from the TARGT pipeline.
    ---------------------------------------------------------------------------
    OUTPUT
        targt_sample_array: The TGP samples in their observed order (array)
        targt_dicc: The TARGT VCF line information and reference allele for
                    every site in the TARGT VCF file (dictionary)
    ###########################################################################
    """
    # Iterate through every line in the TARGT vcf file.
    with open(targt_headerless_vcf, 'rt') as data:
        # Intialize dictionary.
        targt_dicc = {}
        for line in data:
            # Skip all the info lines
            if line.startswith('#'):
                spline = line.split()
                # Save the sample array for the validation step.
                targt_sample_array = np.asarray(spline[9:])
            else:
                spline = line.split()
                # Grab the position
                pos = spline[1]
                # Grab the refernce allele.
                ref = spline[3]
                targt_dicc[int(pos)] = {
                        'site': line,
                        'ref': ref,
                        }
    return targt_sample_array, targt_dicc

def validate_sample_indicies(tgp_sample_array, targt_sample_array):
    """
    ###########################################################################
    INPUT: Sample arrays from the TGP and TARGT VCF files.
    ---------------------------------------------------------------------------
    OUTPUT: A report indicating if the samples are in the correct order in the
            TARGT VCF file.
    ###########################################################################
    """
    # Configure and open the report file.
    report_file = open('targt_sample_validation.txt', 'w')
    # If all the sample IDs are the same and in the correct order tell me
    # things went well and if not tell me something went wrong.
    if (tgp_sample_array == targt_sample_array).all():
        report_file.write('All good Dave!')
    else:
        report_file.write('Something went wrong Dave...')
    report_file.close()
    return

def validate_ref_alleles(tgp_dicc, targt_dicc):
    """
    ###########################################################################
    INPUT: Site dictionaries from the TGP and TARGT VCF files.
    ---------------------------------------------------------------------------
    OUTPUT: A report indicating if the reference allele is concdorant between
            datasets (1) or discordant (2).
    ###########################################################################
    """
    # Configure and open the report file.
    report_file = open('targt_ref_allele_validation.txt', 'w')
    # Using each data set's dicctionary create an array of each position.
    tgp_pos = tgp_dicc.keys()
    targt_pos = targt_dicc.keys()
    tgp_sites = np.asarray(list(tgp_pos))
    targt_sites = np.asarray(list(targt_pos))
    # Determine the sites that occur in both data sets.
    common_sites = np.intersect1d(tgp_sites, targt_sites)
    # Validate that each data set contains the same reference allele.
    for site in np.sort(common_sites, axis=None):
        if tgp_dicc[site]['ref'] == targt_dicc[site]['ref']:
            report_file.write(
                    str(site)\
                    + '\t'\
                    + tgp_dicc[site]['ref']\
                    + '\t'\
                    + targt_dicc[site]['ref']\
                    + '\t'\
                    + '1'\
                    + '\n',
                    )
        else:
            report_file.write(
                    str(site)\
                    + '\t'\
                    + tgp_dicc[site]['ref']\
                    + '\t'\
                    + targt_dicc[site]['ref']\
                    + '\t'\
                    + '0'\
                    + '\n',
                    )
    report_file.close()
    return

def create_new_vcf(header, tgp_dicc, targt_dicc):
    """
    ###########################################################################
    INPUT
        header: Each line of the TGP VCF header (list)
        tgp_dicc: The TGP VCF line information, reference allele, and ancestral
                  call for every site in the TGP VCF file (dictionary)
        targt_dicc: The TARGT VCF line information and reference allele for
                    every site in the TARGT VCF file (dictionary)
    ---------------------------------------------------------------------------
    OUTPUT: Saves an unannotated TGP VCF file with updated genotype calls to
            your working directory.
    ###########################################################################
    """
    # Add the two additional INFO lines necessary at the correct indicies.
    info_1 = '##INFO=<ID=TARGT,Number=0,Type=Flag,Description="indicates that the genotype call was produced by the TARGT pipeline">\n'
    info_2 = '##INFO=<ID=REPLACED_GT,Number=0,Type=Flag,Description="indicates that the TGP genotype call was replaced by TARGT pipeline genotype call">\n'
    info_3 = '##INFO=<ID=NEW_GT,Number=0,Type=Flag,Description="indicates a new genotype call introduced by TARGT pipeline genotype call">\n'
    header.insert(-4, info_1)
    header.insert(-4, info_2)
    header.insert(-4, info_3)
    # Configure and open the new vcf file.
    new_vcf = open('unannotated_replaced_hla_exons.vcf', 'w')
    # Add the TGP header lines to the new vcf file.
    for line in header:
        new_vcf.write(line)
    # Using each data set's dicctionary create an array of each position.
    tgp_pos = tgp_dicc.keys()
    targt_pos = targt_dicc.keys()
    tgp_sites = np.asarray(list(tgp_pos))
    targt_sites = np.asarray(list(targt_pos))
    # Determine which sites are private to each dataset.
    tgp_private_sites = np.setdiff1d(tgp_sites, targt_sites)
    targt_private_sites = np.setdiff1d(targt_sites, tgp_sites)
    # Determine the sites that occur in both data sets.
    common_sites = np.intersect1d(tgp_sites, targt_sites)
    # Combine the three sets to iterate over.
    all_sites = np.concatenate(
            (tgp_private_sites, targt_private_sites, common_sites),
            axis=None,
            )
    # If the site is private to the TGP genotype calls copy over the site info
    # from the TGP vcf files, otherwise use the updated genotype calls from the
    # TARGT pipeline.
    for site in np.sort(all_sites, axis=None):
        if np.any(site == tgp_private_sites):
            new_vcf.write(tgp_dicc[site]['site'])
        else:
            new_vcf.write(targt_dicc[site]['site'])
    new_vcf.close()
    return

def annotate_new_vcf(new_vcf, tgp_dicc, targt_dicc):
    """
    ###########################################################################
    INPUT
        new_vcf: The unannotated TGP VCF file with updated genotype calls (vcf)
        tgp_dicc: The TGP VCF line information, reference allele, and ancestral
                  call for every site in the TGP VCF file (dictionary)
        targt_dicc: The TARGT VCF line information and reference allele for
                    every site in the TARGT VCF file (dictionary)
    ---------------------------------------------------------------------------
    OUTPUT: Saves an annotated TGP VCF file with updated genotype calls to
            your working directory.
    ###########################################################################
    """
    # Using each data set's dicctionary create an array of each position.
    tgp_pos = tgp_dicc.keys()
    targt_pos = targt_dicc.keys()
    tgp_sites = np.asarray(list(tgp_pos))
    targt_sites = np.asarray(list(targt_pos))
    # Determine which sites are private to each dataset.
    tgp_private_sites = np.setdiff1d(tgp_sites, targt_sites)
    targt_private_sites = np.setdiff1d(targt_sites, tgp_sites)
    # Configure and open the new annotated vcf file.
    annotated_vcf = open('annotated_replaced_hla_exons.vcf', 'w')
    # Iterate through every line in the new unannotated vcf file.
    with open(new_vcf, 'rt') as data:
        for line in data:
            # Write the vcf header and first line.
            if line.startswith('##') or line.startswith('#'):
                annotated_vcf.write(line)
            else:
                spline = line.split()
                # Grab the position
                pos = spline[1]
                # Add the ancestral call information to the sites private to
                # the TGP data set.
                if np.any(int(pos) == tgp_private_sites):
                    # Add the TGP ancestral info to the INFO field.
                    if 'AA' in spline[7]:
                        annotated_vcf.write(line)
                    else:
                        spline[7] = spline[7]+';AA=.|||'
                        new_line = '\t'.join(spline)+'\n'
                        annotated_vcf.write(new_line)
                elif np.any(int(pos) == targt_private_sites):
                    # Add the ancestral and call info to the INFO field for all
                    # sites in the TARGT data set.
                    spline[7] = spline[7]+';NEW_GT;AA=.|||'
                    new_line = '\t'.join(spline)+'\n'
                    annotated_vcf.write(new_line)
                else:
                    spline[7] = spline[7]+';REPLACED_GT'+tgp_dicc[int(pos)]['aa']
                    new_line = '\t'.join(spline)+'\n'
                    annotated_vcf.write(new_line)
    annotated_vcf.close()
    return