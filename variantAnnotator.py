''' "variantAnnotator"
    Author: Raul Torres (ratrrs@gmail.com)
    Date: 05-05-2021
    developed with python 3.9 '''

import argparse
import pandas as pd
import annotator as an
from constants import CODONS


# set up arguments so we can pass a VCF file and other required files to the script
parser = argparse.ArgumentParser(description='Annotate variants from a VCF file')
parser.add_argument('-V', '--VCFfile', type=str, metavar='', required=True, help='VCF file to annotate')
parser.add_argument('-C', '--CCDS_coordinates_file', type=str, metavar='', required=True, help='Local CCDS file (uncompressed) of coordinates of cds. \
    (e.g., https://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs37.3/CCDS.current.txt)')
parser.add_argument('-N', '--CCDS_nucleotide_file', type=str, metavar='', required=True, help='Local CCDS file (uncompressed) of cds nucleotide sequences. \
    (e.g., https://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs37.3/CCDS_nucleotide.current.fna.gz)')
parser.add_argument('-P', '--CCDS_protein_file', type=str, metavar='', required=True, help='Local CCDS file (uncompressed) of cds protein sequences. \
    (e.g., https://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs37.3/CCDS_protein.current.faa.gz)')
parser.add_argument('-o', '--out', type=str, metavar='', required=True, help='filename for output tsv file')
args = parser.parse_args()

# assign arguments to variables
VCF_filename = args.VCFfile  # VCF file to annotate
CCDS_coords_filename = args.CCDS_coordinates_file
CCDS_nucleotide_filename = args.CCDS_nucleotide_file
CCDS_protein_filename = args.CCDS_protein_file
out = args.out


# PRE-PROCESSING STEP
# 1) Load VCF
# 2) Load nucleotide positions and sequence for each cds (by ccds_id)
# 3) Load strand types for each ccds_id
# 4) Load protein sequences for each ccds_id

print('pre-processing and loading required files...')

header = an.find_header(VCF_filename)

df_vcf = pd.read_csv(VCF_filename, sep='\t', header=header)

ccds_id_position_and_nucleotide = an.ccds_pos_nuc_sequence(CCDS_coords_filename, CCDS_nucleotide_filename)

ccds_id_strand_type = an.get_strand_type(CCDS_coords_filename)

prot_seq = an.ccds_ids_nuc_prot_sequence(CCDS_protein_filename)


# PROCESS VCF STEP
# 1) we will create a new dataframe with annotations and sequencing
#    stats for each variant in our vcf by making a list of
#    dictionaries via iterating through df_vcf and combining that
#    into a pandas dataframe once we are done iterating.
# 2) we will then query variant frequencies from  ExAC API
#   (http://exac.hms.harvard.edu/) and add this to our
#    newly created pandas dataframe from step 1).

list_of_dicts = []

CHROM = '0'  # 'dummy' chromosome to initiate

for index, row in df_vcf.iterrows():

    if row['#CHROM'] != CHROM:
        CHROM = row['#CHROM']
        print(f'processing annotations for chromosome {CHROM}...')
        exon_gene_ccds_dict = an.chrom_CDS_pos(CHROM, CCDS_coords_filename, CCDS_nucleotide_filename)

    # parse sequencing stats (read depth, read %'s, variant types, etc.)
    # and start dictionary to store annotations and data
    rowdict = an.seq_stats(row)
    rowdict['#CHROM'] = CHROM
    rowdict['POS'] = row.POS
    rowdict['REF'] = row.REF
    rowdict['ALT'] = row.ALT
    rowdict['mutation_effect_type'] = []

    # get mutational effect type of alternate variant(s)
    if row.POS in exon_gene_ccds_dict['gene']:
        if row.POS in exon_gene_ccds_dict['exon']:
            INFO = an.parse_INFO(row.INFO)
            CIGAR = INFO['CIGAR'].split(',')
            ccds_id = exon_gene_ccds_dict['exon'][row.POS]

            # traverse through alternate sequence(s)
            # and record effect type of the variant(s)
            for i, ALTseq in enumerate(row.ALT.split(',')):
                alternate_nuc_seq = an.alt_cds_sequence(CIGAR[i], ccds_id_position_and_nucleotide[ccds_id], ALTseq, row.POS)
                if alternate_nuc_seq == ''.join(ccds_id_position_and_nucleotide[ccds_id].values()):

                    # variant(s) led to no sequence change, so mutation_type is simply 'None(exonic)'
                    rowdict['mutation_effect_type'].append('None(exonic)')
                else:

                    # generate reverse complement and reverse strand if '-' strand
                    if ccds_id_strand_type[ccds_id] == '-':
                        alternate_nuc_seq = an.reverse_complement(alternate_nuc_seq)[::-1]
                    alternate_prot_seq = an.nuc_seq_to_prot_seq(alternate_nuc_seq, CODONS)

                    # NOTE: need to add stop codon to all reference protein sequences since CCDS reference doesn't include
                    reference_prot_seq = prot_seq[ccds_id] + 'O'

                    # get mutational effect of alternate protein sequence
                    rowdict['mutation_effect_type'].append(an.cds_prot_effect_type(reference_prot_seq, alternate_prot_seq))
        else:
            rowdict['mutation_effect_type'].append('intronic')
    else:
        rowdict['mutation_effect_type'].append('intergenic')

    rowdict['mutation_effect_type'] = ';'.join(rowdict['mutation_effect_type'])

    list_of_dicts.append(rowdict)

df = pd.DataFrame(list_of_dicts)

# update df by adding frequencies from ExAC API (http://exac.hms.harvard.edu/)
print('adding variant frequencies from ExAC...')
df = an.request_ExAC_frequency(df)

# reorder columns
column_order = ['#CHROM',
                'POS',
                'REF',
                'ALT',
                'variant_type',
                'mutation_effect_type',
                'total_read_depth',
                'total_reference_reads',
                'percent_reference_reads(%)',
                'total_alternate_reads',
                'percent_alternate_reads(%)',
                'ExAC_allele_frequency']

df = df[column_order]

# substitute ',' for ';' in ALT column so all columns are
# ';'-delmited when multiple alternative variants exist
df['ALT'] = df['ALT'].str.replace(',', ';')

# write out results to tsv
df.to_csv(out, sep='\t', index=False)