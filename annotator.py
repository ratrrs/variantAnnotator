'''variantAnnotator.py required functions'''

import requests
import json
import pandas as pd
import numpy as np


# function to find header index
# of a vcf file
def find_header(filename):
    infile = open(filename, 'r')

    for i, line in enumerate(infile):
        if line[:6] == '#CHROM':
            break

    return(i)


# function to parse metadata lines of VCF
# and return as pandas dataframe with columns
# 'Category', 'ID', 'Number', 'Type', 'Description'
# NOTE: not used directly in variantAnnotator.py but
#       useful for parsing vcf fields
def get_metadata(filename):
    infile = open(filename, 'r')

    # dataframe for storing metadata
    df = pd.DataFrame(columns = ['Category', 'ID', 'Number', 'Type', 'Description'])

    # we will parse out the fields in each
    # line of metadata, store in metadata_dict
    # and then append to df
    for line in infile:
        if line[:6] == '#CHROM':
            break
        else:
            if '=<' in line:
                line = line.strip('#>\n').split('=<')
            else:
                continue
        if line[0] != 'contig':
            metadata_dict = dict()
            metadata_category = line[0]
            metadata_dict['Category'] = metadata_category
            metadata_info = line[1]

            # 'Description' is a weird field in the metadata
            # bc it has quotes so split on it and deal with seperately
            metadata_desc = metadata_info.split(',Description=')[-1]
            metadata_dict['Description'] = metadata_desc 

            # parse out rest of the fields in metadata
            metadata_nondesc = metadata_info.split(',Description=')[0]
            metadata_nondesc = metadata_nondesc.split(',')

            for field in metadata_nondesc:
                key = field.split('=')[0]
                value = field.split('=')[1]
                metadata_dict[key] = value

            df = df.append(metadata_dict, ignore_index=True)

    return(df)


# helper function to get ccds_ids
# that have nucleotide sequence data
def ccds_ids_w_nuc(CCDS_nucleotide_filename):
    infile_nuc = open(CCDS_nucleotide_filename, 'r')
    ccds_ids = set()

    for line in infile_nuc:
        if line[0] == '>':
            line = line.strip('>\n').split('|')
            ccds_ids.add(line[0])

    infile_nuc.close()

    return(ccds_ids)


# chromosome specific function to return a dictionary
# with keys 'gene' and 'exon'
# 'gene': a set of all positions (intron+exon) within all CCDS
#         genes for a given chromosome. NOTE: because based on
#         cds sequence, will miss 5' and 3' UTR exon annotations.
# 'exon': a dictionary where the key is the position for all
#         cds regions within a given chromosome and the
#         value is the ccds_id (which is a unique ID for each gene)
def chrom_CDS_pos(chromosome, CCDS_coords_filename, CCDS_nucleotide_filename):
    df_coords = pd.read_csv(CCDS_coords_filename, sep='\t', header=0)
    df_coords = df_coords[df_coords['#chromosome'] == chromosome]
    ccds_ids = ccds_ids_w_nuc(CCDS_nucleotide_filename)
    df_coords = df_coords[df_coords.ccds_id.isin(ccds_ids)]
    gene_positions = []
    exon_positions = {}

    for index, row in df_coords.iterrows():
        gene_coords_start = int(row.cds_from) + 1  # need to add 1 bc CCDS format is 0-based
        gene_coords_end = int(row.cds_to) + 2  # see above
        gene_positions.extend(range(gene_coords_start, gene_coords_end))
        cds_coords = row.cds_locations
        cds_coords = cds_coords.strip('[]').split(', ')

        for i in cds_coords:
            i = i.strip().split('-')
            start = int(i[0]) + 1  # need to add 1 bc CCDS format is 0-based
            stop = int(i[1]) + 2  # see above

            for pos in range(start, stop):
                exon_positions[pos] = row.ccds_id

    gene_positions = set(gene_positions)
    positions_dict = {'gene': gene_positions,
                      'exon': exon_positions}

    return(positions_dict)


# helper function to get nucleotide/protein sequence
# for each ccds_id. will return a dictionary where
# key is ccds_id and value is nucleotide/protein sequence
def ccds_ids_nuc_prot_sequence(CCDS_nucleotide_protein_filename):
    infile_nuc_prot = open(CCDS_nucleotide_protein_filename, 'r')
    nuc_prot_dict = {}
    nuc_prot_seq = ''  # to initiate, we will remove at the end
    ccds_id = 'NULL'  # we will remove at the end

    for line in infile_nuc_prot:
        if line[0] == '>':
            nuc_prot_dict[ccds_id] = nuc_prot_seq
            nuc_prot_seq = ''
            line = line.strip('>\n').split('|')
            ccds_id = line[0]
        else:
            line = line.strip()
            nuc_prot_seq = nuc_prot_seq + line

    infile_nuc_prot.close()
    nuc_prot_dict[ccds_id] = nuc_prot_seq  # last ccds_id and nucleotide/protein sequence
    nuc_prot_dict.pop('NULL')

    return(nuc_prot_dict)


# helper function to obtain the reverse
# complement of a nucleotide sequence
def reverse_complement(sequence):
    complement_dict = {'A': 'T',
                       'C': 'G',
                       'T': 'A',
                       'G': 'C'}
    
    sequence_complement = \
        ''.join([complement_dict[nucleotide] for nucleotide in sequence])

    return(sequence_complement)


# function to create dictionary mapping ccds_id position
# and nucleotide sequence. will take reverse complement and
# reversed sequence order if cds is on '-' strand.
# returns dictionary of dictionaries, where the primary key is
# ccds_id and the value is another dictionary with position as key and
# nucleotide as value (e.g. {'CCDS_ID': {1: 'A'}})
def ccds_pos_nuc_sequence(CCDS_coords_filename, CCDS_nucleotide_filename):
    df_coords = pd.read_csv(CCDS_coords_filename, sep='\t', header=0)
    ccds_ids_nuc_seq = ccds_ids_nuc_prot_sequence(CCDS_nucleotide_filename)
    df_coords = df_coords[df_coords.ccds_id.isin(ccds_ids_nuc_seq)]
    position_nuc_dict = {}

    for index, row in df_coords.iterrows():
        if row['#chromosome'] == 'Y':  # NOTE: ignore chrY for now
            continue
        else:
            ccds_id = row.ccds_id
            cds_coords = row.cds_locations
            cds_coords = cds_coords.strip('[]').split(', ')
            position_nuc_dict[ccds_id] = {}
            strand = row.cds_strand
            positions = []
            nuc_index = -1

            # if negative strand, need to get reverse
            # complement and then reverse the order
            if strand == '-':
                ccds_ids_nuc_seq[ccds_id] = reverse_complement(ccds_ids_nuc_seq[ccds_id])
                ccds_ids_nuc_seq[ccds_id] = ccds_ids_nuc_seq[ccds_id][::-1]

            for i in cds_coords:
                i = i.strip().split('-')
                start = int(i[0]) + 1  # need to add 1 bc CCDS format is 0-based
                stop = int(i[1]) + 2  # see above
                positions.extend(list(range(start, stop)))

            # double check number positions is same as number of nucleotides
            assert(len(positions) == len(ccds_ids_nuc_seq[ccds_id]))

            for pos in positions:
                nuc_index += 1
                position_nuc_dict[ccds_id][pos] = ccds_ids_nuc_seq[ccds_id][nuc_index]

            # double check keys (positions) are sorted numerically since other operations may assume this
            assert(list(position_nuc_dict[ccds_id].keys()) == sorted(position_nuc_dict[ccds_id].keys()))

    return(position_nuc_dict)


# function to return ending position on reference sequence
# from a sequence alignment given CIGAR string and starting
# reference position.
# NOTE: following notations consume reference positions
# (see https://samtools.github.io/hts-specs/SAMv1.pdf):
#   M: match
#   D: reference deletion
#   N: skipped region
#   =: match
#   X: reference mismatch
def ending_CIGAR_pos(CIGARstring, start_pos):
    poslist = []
    add_pos = ''
    end_pos = start_pos

    for i in CIGARstring:
        if i.isdigit():
            add_pos = add_pos + i
        else:
            if i in 'MDN=X':  # consume reference
                end_pos += int(add_pos)
                add_pos = ''
            else:
                add_pos = ''

    return(end_pos-1)  # subtract 1 to make position closed/inclusive


# get variant type(s) according to CIGAR string.
# variant type options are: insertion, deletion, or substitution
# NOTE: if multiple variant types exist, then list of all
# variant types will be returned.
def CIGAR2variant_type(CIGARstring):
    CIGARvariants = ''.join([i for i in CIGARstring if not i.isdigit()])
    CIGARvariants = ''.join(set(CIGARvariants))
    CIGARvariants = CIGARvariants.replace('M', '')  # ignore matching variants

    variants_dict = {'I': 'insertion',
                     'D': 'deletion',
                     'X': 'substitution'}

    variants_list = [variants_dict[i] for i in CIGARvariants]
    
    if not variants_list:
        return(['None'])
    else:
        return(variants_list)


# helper function to parse fields of INFO
# column of vcf and return as dictionary
def parse_INFO(INFO):
    INFO = INFO.split(';')
    INFOdict = {i.split('=')[0]: i.split('=')[1] for i in INFO}

    return(INFOdict)


# helper function to parse fields of FORMAT
# for each sample in a vcf (by row) and return as
# pd.DataFrame with samples as rows and FORMAT fields as columns
def parse_FORMAT(row):
    dflist = []
    FORMAT_fields = row.FORMAT.split(':')

    for sample_id in row[9:].index:
        sample_data = row[sample_id].split(':')
        sample_dict = {i: j for i,j in zip(FORMAT_fields, sample_data)}
        dflist.append(sample_dict)
    
    df = pd.DataFrame.from_dict(dflist, orient='columns')

    return(df)


# function to return new cds sequence given CIGAR string,
# reference sequence, alternate sequence, and starting position
# NOTE on required arguments:
#   CIGARstring is CIGAR from INFO column of vcf
#   REFsequence is output of ccds_pos_nuc_sequence()
#   ALTsequence is ALT column of vcf
#   start_pos is POS column of vcf
def alt_cds_sequence(CIGARstring, REFsequence, ALTsequence, start_pos):
    end_pos = ending_CIGAR_pos(CIGARstring, start_pos)
    reference_positions = np.asarray(list(REFsequence.keys()))
    reference_sequence = list(REFsequence.values())

    # bc we are changing nucleotide sequence of cds, positions are not consecutive
    # (due to intron removal), therefore adjust ending position of alternate sequence
    # so that it is no greater than last position of current exon
    updated_end_pos = reference_positions[reference_positions <= end_pos].max()
    slice_length = updated_end_pos - start_pos + 1
    relative_pos_start = list(np.where(reference_positions == start_pos))[0][0]
    alternative_sequence = [i for i in ALTsequence]

    if end_pos == updated_end_pos:
        reference_sequence[relative_pos_start:relative_pos_start + slice_length] = alternative_sequence
    else:
        reference_sequence[relative_pos_start:relative_pos_start + slice_length] = alternative_sequence[:slice_length]

    updated_sequence = ''.join(reference_sequence)

    return(updated_sequence)


# function to return strand type ('+' or '-') for each ccds_id
def get_strand_type(CCDS_coords_filename):
    df_coords = pd.read_csv(CCDS_coords_filename, sep='\t', header=0)
    strand_dict = {}

    for index, row in df_coords.iterrows():
        if row['#chromosome'] == 'Y':  # NOTE: ignore chrY for now
            continue
        else:
            ccds_id = row.ccds_id
            strand = row.cds_strand
            strand_dict[ccds_id] = strand

    return(strand_dict)


# function to convert nucleotide sequence to
# AA sequence using a dictionary of codons with mappings
# to AA residues (single letter code)
def nuc_seq_to_prot_seq(nucleotide_sequence, codons_dict):
    prot_seq = ''

    for nuc in range(0, len(nucleotide_sequence), 3):
        codon = nucleotide_sequence[nuc:nuc+3]
        if len(codon) != 3:
            break
        else:
            AA = codons_dict[codon]
            prot_seq += AA

    return(prot_seq)


# obtain what kind of mutation resulted from
# nucleotide change 
# (i.e., nonsense, missense, silent)
def get_mutation_type(AA1, AA2):
    if AA1 != AA2:
        if AA2 == 'O':
            return('nonsense')
        else:
            return('missense')
    else:
        return('silent')


# function to return effect type of an amino acid change to a
# protein sequence (i.e., 'silent', 'missense', 'nonsense', or 'frameshift').
# NOTE: in this function, if we see more than two differences in the protein
# sequence and the first is not a stop codon (i.e., 'nonsense' change), we
# will assign the effect as 'frameshift'
def cds_prot_effect_type(reference_prot_seq, alternate_prot_seq):
    if reference_prot_seq == alternate_prot_seq:
        mutation = 'silent'
    else:
        mutation_types = []

        for ref, alt in zip(reference_prot_seq, alternate_prot_seq):
            if ref == alt:
                continue
            else:
                mutation_types.append(get_mutation_type(ref, alt))

        if mutation_types[0] == 'nonsense':
            mutation = 'nonsense'
        else:
            if len(mutation_types) == 1:
                mutation = mutation_types[0]
            else:
                mutation = 'frameshift'

    return(mutation)


# function to parse following sequence stats from vcf row:
# 1) total read depth
# 2) number of reads supporting reference variant
# 3) number of reads supporting alternate variant(s)
# 4) % of reads supporting reference variant
# 5) % of reads supporting alternate variants(s)
# 6) alternate variant(s) type(s)
# returns above data as a dictionary
def seq_stats(row):
    dictINFO = parse_INFO(row.INFO)
    dfFORMAT = parse_FORMAT(row)
    CIGAR = dictINFO['CIGAR'].split(',')
    variant_type = []

    for i in CIGAR:
        variant_type.append(CIGAR2variant_type(i))
    
    for i, j in enumerate(variant_type):
        variant_type[i] = ','.join(variant_type[i])

    variant_type = ';'.join(variant_type)
    total_reads = float(dictINFO['DP'])
    reference_reads = 0
    alternate_reads = len(row.ALT.split(',')) * [0]

    for index, FORMATrow in dfFORMAT.iterrows():
        reference_reads += int(FORMATrow.RO)
        alternate_reads_list = FORMATrow.AO.split(',')
        for i, j in enumerate(alternate_reads_list):
            alternate_reads[i] += int(j)
    
    reference_reads_percent = reference_reads/total_reads*100
    reference_reads_percent = round(reference_reads_percent, 2)
    alternate_reads_percent = list(np.asarray(alternate_reads)/total_reads*100)

    for i, j in enumerate(alternate_reads_percent):
        alternate_reads_percent[i] = round(alternate_reads_percent[i], 2)

    reference_reads = str(reference_reads)
    reference_reads_percent = str(reference_reads_percent)
    alternate_reads = ';'.join([str(i) for i in alternate_reads])
    alternate_reads_percent = ';'.join([str(i) for i in alternate_reads_percent])
    total_read_depth = dictINFO['DP']

    data_dict = {'variant_type': variant_type,
                 'total_read_depth': total_read_depth,
                 'total_reference_reads': reference_reads,
                 'percent_reference_reads(%)': reference_reads_percent,
                 'total_alternate_reads': alternate_reads,
                 'percent_alternate_reads(%)': alternate_reads_percent
    }

    return(data_dict)


# use python requests module and ExAC API (http://exac.hms.harvard.edu)
# to get frequency of alleles in Exome Aggregation Consortium (ExAC) database
# function takes in a pandas dataframe vcf as input and returns updated dataframe
# with additional ExAC frequencies column
def request_ExAC_frequency(df):
    # ExAC API can take in a bulk list of variants
    positions_variants_list = []

    # generate population bulk list by iterating through vcf
    for index, row in df.iterrows():
        CHROM = row['#CHROM']

        for var in row.ALT.split(','):
            variant_string = f'{CHROM}-{row.POS}-{row.REF}-{var}'
            positions_variants_list.append(variant_string)

    r = requests.post(url='http://exac.hms.harvard.edu/rest/bulk/variant', json=positions_variants_list)
    results = json.loads(r.content)
    list_of_dicts = []

    for index, row in df.iterrows():
        CHROM = row['#CHROM']
        alt_freqs = []

        for var in row.ALT.split(','):
            key = f'{CHROM}-{row.POS}-{row.REF}-{var}'
            if 'allele_freq' in results[key]['variant']:
                freq = results[key]['variant']['allele_freq']
                alt_freqs.append(str(freq))
            else:
                alt_freqs.append('None')

        alt_freqs = ';'.join(alt_freqs)
        rowdict = {'#CHROM': CHROM,
                   'POS': row.POS,
                   'REF': row.REF,
                   'ALT': row.ALT,
                   'ExAC_allele_frequency': alt_freqs
                   }
        list_of_dicts.append(rowdict)

    df_freq = pd.DataFrame(list_of_dicts)
    newdf = pd.merge(df, df_freq,  how='left', left_on=['#CHROM','POS', 'REF', 'ALT'], right_on = ['#CHROM','POS', 'REF', 'ALT'])

    return(newdf)