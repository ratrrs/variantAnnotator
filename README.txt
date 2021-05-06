README

Scripts in this repo, named 'variantAnnotator', allow for parsing of a VCF file and result in a tab-delimited tsv
with variant annotations, sequence metrics, and ExAC variant frequencies (if applicable).

'variantAnnotator' was developed using python3.9.4. To match dependencies, it is suggested you create a python
virtual environment and install requirements.txt. E.g.,

python3 -m venv venv
source venv/bin/activate
pip3 install -r requirements.txt

Usage of 'variantAnnotator' requires downloading necessary reference cds (coding sequence) files from the Consensus
CDS (CCDS) project website (https://www.ncbi.nlm.nih.gov/projects/CCDS). Files must be uncompressed. The required
files are:

1. CCDS.current.txt: coordinates of cds sequences and constituent exons
2. CCDS_nucleotide.current.fna: nucleotide sequences of the cds regions/exons
3. CCDS_protein.current.faa: protein sequence data of cds regions/exons (based on translating
                             CCDS_nucleotide.current.fna).
                             NOTE: stop codons are not included in CCDS_protein.current.faa but are added in when
                             variantAnnotator.py is run.

The most up to date CCDS cds data for build GRCh37 can be found at
https://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs37.3/.

Consuming the above files, variantAnnotator.py will take an input unannotated VCF file and return annotations in a
tab-delimited format.

Example usage:
python3 variantAnnotator.py -V 'example.vcf' -C 'CCDS.current.txt' -N 'CCDS_nucleotide.current.fna' -P 'CCDS_protein.current.faa' -o 'example.tsv'

After running, the output tsv will contain the following columns:

1. #CHROM: chromosome location of variant

2. POS: position location of variant in chromosome

3. REF: reference base

4. ALT: alternate variant base/sequence (can be more than 1 alternate variants total if the site is multi-allelic)
        NOTE: if multiple ALT alleles exist, these will be ';' seperated in the tsv file

5. variant_type: variant type, which will include one or more of the following:
                 1) substitution
                 2) insertion
                 3) deletion
                 NOTE: if multiple variant types exist for a single alternate variant (e.g., if the alternate
                 variant contains both an insertion and a deletion), they will all be included and ',' seperated.
                 additionally, if multiple ALT alleles exist, then these will be ';' seperated in the order
                 of the alleles in the ALT column

6. mutation_effect_type: the effect of a mutation. if a potential functional effect exists (i.e., falls
                         within exon/cds), the following possible annotations are given: 
                         1) frameshift
                         2) nonsense
                         3) missense (i.e., non-synonymous)
                         4) silent (i.e., synonymous)
                         5) None(exonic)
                         else, if the mutation falls outside an exon (i.e., outside a cds),
                         the following possible annotations are given:
                         1) intergenic
                         2) intronic
                         NOTE: only one effect type is given per alternate variant. however, if multiple ALT
                         alleles exist, then the effect types will be ';' seperated in the order of the alleles
                         in the ALT column. the effect type 'None(exonic)' represents examples where the nucleotide
                         sequence itself experienced no change (i.e., ALT == REF)

7. total_read_depth: total sequence read coverage at the site

8. total_reference_reads: total sequence reads at the site belonging to the reference base

9. percent_reference_reads(%): percent of sequence reads as a function of total_read_depth supporting the
                               reference base

10. total_alternate_reads: total sequence reads at the site belonging to the alternate variant(s)
                           NOTE: if multiple ALT alleles exist, then these will be ';' seperated in the order
                           of the alleles in the ALT column

11. percent_alternate_reads(%): percent of sequence reads as a function of total_read_depth supporting the
                                alternate variant(s)
                                NOTE: if multiple ALT alleles exist, then these will be ';' seperated in the order
                                of the alleles in the ALT column

12. ExAC_allele_frequency: frequency of alternate variant(s) from ExAC API (http://exac.hms.harvard.edu/)
                           NOTE: if multiple ALT alleles exist, then these will be ';' seperated in the order
                           of the alleles in the ALT column


NOTES/FUTURE IMPROVEMENTS:
Because of reliance on CCDS cds sequences, annotation of exonic, intronic, and intergenic regions relies strictly
on cds (coding sequence) coordinates. Thus, it is not unlikely that variants technically part of an exon (i.e., 5'
UTR region or 3' UTR region) might be annotated as intergenic (i.e., if up/downstream from start and stop codons).
Future updates may include the addition of 5'UTR and 3'UTR regions using other annotated gene coordinates such as
the 'UCSC Genes' track from UCSC Genome Browser (https://genome.ucsc.edu/cgi-bin/hgTables).

Annotations have been spotchecked throughout for accuracy. However, benchmarking accuracy of annotations from
variantAnnotator should be done with other readily available annotation tools such as Annovar
(https://annovar.openbioinformatics.org/) or SeattleSeq (https://snp.gs.washington.edu/SeattleSeqAnnotation154/) in
order to compare agreement.

Annotations for chromosome Y are currently not accounted for. The reason for this is that the various CCDS input
files give the same ID for both chrX and chrY cds sequences if they are part of the pseudoautosomal region. To
avoid that potential clash, the entirety of chrY is simply skipped if it's in the VCF. Accomadating chrY in future
versions should be trivial to implement.
