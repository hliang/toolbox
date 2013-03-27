
import os
import sys
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO

args = sys.argv
fafile = args[1]
print 'Input file is:', fafile
# fasta_string = open(fafile).read()
# fa_record = SeqIO.read(open(fafile), format="fasta") # get an error if More than one record found in handle
record_iterator = SeqIO.parse(fafile, "fasta")

fasta_string = '''>comp36320_c0_seq1
GGCAGCCATTGAAGCTTGCACGAAGGCGGGCGTCGCCGTCAAGATGGTCACCGGCGACAA
CATCCTCACGGCCCGCGCGATCGCCAAGGAGTGCGGCATCATATCCAGCAACGACCCCAG
CGGCATCGTCATCGAGGGGCACGAGTTCCGCGCCATGTCGCCGGAGCAGCAGCTGGAGAT
>testttt2
CGTGGACAGGATCCGCGTCATGGCGCGGTCCCTGCCGCTGGACAAGCTGGCGCTGGTGCA
GCGGCTGAAGCAGAAGGGGCACGTGGTGGCCGTGACCGGCGACGGCACCAACGACGCGCC
GGCGCTCAAGGAGGCGG'''

for fa_record in record_iterator:
  print '=============', fa_record.id, '==============='
# this is the actual blast search, so it may take a moment
# result_handle = NCBIWWW.qblast('blastx', 'nr', fasta_string, hitlist_size=25, expect=0.001, format_type='XML')
  result_handle = NCBIWWW.qblast('blastx', 'nr', fa_record.seq, hitlist_size=25, expect=0.001, format_type='XML')
# result_handle = NCBIWWW.qblast('blastx', 'nr', fa_record.format("fasta"), hitlist_size=25, expect=0.001, format_type='XML')

# save blast_handle as a file for later use
#  result_handle.seek(0)
  print result_handle.read()

# save_file = open('blast-output.xml', 'w')
# save_file.write(result_handle.read())
# save_file.close()
# result_handle.close()
