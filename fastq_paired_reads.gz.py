##########################################################
# For paired-end sequences, after trimming or filtering, the two ends
# are not "synchronized": some pairs will get one end passed while the
# other end failed.  This is a simple and not-optimized python script
# to re-synchronize the reads.
# modified from scripts written by maubp from SeqAnswers.com

#sys.path.append('/homes/hliang/local/lib/python2.7/site-packages')
import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator #Biopython 1.51 or later
#support for reading gzip file
import gzip

##########################################################
#
# Change the following settings to suit your needs
#

#input_forward_filename = "split_RKQQ_l30_q15_p80_R1.fq"
#input_reverse_filename = "split_RKQQ_l30_q15_p80_R2.fq"
#
#output_paired_forward_filename = "split_RKQQ_l30_q15_p80_paired_R1.fq"
#output_paired_reverse_filename = "split_RKQQ_l30_q15_p80_paired_R2.fq"
#output_orphan_filename = "split_RKQQ_unpaired_orphans.fq"

input_forward_filename = sys.argv[1]
input_reverse_filename = sys.argv[2]

output_paired_forward_filename = "paired."+sys.argv[1]
output_paired_reverse_filename = "paired."+sys.argv[2]
output_orphan_filename = "orphans." + sys.argv[1] + "." + sys.argv[2]

f_suffix = ""
r_suffix = ""

##########################################################

if f_suffix:
    f_suffix_crop = -len(f_suffix)
    def f_name(title):
        """Remove the suffix from a forward read name."""
        name = title.split()[0]
        assert name.endswith(f_suffix), name
        return name[:f_suffix_crop]
else:
    def f_name(title):
        return title.split()[0]

if r_suffix:
    r_suffix_crop = -len(r_suffix)
    def r_name(title):
        """Remove the suffix from a reverse read name."""
        name = title.split()[0]
        assert name.endswith(r_suffix), name
        return name[:r_suffix_crop]
else:
    def r_name(title):
        return title.split()[0]

print "Scaning reverse file to build list of names..."    
reverse_ids = set()
paired_ids = set()
for title, seq, qual in FastqGeneralIterator(gzip.open(input_reverse_filename)):
    reverse_ids.add(r_name(title))

print "Processing forward file..."
forward_handle = gzip.open(output_paired_forward_filename, "wb")
orphan_handle = gzip.open(output_orphan_filename, "wb")
for title, seq, qual in FastqGeneralIterator(gzip.open(input_forward_filename)):
    name = f_name(title)
    if name in reverse_ids:
        #Paired
        paired_ids.add(name)
        reverse_ids.remove(name) #frees a little memory
        forward_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
    else:
        #Orphan
        orphan_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
forward_handle.close()
del reverse_ids #frees memory, although we won't need more now

print "Processing reverse file..."
reverse_handle = gzip.open(output_paired_reverse_filename, "wb")
for title, seq, qual in FastqGeneralIterator(gzip.open(input_reverse_filename)):
    name = r_name(title)
    if name in paired_ids:
        #Paired
        reverse_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
    else:
        #Orphan
        orphan_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
orphan_handle.close()
reverse_handle.close()
print "Done"