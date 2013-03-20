
import os
import sys
import re

#def read_snp():

def read_blast(fname):
    """Read blast output (converted from xml to tab-delimited format),
    and return a dict of dicts of hsp entries."""

    blast_hits = {}
    f = open(fname, 'rU')
    for line in f:
        line = line.rstrip('\n')
        hit_entry = line.split('\t')
        if len(hit_entry) > 2:                # each hit entry has more than 2 fiels
            hit_entry[0] = hit_entry[0].split()[0] # proper hitDef
            hit_entry[3] = hit_entry[3].split()[0] # proper hitDef
            hit_entry[5:7] = [float(x) for x in hit_entry[5:7]] # convert into float point numbers
            hit_entry[7:17] = [int(x) for x in hit_entry[7:17]] # convert into integer numbers
            # hit_entry[11:17] = [float(x) for x in hit_entry[11:17]] # convert into float point numbers
            if hit_entry[0] in blast_hits:
                if hit_entry[3] in blast_hits[hit_entry[0]]:
                    blast_hits[hit_entry[0]][hit_entry[3]].append(hit_entry[5:20])
                else:
                    blast_hits[hit_entry[0]][hit_entry[3]]=[hit_entry[5:20]]
            else:
                blast_hits[hit_entry[0]] = {hit_entry[3]:[hit_entry[5:20]] }
    #print blast_hits
    return blast_hits

# xsltproc --novalid blast2tsv.xsl2 blastout.xml > blastout.tsv
# 0queryDef 1queryLen 2hitId 3hitDef 4hitLen 5Hsp_bit-score 6Hsp_evalue 7Hsp_query-from 8Hsp_query-to 9Hsp_hit-from 10Hsp_hit-to 11Hsp_query-frame 12Hsp_hit-frame 13Hsp_identity 14Hsp_positive 15Hsp_gaps 16Hsp_align-len 17Hsp_qseq 18Hsp_hseq 19Hsp_midline

def get_snps():
    """Return a dict containing seq:{pos1:[Ref,Alt], pos2:[Ref,Alt], pos3:[Ref,Alt]...}
    """
    # snp_pos = {'comp424_c0_seq1':{80:['A', 'T'], 110:['G', 'C']}, 'comp1612_c1_seq1':{100:['A', 'T'], 400:['G', 'C'], 900:['T', 'G'] }, 'comp26268_c0_seq3':{100:['A', 'T'], 900:['T', 'G']}}
    snp_pos = {}
    fname = 'varscan_2.snp.flt20130131.20.4'
    # print 'processing ...', fname
    f = open(fname, 'rU')
    for line in f:
        if not line.startswith('#'): # skip comment line 
            line = line.rstrip('\n')
            snp_entry = line.split('\t')
            # print snp_entry

            # variant in CS (monomorphic)
            var_cs = snp_entry[7] 
            # variant in RWG1 (monomorphic allele which is different from CS allele, or polymorphic site which has both alleles). for the purpose of this program, assign only the unique allele.
            var_rwg1 = snp_entry[2] + snp_entry[3]
            if var_rwg1[0] == snp_entry[11] or var_rwg1[1] == snp_entry[11]:
                var_rwg1 = snp_entry[11]
            #if var_rwg1 == var_cs:
            #    var_rwg1 = snp_entry[3]

            # store in dict
            if snp_entry[0] in snp_pos:
                snp_pos[snp_entry[0]][int(snp_entry[1])] = [var_rwg1, var_cs]
            else:
                snp_pos[snp_entry[0]] = { int(snp_entry[1]): [var_rwg1, var_cs] }
    return snp_pos


def peek_at_neighbor(snp_pos, blast_hits):
    """check the alignment before and after the SNP site,
    see if they have good alignment,
    and also see if the both alleles are present.
    """
    dist = 10
    for query in snp_pos:
        print '============', query, '============'
        query_hits = {}
        if query in blast_hits:                  # query has blast hits
            query_hits = blast_hits[query]       # {subj1:[[hitA], [hitB]], subj2:[[hitC]], subj3:[[hitD],[hitE],[hitF]]}
            for pos in sorted(snp_pos[query]):
                print pos, snp_pos[query][pos]
                # search each alignment, check if the snps pos is inside the alignment
                for subject in query_hits:
                    # print subject
                    for this_hit in query_hits[subject]:
                        q_begin, q_end, s_begin, s_end = this_hit[2], this_hit[3], this_hit[4], this_hit[5] # start and end pos of query and subject in alignment
                        if pos >= q_begin and pos <= q_end:
                            q_aln_seq_w_gap = this_hit[-3]
                            pos_in_aln_wo_gap = 1
                            pos_in_aln_w_gap = 1
                            # print q_aln_seq_w_gap
                            while pos_in_aln_wo_gap < pos - q_begin + 1:
                                # print len(q_aln_seq_w_gap) pos_in_aln_w_gap, pos_in_aln_wo_gap
                                if not q_aln_seq_w_gap[pos_in_aln_w_gap-1] == '-':
                                    pos_in_aln_wo_gap = pos_in_aln_wo_gap + 1
                                pos_in_aln_w_gap = pos_in_aln_w_gap + 1
                                    
                            # print alignment
                            # print '=== hit on', subject, '==='
                            # print '\t'.join([str(q_begin), this_hit[-3] , str(q_end)])
                            # print '\t'.join([' '*len(str(q_begin)), this_hit[-1] , str(q_end)])
                            # print '\t'.join([str(s_begin), this_hit[-2] , str(s_end)])

                            # q_aln_seq_wo_gap
                            # q_aln_seq_wo_gap = re.sub('-', '', q_aln_seq_w_gap)
                            # print '# Q #', q_aln_seq_wo_gap[pos_in_aln_wo_gap-1-dist:pos_in_aln_wo_gap-1], q_aln_seq_wo_gap[pos_in_aln_wo_gap-1], q_aln_seq_wo_gap[pos_in_aln_wo_gap:pos_in_aln_wo_gap+dist]

                            # print '# Q #', this_hit[-3][pos_in_aln-dist:pos_in_aln], this_hit[-3][pos_in_aln], this_hit[-3][pos_in_aln+1:pos_in_aln+dist+1]
                            # print '# | #', this_hit[-1][pos_in_aln-dist:pos_in_aln], this_hit[-1][pos_in_aln], this_hit[-1][pos_in_aln+1:pos_in_aln+dist+1]
                            # print '# S #', this_hit[-2][pos_in_aln-dist:pos_in_aln], this_hit[-2][pos_in_aln], this_hit[-2][pos_in_aln+1:pos_in_aln+dist+1]
                            print '# Q #', this_hit[-3][pos_in_aln_w_gap-1-dist:pos_in_aln_w_gap-1], this_hit[-3][pos_in_aln_w_gap-1], this_hit[-3][pos_in_aln_w_gap:pos_in_aln_w_gap+dist]
                            print '# | #', this_hit[-1][pos_in_aln_w_gap-1-dist:pos_in_aln_w_gap-1], this_hit[-1][pos_in_aln_w_gap-1], this_hit[-1][pos_in_aln_w_gap:pos_in_aln_w_gap+dist]
                            print '# S #', this_hit[-2][pos_in_aln_w_gap-1-dist:pos_in_aln_w_gap-1], this_hit[-2][pos_in_aln_w_gap-1], this_hit[-2][pos_in_aln_w_gap:pos_in_aln_w_gap+dist]
        else:                                # query has no hit
            print query, 'is not in the blast result or it has hit.'
    



def main():
    filename = sys.argv[1]
    blast_hits = read_blast(filename)
    #print read_blast(filename)

    snp_pos = get_snps()
    # print snp_pos
    peek_at_neighbor(snp_pos, blast_hits)


if __name__ == '__main__':
    main()
