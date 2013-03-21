
# python2 good_neighbor.py  blast_xml2tsv  snp_pos
# e.g. python2 good_neighbor.py megablast.varscan_2.snp.flt20130131.20.4.fa_vs_wheat_fs_gDNA.tsv varscan_2.snp.flt20130131.20.4 > identify_false_positive_rwg1_uniq_alelle.txt

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

def get_snps(fname):
    """Return a dict containing seq:{pos1:[Ref,Alt], pos2:[Ref,Alt], pos3:[Ref,Alt]...}
    """
    # snp_pos = {'comp424_c0_seq1':{80:['A', 'T'], 110:['G', 'C']}, 'comp1612_c1_seq1':{100:['A', 'T'], 400:['G', 'C'], 900:['T', 'G'] }, 'comp26268_c0_seq3':{100:['A', 'T'], 900:['T', 'G']}}
    # fname = 'varscan_2.snp.flt20130131.20.4'
    snp_pos = {}
    # print 'processing ...', fname
    f = open(fname, 'rU')
    for line in f:
        if not line.startswith('#'): # skip comment line 
            line = line.rstrip('\n')
            snp_entry = line.split('\t')
            # print snp_entry

            # variant in CS (monomorphic)
            var_cs = snp_entry[2] 
            if snp_entry[6] == '100':
                var_cs = snp_entry[3]

            # varianti(s) in RWG1 (monomorphic allele which is different from CS allele, or polymorphic site which has both alleles).
            var_rwg1 = snp_entry[2] + snp_entry[3]
            # if var_rwg1[0] == snp_entry[11] or var_rwg1[1] == snp_entry[11]: # if this site is monomorphic
            if snp_entry[10] == '0' or snp_entry[10] == '100': # if this site is monomorphic
                var_rwg1 = snp_entry[11]

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
    dist = 20
    for query in snp_pos:
        # print '============', query, '============'
        query_hits = {}
        if query in blast_hits:                  # query has blast hits
            query_hits = blast_hits[query]       # {subj1:[[hitA], [hitB]], subj2:[[hitC]], subj3:[[hitD],[hitE],[hitF]]}
            for pos in sorted(snp_pos[query]):
                # print str(pos) + '\t' + snp_pos[query][pos][0] + '\t' + snp_pos[query][pos][1]

                # search each alignment, check if the snps pos is inside the alignment
                for subject in query_hits:
                    # print subject
                    for this_hit in query_hits[subject]:
                        q_begin, q_end, s_begin, s_end = this_hit[2], this_hit[3], this_hit[4], this_hit[5] # start and end pos of query and subject in alignment
                        if pos >= q_begin and pos <= q_end: # if SNP is located within the alignment
                            # find out how the sequences before and after SNP are aligned
                            q_aln_seq_w_gap = this_hit[-3]
                            pos_in_aln_wo_gap = 1
                            pos_in_aln_w_gap = 1
                            while pos_in_aln_wo_gap < pos - q_begin + 1:
                                if not q_aln_seq_w_gap[pos_in_aln_w_gap-1] == '-':
                                    pos_in_aln_wo_gap = pos_in_aln_wo_gap + 1
                                pos_in_aln_w_gap = pos_in_aln_w_gap + 1
                            if len(this_hit[-1][pos_in_aln_w_gap-1-dist:pos_in_aln_w_gap-1]) == dist and len(this_hit[-1][pos_in_aln_w_gap:pos_in_aln_w_gap+dist]) == dist and this_hit[-1][pos_in_aln_w_gap-1-dist:pos_in_aln_w_gap-1].find(' ') == -1 and this_hit[-1][pos_in_aln_w_gap:pos_in_aln_w_gap+dist].find(' ') == -1 : # if alignments are present in both xx bp before and after SNP, and they are perfect, i.e. 100% match
                                snp_found_by_blast = [this_hit[-3][pos_in_aln_w_gap-1], this_hit[-2][pos_in_aln_w_gap-1]]
                                snp_pos[query][pos].append(snp_found_by_blast)

                snp_pos[query][pos] = uniq_f7_noHash(snp_pos[query][pos])
                # print query, pos, snp_pos[query][pos]
                '''
                # add snps found with blast
                for e in uniq_f7_noHash(snp_pos[query][pos][2:]):
                    if len(e):  # not empty
                        print 'blast:\t'+'\t'.join(e)
                # print '-------'
                '''
        else:                                # query has no hit
            for pos in sorted(snp_pos[query]):
                # print query, pos, snp_pos[query][pos] '# query is not in the blast result or it has hit.'
                pass
    
    goodneighbor = snp_pos
    return goodneighbor
'''
dict of dicts of lists
comp50115_c0_seq2 393 ['GC', 'G', ['G', 'G']]
comp56990_c0_seq11 240 ['CT', 'C', ['C', 'C']]
comp50115_c0_seq1 511 ['GC', 'G', ['G', 'G']]
comp51266_c1_seq3 711 ['CT', 'T', ['C', 'C'], ['C', 'T']]
comp49596_c0_seq2 394 ['AC', 'C', ['A', 'C']]
comp49596_c0_seq2 515 ['AG', 'A', ['A', 'A']]
comp49596_c0_seq1 451 ['AC', 'C', ['A', 'C']]
comp49596_c0_seq1 572 ['AG', 'A', ['A', 'A']]
comp60624_c0_seq1 541 ['CT', 'T', ['C', 'T'], ['C', 'C']]
comp59948_c0_seq3 464 ['AG', 'A', ['A', 'G'], ['A', 'A']]
comp59948_c0_seq2 599 ['AG', 'A', ['A', 'G'], ['A', 'A']]
comp42124_c0_seq1 146 ['AG', 'A', ['A', 'A']]
comp53773_c0_seq11 314 ['TC', 'T', ['T', 'T']]
comp53773_c0_seq15 246 ['TC', 'T', ['T', 'T']]
'''

def kill_gn(snp_pos_pm):
    '''
    if a rwg1-unique allele is found in blast result, with perfect alignment itself and before and after it,
    that indicates that the rwg1-unique allele is present in CS, it's a false positive unique allele.
    we should be interested in rwg1-unique alleles.
    '''
    print '#seqId pos rwg1_RNA cs_RNA cs_gDNA'
    for seqId in sorted(snp_pos_pm.keys()):
        for pos in sorted(snp_pos_pm[seqId].keys()):
            if len(snp_pos_pm[seqId][pos][2:]) == 0: # if there is NO alignment found by blast
                print seqId, pos, snp_pos_pm[seqId][pos][0], snp_pos_pm[seqId][pos][1]
            else: # there are alignments found by blast
                alleles_in_cs_gDNA = ''
                for cs_allele in snp_pos_pm[seqId][pos][2:]:
                    alleles_in_cs_gDNA = alleles_in_cs_gDNA + cs_allele[1]
                uniq_rwg1_allele = re.sub( snp_pos_pm[seqId][pos][1] , '' , snp_pos_pm[seqId][pos][0])  # cs allele (monomorphic, found in RNA-seq), empty string(remove), rwg1 allele(s)
                if alleles_in_cs_gDNA.find(uniq_rwg1_allele) >= 0: # uniq_rwg1_allele found in blast alignment, meaning that it's actually NOT unique.
                    print seqId, pos, snp_pos_pm[seqId][pos][0], snp_pos_pm[seqId][pos][1], alleles_in_cs_gDNA, '#NotUniq'
                else:
                    print seqId, pos, snp_pos_pm[seqId][pos][0], snp_pos_pm[seqId][pos][1], alleles_in_cs_gDNA




def uniq_f6(seq):
    """input a list,
    return a list without duplicates"""
    # Not order preserving
    keys = {}
    for e in seq:
        keys[e] = 1
    return keys.keys()

def uniq_f7(seq):
    """input a list,
    return a list without duplicates"""
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if x not in seen and not seen_add(x)]

def uniq_f7_noHash(seq):
    """input a list,
    return a list without duplicates"""
    seen = set()
    return [ x for x in seq if str( x ) not in seen and not seen.add( str( x ) )]

def main():
    blast_tsv_fname = sys.argv[1]
    varscan_snp_pos_fname = sys.argv[2]

    blast_hits = read_blast(blast_tsv_fname)

    snp_pos = get_snps(varscan_snp_pos_fname)

    snp_pos_pm = peek_at_neighbor(snp_pos, blast_hits) # pm = perfect match
    kill_gn(snp_pos_pm)


if __name__ == '__main__':
    main()

