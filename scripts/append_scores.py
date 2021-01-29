"""
Appends Doench scores to the Guidescan databases in BAM
format.
"""

import pickle
import argparse
import os
import pysam
import numpy as np
import binascii
import sys

from Bio import SeqIO

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(SCRIPT_DIR + '/Rule_Set_2_scoring/analysis')

import model_comparison

###############
## CFD Score ##
###############

def calc_cfd(wt,sg,pam,mm_scores,pam_scores):
    score = 1
    sg = sg.replace('T','U')
    wt = wt.replace('T','U')
    s_list = list(sg)
    wt_list = list(wt)
    for i,sl in enumerate(s_list):
        if wt_list[i] == sl:
            score*=1
        else:
            try:
                key = 'r'+wt_list[i]+':d'+revcom(sl)+','+str(i+1)
                score*= mm_scores[key]
            except KeyError:
                continue
    score*=pam_scores[pam]
    return (score)

def get_mm_pam_scores(mms,pams):
    try:
        mm_scores = pickle.load(open(mms,'rb'))
        pam_scores = pickle.load(open(pams,'rb'))
        return (mm_scores,pam_scores)
    except:
        raise Exception("Could not find file with mismatch scores or PAM scores")

def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)

PAM_PKL = SCRIPT_DIR + "/CFD_scoring/pam_scores.pkl"
MM_PKL = SCRIPT_DIR + "/CFD_scoring/mismatch_score.pkl"
mm_scores, pam_scores = get_mm_pam_scores(MM_PKL, PAM_PKL)

######################
## Rule Set 2 Score ##
######################

def load_rs2_model(model_file):
    try:
        with open(model_file, 'rb') as f:
            model = pickle.load(f)    
            return model
    except Exception as e:
        print(e)
        raise Exception("could not find model stored to file %s" % model_file)

RULE_SET_2_PKL = SCRIPT_DIR + "/Rule_Set_2_scoring/saved_models/V3_model_nopos.pickle"
rule_set_2_model = load_rs2_model(RULE_SET_2_PKL)

##########################
## Guidescan Processing ##
##########################

def repeat(item):
    while True:
        yield item

def hex_to_array(hexstr):
    return np.fromstring(binascii.unhexlify(hexstr), dtype=int)

def hex_to_offtargetinfo(hexstr, delim):
    mainarr = hex_to_array(hexstr)
    index = np.where(mainarr == delim)[0]
    slices = zip(index, np.append(index[1:], [len(mainarr)]))
    out = []
    for start, end in slices:
        out += zip(repeat(mainarr[end - 1]), mainarr[start + 1:end - 1])
    return out

def map_int_to_coord(x, genome, onebased=False):
    strand = '+' if x > 0 else '-'
    x = abs(x)
    x -= 1
    i = 0
    while genome[i]['LN'] <= x:
        x -= genome[i]['LN']
        i += 1
    chrom = genome[i]['SN']
    coord = x
    if onebased:
        coord += 1
    t = (chrom, coord, strand)
    return t

def get_length(genome):
    return sum(p['LN'] for p in genome)

def get_nonexist_int_coord(genome):
    return -(get_length(genome) + 1)

def map_coord_to_sequence(fasta_record_dict, length, chr, pos, strand):
    if strand == '+':
        pos_start = pos - 1 - length
        pos_end = pos - 1
    else:
        pos_start = pos + 4
        pos_end = pos + 4 + length
    return fasta_record_dict[chr].seq[pos_start:pos_end].upper()

def map_coord_to_30nt_context(fasta_record_dict, chr, start, end, antisense):
    if not antisense:
        pos_start = start - 4
        pos_end = end + 3
        return fasta_record_dict[chr].seq[pos_start:pos_end].upper()
    else:
        pos_start = start - 3
        pos_end = end + 4
        return revcom(fasta_record_dict[chr].seq[pos_start:pos_end].upper())

def offtarget_hex_to_cfd_score(fasta_record_dict, genome, delim, ots_hex, sgrna):
    cfd = 0 
    for distance, offtarget_pos in hex_to_offtargetinfo(ots_hex, delim=delim):
        chr, pos, strand = map_int_to_coord(offtarget_pos, genome)
        seq = str(map_coord_to_sequence(fasta_record_dict, len(sgrna),
                                        chr, pos, strand))
        seq = revcom(seq) if strand == '-' else seq 
        cfd += calc_cfd(sgrna, seq, "GG", mm_scores, pam_scores)
    return 1 / (1 + cfd)

##########################
##     Main Method      ##
##########################

def argument_parser():
    parser = argparse.ArgumentParser(description=(
        'Appends Doench et al. scores to the Guidescan database, '
        'writing the database (with scores) as SAM file to stdout.'
    ))

    parser.add_argument('database',
        type=str,
        help='Guidescan database file (SAM/BAM format acceptable)')
    parser.add_argument('ref_seq',
        type=str,
        help='Reference sequence used to construct Guidescan DB')

    return parser

if __name__ == "__main__":
    args = argument_parser().parse_args()

    guidescan_db = pysam.AlignmentFile(args.database, "rb")
    fasta_record_dict = SeqIO.to_dict(SeqIO.parse(args.ref_seq, "fasta"))

    genome = guidescan_db.header['SQ']
    delim = get_nonexist_int_coord(genome)

    def compute_cfd(guide_record):
        if not guide_record.has_tag("of"):
            return 1

        ots_hex = guide_record.get_tag("of")
        sgrna   = guide_record.query_name[:20]
        cfd = offtarget_hex_to_cfd_score(fasta_record_dict, genome, delim, ots_hex, sgrna)

        return cfd

    def compute_rs2(guide_record):
        chr, start, end, antisense = (guide_record.reference_name, guide_record.reference_start,
                                      guide_record.reference_end, guide_record.is_reverse)
        seq = map_coord_to_30nt_context(fasta_record_dict, chr, start, end, antisense)
        seq = ''.join([nuc if nuc != 'N' else 'A' for nuc in list(seq)])
        rs2 = model_comparison.predict(seq, -1, -1, model=rule_set_2_model)
        return rs2

    outfile = pysam.AlignmentFile("-", "w", template=guidescan_db)
    for guide_record in guidescan_db:
        cfd = compute_cfd(guide_record)
        rs2 = compute_rs2(guide_record)
        
        guide_record.set_tag("cs", cfd)
        guide_record.set_tag("ds", rs2)
        outfile.write(guide_record)
