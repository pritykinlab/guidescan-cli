import pickle
import os
import numpy as np
import argparse
import pysam
import sys
import binascii

from Bio import SeqIO

from functools import reduce

def repeat(item):
    while True:
        yield item

def get_length(genome):
    return sum(p['LN'] for p in genome)

def get_nonexist_int_coord(genome):
    return -(get_length(genome) + 1)

def ilen(iterable):
    return reduce(lambda sum, _: sum + 1, iterable, 0)

def hex_to_array(hexstr):
    return np.frombuffer(binascii.unhexlify(hexstr), dtype=int)

def hex_to_offtargetinfo(hexstr, delim):
    mainarr = hex_to_array(hexstr)
    index = np.where(mainarr == delim)[0]
    slices = list(zip(np.insert(index, 0, -1), index))
    out = []
    for start, end in slices:
        out += zip(repeat(mainarr[end - 1]), mainarr[start + 1:end - 1])
    return out

def map_int_to_coord(x, genome, onebased=False):
    strand = '+' if x > 0 else '-'
    x = abs(x)
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

def map_coord_to_sequence(fasta_record_dict, sgrna, chr, pos, strand):
    if strand == '+':
        pos_start = pos + 1 - len(sgrna)
        pos_end = pos + 1
    else:
        pos_start = pos
        pos_end = pos + len(sgrna)
    return str(fasta_record_dict[chr].seq[pos_start:pos_end].upper()), pos_start, pos_end

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

###############
## CFD Score ##
###############

def calc_cfd_e(sg,wt,pam,mm_scores,pam_scores):
    score = 1
    wt = wt.replace('T','U')
    sg = sg.replace('T','U')
    wt_list = list(wt)
    sg_list = list(sg)
    for i,wl in enumerate(wt_list):
        if sg_list[i] == wl:
            score*=1
        else:
            try:
                key = 'r'+sg_list[i]+':d'+revcom(wl)+','+str(i+1)
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
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A', 'N': 'N'}
    letters = list(s[::-1])
    letters = [(basecomp[base] if base in basecomp else base) for base in letters]
    return ''.join(letters)

PAM_PKL = SCRIPT_DIR + "/cfd/pam_scores.pkl"
MM_PKL = SCRIPT_DIR + "/cfd/mismatch_score.pkl"
mm_scores, pam_scores = get_mm_pam_scores(MM_PKL, PAM_PKL)

# Curry last two arguments for ease of use
calc_cfd = lambda sg, wt, pam: calc_cfd_e(sg, wt, pam, mm_scores, pam_scores)

def decode_off_targets(sam_record, genome, delim, fasta_record_dict):
    if not sam_record.has_tag('of'):
        return

    ots = sam_record.get_tag('of')
    ots = hex_to_offtargetinfo(ots, delim)

    count = 0
    for distance, pos in ots:
        chrm, pos, strand = map_int_to_coord(pos, genome)

        sgrna = sam_record.query_sequence
        if sam_record.is_reverse:
            sgrna = revcom(sgrna)

        # Note - offtarget_pos_start is always the lower in value of {offtarget_pos_start, offtarget_pos_end},
        # regardless of strand
        offtarget, offtarget_pos_start, offtarget_pos_end = map_coord_to_sequence(fasta_record_dict, sgrna, chrm, pos, strand)

        if len(offtarget) == 23:
            seq = revcom(offtarget) if strand == '-' else offtarget
            cfd = calc_cfd(sgrna, seq[:20], seq[21:23])
        else:
            cfd = None

        offtarget_map = {
            "identifier": sam_record.query_name,
            "distance": distance,
            "chr": chrm,
            # 0-based indexing is used by pysam and inside this script;
            # For reporting, purposes, switch to 1-based indexing to avoid confusion with samtools interoperability
            "pos": offtarget_pos_start + 1,
            "sense": strand,
            "offtarget": revcom(offtarget) if strand == '-' else offtarget,
            "cfd": cfd
        }

        yield offtarget_map
        count += 1

def load_guide_db(samfile):
    sam_db = pysam.AlignmentFile(samfile, "rb")
    genome = sam_db.header['SQ']
    delim = get_nonexist_int_coord(genome)
    return sam_db, delim, genome

def output_complete(offtargets):
    for i, offtarget in enumerate(offtargets):
        print(','.join(map(str, [
            offtarget['identifier'], i, offtarget['offtarget'], offtarget['chr'],
            offtarget['pos'], offtarget['sense'], offtarget['distance'],
            offtarget['cfd'] or ''
        ])))

def output_succinct(record, offtargets):
    identifier = record.query_name
    sequence = record.query_sequence
    chrm = record.reference_name
    position = record.reference_start
    sense = '-' if record.is_reverse else '+'

    match_counts = [0, 0, 0, 0]
    cfd_sum = None
    if offtargets:
        if all(offtarget['cfd'] is not None for offtarget in offtargets):
            cfd_sum =  sum((offtarget['cfd'] for offtarget in offtargets))

        flag = False

        for offtarget in offtargets:
            match_counts[offtarget['distance']] += 1

            if (offtarget['distance'] == 0 and flag == False and offtarget['cfd'] is not None and cfd_sum is not None):
                cfd_sum -= offtarget['cfd']
                flag = True

    if cfd_sum:
        specificity = 1 / (1 + cfd_sum)
        if cfd_sum == 0:
            specificity = 1
    else:
        specificity = ''

    print(','.join(list(map(str, [identifier, sequence, chrm, position, sense,
                                  match_counts[0], match_counts[1], match_counts[2],
                                  match_counts[3], specificity]))))

def parse_args():
    p = argparse.ArgumentParser(
        "Decodes Guidescan2 off-target information in hex-encoded SAM/BAM format."
    )

    p.add_argument(
        "grna_database",
        help="SAM/BAM file containing Guidescan2 processed gRNAs."
    )

    p.add_argument(
        "fasta_file",
        help="FASTA file for resolving off-target sequences."
    )

    p.add_argument(
        "--mode",
        help="Succinct or complete off-target information.",
        choices=['succinct', 'complete'],
        default='succinct'
    )

    return p.parse_args()

SUCCINCT_HEADER = ('id,sequence,chromosome,position,sense,'
                   'distance_0_matches,distance_1_matches,'
                   'distance_2_matches,distance_3_matches,'
                   'specificity')
COMPLETE_HEADER = ('id,match_number,sequence,chromosome,position,sense,distance,cfd')

if __name__ == "__main__":
    args = parse_args()

    sam_db, delim, genome = load_guide_db(args.grna_database)
    fasta_dict = SeqIO.to_dict(SeqIO.parse(args.fasta_file, "fasta"))
    decode_ot = lambda record: decode_off_targets(record, genome, delim, fasta_dict)

    if args.mode == 'succinct':
        print(SUCCINCT_HEADER)
        for sam_record in sam_db:
            output_succinct(sam_record, list(decode_ot(sam_record)))
    elif args.mode == 'complete':
        print(COMPLETE_HEADER)
        for sam_record in sam_db:
            output_complete(decode_ot(sam_record))