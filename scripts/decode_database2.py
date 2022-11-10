import pickle
import os
import numpy as np
import argparse
import pysam
import sys
import binascii

from Bio import SeqIO

from functools import reduce

FIX_LEGACY_BUG = True
FIX_ANOTHER_LEGACY_BUG = True

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
    x -= 1  # WHY DO WE DO THIS??
    i = 0
    while genome[i]['LN'] <= x:
        x -= genome[i]['LN']
        i += 1
    chrom = genome[i]['SN']
    coord = x
    if onebased:
        coord += 1
    t = [chrom, coord, strand]

    if FIX_ANOTHER_LEGACY_BUG:
        # Later on we're going to be incrementing the coord by 1,
        # so if we currently happen to be exactly at the last position
        # we WILL (conceptually) end up on the next chromosome at the first position
        # Do this explicitly here, since the rest of the code is not equipped to handle this case
        if coord+1 == genome[i]['LN']:
            t[0] = genome[i+1]['SN']
            t[1] = 0

    return tuple(t)

def map_coord_to_sequence(fasta_record_dict, sgrna, chr, pos, strand):
    if strand == '+':
        pos_start = pos + 2 - len(sgrna)
        pos_end = pos + 2
    else:
        pos_start = pos + 1
        pos_end = pos + 1 + len(sgrna)
    return fasta_record_dict[chr].seq[pos_start:pos_end].upper()

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
        import json
        scores = {'mm': mm_scores, 'pam': pam_scores}
        with open('scores.json', 'w') as f:
            j = json.dumps(scores, indent=4)
            print(j, file=f)
        return (mm_scores,pam_scores)
    except:
        raise Exception("Could not find file with mismatch scores or PAM scores")

def com(s, ignore_case=False):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A', 'N': 'N'}
    if ignore_case:
        s = s.upper()
    letters = list(s)
    letters = [(basecomp[base] if base in basecomp else base) for base in letters]
    return ''.join(letters)

def revcom(s, ignore_case=False):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A', 'N': 'N'}
    letters = list(s[::-1])
    if ignore_case:
        letters = [l.upper() for l in letters]
    letters = [(basecomp[base] if base in basecomp else base) for base in letters]
    return ''.join(letters)

def rev(s, ignore_case=False):
    if ignore_case:
        s = s.upper()
    return s[::-1]

def same(s, ignore_case=False):
    if ignore_case:
        s = s.upper()
    return s

PAM_PKL = SCRIPT_DIR + "/cfd/pam_scores.pkl"
MM_PKL = SCRIPT_DIR + "/cfd/mismatch_score.pkl"
mm_scores, pam_scores = get_mm_pam_scores(MM_PKL, PAM_PKL)

# Curry last two arguments for ease of use
calc_cfd = lambda sg, wt, pam: calc_cfd_e(sg, wt, pam, mm_scores, pam_scores)

def my_distance(seq, sgrna, strand):
    assert len(seq) == 23
    assert len(sgrna) == 23

    distance = 0
    if strand == '+':
        seq = seq[:-3]  # remove PAM at end
        sgrna = sgrna[:-3]  # remove PAM at end
    else:
        seq = seq[:-3]
        sgrna = sgrna[:-3]

    for s1, s2 in zip(seq[:20], sgrna[:20]):
        if s1 != s2:
            distance += 1

    return distance


def my_cfd(seq, sgrna, strand, mm_scores, pam_scores, buggy=False):

    # PAM is always at the end of sgrna, regardless of strand
    pam = sgrna[-2:]  # Why is this 2 letters instead of 3?

    if strand == '-':
        if buggy:
            sgrna = revcom(sgrna)

    score = 1
    sgrna = sgrna.replace('T', 'U')
    seq = seq.replace('T', 'U')

    for i, (_sgrna, _seq) in enumerate(zip(sgrna, seq), start=1):
        if _sgrna != _seq:
            key = f'r{_sgrna}:d{com(_seq.upper())},{i}'
            score *= mm_scores.get(key, 1)

    score *= pam_scores[pam]

    return score

def seq_is_same(theirs, mine, strand):
    if len(theirs) != len(mine):
        raise AssertionError
    if strand == '+':
        assert theirs == same(mine, ignore_case=True)
    else:
        assert theirs == revcom(mine, ignore_case=True)
    return True


def decode_off_targets(sam_record, genome, delim, fasta_record_dict):
    if not sam_record.has_tag('of'):
        return

    ots = sam_record.get_tag('of')
    ots = hex_to_offtargetinfo(ots, delim)
    if sam_record.has_tag('cf'):
        cfd_scores = np.frombuffer(binascii.unhexlify(sam_record.get_tag('cf')), dtype=np.float32)
    else:
        cfd_scores = None

    temp_chr = sam_record.get_tag('ZC').split('|') if sam_record.has_tag('ZC') else None
    if sam_record.has_tag('ZP'):
        temp_pos = sam_record.get_tag('ZP').split('|')
        temp_pos = [int(p) for p in temp_pos if p != '']
    else:
        temp_pos = None
    temp_strand = list(sam_record.get_tag('ZS')) if sam_record.has_tag('ZS') else None
    temp_seq_all = sam_record.get_tag('ZQ').split('|') if sam_record.has_tag('ZQ') else None
    temp_scores = [float(s) for s in sam_record.get_tag('ZN').split('|')[:-1]] if sam_record.has_tag('ZN') else None

    count = 0
    for i, (distance, pos) in enumerate(ots):
        original_pos = pos
        chrm, pos, strand = map_int_to_coord(pos, genome)

        sgrna = sam_record.query_sequence
        if sam_record.is_reverse:
            sgrna = revcom(sgrna)

        if temp_strand is not None:
            assert strand == temp_strand[i]
        if temp_chr is not None:
            assert chrm == temp_chr[i]
        if temp_pos is not None:
            assert pos == temp_pos[i]

        offtarget = str(map_coord_to_sequence(fasta_record_dict, sgrna, chrm, pos, strand))
        if temp_seq_all is not None:
            temp_seq = temp_seq_all[i]
        else:
            temp_seq = None

        cfd = None
        if len(offtarget) == 23:
            if temp_seq is not None:
                assert seq_is_same(offtarget, temp_seq, strand)
                assert distance == my_distance(temp_seq, sgrna, strand)

            seq_ = revcom(offtarget) if strand == '-' else offtarget
            if FIX_LEGACY_BUG:
                sgrna_ = sgrna
            else:
                sgrna_ = revcom(sgrna) if strand == '-' else sgrna

            cfd = calc_cfd(sgrna_, seq_[:20], seq_[21:23])

            if temp_seq is not None:
                my_py_cfd_incorrect = my_cfd(temp_seq, sgrna, strand, mm_scores, pam_scores, buggy=True)
                my_py_cfd_correct = my_cfd(temp_seq, sgrna, strand, mm_scores, pam_scores)

                if FIX_LEGACY_BUG:
                    assert abs(cfd - my_py_cfd_correct) < 1e-6
                else:
                    assert abs(cfd - my_py_cfd_incorrect) < 1e-6

                assert abs(my_py_cfd_correct - temp_scores[i]) < 1e-6

                if cfd_scores is not None:
                    assert abs(my_py_cfd_correct - cfd_scores[i]) < 1e-6
        else:
            _ = map_int_to_coord(original_pos, genome)
            offtarget2 = str(map_coord_to_sequence(fasta_record_dict, sgrna, chrm, pos, strand))
            raise RuntimeError('should never happen!')

        offtarget_map = {
            "identifier": sam_record.query_name,
            "distance": distance,
            "chr": chrm,
            "pos": pos,
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
                assert offtarget['cfd'] == 1
                cfd_sum -= offtarget['cfd']
                flag = True

    if cfd_sum:
        specificity = 1 / (1 + cfd_sum)
        if record.has_tag('sp') and FIX_LEGACY_BUG:
            if abs(record.get_tag('sp') - specificity) > 1e-5:
                raise AssertionError
        if cfd_sum == 0:
            specificity = 1
            if record.has_tag('sp') and FIX_LEGACY_BUG:
                assert record.get_tag('sp') == 1.
    else:
        specificity = ''
        if record.has_tag('sp') and FIX_LEGACY_BUG:
            if record.get_tag('sp') != 1.:
                raise AssertionError

    assert match_counts[0] == 1

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
        for sam_i, sam_record in enumerate(sam_db):
            output_succinct(sam_record, list(decode_ot(sam_record)))
    elif args.mode == 'complete':
        print(COMPLETE_HEADER)
        for sam_i, sam_record in enumerate(sam_db):
            output_complete(decode_ot(sam_record))
