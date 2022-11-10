from Bio import SeqIO
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generates a set of kmers for processing by Guidescan2."
    )

    parser.add_argument(
        "fasta",
        type=str,
        help="FASTA file to use as a reference for kmer generation."
    )

    parser.add_argument(
        "--pam",
        help="Protospacer adjacent motif to match.",
        default="NGG"
    )

    parser.add_argument(
        "--kmer-length",
        help="Length of kmers to generate.",
        type=int,
        default=20
    )

    parser.add_argument(
        "--min-chr-length",
        help="Minimum chromosome length to consider for kmer generation.",
        type=int,
        default=10_000
    )

    parser.add_argument(
        "--prefix",
        help="Prefix to use for kmer identifiers.",
        default=""
    )

    parser.add_argument(
        "--start",
        help="Match PAM at start of kmer instead at end (default).",
        action='store_true'
    )

    return parser.parse_args()

NUCS = list("ACTG")
NUC_MAP = {"A": "T", "T": "A", "C": "G", "G": "C"}

def revcom(dna):
    return "".join(list(map(lambda n: NUC_MAP[n], list(dna)))[::-1])

def generate_pam_set(pam):
    pam_stack = [pam]

    while any(["N" in pam for pam in pam_stack]):
        pam = pam_stack.pop(0)

        if "N" not in pam:
            pam_stack.append(pam)
            continue 

        for nuc in NUCS:
            pam_stack.append(pam.replace("N", nuc, 1))

    return pam_stack

def find_kmers(pam, k, chrm, forward=True, end=True):
    index = 0

    while True:
        index = chrm.find(pam, index)

        if index == -1:
            break 

        if end:
            if forward:
                kmer = chrm[index - k:index]
                position = index - k
            else:
                kmer = chrm[index + len(pam):index + k + len(pam)]
                position = index
        else:
            if forward:
                kmer = chrm[index + len(pam):index + k + len(pam)]
                position = index
            else:
                kmer = chrm[index - k:index]
                position = index - k

        index += 1

        if position < 0:
            continue

        # Return the 1-indexed position to caller
        yield kmer.upper(), position+1

def find_all_kmers(pam, k, chrm, end=True):
    pam_set = generate_pam_set(pam)
    rev_pam_set = list(map(revcom, pam_set))

    for p in pam_set:
        for kmer, pos in find_kmers(p, k, chrm, end=end):
            yield {"sequence" : kmer, "position" : pos, "pam" : pam, "sense": "+"}

    for p in rev_pam_set:
        for kmer, pos in find_kmers(p, k, chrm, forward=False, end=end):
            yield {"sequence" : revcom(kmer), "position" : pos,
                   "pam" : pam, "sense": "-"}

def output_kmer(prefix, chrm_name, kmer):
    identifier = f"{prefix}{chrm_name}:{kmer['position']}:{kmer['sense']}"
    row = [identifier, str(kmer['sequence']),
           kmer['pam'], chrm_name,
           str(kmer['position']), kmer['sense']]
    print(",".join(row))

if __name__ == "__main__":
    args = parse_arguments()

    print("id,sequence,pam,chromosome,position,sense")
    for record in SeqIO.parse(args.fasta, "fasta"):
        if len(record) < args.min_chr_length:
            continue

        for kmer in find_all_kmers(args.pam, args.kmer_length, record.seq, end=not args.start):
            output_kmer(args.prefix, record.name, kmer)
