import sys

if len(sys.argv) < 3:
    print("usage: python chr2acc.py [chr2acc.txt] [kmers.csv]")
    sys.exit(1)

chr2acc = {}
with open(sys.argv[1]) as f:
    next(f)
    for l in f:
        words = l.split() 
        chr2acc[words[0]] = words[1]

with open(sys.argv[2]) as f:
    print(next(f), end='')
    for l in f:
        words = l.split(',')
        words[3] = chr2acc[words[3][3:]] # strips 'chr' prefix
        print(','.join(words), end='')
