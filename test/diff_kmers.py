import sys

if __name__ == "__main__":
    f1 = open(sys.argv[1])
    f2 = open(sys.argv[2])

    k = 23

    lines1 = f1.readlines()
    lines1 = [line[0:k] for line in lines1]

    lines2 = f2.readlines()
    lines2 = [line[0:k] for line in lines2]

    diff1 = set(lines1).difference(set(lines2))
    diff2 = set(lines2).difference(set(lines1))

    print(f"Only in {sys.argv[1]}:\n{list(diff1)[:10]}")
    print(f"Only in {sys.argv[2]}:\n{list(diff2)[:10]}")
    print(f"Number different: {len(diff1)}")
    print(f"Number different: {len(diff2)}")
