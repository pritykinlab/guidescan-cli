# guidescan-cli

Guidescan-cli is a tool for memory efficient gRNA database generation.
It takes as input a genome along with a targetable space specified as
a PAM and outputs a database of gRNAs in SAM format that have few
off-targets.

# Installation 

Installation is rather straightforward once you have
[CMake](https://cmake.org/) installed. Currently, the software has
been tested on Ubuntu 19.04 and CentOS 7, but it should work on
other systems.

To build, run, 

``` shell
$ mkdir build; cd build
$ cmake -DCMAKE_BUILD_TYPE=RELEASE ..
$ make
```

from the root directory of the project. The output will be a binary in
build/bin. If you want the executable to be statically linked for
deployment onto a compute cluster, instead run,

``` shell
$ cmake -DCMAKE_BUILD_TYPE=RELEASE -DLINK_STATICALLY=ON ..
$ make
```

The dependencies are listed below:
- [CMake](https://cmake.org/) version >= 3.1.0
- C++ compiler supports C++11 standard
- pthreads

# Usage

To use the tool simply run,

```
$ guidescan --help
Guidescan all-in-one interface.

Usage: guidescan [OPTIONS] SUBCOMMAND

Options:
  -h,--help                   Print this help message and exit

Subcommands:
  build                       Builds a gRNA database over the given genome.
  kmers                       Generates a list of kmers for a specific PAM written and writes them to stdout.
```

There are two subcommands `build` and `kmers`, each with their own
purpose. 

## Build

The subcommand `build` constructs a gRNA database over a given genome,
input in FASTA file format. By default, it lets all kmers matching the
input PAM be potential guideRNAs, and prunes those that have bad
off-targets. However, you can explicitly specify the kmers to check
with the option `--kmers-file`.

The usage for the subcommand is as follows:

``` shell
$ guidescan build --help
Builds a gRNA database over the given genome.
Usage: guidescan build [OPTIONS] genome

Positionals:
  genome TEXT:FILE REQUIRED   Genome in FASTA format

Options:
  -h,--help                   Print this help message and exit
  -k,--kmer-length UINT=20    Length of kmers excluding the PAM
  -n,--threads UINT=8         Number of threads to parallelize over
  -p,--pam TEXT=NGG           Main PAM to generate gRNAs and find off-targets
  -a,--alt-pam TEXT=[NAG] ... Alternative PAMs used to find off-targets
  -m,--mismatches UINT=3      Number of mismatches to allow when finding off-targets
  -f,--kmers-file TEXT:FILE   File containing kmers to build gRNA database over, if not specified, will generate the database over all kmers with the given PAM
  -o,--output TEXT REQUIRED   Output database file.
```

## Kmers

The subcommand `kmers` finds all kmers matching a given PAM in the
genome and outputs them to a file. This file, or a subset of it, can
then be input to the `build` command. This is particularly useful if
you want to do some special processing on the set of kmers used to
build the gRNA database or analyze a small set of kmers instead of the
whole genome.

As an example, suppose one wanted to parallelize gRNA database
construction over a bunch of nodes (as it already runs in parallel on one
node). 

One could run,

``` shell
$ guidescan kmers hg38.fasta -o hg38.kmers
$ split -d -l 100000 hg38.kmers hg38.kmers_
```

to split the set of kmers into seperate kmers files of 100,000 lines long.

Then, one could run

``` shell
$ guidescan build hg38.fasta --kmers-file hg38.kmers_XX
```

seperately on each node, and merge the databases output with `cat` after
stripping the headers. Of course, more complex tasks could be performed as well.
The format for the kmers file is extremely simple, so that it can be easily
processed and generated with familiar tools.

The usage for the subcommand is as follows:

``` shell
$ guidescan kmers --help
Generates a list of kmers for a specific PAM written and writes them to stdout.
Usage: ./bin/guidescan kmers [OPTIONS] genome

Positionals:
  genome TEXT:FILE REQUIRED   Genome in FASTA format

Options:
  -h,--help                   Print this help message and exit
  -k,--kmer-length UINT=20    Length of kmers excluding the PAM
  -p,--pam TEXT=NGG           PAM to generate kmers for
  -o,--output TEXT REQUIRED   Output kmers file.
```
