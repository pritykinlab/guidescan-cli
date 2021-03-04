# guidescan-cli

Guidescan-cli is a tool for memory efficient gRNA off-target
enumeration. 

# Installation 

## Singularity Container

There is a [Singularity](https://sylabs.io/docs/) container that
bundles all dependencies for the application together. This is the
easiest and most portable method to set up. It is **highly recommended**
way to build the application this way. Additionally, the auxillary
scripts (not necessary for most purposes) are only supported from
within the the container environment.

You can either build the container image from scratch or download it
here. To build the container from scratch, simply run, 

```shell
$ sudo singularity build guidescan.sif Singularity
```

from the project directory `container/`.

This take a few minutes and result in a container file
`guidescan.sif`. Though root access is required to build the
container, it is not required for execution.

Once the container is built (or downloaded) you can run the
`guidescan` tool with,

```
singularity exec guidescan.sif guidescan [ARGS]
```

refer to usage section for a description of the arguments.

## Manual Build

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

# Standard Usage

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
  http-server                 Starts a local HTTP server to receive gRNA processing requests.
```

There are three subcommands `build`, `kmers`, and `http-server`. For
the majority of use-cases, the first two commands are the most useful.

## Build

The subcommand `build` constructs a gRNA database over a given genome,
input in FASTA file format. By default, it lets all kmers matching the
input PAM be potential guideRNAs, and prunes those that have bad
off-targets. However, you can explicitly specify the kmers to check
with the option `--kmers-file`. Because of the time taken to build a
genome-wide database, it is **almost always** desirable to first
generate the kmers prior to building the database. 

The usage for the subcommand is as follows:

``` shell
$ guidescan build --help
Builds a gRNA database over the given genome.
Usage: ./bin/guidescan build [OPTIONS] genome

Positionals:
  genome TEXT:FILE REQUIRED   Genome in FASTA format

Options:
  -h,--help                   Print this help message and exit
  --min-chr-length UINT=1000  Minimum length of chromosomes to consider for gRNAs
  -k,--kmer-length UINT=20    Length of kmers excluding the PAM
  -n,--threads UINT=8         Number of threads to parallelize over
  -p,--pam TEXT=NGG           Main PAM to generate gRNAs and find off-targets
  -a,--alt-pam TEXT=[NAG] ... Alternative PAMs used to find off-targets
  -m,--mismatches UINT=3      Number of mismatches to allow when finding off-targets
  -t,--threshold UINT=1       Filters gRNAs with off-targets at a distance at or below this threshold
  -f,--kmers-file TEXT:FILE   File containing kmers to build gRNA database over, if not specified, will generate the database over all kmers with the given PAM
  -o,--output TEXT REQUIRED   Output database file.
```

which again, can (and should) be executed as,

```shell
$ singularity exec guidescan.sif guidescan --help
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
$ singularity exec guidescan.sif guidescan kmers hg38.fasta -o hg38.kmers
$ split -d -l 100000 hg38.kmers hg38.kmers_
```

to split the set of kmers into seperate kmers files of 100,000 lines long.

Then, one could run

``` shell
$ singularity exec guidescan.sif guidescan build hg38.fasta --kmers-file hg38.kmers_XX
```

seperately on each node, and merge the databases output with `cat` after
stripping the headers. Of course, more complex tasks could be performed as well.
The format for the kmers file is extremely simple, so that it can be easily
processed and generated with familiar tools.

The usage for the subcommand is as follows:

``` shell
$ singularity exec guidescan.sif guidescan kmers --help
Generates a list of kmers for a specific PAM written and writes them to stdout.
Usage: /usr/bin/guidescan kmers [OPTIONS] genome

Positionals:
  genome TEXT:FILE REQUIRED   Genome in FASTA format

Options:
  -h,--help                   Print this help message and exit
  --min-chr-length UINT=1000  Minimum length of chromosone for kmers to be included in output
  -k,--kmer-length UINT=20    Length of kmers excluding the PAM
  -p,--pam TEXT=NGG           PAM to generate kmers for
  -o,--output TEXT REQUIRED   Output kmers file.
```

# Auxillary Scripts

There are two auxillary scripts located under `scripts/`. These are
not strictly necessary, but were useful for some of the analysis that
we performed. Because of their niche use case and the tricky
dependencies required, we only support them from inside the
Singularity container. However, with a little bit of work to set up
the environment, they can be run standalone.

## append_scores script

This script takes the SAM file output from `guidescan build` and
appends two sets of scores in the tags field of the SAM file, printout
the output to stdout. These scores are our *specificity* score derived
from Doench's CFD score and Doench's *Rule Set 2 Score*, which we name
*cutting efficiency*. In the resultant file,

- Tag `cs` corresponds to *specificity*
- Tag `ds` corresponds to *cutting efficiency*

To run the script simply execute,

```shell
singularity exec guidescan.sif append_scores [ORGANISM_FASTA] [GUIDESCAN_DB_SAM]
```

## build_dbs script

This script generates guidescan databases on the distributed
architecture *Platform Load Sharing Facility* also known as LSF. It
takes advantage of the highly parallel architecture to be able to
perform a genome wide database construction in a reasonable amount of
time. Unfortunately, this script is highly specific to the LSF
architecture, but the procedure taken by the script can be performed
manually.

The script performs the following steps:

1. Constructs genome index 
2. Constructs and randomizes KMER file
3. Splits KMER file into N parts
4. Builds databases in parallel for each of the N parts
5. Appends scores in parallel to N constructed databases
6. Merges databases together into a single SAM file

As an example, to run the script for 4 days on *N* nodes each
utilizing *k* cores execute,

```shell
singularity exec guidescan.sif build_dbs --num_jobs N --num_cores k --time 96:00 [ORGANISM_FASTA]
```

The BSUB files will be output to a directory `bsub_files` and the
results will be output to a directory `results`. The final output is a
genome wide guidescan database in SAM format named
`results/guide_db.sam`. Another file named `state_[JOBID].txt` will
also be generated and contain the status of the job.
