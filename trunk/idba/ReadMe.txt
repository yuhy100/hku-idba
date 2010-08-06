Installation Guide

Exract the package, then use make to compile the source code.
$ tar -xvzf idba-0.1.tar.gz
$ make


Usage

The tool kit contain IDBA (single end) and IDBA (pair end).

usage: script/runIdba --read read-file --prefix prefix --ref ref-file
       [--build] [--assemble] [--align] [--mate] [--validate]
       [--prefixLength p] [--kmer k] [--maxK k2]
       [--minCount m] [--minContig l] [--cover c]

read-file:  The input read file.
prefix:     The path prefix for intermediate files and output file.
ref-file:   The reference file.
build:      Build hash table.
assemble:   Perform single end IDBA.
align:      Align the reads to contigs for next step.
mate:       Perform pair end IDBA.
validate:   Validate the contigs by aligning them onto the reference.
p:          The length of prefix that is used for building hash table.
k:          The minimum k IDBA start with.
k2:         The maximum k IDBA end with.
m:          The filtering threshold.
l:          The minimum length of contigs that will be outputed.
c:          The minimum average coverage of contigs.

The output contigs will be found in file prefix.contig and prefix.contig-mate.
The former is single end output, the latter is pair end output.
runIdba script can perform five tasks: build, assemble, align, mate, validate.
You can specify at one or more task in one run.

Example (using the sample input)
$ script/runIdba --read Lacto.reads-mate-30-0.01-75 --prefix lacto --ref Lacto.fasta --build --assemble --align --mate --validate

Comment
IDBA(pair end) requires pair end reads stored in single file and a pair of
reads is in consecutive two lines. If not, please use mergeReads to merge two
read files to single file.
$ bin/mergeReads read-file1 read-file2 merge-read-file

If the data is in fastq format, please use fq2fa tool to do conversion first.
$ bin/fq2fa fq-file fa-file

Reads with same length are prefered. normReads tool can help to normlize reads.
$ bin/normReads read-file output-read-file --length l --mate

