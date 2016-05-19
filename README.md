# cbiC1

## Tools for the Computational Biology Institute at The George Washington Univeristy

### Using cbiC1

The cbiC1 module adds a few useful utilities and wrapper scripts to your path.
Most of these are designed to work on a pipe (using `stdin` and `stdout`) so they can be
easily incorporated into workflows without creating temporary files. In order to load the
module you will also need to use the CBI modules on C1

```bash
# Add the CBI modules to your search path
module use /groups/cbi/shared/modulefiles 
# Load the cbiC1 module
module load cbiC1
```

### Utilities

#### `fastqCount`

**Count number of reads and bases in fastq files.** Finds files that end in .fastq or .fq,
then compute read and base counts for each file. Uses linux utilities (sed and wc), so it
is reasonably fast. Run without arguments to find files in working directory, or provide
an expandable path as command line argument.

###### Example:
```bash
# Count files in current working directory
fastqCount

# Count files in another directory
fastqCount readdir/*
```

#### `inlineFastq2Fasta`

**Convert fastq file on `stdin` to fasta file on `stdout`.** Easy enough.

###### Example:
```bash
cat reads.fastq | inlineFastq2Fasta > reads.fasta
```

#### `inlineFasta2Fastq`

**Convert fasta file on `stdin` to fastq file on `stdout`.** All bases receive quality score of 40 (I).

###### Example:
```bash
cat reads.fasta | inlineFasta2Fastq > reads.fastq
```

#### `deinterleave_fasta`

**Deinterleaves a fasta file of paired reads into two fasta files.** Interleaved file should be
passed on `stdin`, provide output file names as command line arguments. Uses linux utilities 
so it is reasonably fast. Output files are compressed if 3rd command line argument
is "compress" (requires pigz).

###### Example:
```bash
cat interleaved.fasta | deinterleave_fasta fwd.fasta rev.fasta
```

#### `deinterleave_fastq`

**Deinterleaves a fastq file of paired reads into two fastq files.** Interleaved file should be
passed on `stdin`, provide output file names as command line arguments. Uses linux utilities,
can deinterleave 100 million paired reads (200 million total reads; a 43Gbyte file), in 
memory (/dev/shm), in 4m15s (255s). Output files are compressed if 3rd command line argument
is "compress" (requires pigz).

###### Example:
```bash
cat interleaved.fastq | deinterleave_fastq fwd.fastq rev.fastq
```

#### `cbi_guess_encoding`

**Guess encoding (Phred+33 or Phred+64) of fastq file.** Calculates the range of quality scores
on the first 10K reads and returns the best guess. Fastq can be passed on `stdin` and encoding
returned on `stdout`

###### Example:
```bash
cbi_guess_encoding < reads.fastq
# Phred+33
```

#### `cbi_sra_metadata`

**Retrieve metadata for all runs in SRA project.** Uses Entrez to find all runs for
a given SRA project accession (begins with "SRP") and retrieve XML metadata for all
runs. In addition to the XML files, generates a sample matrix with attribute values for
each sample and a list of URLs for downloading the SRA files.

###### Example:
```bash
# Save data for SRA runs in "metadata" directory
cbi_sra_metadata SRP055981

# Download all the SRA files using wget
wget -i sample_urls.txt
```

### Wrappers

#### `cbi_fastq_filter`

Need to write description

#### `cbi_prinseq_wrapper`

**Creates SLURM jobs for running [PRINSEQ](http://prinseq.sourceforge.net/) in parallel.**
Sequence read files can be discovered (in current directory with *.fq or *.fastq extension) or
specified using a file of file names; paired files are automatically detected. Several processes
are spawned for each SLURM job for multi-core utilization, and may be split into several jobs
depending on the number of inputs.

For detailed information:
```bash
cbi_prinseq_wrapper -h
```

###### Example:
```bash
# Find all *.fq and *.fastq in the current directory, identify
# possible paired files, and trim/filter using Illumina presets.

cbi_prinseq_wrapper --illumina

# For the files listed in "samples.fofn", trim and filter using
# IonTorrent presets, assuming all are unpaired

cbi_prinseq_wrapper --fofn samples.fofn --single --iontorrent
```

#### `cbi_fastq_dump_wrapper`

**Creates SLURM jobs for running [fastq-dump](http://www.ncbi.nlm.nih.gov/books/NBK158900/) in parallel.**
SRA files can be discovered (in current directory with *.sra extension) or
specified using a file of file names. Files are extracted to FASTQ, mate pairs are split 
into seperate files, and output is gzipped. This wrapper replaces `extractSRA.sh`.

###### Example:
```bash
# Extract all *.sra files in current directory
cbi_fastq_dump_wrapper
# [--- job: fastq-dump01 ---] Submitted batch job 1789048. Running 8 fastq-dump processes.
```
