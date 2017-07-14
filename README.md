# peak-type-summary
ChIPSeeker peaks plotting and summarizing script

R script to run the ChIPSeeker peak annotation and summary plotting pipeline on .bed files. 

# Usage

```bash
./peak-summary.R /path/to/Sample1_peaks.bed /path/to/Sample2_peaks.bed 
```

## Options

- `-d`, `--dir`: Dir mode; treat input items as directories in which to search for .bed files
- `-id-dirname`: Use the directory name of each input file as the sample ID, instead of the filename
- `--out-dir`: Parent directory to save the output to
- `--tss-dist`: TSS region distance

# Examples

Annotate & plot .bed files

```bash
$ ./peak-summary.R example-data/Sample1.bed example-data/Sample2.bed
```

Annotate all .bed files in a directory, get sample IDs from the dir name

```bash
$ ./peak-summary.R example-data/Sample3 --dir
```

# Software
- Tested with R version 3.2.3 and 3.3.0, with the following packages:
  - `ChIPseeker_1.6.7`
  - `clusterProfiler_2.4.3`
  - `TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2`
  - `optparse_1.3.2 `
