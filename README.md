# peak-type-summary
ChIPSeeker peaks plotting and summarizing script

R script to run the ChIPSeeker peak annotation and summary plotting pipeline on .bed files. 

__[ Example output files can be seen [here](https://github.com/stevekm/peak-type-summary/tree/d443be906c8a6bfc52c336b31761ef7230a64328/example-output/example-data) ]__

__[ Example report output can be seen [here](https://cdn.rawgit.com/stevekm/peak-type-summary/7285ef4bba3622b1f7edd02161e3b8ee8d895a75/peak-type-summary-report.html)]__

# Usage

To run the basic script on selected files:

```bash
./peak-summary.R /path/to/Sample1_peaks.bed /path/to/Sample2_peaks.bed 
```

To compile a report on a directory of files:

```bash
./compile-peak-type-summary-report.R -d example-data
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

Annotate all .bed files in a directory, and save results to a different output directory

```bash
$ ./peak-summary.R example-data/ -d --out-dir example-output
```

# Software
- Tested with R version 3.2.3 and 3.3.0, with the following packages:
  - `ChIPseeker_1.6.7`
  - `clusterProfiler_2.4.3`
  - `TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2`
  - `optparse_1.3.2 `
