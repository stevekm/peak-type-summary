#!/usr/bin/env Rscript

## USAGE: compile_RMD_report.R report.Rmd
## EXAMPLE: ./compile-report.R -n "peaks.by_sample.macs_broad" --height 42 --width 12 results/peaks.by_sample.macs_broad
## Requires pandoc version 1.13+
# module load pandoc/1.13.1

message("
        Remember to have pandoc/1.13.1+ loaded!
        module load pandoc/1.13.1
        ")

# ~~~~~ PACKAGES ~~~~~ # 
library("optparse")
library("tools")

# ~~~~~ SCRIPT ARGS ~~~~~ # 
option_list <- list(
    make_option(c("-n", "--name"), type="character", default=FALSE,
                dest="report_name", help="A different output name to use for the report file (excluding file extension)"),
    make_option(c("-d", "--dir"), action="store_true", default=FALSE,
                dest="dir_mode", help="Treat input items as directories to be searched for .bed files"),
    make_option(c("--id-dirname"), action="store_true", default=FALSE,
                dest="id_dirname", help="Take the sample ID from the .bed file's dirname, not its basename"),
    make_option(c("--out-dir"), type="character", default=FALSE,
                dest="out_dir", help="Path to the parent output directory. Will be created and appended to the input item's file path"),
    make_option(c("--tss-dist"), type="numeric", default=3000,
                dest = "tss_dist", help="TSS distance to use [default %default]",
                metavar="tss-dist"),
    make_option(c("--height"), type="numeric", default=8,
                dest = "plot_height", help="Height for boxplot [default %default]",
                metavar="plot_height"),
    make_option(c("--width"), type="numeric", default=8,
                dest = "plot_width", help="Width for boxplot [default %default]",
                metavar="plot_width")
)
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments = TRUE)

report_name <- opt$options$report_name
dir_mode <- opt$options$dir_mode 
out_dir <- opt$options$out_dir
tss_dist <- opt$options$tss_dist
id_dirname <- opt$options$id_dirname
plot_height <- opt$options$plot_height
plot_width <- opt$options$plot_width

input_items <- opt$args
input_dir <- opt$args[1] # input directory
input_items <- opt$args

# get script dir
args <- commandArgs(trailingOnly = F)  
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
save.image(file.path(scriptPath, "loaded_args.Rdata"))

#  ~~~~~ COMPILE ~~~~~ #
# the default report filename
Rmdfile <- "peak-type-summary-report.Rmd"

# rename the report file if arg was passed
if(report_name != FALSE){
    old_ext <- file_ext(Rmdfile)
    new_Rmdfile <- sprintf('%s.%s', report_name, old_ext)
    file.copy(from = Rmdfile, to = new_Rmdfile, overwrite = TRUE)
    Rmdfile <- new_Rmdfile
}

# compile the report with params passed
rmarkdown::render(input = Rmdfile, params = list(input_dir = input_dir,
                                                 input_items = input_items,
                                                 dir_mode = dir_mode, 
                                                 out_dir = out_dir,
                                                 tss_dist = tss_dist,
                                                 plot_height = plot_height,
                                                 plot_width = plot_width, 
                                                 is_report = TRUE))
