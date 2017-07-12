#!/usr/bin/env Rscript

# Run the ChIPSeeker pipeline on .BED files


# ~~~~~ PACKAGES ~~~~~ # 
library("optparse")
library("tools")



# ~~~~~ FUNCTIONS ~~~~~ # 
msprintf <- function(fmt, ...) {
    message(sprintf(fmt, ...))
}


make_filename <- function (input_file, new_ext) {
    # Convert '/path/to/file.bed' to '/path/to/file_annotations.tsv'
    old_ext <- file_ext(input_file)
    filename_base <- gsub(pattern = sprintf('.%s$', old_ext), replacement = '', x = basename(input_file))
    filename_new <- sprintf('%s.%s', filename_base, new_ext)
    return(file.path(dirname(input_file), filename_new))
}

check_numlines <- function(input_file, min_value = 0) {
    # make sure a file has >0 lines
    has_enough_lines <- FALSE
    if (length(readLines(input_file)) > min_value) has_enough_lines <- TRUE
    return(has_enough_lines)
}


validate_file <- function(input_file) {
    # make sure that all files are .bed, and that they have >0 lines
    # validation passes if all files are .bed
    all_exist <- all(file.exists(input_file))
    if ( ! isTRUE(all_exist)) {
        msprintf("WARNING: Input file do not exist:\n%s\nFile will not be processed\n\n", input_file)
        return(FALSE)
    }
    all_bed <- all(grepl(pattern = '*.bed$', x = basename(input_file)))
    if ( ! isTRUE(all_bed)) {
        msprintf("WARNING: Input file is not .bed:\n%s\nFile will not be processed\n\n", input_file)
        return(FALSE)
    }
    all_min_linenum <- all(sapply(input_file, check_numlines))
    if ( ! isTRUE(all_min_linenum)) {
        msprintf("WARNING: Input file does not have enough lines:\n%s\nFile will not be processed\n\n", input_file)
        return(FALSE)
    }
    return(TRUE)
}



find_all_beds <- function (input_dirs) {
    # find all .bed files in the supplied dirs
    return(dir(input_dirs, pattern = '.bed', full.names = TRUE, recursive = TRUE))
}


get_sampleID <- function(input_file, id_dirname = FALSE){
    # get the sample ID for a file
    # right now just use the basename but maybe some day do something fancier here
    sampleID <- basename(input_file)
    if(isTRUE(id_dirname)) sampleID <- basename(dirname(input_file))
    return(sampleID)
}

get_sample_outdir <- function(parent_outdir, sampleID, create = TRUE){
    # make a path for the sample's output directory
    output_path <- file.path(parent_outdir, sampleID)
    if(isTRUE(create)) dir.create(output_path, recursive = TRUE)
    return(output_path)
}

summarize_beds <- function(bed_files, tss_dist, id_dirname = FALSE, annoDb = "org.Hs.eg.db") {
    # run the ChIPSeeker pipeline on all the .bed files
    
    
    # ~~~~~ VALIDATION ~~~~~ # 
    # check to make sure at least one files has >0 lines before we try to load data, because it takes a while to load
    any_min_linenum <- any(sapply(names(bed_files), check_numlines))
    if ( ! isTRUE(any_min_linenum)) {
        msprintf("ERROR: No input files have enough lines to be processed\nExiting...\n\n")
        quit()
    }
    
    # ~~~~~ LOAD DATA ~~~~~ # 
    message("\nLoading packages and data...\n")
    # source("http://bioconductor.org/biocLite.R")
    # biocLite("ChIPseeker")
    suppressPackageStartupMessages(library("ChIPseeker"))
    suppressPackageStartupMessages(library("clusterProfiler"))
    suppressPackageStartupMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene"))
    txdb <- get("TxDb.Hsapiens.UCSC.hg19.knownGene")
    
    
    # ~~~~~ RUN ~~~~~ # 
    # iterate over bed files
    msprintf('\n------------------------------\n')
    msprintf('\n------------------------------\n')
    for(i in seq_along(bed_files)){
        bed_file <- names(bed_files[i])
        process_file <- bed_files[i] # TRUE or FALSE
        
        
        
        promoter_dist <- tss_dist
        
        msprintf("Input File:\n%s\n\n\nFile will be processed:\n%s\n\n", bed_file, process_file)
        if(isTRUE(as.logical(process_file))){
            sampleID <- get_sampleID(input_file = bed_file, id_dirname = id_dirname)
            output_directory <- get_sample_outdir(parent_outdir = global_parent_outdir, sampleID = sampleID)
            
            msprintf("Sample ID:\n%s\n\n\nOutput directory:\n%s\n\n", sampleID, output_directory)
            
            msprintf("Reading peaks file...\n\n")
            peak <- readPeakFile(bed_file)
            
            
            msprintf("Making Chrom Coverages plot...\n\n")
            peaks_coverage_plot_file <- file.path(output_directory, sprintf("%s_peaks-coverage.pdf", sampleID))
            sample_title <- paste0(sampleID, " ChIP Peaks over Chromosomes")
            pdf(file = peaks_coverage_plot_file)
            covplot(peak, weightCol="V5", title = sample_title) # title = "ChIP Peaks over Chromosomes"
            dev.off()
            
            msprintf("Getting peak annotations...\n\n")
            peakAnno <- annotatePeak(peak, tssRegion = c(-promoter_dist, promoter_dist), 
                                     TxDb = txdb, 
                                     annoDb = annoDb)
            
            msprintf("Saving tables...\n\n")
            peak_anno_table_file <- file.path(output_directory, sprintf("%s_peak_anno.tsv", sampleID))
            write.table(peakAnno, quote=FALSE, sep="\t", row.names =FALSE, file=peak_anno_table_file)
            
            peak_anno_stats_file <- file.path(output_directory, sprintf("%s_peak_anno_stats.tsv", sampleID))
            write.table(peakAnno@annoStat, quote=FALSE, sep="\t", row.names =FALSE, file=peak_anno_stats_file)
            
            tss_dist_file <- file.path(output_directory, sprintf("%s_tss_distance.txt", sampleID))
            cat(as.character(promoter_dist), file = tss_dist_file)
            
            
            msprintf("Making Peak Anno pie chart...\n\n")
            anno_piechart_plot_file <- file.path(output_directory, sprintf("%s_anno-piechart.pdf", sampleID))
            sample_title <- paste0("\n\n", sampleID, " Peak Types")
            pdf(file = anno_piechart_plot_file, height = 8, width = 8)
            plotAnnoPie(peakAnno, main = sample_title)
            dev.off()
            
            
            msprintf("Making Upset plot...\n\n")
            upset_plot_file <- file.path(output_directory, sprintf("%s_upsetplot.pdf", sampleID))
            sample_title <- paste0(sampleID, " Peak Overlaps")
            pdf(file = upset_plot_file, width = 9, height = 4.5, onefile = F)
            upsetplot(peakAnno, vennpie=TRUE) 
            text(x = 0, y = 1, sample_title) # add a title
            dev.off()
            
            
            
            
            
        }
        msprintf('\n------------------------------\n')
    }
}



# ~~~~~ SCRIPT ARGS ~~~~~ # 
option_list <- list(
    make_option(c("-d", "--dir"), action="store_true", default=FALSE,
                dest="dir_mode", help="Directories with peaks to annotate"),
    make_option(c("--id-dirname"), action="store_true", default=FALSE,
                dest="id_dirname", help="Take the sample ID from the file's dirname, not its basename"),
    make_option(c("--out-dir"), type="character", default=FALSE,
                dest="out_dir", help="Path to the parent output directory. Defaults to the current directory"),
    make_option(c("--tss-dist"), type="numeric", default=3000,
                dest = "tss_dist", help="TSS distance to use [default %default]",
                metavar="tss-dist")
)
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments = TRUE)

dir_mode <- opt$options$dir_mode 
out_dir <- opt$options$out_dir
tss_dist <- opt$options$tss_dist
id_dirname <- opt$options$id_dirname

input_items <- opt$args

# get script dir
args <- commandArgs(trailingOnly = F)  
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
save.image(file.path(scriptPath, "loaded_args.Rdata"))




# ~~~~~ RUN ~~~~~ #
# default output dir
global_parent_outdir <- getwd()
if(out_dir != FALSE) global_parent_outdir <- out_dir 

if (isTRUE(dir_mode)) input_items <- find_all_beds(input_items)

validated_items <- sapply(input_items, validate_file)

summarize_beds(bed_files = validated_items, tss_dist = tss_dist, id_dirname = id_dirname)