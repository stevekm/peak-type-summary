#!/usr/bin/env Rscript

# Run the ChIPSeeker pipeline on .BED files


# ~~~~~ PACKAGES ~~~~~ # 
library("optparse")
library("tools")
library("grid")
library("gridExtra")

timestamp <- format(Sys.time(), "%Y-%m-%d-%H-%M-%S")
logfile <- file.path(".", sprintf("report_log.%s.txt", timestamp))

# default value, overwrite in report
is_report <- FALSE

# ~~~~~ FUNCTIONS ~~~~~ # 
tsprintf <- function(fmt, ...){
    # print a formatted message with timestamp
    # base message
    m <- sprintf(fmt, ...)
    # message with timestamp
    tm <- sprintf('[%s] %s', format(Sys.time(), "%H:%M:%S"), m)
    # emit message
    message(tm)
    # add to log file
    if(isTRUE(is_report)) cat(sprintf("%s\n", tm), file = logfile, append = TRUE)
}

msprintf <- function(fmt, ...) {
    message(sprintf(fmt, ...))
}

mycat <- function(text){
    # print formatted text in Rmd
    cat(gsub(pattern = "\n", replacement = "  \n", x = text))
}

make_filename <- function (input_file, new_ext, out_dir = FALSE) {
    # Convert '/path/to/file.bed' to '/path/to/file_annotations.tsv'
    old_ext <- file_ext(input_file)
    filename_base <- gsub(pattern = sprintf('.%s$', old_ext), replacement = '', x = basename(input_file))
    filename_new <- sprintf('%s.%s', filename_base, new_ext)
    new_path <- file.path(dirname(input_file), filename_new)
    if(out_dir != FALSE){
        new_path <- file.path(out_dir, new_path)
        dir.create(path = dirname(new_path), recursive = TRUE, showWarnings = FALSE)
    }
    return(new_path)
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
        tsprintf("WARNING: Input file do not exist:\n%s\nFile will not be processed\n\n", input_file)
        return(FALSE)
    }
    all_bed <- all(grepl(pattern = '*.bed$', x = basename(input_file)))
    if ( ! isTRUE(all_bed)) {
        tsprintf("WARNING: Input file is not .bed:\n%s\nFile will not be processed\n\n", input_file)
        return(FALSE)
    }
    all_min_linenum <- all(sapply(input_file, check_numlines))
    if ( ! isTRUE(all_min_linenum)) {
        tsprintf("WARNING: Input file does not have enough lines:\n%s\nFile will not be processed\n\n", input_file)
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




chipseeker_pipeline <- function(bed_file, sampleID, tss_dist, txdb, out_dir = FALSE, annoDb = "org.Hs.eg.db"){
    # the pipeline for ChIPSeeker peak annotations and plots
    
    # list to hold output data
    output_list <- list()
    output_list[["plots"]] <- list()
    
    # read in the peaks
    tsprintf("Reading peaks file...\n\n")
    peak <- readPeakFile(bed_file)
    
    
    # make coverage plot
    peaks_coverage_plot_file <- make_filename(input_file = bed_file, new_ext = 'coverage.pdf', out_dir = out_dir)
    tsprintf("Making Chrom Coverages plot:\n%s\n\n", peaks_coverage_plot_file)
    sample_title <- paste0(sampleID, " ChIP Peaks over Chromosomes")
    pdf(file = peaks_coverage_plot_file)
    
    output_list[["plots"]][["covplot"]] <- covplot(peak, weightCol = "V5", title = sample_title)
    print(output_list[["plots"]][["covplot"]]) # title = "ChIP Peaks over Chromosomes"
    
    dev.off()
    
    # make annotations
    tsprintf("Getting peak annotations...\n\n")
    peakAnno <- annotatePeak(peak, tssRegion = c(-tss_dist, tss_dist), 
                             TxDb = txdb, 
                             annoDb = annoDb)
    
    peak_anno_table_file <- make_filename(input_file = bed_file, new_ext = 'peak_anno.tsv', out_dir = out_dir)
    tsprintf("Saving table:\n%s\n\n", peak_anno_table_file)
    output_list[["peakAnno"]] <- peakAnno
    write.table(output_list[["peakAnno"]], quote=FALSE, sep="\t", row.names =FALSE, file=peak_anno_table_file)
    
    peak_anno_stats_file <- make_filename(input_file = bed_file, new_ext = 'peak_anno_stats.tsv', out_dir = out_dir)
    tsprintf("Saving table:\n%s\n\n", peak_anno_stats_file)
    output_list[["annoStat"]] <- peakAnno@annoStat
    write.table(output_list[["annoStat"]], quote=FALSE, sep="\t", row.names =FALSE, file=peak_anno_stats_file)
    
    tss_dist_file <- make_filename(input_file = bed_file, new_ext = 'tss_distance.txt', out_dir = out_dir)
    tsprintf("Saving table:\n%s\n\n", tss_dist_file)
    output_list[["tss_dist"]] <- tss_dist
    cat(as.character(output_list[["tss_dist"]]), file = tss_dist_file)
    
    
    # make annot pie chart
    anno_piechart_plot_file <- make_filename(input_file = bed_file, new_ext = 'anno-piechart.pdf', out_dir = out_dir)
    tsprintf("Making Peak Anno pie chart:\n%s\n\n", anno_piechart_plot_file)
    sample_title <- paste0("\n\n", sampleID, " Peak Types")
    
    output_list[["plots"]][["plotAnnoPie_sample_title"]] <- sample_title
    
    pdf(file = anno_piechart_plot_file, height = 8, width = 8)
    output_list[["plots"]][["plotAnnoPie"]] <- plotAnnoPie(peakAnno, main = sample_title)
    print(output_list[["plots"]][["plotAnnoPie"]])
    dev.off()
    
    
    # make UpSet plot
    tsprintf("Making Upset plot...\n\n")
    # upset_plot_file <- file.path(output_directory, sprintf("%s_upsetplot.pdf", sampleID))
    upset_plot_file <- make_filename(input_file = bed_file, new_ext = 'upsetplot.pdf', out_dir = out_dir)
    sample_title <- paste0(sampleID, " Peak Overlaps")
    output_list[["plots"]][["upsetplot_sample_title"]] <- sample_title
    pdf(file = upset_plot_file, width = 9, height = 4.5, onefile = F)
    output_list[["plots"]][["upsetplot"]] <- upsetplot(peakAnno, vennpie=TRUE)
    print(output_list[["plots"]][["upsetplot"]])
    text(x = 0, y = 1, sample_title) # add a title
    dev.off()
    
    return(output_list)
}


summarize_beds <- function(bed_files, tss_dist, id_dirname = FALSE, out_dir = FALSE, txdb_file = "TxDb.Hsapiens.UCSC.hg19.knownGene.Rdata") {
    # run the ChIPSeeker pipeline on all the .bed files
    
    
    # ~~~~~ VALIDATION ~~~~~ # 
    # check to make sure at least one files has >0 lines before we try to load data, because it takes a while to load
    any_min_linenum <- any(sapply(names(bed_files), check_numlines))
    if ( ! isTRUE(any_min_linenum)) {
        tsprintf("ERROR: No input files have enough lines to be processed\nExiting...\n\n")
        quit()
    }
    
    # ~~~~~ LOAD DATA ~~~~~ # 
    message("\nLoading packages and data...\n")
    # source("http://bioconductor.org/biocLite.R")
    # biocLite("ChIPseeker")
    suppressPackageStartupMessages(library("ChIPseeker"))
    suppressPackageStartupMessages(library("clusterProfiler"))
    suppressPackageStartupMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene"))
    
    # load data from saved Rdata file if it exists
    if(file.exists(txdb_file)){
        # txdb <- load(txdb_file)
        # TODO: fix this, does not load correctly
        txdb <- get("TxDb.Hsapiens.UCSC.hg19.knownGene")
    } else {
        txdb <- get("TxDb.Hsapiens.UCSC.hg19.knownGene")
    }
    
    
    
    # ~~~~~ RUN ~~~~~ # 
    # list to hold the output data
    output_list <- list()
    
    # iterate over bed files
    tsprintf('\n------------------------------\n')
    tsprintf('\n------------------------------\n')
    for(i in seq_along(bed_files)){
        bed_file <- names(bed_files[i])
        process_file <- bed_files[i] # TRUE or FALSE
        
        files_errors <- character()
        files_warnings <- character()
        
        tsprintf("Input File:\n%s\n\n\nFile will be processed:\n%s\n\n", bed_file, process_file)
        if(isTRUE(as.logical(process_file))){
            
            sampleID <- get_sampleID(input_file = bed_file, id_dirname = id_dirname)
            
            tsprintf("Sample ID:\n%s\n\n\n", sampleID)
            
            result <- tryCatch(
                {
                    tsprintf("Running ChIPSeeker pipeline for sample %s, file:\n%s\n\n", sampleID, bed_file)
                    chipseeker_pipeline_output <- chipseeker_pipeline(bed_file = bed_file, 
                                                                      sampleID = sampleID, 
                                                                      tss_dist = tss_dist, 
                                                                      txdb = txdb, 
                                                                      out_dir = out_dir)
                    output_list[[sampleID]] <- list()
                    output_list[[sampleID]][["pipeline_output"]] <- chipseeker_pipeline_output
                    output_list[[sampleID]][["bed_file"]] <- bed_file
                    output_list[[sampleID]][["process_file"]] <- process_file
                    
                },
                error = function(cond) {
                    tsprintf("An error occured while running ChIPSeeker pipeline for sample %s, file:\n%s\n\n", sampleID, bed_file)
                    message("Original error message:")
                    message(cond)
                    return("error")
                },
                warning = function(cond) {
                    tsprintf("An warning occured while running ChIPSeeker pipeline for sample %s, file:\n%s\n\n", sampleID, bed_file)
                    message("Original warning message:")
                    message(cond)
                    return("warning")
                },
                finally={
                    tsprintf("Finished running ChIPSeeker pipeline for sample %s, file:\n%s\n\n", sampleID, bed_file)
                }
            )
            
            if(result == "error") files_errors <- c(files_errors, bed_file)
            if(result == "warning") files_warnings <- c(files_warnings, bed_file)
        }
        tsprintf('\n------------------------------\n')
    }
    
    tsprintf('The following files had errors:\n')
    tsprintf('%s\n', files_errors)
    
    tsprintf('The following files had warnings:\n')
    tsprintf('%s\n', files_warnings)
    
    cat(files_errors, file = "file_errors.txt", append = TRUE)
    cat(files_warnings, file = "file_warnings.txt", append = TRUE)
    
    return(output_list)
}


sysinfo <- function(){
    # print custom information about the system
    
    # check if 'mycat' is loaded in case I copy/pasted this from elsewhere
    if( ! exists('mycat')) mycat <- function(text){cat(gsub(pattern = "\n", replacement = "  \n", x = text))}
    
    # system info for use on Linux with GNU tools installed
    # mycat(sprintf("System:\n%s\n%s", system("hostname", intern = TRUE), system("uname -srv", intern = TRUE)))
    # mycat(sprintf("System user:\n%s", system("whoami", intern = TRUE)))
    
    # dir
    # mycat(sprintf("System location:\n%s", system('pwd', intern = T, ignore.stderr = TRUE)))
    mycat(sprintf("System location:\n%s", getwd()))
    
    
    # repo info
    mycat(sprintf("Git Remote:\n%s\n", system('git remote -v', intern=T)))
    mycat(sprintf("Git branch and commit\n%s", system('printf "%s: %s" "$(git rev-parse --abbrev-ref HEAD)" "$(git rev-parse HEAD)"', 
                                                      intern = TRUE,  ignore.stderr = TRUE)))
    
    # date time
    mycat(sprintf("Time and Date of report creation:\n%s", system("date", intern = TRUE)))
    
    # R system info, packages, etc
    print(sessionInfo())
    print(Sys.info())
}
