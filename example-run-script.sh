#!/bin/bash

# example script for running the peak-summary.R script against a slightly more complicated input file directory structure

project_parent_dir="/ifs/home/kellys04/projects/SmithLab-split2"
sub_dirs="SmithLab_CTCF SmithLab_CTCF_no_input SmithLab_H3K27AC SmithLab_H3K27ME3 SmithLab_H3K9AC"
sub_paths="peaks.by_sample.macs_broad peaks.by_group.macs_broad"

printf "%s" "$(date +"%Y-%m-%d-%H-%M-%S")" >> "start.txt"

for sub_dir in $sub_dirs; do
    printf "${sub_dir}\n\n"
    sub_dir_name="$(basename $sub_dir)"
    for sub_path in $sub_paths; do
        (
            printf "${sub_path}\n\n"
        {
            find "${project_parent_dir}/${sub_dir}/pipeline/peaks/results" -path "*${sub_path}*" -type f -name "peaks.bed" | xargs ./peak-summary.R --id-dirname --out-dir "${sub_dir_name}/${sub_path}"
        } | tee log.txt
        printf "%s\t%s\t%s\n" "$(date +"%Y-%m-%d-%H-%M-%S")" "${sub_dir}" "${sub_path}" >> "end.txt"
        
        ) &
        
    done
done


# Copy over the results
#  rsync --dry-run -vhtrP SmithLab_* results_dir/

# Aggregate Plots
# 
# project_dir="/ifs/home/kellys04/projects/ChIpSeq_2017-12-31/project_notes/peak_annotation_stats"
# peak_summary_dir="${project_dir}/peaks_summaries"
# 
# ls -1 "$peak_summary_dir" | while read summary_dir; do
# (
# echo "$summary_dir"
# echo "${peak_summary_dir}/${summary_dir}"
# 
# upset_pdfs="$(find "${peak_summary_dir}/${summary_dir}" -name "upsetplot.pdf")"
# # echo $upset_pdfs
# anno_piechart_pdfs="$(find "${peak_summary_dir}/${summary_dir}" -name "anno-piechart.pdf")"
# peaks_coverage_pdfs="$(find "${peak_summary_dir}/${summary_dir}" -name "peaks-coverage.pdf")"
# 
# 
# gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=${project_dir}/${summary_dir}_upsetplot.pdf $upset_pdfs
# gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=${project_dir}/${summary_dir}_anno_piechart.pdf $anno_piechart_pdfs
# gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=${project_dir}/${summary_dir}_peaks_coverage.pdf $peaks_coverage_pdfs
# ) & 
# done
# 
# # FILES="$(find /path/to/dir/ -type f -name "*.pdf" | sort)"
# # gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=merged_output.pdf $FILES
