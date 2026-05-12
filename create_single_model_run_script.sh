#!/bin/bash

combine_r_files() {
    local folder_name="${1:-R/}"
    local output_file="${2:-all_funcs_combined.R}"
    local model_run_file="${3:-runModelRCS.R}"

    > "$output_file"

    for file in "$folder_name"*; do
        if [[ -f "$file" ]]; then
            echo -e "\n\n#### Current file: $file\n" >> "$output_file"
            cat "$file" >> "$output_file"
        fi
    done

    echo -e "\n\n#### Current file: $model_run_file\n" >> "$output_file"
    cat "$model_run_file" >> "$output_file"
}

output_file_name="$1"
model_run_file_name="$2"

if [[ -n "$output_file_name" ]]; then
    if [[ -n "$model_run_file_name" ]]; then
        combine_r_files "R/" "$output_file_name" "$model_run_file_name"
    else
        combine_r_files "R/" "$output_file_name"
    fi
elif [[ -n "$model_run_file_name" ]]; then
    combine_r_files "R/" "" "$model_run_file_name"
else
    combine_r_files
fi