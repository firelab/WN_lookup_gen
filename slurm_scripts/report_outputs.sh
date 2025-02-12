#!/bin/bash

base_dir="/mnt/fsim/windninja/CONUS"
required_dirs=("mass" "momentum")
direction_folders=("0-0-deg" "22-5-deg" "45-0-deg" "67-5-deg" "90-0-deg" "112-5-deg" "135-0-deg" "157-5-deg" "180-0-deg" "202-5-deg" "225-0-deg" "247-5-deg" "270-0-deg" "292-5-deg" "315-0-deg" "337-5-deg")
required_files=("*_ang.asc" "*_ang.prj" "*_cld.asc" "*_cld.prj" "*_vel.asc" "*_vel.prj")  # Excluding cli.cfg

valid_folders=()
unfinished_folders=()
issue_folders=()

declare -A missing_files_count_per_folder
declare -A convergence_count_per_folder
declare -A convergence_wind_dirs
declare -A wind_dir_failure_count  # Map to count failures per wind direction

total_failed_simulations=0  # Sum of all failed simulations

for folder in "$base_dir"/*; do
    folder_name=$(basename "$folder")
    dem0_path="$folder/dems_folder/dem0"

    # Initialize counters
    missing_count=0
    convergence_count=0
    converge_dirs=()

    # Check if dem0 exists
    if [ ! -d "$dem0_path" ]; then
        unfinished_folders+=("$folder_name")
        continue
    fi

    # Check mass & momentum exist
    if [ ! -d "$dem0_path/mass" ] || [ ! -d "$dem0_path/momentum" ]; then
        unfinished_folders+=("$folder_name")
        continue
    fi

    # Check for missing files
    for req_dir in "${required_dirs[@]}"; do
        for dir_folder in "${direction_folders[@]}"; do
            dir_path="$dem0_path/$req_dir/$dir_folder"
            if [ ! -d "$dir_path" ]; then
                ((missing_count += 6))  # Each missing direction contributes 6 missing files
            else
                for file in "${required_files[@]}"; do
                    if ! ls "$dir_path"/$file &>/dev/null; then
                        ((missing_count++))  # Increment integer safely
                    fi
                done
            fi
        done
    done

    # Count convergence errors from log and extract wind directions
    log_file="$folder/simulation.log"
    if [ -f "$log_file" ]; then
        convergence_count=$(grep -c "The flow solution did not converge" "$log_file")

        # Extract wind directions just before an exception
        prev_wind_dir=""
        while IFS= read -r line; do
            if [[ "$line" =~ running\ \"momentum\"\ simulation\ for\ dem\ \"dem0\"\ windDir\ \"([0-9]+\.[0-9]+)\ deg\" ]]; then
                prev_wind_dir="${BASH_REMATCH[1]}"
            elif [[ "$line" =~ Exception\ caught:\ The\ flow\ solution\ did\ not\ converge ]]; then
                converge_dirs+=("$prev_wind_dir")

                # Increment count for this wind direction in the failure map
                ((wind_dir_failure_count["$prev_wind_dir"]++))
            fi
        done < "$log_file"
    fi

    # Ensure variables remain integers
    missing_count=${missing_count:-0}
    convergence_count=${convergence_count:-0}

    # Remove duplicate wind directions (halve them)
    unique_converge_dirs=($(echo "${converge_dirs[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

    # Sum up total failed simulations
    total_failed_simulations=$((total_failed_simulations + convergence_count / 2))

    # Classify folder
    if [[ "$missing_count" -eq 0 && "$convergence_count" -eq 0 ]]; then
        valid_folders+=("$folder_name")
    else
        issue_folders+=("$folder_name")
        missing_files_count_per_folder["$folder_name"]=$missing_count
        convergence_count_per_folder["$folder_name"]=$((convergence_count / 2))  # Halve the count
        convergence_wind_dirs["$folder_name"]="${unique_converge_dirs[*]}"  # Unique wind directions
    fi
done

# Summary Report
echo "--------------------------------"
echo "Valid Folders: ${#valid_folders[@]}"
echo "Unfinished Folders: ${#unfinished_folders[@]}"
echo "Issue Detected Folders: ${#issue_folders[@]}"
echo "Total Failed Simulations: $total_failed_simulations"
echo "--------------------------------"

# Detailed Report for Issue Folders
for folder in "${issue_folders[@]}"; do
    echo "$folder:"
    echo "   Missing Files: ${missing_files_count_per_folder[$folder]}"
    echo "   Converge Error Count: ${convergence_count_per_folder[$folder]}"
    echo "   Converge Error Wind Directions: ${convergence_wind_dirs[$folder]}"
    echo ""
done

# Mapping of issue folders and failed wind directions
echo "Failed Simulations Map:"
echo "{"
for folder in "${issue_folders[@]}"; do
    echo "  \"$folder\": [ ${convergence_wind_dirs[$folder]} ],"
done
echo "}"

# Space-separated list of folders with failed simulations
echo ""
echo "Folders with failed simulations (space-separated):"
echo "${issue_folders[*]}"

# WindNinja Simulation Convergence Error Analysis (Count of Each Wind Direction Failure)
echo "WindNinja Simulation Convergence Error Analysis:"
echo "{"
for dir in "${!wind_dir_failure_count[@]}"; do
    echo "  \"$dir\": ${wind_dir_failure_count[$dir]},"
done
echo "}"
