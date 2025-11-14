#!/bin/bash

# Function to print usage information if the user runs the script incorrectly
usage() {
    echo "Usage: $0 -i <input_file> -s <success_output_file> -f <failure_output_file> -p <max_percentage_threshold> -t <sum_threshold>"
    exit 1
}

# Default values for thresholds
max_percentage_threshold=95  # Default max percentage threshold is 95%
sum_threshold=100            # Default sum threshold is 100

# Parse command-line arguments using getopt
while getopts ":i:s:f:p:t:" opt; do
  case $opt in
    i) input_file="$OPTARG" ;;    # Input file to process
    s) success_file="$OPTARG" ;;  # Output file for entries that pass the threshold
    f) failure_file="$OPTARG" ;;  # Output file for entries that fail the threshold
    p) max_percentage_threshold="$OPTARG" ;;  # Max percentage threshold
    t) sum_threshold="$OPTARG" ;;  # Sum threshold
    *) usage ;;  # Show usage message if the argument is invalid
  esac
done

# Check if all required arguments are provided (input and output files)
if [ -z "$input_file" ] || [ -z "$success_file" ] || [ -z "$failure_file" ]; then
    usage  # Exit and show usage if any file argument is missing
fi

awk -F'\t' -v max_threshold="$max_percentage_threshold" -v sum_threshold="$sum_threshold" '
BEGIN {
    FS="\t"; OFS="\t"
}
{
    # Split the second column into items
    n = split($2, items, ",")
    
    # Initialize variables for total sum, max value, and max item
    total_sum = 0
    max_value = 0
    max_item = ""
    other_items = ""
    
    # Process each item in the second column
    for (i = 1; i <= n; i++) {
        # Extract the value from the item (after the last underscore)
        split(items[i], parts, "_")
        value = parts[length(parts)]
        
        # Update total sum
        total_sum += value
        
        # Track the max value and corresponding item
        if (value > max_value) {
            max_value = value
            max_item = items[i]
        }
    }
    
    # After identifying max, process the items to separate them
    for (i = 1; i <= n; i++) {
        if (items[i] != max_item) {
            if (other_items != "") {
                other_items = other_items "," items[i]
            } else {
                other_items = items[i]
            }
        }
    }

    # Calculate the max percentage
    max_percentage = (max_value / total_sum) * 100

    # Check if both thresholds are met
    if (max_percentage >= max_threshold && total_sum >= sum_threshold) {
        print $1, max_item, other_items > "'$success_file'"
    } else {
        print $1, max_item, other_items > "'$failure_file'"
    }
}
' "$input_file"