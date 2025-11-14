#!/bin/bash

# Function to display usage instructions
usage() {
    echo "Usage: $0 --length <match_length> --input <input.bam> --output <output.csv>"
    exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --length)
            length="$2"
            shift 2
            ;;
        --input)
            input="$2"
            shift 2
            ;;
        --output)
            output="$2"
            shift 2
            ;;
        *)
            echo "Error: Invalid argument '$1'"
            usage
            ;;
    esac
done

# Validate required arguments
if [[ -z "$length" || -z "$input" || -z "$output" ]]; then
    echo "Error: Missing required arguments."
    usage
fi

# Check if length is a positive integer
if ! [[ "$length" =~ ^[0-9]+$ ]]; then
    echo "Error: Length must be a positive integer."
    exit 1
fi

# Check if input file exists
if [ ! -f "$input" ]; then
    echo "Error: Input file '$input' does not exist."
    exit 1
fi

# Use samtools and awk to process the input file
samtools view "$input" | awk -v FS="\t" -v OFS="," -v m="${length}M" '
{
    # Check if the last field matches the pattern
    if ($NF ~ /^XA:Z:/) {
        split($NF, arr, ";")
        
        # Process each item in the array
        for (i = 1; i <= length(arr); i++) {
            if (arr[i] == "") continue
            
            split(arr[i], fields, ",")
            
            # Extract the relevant fields
            s1 = substr(fields[1], 6)
            s2 = fields[2]
            s3 = fields[3]
            s4 = fields[4]
            
            if (s3 == m && s4 <= substr($12, 6, length($12)) && s2 == "+1" && $4 != 1) {
                $3 = s1
                $12 = "NM:i:" s4
                $2 = 0
                $4 = 1
            }
        }
        print $1, $2, $3, $6, $12, $4
    }
    else {
        print $1, $2, $3, $6, $12, $4
    }
}' |sort -k1,1 > "$output"
