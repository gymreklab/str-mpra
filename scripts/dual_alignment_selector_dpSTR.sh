#!/bin/bash

# Function to print usage information
usage() {
    echo "Usage: $0 -m <mode> -l <length> <input_file> <file1> <file2>"
    echo "  -m <mode>       Mode of operation: 'indel' or 'match'"
    echo "  -l <length>     Length value to be used in the script"
    echo "  <input_file>    Path to the input file"
    echo "  <file1>         Path to the output file for matching records"
    echo "  <file2>         Path to the output file for non-matching records"
    exit 1
}

# Parse command-line arguments
while getopts "m:l:" opt; do
    case $opt in
        m) mode=$OPTARG ;;
        l) length=$OPTARG ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))

# Check if the correct number of arguments are provided
if [ "$#" -ne 3 ]; then
    usage
fi

# Assign arguments to variables
input_file="$1"
file1="$2"
file2="$3"

# Validate the mode
if [[ "$mode" != "indel" && "$mode" != "match" ]]; then
    echo "Invalid mode specified. Use 'indel' or 'match'."
    exit 1
fi

# Validate length
if ! [[ "$length" =~ ^[0-9]+$ ]]; then
    echo "Length must be a non-negative integer."
    exit 1
fi

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "Input file $input_file does not exist."
    exit 1
fi

# Process the input file based on the mode
if [ "$mode" == "match" ]; then
    awk -F "," -v OFS="\t" -v file1="$file1" -v file2="$file2" -v m="${length}M" '
    {
        len3 = length($3)
        len9 = length($9)
        
        tag5 = (match($5, /^NM/) ? substr($5, 1, 4) : "")
        s5 = (match($5, /^NM/) ? substr($5, 6, 6) : "")
        tag11 = (match($11, /^NM/) ? substr($11, 1, 4) : "")
        s11 = (match($11, /^NM/) ? substr($11, 6, 6) : "")

        substr3 = substr($3, 1, len3 - 2)
        substr9 = substr($9, 1, len9 - 2)
        
        # The order of the input file
        # QNAME, FLAG, RNAME, POS, CIGAR, TAGs(last column for NM:i:n which denote the editing distance)

        if (substr3 == substr9) { # exclude the NA case
            if (($2 == $8 && $3 == $9 && $4 == $10 && $5 == $11 && $4 == m && $6 == 1 && $12 == 1) || 
                ($4 == m && $2 == 0 && $8 == 16 && tag5 == "NM:i" && $6 == 1 && tag11 == "NM:i" && s5 <= 0)) { 
                # when the first alignment is better than the second on CIGAR, FLAG, Editing distance, position
                print $1, $2, $3, $4, $5, $6 > file1
            } else if ($10 == m && $8 == 0 && $2 == 16 && tag11 == "NM:i" && $12 == 1 && tag5 == "NM:i" && s11 <= 0) { # change s11 <= 0 to s11 <= s5 if want to relax
                # when the second alignment is better than the second on CIGAR, FLAG, Editing distance, position
                print $7, $8, $9, $10, $11, $12 > file1
            } 
        } else if (substr3 != substr9) {
            if ($4 == m && $2 == 0 && $6 == 1) {
                if ($10 != m) { # NA case is saved here
                    print $1, $2, $3, $4, $5, $6 > file2
                } else if (tag11 == "NM:i" && s5 < s11 && s5 == 0) { 
                    print $1, $2, $3, $4, $5, $6 > file2
                }
            } else if ($10 == m && $8 == 0 && $12 == 1) {
                if ($4 != m) { # NA case is saved here
                    print $7, $8, $9, $10, $11, $12 > file2
                } else if (tag5 == "NM:i" && s11 < s5 && s11 == 0) { 
                    print $7, $8, $9, $10, $11, $12 > file2
                }
            }
        }
    }' "$input_file"

elif [ "$mode" == "indel" ]; then
    awk -F "," -v OFS="\t" -v file1="$file1" -v file2="$file2" -v m="${length}M" '
    {
        len3 = length($3)
        len9 = length($9)

        tag5 = (match($5, /^NM/) ? substr($5, 1, 4) : "")
        s5 = (match($5, /^NM/) ? substr($5, 6, 6) : "")
        tag11 = (match($11, /^NM/) ? substr($11, 1, 4) : "")
        s11 = (match($11, /^NM/) ? substr($11, 6, 6) : "")

        substr3 = substr($3, 1, len3 - 2)
        substr9 = substr($9, 1, len9 - 2)

        if (substr3 == substr9) { # exclude the NA case
            if ($4 != m && $2 == 0 && tag5 == "NM:i" && $6 == 1) { # when the first entry editing distance is better, started at 1 and has paired flag with second
                if ($8 == 16 && tag11 == "NM:i" && s5 <= s11) {
                    print $1, $2, $3, $4, $5, $6 > file1
                }
            } else if ($10 != m && $8 == 0 && tag11 == "NM:i" && $12 == 1) { # when the second entry editing distance is better, started at 1 and has paired flag with second
                if ($2 == 16 && tag5 == "NM:i" && s11 <= s5) {
                    print $7, $8, $9, $10, $11, $12 > file1
                }
            }
        } else if (substr3 != substr9) {
            if ($4 != m && $2 == 0 && $6 == 1 && tag5 == "NM:i" && $11 == "NA") {
                print $1, $2, $3, $4, $5, $6 > file2
            } else if ($10 != m && $8 == 0 && $12 == 1 && tag11 == "NM:i" && $5 == "NA") {
                print $7, $8, $9, $10, $11, $12 > file2
            }
        }
    }' "$input_file"
fi

echo "Processing completed."
