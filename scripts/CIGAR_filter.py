# cigar_filter.py

from pysam import AlignmentFile
from cigar import Cigar
import argparse

MIN_START_MATCH = 25
MIN_END_MATCH = 20
MAX_INDEL_ALLOWANCE = 2

def cigar_filter(cigar_string, read_length, end_threshold):
    """
    Check the CIGAR string for validity.

    Parameters:
      cigar_string (str): The CIGAR string from column 6 of a SAM file.
      read_length (int): The expected read length.
      end_threshold (int): The threshold for the end condition.

    Returns:
      bool: True if the CIGAR string is valid; False otherwise.
    """
    
    MIN_START_MATCH = 25
    MIN_END_MATCH = 20
    MAX_INDEL_ALLOWANCE = 2
    
    if cigar_string is None or cigar_string == "*":
        return False
    
    parse_cigar = list(Cigar(cigar_string).items())
    # print(len(parse_cigar))
    # print(parse_cigar)
    
    if len(parse_cigar) > 4:
        return False
    
    if len(parse_cigar) == 1:
        if parse_cigar[0][1] == "M" and parse_cigar[0][0] < read_length:
            return True
        else:
            return False
    
    elif len(parse_cigar) == 2:
        if parse_cigar[0][1] == "M" and parse_cigar[1][0] <= end_threshold:
            return True
        elif parse_cigar[1][1] == "M" and parse_cigar[0][0] <= end_threshold:
            return True
        else:
            return False
    
    elif len(parse_cigar) == 3:
        start = parse_cigar[0]
        end = parse_cigar[2]
        mid = parse_cigar[1]

        start_condition = start[1] == "M" and start[0] >= MIN_START_MATCH and start[0] <= read_length
        end_condition = end[1] == "M" and end[0] >= MIN_END_MATCH and end[0] <= read_length
        mid_condition = mid[1] in ["I", "D"] and mid[0] <= MAX_INDEL_ALLOWANCE
        
        if (start_condition and end_condition and mid_condition):
            return True
        elif mid[1] == "M" and start[0] <= end_threshold and end[0] <= end_threshold:
            return True
        else:
            return False
        
    elif len(parse_cigar) == 4:
        start = parse_cigar[0]
        end = parse_cigar[3]
        mid1 = parse_cigar[1]
        mid2 = parse_cigar[2]

        if start[1] != "M" and end[1] == "M" and mid1[1] == "M" and mid2[1] in ["I", "D"]:
            if start[0] < end_threshold:
                return True
            else:
                return False

        elif start[1] == "M" and end[1] != "M" and mid2[1] == "M" and mid1[1] in ["I", "D"]:
            if end[0] < end_threshold:
                return True
            else:
                return False
        else:
            return False
        
    return False 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter BAM file records based on CIGAR strings.")
    parser.add_argument("-i", "--input", help="Input BAM file path", required=True)
    parser.add_argument("-o", "--output", help="Output BAM file path", required=True)
    parser.add_argument("-l", "--readLength", help="Read length", default=500, type=int)
    parser.add_argument("-t", "--lengthThreshold", help="End mismatch threshold", default=10, type=int)
    args = parser.parse_args()

    input_bam_file = args.input
    output_bam_file = args.output
    length = args.readLength
    end_threshold = args.lengthThreshold  

    with AlignmentFile(input_bam_file, "rb") as input_bam, AlignmentFile(output_bam_file, "wb", header=input_bam.header) as output_bam:
        for alignment in input_bam:
            if alignment.cigarstring is None:
                break
            if cigar_filter(alignment.cigarstring, length, end_threshold):
                output_bam.write(alignment)







