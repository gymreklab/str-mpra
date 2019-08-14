import sys
import random


oligos_file = open(sys.argv[1], 'r')
repeat_file = open('/storage/mlamkin/projects/eSTR-data/repeats-6-12.txt', 'r')
output = open(sys.argv[2], 'a')
full_output = open(sys.argv[3], 'a')
neg_cntrls = []
eSTRs = []
fun = []
repeats = []

# remove header
oligos_file.readline()

# grab all data
for line in oligos_file:
    oligo_data = line.rstrip('\n').split('\t') 
    if oligo_data[0] == 'Negative_Control':
        neg_cntrls.append(oligo_data)
    elif oligo_data[0] == 'eSTR':
        eSTRs.append(oligo_data)
    else:
        fun.append(oligo_data)

# repeats to search for
for line in repeat_file:
    unit = line.rstrip('\n').split('\t')
    repeats.append([unit[0], unit[1].replace('\r', '')])

# Remove repeats
stop_filter = False
for repeat in repeats:
    if stop_filter: break
    to_remove = []
    cur_remove = 0
    for eSTR in eSTRs:
        if not eSTR[5].replace(' ', '').find(repeat[1]) == -1:
            to_remove.append(eSTR)
            cur_remove += 1
            if int(repeat[0]) - cur_remove == 3:
                break

    # remove all eSTRs 
    for repeat_seq in to_remove:
        eSTRs.remove(repeat_seq)
        if len(eSTRs) + len(fun) + len(neg_cntrls) == 7500:
            stop_filter = True
            break

all_oligos = eSTRs + neg_cntrls + fun
random.shuffle(all_oligos)


full_output.write("Oligo_Type\tChrom\tPos\tGene\tRepeat_Number\tFoward_Primer_1 Variant_Sequence Restriction_Enzyme_1 Filler_Sequence Restriction_Enzyme_2 Tag Reverse_Primer_1\n")
for oligo in all_oligos:
    full_output.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(oligo[0], oligo[1], oligo[2], oligo[3], oligo[4], oligo[5]))
    output.write(oligo[5].replace(' ', '')+'\n')
