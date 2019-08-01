import sys
from oligos import hamming_distance


KpnI = "GGTACC"
XbaI = "TCTAGA"

# filter once all oligos have been generated
def tag_generator(tag_file, num_tags):
    tags = []
    filter_suffix = ['TCT', 'TCA', 'TCC','TCG', 'TAT', 'TGT', 'TTT', 'ACT', 'CCT', 'GCT']
    for line in tag_file:
        tags.append(line.rstrip('\n'))
   
    oligo_tags = [tags[0]]
    
    # iterate through all current tags and grab a list of num_tags
    for tag in tags[1:]:
        if len(oligo_tags) == num_tags: break
        grab_tag = True
        for check_tag in oligo_tags:
            if not grab_tag: break
            # ensure tag is distinct from other tags
            if hamming_distance(tag, check_tag) < 2:
                grab_tag = False
            
            # ensure tag does not contain restriction enzymes
            for motif in [XbaI, KpnI]:
                if not grab_tag: break
                # check tag itself
                for i in range(len(tag)-len(motif)+1):
                    if hamming_distance(motif, tag[i:i+len(motif)]) < 2:
                        grab_tag = False
                        break
                # check context of tag
                for prefix in [XbaI[:5], KpnI[:5]]:
                    if hamming_distance(tag[-5:], prefix) < 2:
                        grab_tag = False

        if tag[-3:] in filter_suffix or tag[:4] == 'TACC':
            grab_tag = False
        if grab_tag:
            oligo_tags.append(tag)
 
    return oligo_tags


tag_file = open(sys.argv[1], 'r')
output_file = open(sys.argv[2], 'a')
oligo_tags = tag_generator(tag_file, 8000)
for tag in oligo_tags:
    output_file.write("%s\n"%tag)

