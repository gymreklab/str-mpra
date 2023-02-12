# MPRA design

We designed arrays with 100K sequences for both mm10 and hg38 genomes.

Array design files: https://drive.google.com/drive/folders/1ad2QJ2hmWKugiGwwPil3pdxGk1pr4tIz 
See code: https://github.com/mikmaksi/111821_mouse_human_array_design

v4_121321_30bp: v4 array design with 30bp cutoff for homopolymers 

Workflow:
1. Transcripts from EMBL -> protein coding only -> keep transcript with best support level or longest
2. STRs from HipSTR reference -> remove any with "/"
3. Make list of TSS/STR, where STRs are within 3400 bp of TSS for hg38 and 2700bp of TSS for mm10; large STR 
 allele (ref + 5) has to be smaller than 60bp
5. Get genomic context for every STR such that large allele + context = 168bp
6. Generate small, medium and large alleles and add filler sequence; censor smallest sequence at zero if central allele - 5 is less than zero
7. Add other required elements
8. Remove probes where BsaI recognition site is in the genomic context 
