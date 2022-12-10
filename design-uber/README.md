# Generate uber array

## STR set

The set of STRs to consider is in: `uber_str_set.bed`. It contains top STRs from round 1:

* The top 200 STRs (ranked by p-value) with positive correlation between rpt. num and expr
* The top 100 with negative correlation

We are aiming for a total of 60,000 probe sequences.

## Array design

Each oligo consists of: 
FIVE_PRIME_ADAPT, vreg, GIBSON_ASISI, filler_seq, BSAI_RECOG, GIBSON_BSAI_CUT

Where vreg consists of the STR + genomic flanking regions on either side.

For each input STR we consider:
* Multiple different repeat lengths. The exact number is determined by the reference repeat unit length. 
* Multiple different alternate motifs and their reverse complements, including random and random with matched GC%
* For each full variable region, we consider both the forward and reverse complement

Multiple checks are implemented:
* Remove regions containing cloning cut sites
* Remove redundant probes
* Remove probes with a repeat region longer than a pre-set maximum

## Extract the old probes for these sequences

```
cat uber_str_set.bed | cut -f 5 | sort > uber_str_ids.txt
cat array_probes.tsv| grep hg38 | awk '{print $2 "\t" $2"_"$3"_"$5 "\t" $6}' | \
	grep -w -f uber_str_ids.txt | cut -f 1,3 > round1_probe_seqs.tab
```

## Generating the probes

To run:

```
./uber-array-design.py \
	--reffa GRCh38_full_analysis_set_plus_decoy_hla.fa \
	--rptsbed uber_str_set.bed \
	--seed 12345 \
	--add-probes round1_probe_seqs.tab \
	--out uber-v2
```
