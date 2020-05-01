#!/bin/bash

GTF=Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf

## ---- prepare DEXSeq ---- 
## flatten gennome annotation file
GFF=simulation/dexseq/GRCh37.71.dexseq.gff
python2 dexseq/dexseq_prepare_annotation.py -r no $GTF $GFF

## count from bam
for n in 1 2 3 4 5 6
do
	python dexseq/dexseq_count.py -p yes -r pos -s no -f bam $GFF \
	simulation/bam/Hs_sample${n}.bam \
	simulation/dexseq/Hs_sample${n}.txt
done

## ---- rMATS ---- 

DIR=simulation/rmats
rm -f $DIR/rmats.log
rmats.py --b1 $DIR/bam_I.txt --b2 $DIR/bam_II.txt --gtf $GTF --od $DIR -t paired --readLength 101 --cstat 0.0001 --libType fr-unstranded  --nthread 6 &> $DIR/rmats.log

## ---- MAJIQ ---- 

GFF3=Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gff3
DIR=simulation/majiq

## prepare annotation file (GFF3)
perl $DIR/gtf2gff3.pl $GTF > $GFF3

## sort & index & visualize bam files
for i in {1..6}
do
  samtools sort -@ 6 simulation/bam/Hs_sample${i}.bam -o simulation/bam/Hs_sample${i}.sorted.bam
  samtools index -@ 6 simulation/bam/Hs_sample${i}.sorted.bam
  bamCoverage -b simulation/bam/Hs_sample${i}.sorted.bam -o simulation/bam/Hs_sample${i}.bw -bs 1 -p 16 --effectiveGenomeSize 2776919808 --normalizeUsing RPKM
done

majiq build $GFF3 -c $DIR/config -j 6 -o $DIR/build

majiq deltapsi -grp1 $DIR/build/Hs_sample1.sorted.majiq $DIR/build/Hs_sample2.sorted.majiq $DIR/build/Hs_sample3.sorted.majiq -grp2 $DIR/build/Hs_sample4.sorted.majiq $DIR/build/Hs_sample5.sorted.majiq $DIR/build/Hs_sample6.sorted.majiq -j 6 -o $DIR/dpsi -n I II

voila tsv $DIR/build/splicegraph.sql $DIR/dpsi/I_II.deltapsi.voila --threshold 0.1 -f $DIR/I_II.deltapsi.1.tsv 
voila tsv $DIR/build/splicegraph.sql $DIR/dpsi/I_II.deltapsi.voila -f $DIR/I_II.deltapsi.2.tsv 
voila tsv $DIR/build/splicegraph.sql $DIR/dpsi/I_II.deltapsi.voila --show-all -f $DIR/I_II.deltapsi.all.tsv 

