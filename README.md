# Alignment of the reads with respect to the reference genome
```
# Indexing in order to speed up the computation of the alignment
./minimap2 -x map-ont -d ecoli_ont.mmi ecoli.fasta
# The alignment itself
./minimap2 -ax map-ont ecoli_ont.mmi ecoli_simulated_reads.fasta > ecoli.sam

# Converting sam file to the bam file in order to use the BamTools library
samtools view -bS ecoli.sam > ecoli.bam

#Sorting the bam file in order for the proper functioning of the manual algorithm and FreeBayes
samtools sort -o ecoli_sorted.bam ecoli.bam

#Indexing for the sorted bam file and reference genome (required for the correct work of FreeBayes(
samtools index ecoli_sorted.bam
samtools faidx ecoli.fasta

#At least 5 reads to cover the position and 50%+ must have the same alternative
./freebayes -f ecoli.fasta --min-coverage 5 --min-alternate-fraction 0.5 ecoli_sorted.bam > ecoli_sorted.vcf

```
**ecoli.fasta** - reference genome
