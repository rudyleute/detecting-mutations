# Alignment of the reads with respect to the reference genome
```
# Indexing in order to speed up the computation of the alignment
./minimap2 -x map-ont -d ecoli_ont.mmi ecoli.fasta
# The alignment itself
./minimap2 -ax map-ont ecoli_ont.mmi ecoli_simulated_reads.fasta > ecoli.sam
```
**ecoli.fasta** - reference genome
