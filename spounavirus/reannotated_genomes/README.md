# GFF Files

The gff files are the reannotated files provided by Evelien. This is how I made them usable

Convert all the GFF files to Genbank files:

```
for i in *.gff; 
	do o=$(echo $i | sed -e 's/gff/gbk/'); 
	seqret -sequence $i -outseq $o -osformat2 genbank -sformat1 gff -feature 1;
done
```



