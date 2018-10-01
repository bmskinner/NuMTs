#!/bin/bash

for($f in file){
	samtools view -h -F 14 SRR5337620.bam | awk '($3!=$7 && $7!="=")' > SRR5337620.filt.sam
	awk '{$3=="MT" || $7=="MT"}' SRR5337620.filt.sam > SRR5337620.filt.MT.sam
	samtools view -b SRR5337620.filt.MT.sam > SRR5337620.filt.MT.bam
	samtools reheader SRR5337620.filt.sam SRR5337620.filt.MT.bam > SRR5337620.filt.MT.reheader.bam
}