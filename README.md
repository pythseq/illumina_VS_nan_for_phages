# ILMN-vs-PACB-Helgoland

This repository collects all code and extended methodology for the manuscript "" currently under revision.

Raw reads, assemblies and MAGs are available in 
All MAGs are locally available [here](/MAGs-ILMN-PACB.tar.gz). 
Table of contents<br>
[1. Raw read analyses](#1-raw-read-analyses)<br>
[2. Assemblies](#2-assemblies-and-recovery-of-metagenome-assembled-genomes)<br>
[3. MAG detection comparison between Illumina and PacBio](#3-mag-detection-and-comparison-between-illumina-and-pacbio)<br>

## 1 Raw read analyses

The following table collects the statistics for unassembled reads.

| Library | Sampling date| Technology |# raw reads | # trimmed reads | avg. length |
|---|:-:|:-:|:-:|:-:|:-:|
| 20200310  | 10.3.2020  | Illumina HiSeq2500  | 158,123,188 |135,156,020 |  240.4 |
| 20200414  | 14.4.2020  | Illumina HiSeq2500  | 174,007,074 |151,208,718 | 242.0  |
|  20200430 |  30.4.2020 | Illumina HiSeq2500  | 177,605,098 |154,874,006 | 241.0  |
|  20200506 | 6.5.2020  |  Illumina HiSeq2500 | 179,184,628 |157,094,252 | 241.1 |
| B  | 10.3.2020  | PacBio Sequel 2  | -  | 1,506,023  | 8951.2|
|  I |  14.4.2020 |  PacBio Sequel 2  | -  |  1,932,822 |6379.1 |
|  Q |  30.4.2020 |  PacBio Sequel 2  | -  | 1,770,521  |3748 |
|  U |  6.5.2020  |  PacBio Sequel 2  | -  | 3,040,187  |5405.9 |


## 2 Assemblies and recovery of metagenome-assembled genomes

For Illumina reads, the assemblies were done using metaSPADES 3.14.1 (options: meta, k-mers 21, 33, 55, 77, 99, 127, and read-error correction enabled), read mapping using BBmap v37.02 (options: minid=0.99, idfilter=0.97), and binning using a combination of CONCOCT v.1.1.0, metabat2 v2.12.1, binsanity v0.4.4, and maxbin2 v2.2.7 and DASTool v1.1.2 as [previously described](https://doi.org/10.1016/j.syapm.2018.08.007). 

For PacBio, assemblies were obtained using Flye. v2.8-b1674 (options: --meta, --pacbio-hifi). Read mapping was performed using pbmm2 v1.4.0 (options: -HIFI, -x 97) and binning was done using metabat2 v2.12.1.

MAGs from both technologies with quality values above 50 (quality = Completion -5xContamination >=50) later refined in anvi'o v6.2.

### Mapping of Illumina short-reads to unassembled PacBio reads

To determine the degree of overlapping between the sequences from PACB and ILMN, I mapped the short-reads against the long-reads and then counted the number of reads mapped (summary of the SLURM script below):

```
#SBATCH --array=310,414,430,506

bowtie2  --reorder --no-unal -p 42 -x 20200$SLURM_ARRAY_TASK_ID.bwt2.db -1 20200$SLURM_ARRAY_TASK_ID.clean.R1.fastq.gz -2 20200$SLURM_ARRAY_TASK_ID.clean.R2.fastq.gz  > 20200$SLURM_ARRAY_TASK_ID-ILMNreads-against-PACBreads-as-DB.sam
samtools view -@ 10 -bS 20200$SLURM_ARRAY_TASK_ID-ILMNreads-against-PACBreads-as-DB.sam -o 20200$SLURM_ARRAY_TASK_ID-ILMNreads-against-PACBreads-as-DB.bam;
samtools sort -@ 10 20200$SLURM_ARRAY_TASK_ID-ILMNreads-against-PACBreads-as-DB.bam -o 20200$SLURM_ARRAY_TASK_ID-ILMNreads-against-PACBreads-as-DB.sorted.bam;

```
Then to count the short-reads mapped:

```
for i in 310 414 430 506; do samtools view -@ 20  -c -F 260 20200${i}-ILMNreads-against-PACBreads-as-DB.sorted.bam; done
```
The -F 260 is to exclude unmmapped and not primarily mapped reads (Flags 4 + 256 =260)

To check how many PacBio reads were not covered, I used:

```
for i in 310 414 430 506; do samtools coverage 20200${i}-ILMNreads-against-PACBreads-as-DB.sorted.bam |awk -F'\t' '{if($4==0) print $0}'|wc -l ; done

```
### Maping of reads to contigs derived from short- and long-read metagenomes

We first mapped the SRs to the contigs >=500 bp and generated sorted bam files that were used many tests:
```
#SBATCH --array=310,414,430,506
bowtie2-build --threads 24 20200${SLURM_ARRAY_TASK_ID}_Illumina-500bp.fasta  20200$SLURM_ARRAY_TASK_ID.bwt2.db;
bowtie2 --reorder --no-unal -p 24 -x 20200$SLURM_ARRAY_TASK_ID.bwt2.db -1 20200$SLURM_ARRAY_TASK_ID.clean.R1.fastq.gz -2 20200$SLURM_ARRAY_TASK_ID.clean.R2.fastq.gz  > 20200$SLURM_ARRAY_TASK_ID-ILMN.contigs.sam
samtools view -bS -F 260  20200$SLURM_ARRAY_TASK_ID-ILMN.contigs.sam -o 20200$SLURM_ARRAY_TASK_ID-ILMN.contigs.bam;
samtools sort 20200$SLURM_ARRAY_TASK_ID-ILMN.contigs.bam -o 20200$SLURM_ARRAY_TASK_ID-ILMN.contigs.sorted.bam
```
We also used them to count the number of reads mapped back to contigs:

```
for i in 310 414 430 506; do 
	samtools view -@ 20  -c -F 260 20200${i}-ILMN.contigs.sorted.bam; 
done
```

For LR assemblies, we first filtered all contigs >= 500 bp and then used them as a reference for mapping using `pbmm2` as described above without using a identity cutoff:

```
for MG in 310 414 430 506; do 
	pbmm2 align -N 1 20200${MG}_MetaFlye-500bp.fasta  20200${MG}.PacBio.fasta 20200${MG}.PacBio.sorted.noid.bam --sort --preset CCS  -j 50
	samtools view -@ 20  -c -F 260 20200${MG}.PacBio.sorted.noid.bam; 
done
```

### Mapping of reads to all MAG contigs

In this case I used the concatenad files from Illumina and PacBio MAGs. I used these files then to map all reads without using an identity threshold.

For PacBio:

```
for MG in 310 414 430 506; do
	pbmm2 align -N 1 all-PACB-20200${MG}.fa 20200${MG}.PacBio.fasta PACB-20200${MG}.sorted.bam --sort --preset HIFI 
	samtools view -c -F 260 PACB-20200${i}.sorted.bam
done
```
For Illumina:
```
for MG in 310 414 430 506;
	bowtie2 --reorder --no-unal -p 24 -x 20200${MG}.bwt2.db -1 20200${MG}.clean.R1.fastq.gz -2 20200${MG}.clean.R2.fastq.gz  > 20200${MG}-ILMN.contigs.sam;
	samtools view -c -F 260 20200${MG}-ILMN.contigs.sam
done
```

### Abundance of MAGs

#### Filtering mappings using an identity threshold

The general strategy was to first map either SRs or LRs to all MAGs (i.e., a concatenated file), then separate the resulting SAM into individual files (one per MAG). Using these individual SAM files, then I calaculated abundance metrics (e.g., sequencing depth, breadth, etc).

For SRs, the mapping of reads were done just as the section above but `sam` files were further filtered using a 97% identity threhold with the `sam.filter.rb` scripts from the [enveomics collection](https://github.com/lmrodriguezr/enveomics). Here I show the 2020.03.10 sample as example:
```
#SBATCH --array=1-82
DATE=20200310;
TECH=ILMN;
echo MAG $SLURM_ARRAY_TASK_ID;	
sam.filter.rb -i 97 -m ../all-${TECH}-${DATE}.sam -o temp.${TECH}-${DATE}-m${SLURM_ARRAY_TASK_ID}.sam -g ../${TECH}-${DATE}-m${SLURM_ARRAY_TASK_ID}.fa	
grep -F m${SLURM_ARRAY_TASK_ID}_ temp.${TECH}-${DATE}-m${SLURM_ARRAY_TASK_ID}.sam > ${TECH}-${DATE}-m${SLURM_ARRAY_TASK_ID}.sam;
rm  temp.${TECH}-${DATE}-m${SLURM_ARRAY_TASK_ID}.sam;
samtools view -F 4 -bS ${TECH}-${DATE}-m${SLURM_ARRAY_TASK_ID}.sam -o ${TECH}-${DATE}-m${SLURM_ARRAY_TASK_ID}.bam;
samtools sort ${TECH}-${DATE}-m${SLURM_ARRAY_TASK_ID}.bam -o ${TECH}-${DATE}-m${SLURM_ARRAY_TASK_ID}.sorted.bam
bedtools genomecov -ibam ${TECH}-${DATE}-m${SLURM_ARRAY_TASK_ID}.sorted.bam -bga > ${TECH}-${DATE}-m${SLURM_ARRAY_TASK_ID}.sorted.bam.bg;	
echo -e ${TECH}-${DATE}-m${SLURM_ARRAY_TASK_ID}"\t"$(BedGraph.tad.rb -i ${TECH}-${DATE}-m${SLURM_ARRAY_TASK_ID}.sorted.bam.bg -r 0.8) >>${DATE}-${TECH}-tad80.tsv;
done

```
For LRs, a similar strategy was followed but here we don't truncate the sequencing depth (i.e., -r 1 in the `BedGraph.tad.rb` script) since the values are much lower compared to those we can get from SRs.

```
DATE=20200310;
TECH=PACB;
for MAG in {1..83}; do
	echo $MAG;	
	sam.filter.rb -t 5 -i 97 -m ${TECH}-${DATE}.MD.sam -o temp.${TECH}-${DATE}-m${MAG}.sam -g ${TECH}-${DATE}-m${MAG}.fa	
	grep -F m${MAG}_ temp.${TECH}-${DATE}-m${MAG}.sam > ${TECH}-${DATE}-m${MAG}.sam;
	rm  temp.${TECH}-${DATE}-m${MAG}.sam;
	samtools view -F 4 -bS ${TECH}-${DATE}-m${MAG}.sam -o ${TECH}-${DATE}-m${MAG}.bam;
	samtools sort ${TECH}-${DATE}-m${MAG}.bam -o ${TECH}-${DATE}-m${MAG}.sorted.bam
	bedtools genomecov -ibam ${TECH}-${DATE}-m${MAG}.sorted.bam -bga > ${TECH}-${DATE}-m${MAG}.sorted.bam.bg;	
	echo -e ${TECH}-${DATE}-m${MAG}"\t"$(BedGraph.tad.rb -i ${TECH}-${DATE}-m${MAG}.sorted.bam.bg -r 1) >>${DATE}-${TECH}-tad100.tsv;
done


```


## 3 MAG detection and comparison between Illumina and Pacbio
### 3.1 MAG taxonomy
First, we get the taxonomy for all MAGs.

All MAGs:
```
gtdbtk classify_wf -x fa --cpus 48 --genome_dir 02.gtdb-tk --out_dir 02.gtdb-tk/gtdb-tk-595-MAGs/
```
For Illumina MAGs:
```
gtdbtk classify_wf -x fa --cpus 20 --genome_dir 01.db/01.Illumina/00.all/ --out_dir 02.gtdb-tk/gtdb-tk-ILMN-341-MAGs/
```
For PacBio MAGs:

```
gtdbtk classify_wf -x fa --cpus 20 --genome_dir 01.db/02.PacBio/00.all/ --out_dir 02.gtdb-tk/gtdb-tk-PACB-254-MAGs/
```

### 3.2 ANIs

I used fastANI (v1.3.2) to determine ANIs between technologies and sampling dates:

```
fastANI --ql ILMN.list --rl PACB.list -t 32 -o 20200310.fastANI133.tsv
```
Then I got the pairs that matched between technologies as follows:
```
cat 20200310.fastANI133.tsv |awk '{if($3>=99) print $1"\t"$2}'|sed 's/.fa//g'|tr "\t" "-"|awk -F'-' '{print $1"-"$2"-"$3"\t"$4"-"$5"-"$6"\t"$1"-"$4"-"$2"-"$3"-"$6}' >ILMN-PACB-renamed.20200310.tsv
```

Not found in previous

```
comm -3  <(sort ILMN.20200310.list) <(sort ILMN-and-PACB.20200310.list)

```

