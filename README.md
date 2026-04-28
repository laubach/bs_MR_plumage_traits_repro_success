# bs_MR_plumage_traits_repro_success

# Assessing the causal effect of plumage traits on reproduction using Mendelian randomization

## Genomics pipeline

Summary...
Lists and reference files can be found in various analysis directories, as can Shell and Python scripts used in elements of the analysis workflow. R scripts are in the R directory. Note that you may need to adjust the organization of path/file locations to suit your environment. This workflow assumes that you return to the main working directory ~/hirundo_speciation_genomics after completing major analysis sections. Adjust the path before the main working directory as necessary. Walkthrough of software installations will take place in ~/hirundo_speciation_genomics/tmp/.

## Contents

1. [Software and dependencies](#software-and-dependencies)
2. [Setup environment](#setup-environment)
3. [Genome data processing](#genome-data-processing)
4. [Read mapping](#read-mapping)
5. [Variant calling](#variant-calling)
6. [Filtering for relatedness](#filtering-for-relatedness)
7. [Calculate relatedness](#calculate-relatedness)
8. [Filtering for imputation](#filtering-for-imputation)
9. [Process sex chromosomes for imputation](#process-sex-chromosomes-for-imputation)
10. [Repeat masking for imputation](#repeat-masking-for-imputation)
11. [Additional processing and filtering for imputation](#additional-processing- and-filtering-for-imputation)
12. [Impute and call genotypes](#impute-and-call-genotypes)



## 1. Software and dependencies

The software and their dependencies used in this genomics pipleine are listed here and assumed to be on the user path.

wget -O picard.jar "https://github.com/broadinstitute/picard/releases/download/2.8.1/picard.jar"

check it is working
java -jar picard.jar -h

## 2. Setup environment

Create subdirectories for storing data, logs, and list files. Raw read data are available from [NCBI SRA](https://www.ncbi...), and should be placed in the `fastq` subdirectory.

```
mkdir fastq
mkdir fastq_filtered
mkdir bam
mkdir gvcf
mkdir vcf
mkdir vcf_imputed
mkdir genome_lists
mkdir phenotype_lists
mkdir input
mkdir log
mkdir relate
```

Note: The intermediate data files generated and used in this workflow are stored in local multiple directories outside the main working directory `hirundo_plumage_MR` becuase they are used in multiple projects. Reproducing this work flow from a single parent directory is described in this `README.md`, and it is noted where mapping to other subdirectories is done. To reproduce these analyses, first, create the above subdirectories in the working directory. Second, copy scripts from the `scripts` directory to the working directory. 



## 3. Genome data processing

### 3.1 Filter raw read data with Trimmomatic

Use `Trimmomatic` to filter and quality trim raw read data, using the settings:
* Remove 5' end bases if quality is below 20
* Remove 3' end bases if quality is below 20
* Minimum read length = 32
* Remove reads if average quality is < 30


### 3.2 generate sample list

```
ls -lh ./fastq/ | tail -n+2 | cut -d" " -f9 > ./genome_lists/sample.ID.list
```


### 3.3 Concatonate files and run Trimmomatic on samples in `./genom_lists/sample.ID.list`

```
nohup sh cat_trimmomatic.sh ./genom_lists/sample.ID.list > ./log/runTrimmomatic.log &
```


### 3.4 Remove unpaired filtered reads
```
cd fastq_filtered
rm ./*U.trim.fq.gz
cd ..
```

Analysis ready paired read data are in `fastq_filtered` subdirectory.



## 4. Map reads using BWA

### 4.1 Index the reference genome

```
bwa index Hirundo_rustica_bHirRus1.final.fasta
```


### 4.2  Map reads for each sample to the indexed reference genome using BWA `mem`.

Note: Need to provide this shell script with a sample list in `./genom_lists/sample.ID.list`. The script `runBWAmem.sh` runs BWA `mem` and also indexes using samtools.

```
nohup sh runBWAmem.sh ./genom_lists/sample.ID.list > ./log/runBWAmem.log &
```

Analysis-ready bam files are in `./bam`.



## 5. Variant calling

### 5.1 Individual variant calling
Call individual raw variants using GATK `HaplotypeCaller`.

Run runGATKhaplotypecaller.sh to run `HaplotypeCaller` on samples in `./genom_lists/sample.ID.list`.

```
sh runGATKhaplotypecaller.sh ./genom_lists/sample.ID.list > ./log/runGATKhaplotypecaller.log
```

The output genomic VCFs for each individual are in `./gvcf/`.

Note: This takes a long time to run on all samples in sequence (15 hours per sample). We can split out the job to parallel processes by breaking up the list of samples and running the script on each list separately (constrained by CPU threads).


### 5.2 Index all input gvcf files

* Make a path list for all inpubut gvcf files. Format `./genom_lists/sample.ID.list`, by adding the path names to all input gvcf files using the `write_path.sh` shell script 

```
sh write_path.sh ./genom_lists/sample.ID.list > ./genom_lists/sample.ID.list.gvcf.paths
```

* Index all input gvcf files. Run `runTabix.sh` to index all input gvcf files providing the `./genom_lists/sample.ID.list.gvcf.paths` list.

```
sh tabix.sh ./genom_lists/sample.ID.list.gvcf.paths > ./log/tabix.log
```

Note: If making a list of "old" vcf.gz files to combine with the samples in `./genom_lists/sample.ID.list`. Use grep to find "CO_[wildcard].vcf.gz". Because there are vcf.gz, vcf.gz.tbi, and vcf.idx files in a subdirectory called `gvcf_older`, and because we need to remove the duplicate matches from the files, then use sort -u and only keep unique lines. Concatenate with the new list before running the above scripts in section 5.2.
```
ls -lh gvcf_older/ | grep -ho '\S*CO_\S*.vcf.gz' | sort -u > sample.ID.list.old
```


### 5.3 Cohort variant calling

* Format intervals files with lists of scaffolds to run parallel analyses on in ./genome_lists/scaffold.list1.intervals, etc. (13 files total). These are split out based on the largest scaffold, which represents ~14% of the genome.

* Run GATK `GenotypeGVCFs` to call cohort variants on the scaffold intervals.

```
nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -stand_call_conf 15 -L ./genome_lists/scaffold.list1.intervals -V ./sample.all.gvcf.batch2.list -o ./vcf/hirundo_rustica.raw.scaffold.list1.batch2.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list1.vcf.recalc.batch2.log &

nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -stand_call_conf 15 -L ./genome_lists/scaffold.list2.intervals -V ./sample.all.gvcf.batch2.list -o ./vcf/hirundo_rustica.raw.scaffold.list2.batch2.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list2.vcf.recalc.batch2.log &

nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -stand_call_conf 15 -L ./genome_lists/scaffold.list3.intervals -V ./sample.all.gvcf.batch2.list -o ./vcf/hirundo_rustica.raw.scaffold.list3.batch2.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list3.vcf.recalc.batch2.log &

nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -stand_call_conf 15 -L ./genome_lists/scaffold.list4.intervals -V ./sample.all.gvcf.batch2.list -o ./vcf/hirundo_rustica.raw.scaffold.list4.batch2.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list4.vcf.recalc.batch2.log &

nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -stand_call_conf 15 -L ./genome_lists/scaffold.list5.intervals -V ./sample.all.gvcf.batch2.list -o ./vcf/hirundo_rustica.raw.scaffold.list5.batch2.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list5.vcf.recalc.batch2.log &

nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -stand_call_conf 15 -L ./genome_lists/scaffold.list6.intervals -V ./sample.all.gvcf.batch2.list -o ./vcf/hirundo_rustica.raw.scaffold.list6.batch2.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list6.vcf.recalc.batch2.log &

nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -stand_call_conf 15 -L ./genome_lists/scaffold.list7.intervals -V ./sample.all.gvcf.batch2.list -o ./vcf/hirundo_rustica.raw.scaffold.list7.batch2.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list7.vcf.recalc.batch2.log &

nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -stand_call_conf 15 -L ./genome_lists/scaffold.list8.intervals -V ./sample.all.gvcf.batch2.list -o ./vcf/hirundo_rustica.raw.scaffold.list8.batch2.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list8.vcf.recalc.batch2.log &

nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -stand_call_conf 15 -L ./log/scaffold.list9.intervals -V ./sample.all.gvcf.batch2.list -o ./vcf/hirundo_rustica.raw.scaffold.list9.batch2.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list9.vcf.recalc.batch2.log &

nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -stand_call_conf 15 -L ./genome_lists/scaffold.list10.intervals -V ./sample.all.gvcf.batch2.list -o ./vcf/hirundo_rustica.raw.scaffold.list10.batch2.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list10.vcf.recalc.batch2.log &

nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -stand_call_conf 15 -L ./genome_lists/scaffold.list11.intervals -V ./sample.all.gvcf.batch2.list -o ./vcf/hirundo_rustica.raw.scaffold.list11.batch2.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list11.vcf.recalc.batch2.log &

nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -stand_call_conf 15 -L ./genome_lists/scaffold.list12.intervals -V ./sample.all.gvcf.batch2.list -o ./vcf/hirundo_rustica.raw.scaffold.list12.batch2.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list12.vcf.recalc.batch2.log &

nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -stand_call_conf 15 -L ./genome_lists/scaffold.list13.intervals -V ./sample.all.gvcf.batch2.list -o ./vcf/hirundo_rustica.raw.scaffold.list13.batch2.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list13.vcf.recalc.batch2.log &
```

Note: to avoid duplicating large vcf files, these cohort calls are run in a shared directory `/data2/hirundo_fertilization_network_2023`. The path list `sample.all.gvcf.batch2.vcf.gz` is the same as `/genom_lists/sample.ID.list.gvcf.paths`, and the former properly maps to the correct subdirectories containing all 1034 indexed input gvcf files. If all work is being done in a single parent directory modify the code for scaffolds 1-13.
`nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -stand_call_conf 15 -L ./genome_lists/scaffold.list1.intervals -V ./genome_lists/sample.ID.list.gvcf.paths -o ./vcf/hirundo_rustica.raw.scaffold.list1.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list1.vcf.log &`

* Check on running progress.
```
for i in ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list*.log; do echo $i; grep 'INFO' $i | tail -n 1; done
```

## 6. Filtering for relatedness

### 6.1 Run Picard `MergeVCFs` to merge the VCF files.

* Make a formated list of paths to scaffolds that do not contain sex chromosomes (scaffold lists 1-11) to be merged for relatedness calculations `vcf.relatedness.path.list` using nano.
Note: In the shared directory use `vcf.recalc.batch2.list`.

* Run Picard `MergeVCFs` to merge the VCF files.

```
nohup java -jar picard.jar MergeVcfs -I vcf.recalc.batch2.list -O ./vcf/hirundo_rustica.recalc.batch2.vcf.gz > ./log/mergeVcfs.recalc.batch2.picard.log &
```

Note: In the shared directory `/data2/hirundo_fertilization_network_2023` use `vcf.relatedness.path.list`. If all work is being done in a single parent directory modify the code.
`nohup java -jar picard.jar MergeVcfs -I vcf.relatedness.path.list -O ./vcf/hirundo_rustica.merge.relatedness.vcf.gz > ./log/mergeVcfs.relatedness.picard.log &`. Also might need to select correct java version. Check current `java -version` and run `sudo update-alternatives --config java`, and enter the number to choose java-8.


### 6.2 Mask repeats, indels, and non-bialleic SNPs

Use bcftools to mask repeats, indels and non-bialleic SNPs. Usage: bcftools view --threads 24 - T <refgenome> <fileoutname> <vcf.gzfile>

```
nohup bcftools view --threads 24 -T ^/data3/hirundo/genome_annotation/repeatmasker/Hirundo_rustica_bHirRus1.final.fasta.repeat.bed -m2 -M2 -U -v snps -O z -o ./vcf/hirundo_rustica.recalc.batch2.mask.vcf.gz ./vcf/hirundo_rustica.recalc.batch2.vcf.gz &
```

Produces this file `hirundo_rustica.recalc.batch2.mask.vcf.gz` in the `vcf` directory which includes 1034 samples
 

Note: In the shared directory `/data2/hirundo_fertilization_network_2023` use `hirundo_rustica.recalc.batch2.vcf.gz` as the input file. If all work is being done in a single parent directory modify the code.
`nohup bcftools view --threads 24 -T ^/data3/hirundo/genome_annotation/repeatmasker/Hirundo_rustica_bHirRus1.final.fasta.repeat.bed -m2 -M2 -U -v snps -O z -o ./vcf/hirundo_rustica.merge.relatedness.mask.vcf.gz ./vcf/hirundo_rustica.merge.relatedness.vcf.gz &`


### 6.3 Filer MAF >= 0.05

Use bcftools filter by MAF >= 0.05. Usage: bcftools view --threads 16 - S <samplelist> <fileoutname> <target_vcf.gzfile>

```
nohup bcftools view --threads 16 -S list.CO-2021-2023.no.rpts.indiv.txt -i 'MAF >= 0.05' -m2 -M2 -U -v snps -O z -o ./vcf/hirundo_rustica.recalc.mask.CO-2021-2023.maf05.vcf.gz ./vcf/hirundo_rustica.recalc.batch2.mask.vcf.gz &

```

Note: In the shared directory `/data2/hirundo_fertilization_network_2023` use `list.CO-2021-2023.no.rpts.indiv.txt` as the sample list. If all work is being done in a single parent directory modify the code.
`nohup bcftools view --threads 16 -S sample.ID.list -i 'MAF >= 0.05' -m2 -M2 -U -v snps -O z -o ./vcf/hirundo_rustica.merge.relatedness.mask.maf05.vcf.gz ./vcf/hirundo_rustica.merge.relatedness.mask.vcf.gz &`


### 6.4 Thin SNPs to prune for LD

Use VCF tools to thin SNPs and prune for LD. Usage: vcftools --gzvcf <target_vcf.gzfile> <fileoutname> 

```
vcftools --gzvcf ./vcf/hirundo_rustica.merge.relatedness.mask.maf05.vcf.gz --thin 10000 --recode --stdout | bgzip -c > ./vcf/hirundo_rustica.merge.relatedness.mask.maf05.thin10k.vcf.gz
```


### 6.5 Make formatted SNP positions file
    	
Use bcftools to generate a formatted SNP positions file. Run from `relate` subdirectory. Usage: bcftools query -f <target_thin10k.vcf.gzfile> <fileoutname> 

```    	
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ./vcf/hirundo_rustica.merge.relatedness.mask.maf05.thin10k.vcf.gz > 
hirundo_rustica.relatedness.mask.maf05.thin10.snp.list
```



## 7. Calculate relatedness

### 7.1 Call genotype likelihoods 

* Format list of paths to all individual bam files using `touch`. 
```
touch indiv.CO-2022.all-bam; for i in `cat list.CO-2022.indv`; do echo /media/safranlab/DataStore1/data/hirundo_lcwgs_2022/bam/${i}.bam >> indiv.CO-2022.all-bam; done
```

Note: In the shared directory `/data2/hirundo_fertilization_network_2023` use `list.CO-2021-2023.no.rpts.indiv.txt` as the sample list. If all work is being done in a single parent directory modify the code.
`touch sample.ID.list.bam.paths; for i in `cat sample.ID.list`; do echo ./bam/${i}.bam >> sample.ID.list.bam.paths; done`. 

* Use cat to combine the list of file path names to BAM files. Note bam path file lists were already created for 2021 `indiv.CO-2021.all-bam ` and 2023 `indiv.CO-2023batch2.all-bam` birds, but skip this step if all files are stored in a single parent directory. 

```
cat indiv.CO-2021.all-bam indiv.CO-2022.all-bam indiv.CO-2023batch2.all-bam > indiv.CO-2021-2023.all-bam
```

Note: BAM files are repeated for several birds that were sampled in multiple years. Remove the older bam files from the `indiv.CO-2021-2023.all-bam` list. 

Note: If bam files not indexed, then loop through all .bam files in the `./bam/` directory and index them.

```
for bam_file in *.bam; do
  echo "Indexing: $bam_file"
  samtools index "$bam_file"
done
```


* Run `SNPbam2vcf.py` to call gentoype likelihoods of selected SNPs. Usage: ./SNPbam2vcf.py <bamlist> <fileoutname> <target_SNP_file>. This took about 24 hrs to run. 

```
sudo  ./SNPbam2vcf.py sample.ID.list.bam.paths ./vcf/hirundo_rustica.merge.relatedness.mask.maf05.thin10k.rawGL.vcf.gz hirundo_rustica.relatedness.mask.maf05.thin10.snp.list
```


### 7.2 Generate pairwise relatedness - lcMLkin
Use `lcmlkin` on recalculated genotype likelihood data to calculate pairwise relatedness

```			
./lcmlkin -i ./vcf/hirundo_rustica.merge.relatedness.mask.maf05.thin10k.rawGL.vcf.gz -o ./relate/hirundo_rustica_pw_relate -g all -t 12
```



## 8. Filtering for imputation

### 8.1 Hard filtering of variants

* Run GATK `VariantFiltration` to retain at set of high-quality SNPs from the raw variant calls, stored in 13 separate interval files from step 5.3 above. Set hard filters on gentoypes in parallel to reduce runtime. Note: Hard filtering follows [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants)

```
./gatk-4.0.8.1/gatk VariantFiltration -V /data2/hirundo_fertilization_network_2023/vcf/hirundo_rustica.raw.scaffold.list1.batch2.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica.allsites.HardFilter.scaffold.list1.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list1.log

./gatk-4.0.8.1/gatk VariantFiltration -V /data2/hirundo_fertilization_network_2023/vcf/hirundo_rustica.raw.scaffold.list2.batch2.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica.allsites.HardFilter.scaffold.list2.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list2.log

./gatk-4.0.8.1/gatk VariantFiltration -V /data2/hirundo_fertilization_network_2023/vcf/hirundo_rustica.raw.scaffold.list3.batch2.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica.allsites.HardFilter.scaffold.list3.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list3.log

./gatk-4.0.8.1/gatk VariantFiltration -V /data2/hirundo_fertilization_network_2023/vcf/hirundo_rustica.raw.scaffold.list4.batch2.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica.allsites.HardFilter.scaffold.list4.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list4.log

./gatk-4.0.8.1/gatk VariantFiltration -V /data2/hirundo_fertilization_network_2023/vcf/hirundo_rustica.raw.scaffold.list5.batch2.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica.allsites.HardFilter.scaffold.list5.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list5.log

./gatk-4.0.8.1/gatk VariantFiltration -V /data2/hirundo_fertilization_network_2023/vcf/hirundo_rustica.raw.scaffold.list6.batch2.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica.allsites.HardFilter.scaffold.list6.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list6.log

./gatk-4.0.8.1/gatk VariantFiltration -V /data2/hirundo_fertilization_network_2023/vcf/hirundo_rustica.raw.scaffold.list7.batch2.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica.allsites.HardFilter.scaffold.list7.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list7.log

./gatk-4.0.8.1/gatk VariantFiltration -V /data2/hirundo_fertilization_network_2023/vcf/hirundo_rustica.raw.scaffold.list8.batch2.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica.allsites.HardFilter.scaffold.list8.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list8.log

./gatk-4.0.8.1/gatk VariantFiltration -V /data2/hirundo_fertilization_network_2023/vcf/hirundo_rustica.raw.scaffold.list9.batch2.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica.allsites.HardFilter.scaffold.list9.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list9.log

./gatk-4.0.8.1/gatk VariantFiltration -V /data2/hirundo_fertilization_network_2023/vcf/hirundo_rustica.raw.scaffold.list10.batch2.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica.allsites.HardFilter.scaffold.list10.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list10.log

./gatk-4.0.8.1/gatk VariantFiltration -V /data2/hirundo_fertilization_network_2023/vcf/hirundo_rustica.raw.scaffold.list11.batch2.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica.allsites.HardFilter.scaffold.list11.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list11.log

./gatk-4.0.8.1/gatk VariantFiltration -V /data2/hirundo_fertilization_network_2023/vcf/hirundo_rustica.raw.scaffold.list12.batch2.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica.allsites.HardFilter.scaffold.list12.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list12.log

./gatk-4.0.8.1/gatk VariantFiltration -V /data2/hirundo_fertilization_network_2023/vcf/hirundo_rustica.raw.scaffold.list13.batch2.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica.allsites.HardFilter.scaffold.list13.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list13.log

```


### 8.2 Run Picard `MergeVCFs` to merge the VCF files.

* Make a formatted list of paths to hard filtered vcf files (scaffold lists 1-13) to be merged `vcf.HardFilter.path.list` using nano.

* Select the correct version of java to work with Picard (java-17-openjdk-amd64)
sudo update-alternatives --config java

* Run Picard `MergeVCFs` to merge the VCF files.

```
java -jar picard.jar MergeVcfs -I ./genome_lists/vcf.HardFilter.path.list -O ./vcf/hirundo_rustica.merge.HardFilter.vcf.gz > ./log/mergeVcfs.HardFilter.picard.log &
```


### 8.3 Mask indels and hard filter genotypes

Use bcftools to mask indels, and exclude/filter variants where quality by depth <= 2, Fisher Strand value >= 60 etc. Then index the recoded vcf file. Usage: bcftools view --threads 20 - T  <fileoutname> <vcf.gzfile>

```
bcftools filter --threads 20 -e 'TYPE="indel" || FILTER="QD2" || FILTER="FS60" || FILTER="MQ40" || FILTER="MQRankSum-12.5" || FILTER="ReadPosRankSum-8"' --set-GTs . -O z -o ./vcf/hirundo_rustica.merge.HardFilter.recode.vcf.gz ./vcf/hirundo_rustica.merge.HardFilter.vcf.gz

bcftools index --threads 48 -t ./vcf/hirundo_rustica.merge.HardFilter.recode.vcf.gz 
```


### 8.4 Filter samples with high missing data

Query a subset of the filtered VCF to quantify proportion of missing genotypes across samples.


```
bcftools view --threads 16 -r NC_053451.1 -m2 -M2 -U -v snps -O z -o ./vcf/hirundo_rustica.merge.HardFilter.recode.snps.chr1.vcf.gz ./vcf/hirundo_rustica.merge.HardFilter.recode.vcf.gz > ./log/hirundo_rustica.merge.HardFilter.recode.snps.chr1.log

vcftools --gzvcf ./vcf/hirundo_rustica.merge.HardFilter.recode.snps.chr1.vcf.gz --missing-indv --out ./log/hirundo_rustica.merge.HardFilter.recode.snps.chr1

rm ./vcf/hirundo_rustica.merge.HardFilter.recode.snps.chr1.vcf.gz
```

Based on this analysis, samples CO_E101 and CO_E087 had > 50% missing data.

* Make a formatted list of samples to be removed due to high missingness `./genome_lists/sample.ID.remove.list` using nano.

* Remove samples in list, and then re-index the merged vcf file.
```
bcftools view --threads 16 -S ^./genome_lists/sample.ID.remove.list -O z -o ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.vcf.gz ./vcf/hirundo_rustica.merge.HardFilter.recode.vcf.gz 
tabix -p vcf ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.vcf.gz
```



## 9. Process sex chromosomes for imputation
Females are hemizygous for the Z chromosome, and so cannot have heterozygous sites on Z-linked scaffolds. We'll extract SNPs for Z-linked scaffolds, then identify SNPs in which any known females have heterozygous genotypes, then mask these sites in all individuals. We also want to remove males from the W chromosome, as they don't have one, and also to mask heterozygous sites on the W in females, since these can't exist.

We'll begin by extracting VCF tables for chromosome-assigned autosome, Z, and W scaffolds, respectively. From here forward we'll also omit un-assigned scaffolds from analysis.

### 9.1 Extract VCF tables for chromosome-assigned autosome, Z, and W scaffolds, respectively

* Format autosome, Z chromosome, and W chromosome-assigned scaffold BED files.
`
./genome_lists/scaffold.assigned-autosome.bed
./genome_lists/scaffold.assigned-chrZ.bed
./genome_lists/scaffold.assigned-chrW.bed
` 
Note: After this step, un-assigned scaffolds are omitted from downstream analysis.


* Extract VCFs based on BED files.
```
bcftools view --threads 16 -R ./genome_lists/scaffold.assigned-autosome.bed -O z -o ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.auto.vcf.gz ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.vcf.gz

bcftools view --threads 16 -R ./genome_lists/scaffold.assigned-chrZ.bed -O z -o ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.chrZ.vcf.gz ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.vcf.gz

bcftools view --threads 16 -R ./genome_lists/scaffold.assigned-chrW.bed -O z -o ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.chrW.vcf.gz ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.vcf.gz
```

* Rename and index autosome all-sites VCF.
```
mv ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.auto.vcf.gz ./vcf/hirundo_rustica.final.auto.vcf.gz
tabix -p vcf ./vcf/hirundo_rustica.final.auto.vcf.gz
```


### 9.2 Mask female heterozygous sites Z-chromosome

* Extract biallelic SNPs on Z-linked scaffolds.
```
bcftools view --threads 16 -m2 -M2 -U -v snps -O z -o ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.chrZ.snps.vcf.gz ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.chrZ.vcf.gz > ./log/hirundo_rustica.merge.HardFilter.recode.remove.chrZ.snps.vcf.gz
```

* Format `./genome_lists/sample.ID.female.list` using nano.

* Run `ZchrFemaleHeterozygous.py` to identify and write positions of the Z chromosome with heterozygous genotype calls in any females.
```
gunzip ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.chrZ.snps.vcf.gz 
python ZchrFemaleHeterozygous.py ./genome_lists/sample.ID.female.list ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.chrZ.snps.vcf ./log/hirundo_rustica.merge.HardFilter.recode.remove.chrZ.snps.FemaleZhetSites.txt
```

* Convert output to BED format and index it as a GATK feature file.
```
awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2}' ./log/hirundo_rustica.merge.HardFilter.recode.remove.chrZ.snps.FemaleZhetSites.txt > ./log/hirundo_rustica.merge.HardFilter.recode.remove.chrZ.snps.FemaleZhetSites.bed
./gatk-4.0.8.1/gatk IndexFeatureFile --feature-file ./log/hirundo_rustica.merge.HardFilter.recode.remove.chrZ.snps.FemaleZhetSites.bed
```

* Select the correct version of java to work with GATK (java-8-openjdk-amd64)
sudo update-alternatives --config java

* Run GATK `VariantFiltration` to annotate female heterozygous sites and mask with `bcftools`.
```
tabix -p vcf ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.chrZ.vcf.gz

./gatk-4.0.8.1/gatk VariantFiltration -V ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.chrZ.vcf.gz --mask ./log/hirundo_rustica.merge.HardFilter.recode.remove.chrZ.snps.FemaleZhetSites.bed --mask-name ZHET -O ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.chrZ.filter.vcf.gz > ./log/GATK_VariantFiltration_hirundo_rustica.merge.HardFilter.recode.remove.chrZ.filter.log

bcftools filter --threads 16 -e 'FILTER="ZHET"' --set-GTs . -O z -o ./vcf/hirundo_rustica.final.chrZ.vcf.gz ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.chrZ.filter.vcf.gz

tabix -p vcf ./vcf/hirundo_rustica.final.chrZ.vcf.gz
```


### 9.3 Mask female heterozygous sites W-chromosome

* Remove males from W chromosome VCF.
```
bcftools view --threads 16 -S ./genome_lists/sample.ID.female.list -O z -o ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.chrW.female.vcf.gz ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.chrW.vcf.gz
```

* Extract W-linked SNPs.
```
bcftools view --threads 16 -m2 -M2 -U -v snps -O v -o ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.chrW.female.snps.vcf ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.chrW.female.vcf.gz
```

* Identify heterozygous sites on the W chromosome.
```
python ZchrFemaleHeterozygous.py ./genome_lists/sample.ID.female.list ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.chrW.female.snps.vcf  ./log/hirundo_rustica.merge.HardFilter.recode.remove.chrW.female.snps.vcf.FemaleHetSites.txt
```

* Convert output to BED format and index it as a GATK feature file.
```
awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2}' ./log/hirundo_rustica.merge.HardFilter.recode.remove.chrW.female.snps.vcf.FemaleHetSites.txt > ./log/hirundo_rustica.merge.HardFilter.recode.remove.chrW.female.snps.vcf.FemaleHetSites.bed
./gatk-4.0.8.1/gatk IndexFeatureFile --feature-file ./log/hirundo_rustica.merge.HardFilter.recode.remove.chrW.female.snps.vcf.FemaleHetSites.bed
```

* Select the correct version of java to work with GATK (java-8-openjdk-amd64)
sudo update-alternatives --config java

* Run GATK `VariantFiltration` to annotate female heterozygous sites and mask with `bcftools`.
```
tabix -p vcf ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.chrW.female.vcf.gz

./gatk-4.0.8.1/gatk VariantFiltration -V ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.chrW.female.vcf.gz --mask ./log/hirundo_rustica.merge.HardFilter.recode.remove.chrW.female.snps.vcf.FemaleHetSites.bed --mask-name WHET -O ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.chrW.female.filter.vcf.gz > ./log/GATK_VariantFiltration_hirundo_rustica.merge.HardFilter.recode.remove.chrW.female.filter.log 

bcftools filter --threads 16 -e 'FILTER="WHET"' --set-GTs . -O z -o ./vcf/hirundo_rustica.final.chrW.vcf.gz ./vcf/hirundo_rustica.merge.HardFilter.recode.remove.chrW.female.filter.vcf.gz 

tabix -p vcf ./vcf/hirundo_rustica.final.chrW.vcf.gz 
```



## 10. Repeat masking for imputation

We want to avoid analyzing genotypes in repetitive regions of the genome. We'll generate a repeat annotation for filtering using `repeatmasker`. The steps to annotate repeats using RepeatMaster are detailed in the Schield et al 2025.

### 10.1. Index repeat BED file.

* Run GATK to index feature file 
```
./gatk-4.0.8.1/gatk IndexFeatureFile --feature-file ./genome_annotation/Hirundo_rustica_bHirRus1.final.fasta.repeat.bed
```


### 10.2 Annotate repeats

* Annotate repeats using GATK `VariantFiltration`.
```
./gatk-4.0.8.1/gatk VariantFiltration -V ./vcf/hirundo_rustica.final.auto.vcf.gz --mask ./genome_annotation/Hirundo_rustica_bHirRus1.final.fasta.repeat.bed --mask-name REP -O ./vcf/hirundo_rustica.final.auto.tmp-rep.vcf.gz > ./log/GATK_VariantFiltration_repeat-mask.hirundo_rustica.auto.log

./gatk-4.0.8.1/gatk VariantFiltration -V ./vcf/hirundo_rustica.final.chrZ.vcf.gz --mask ./genome_annotation/Hirundo_rustica_bHirRus1.final.fasta.repeat.bed --mask-name REP -O ./vcf/hirundo_rustica.final.chrZ.tmp-rep.vcf.gz > ./log/GATK_VariantFiltration_repeat-mask.hirundo_rustica.chrZ.log

./gatk-4.0.8.1/gatk VariantFiltration -V ./vcf/hirundo_rustica.final.chrW.vcf.gz --mask ./genome_annotation/Hirundo_rustica_bHirRus1.final.fasta.repeat.bed --mask-name REP -O ./vcf/hirundo_rustica.final.chrW.tmp-rep.vcf.gz > ./log/GATK_VariantFiltration_repeat-mask.hirundo_rustica.chrW.log
```

### 10.3 Recode repeats as missing genotypes and index VCFs.

```
bcftools filter --threads 24 -e 'FILTER="REP"' --set-GTs . -O z -o ./vcf/hirundo_rustica.final.auto.vcf.gz  ./vcf/hirundo_rustica.final.auto.tmp-rep.vcf.gz

bcftools filter --threads 24 -e 'FILTER="REP"' --set-GTs . -O z -o ./vcf/hirundo_rustica.final.chrZ.vcf.gz ./vcf/hirundo_rustica.final.chrZ.tmp-rep.vcf.gz

bcftools filter --threads 24 -e 'FILTER="REP"' --set-GTs . -O z -o ./vcf/hirundo_rustica.final.chrW.vcf.gz ./vcf/hirundo_rustica.final.chrW.tmp-rep.vcf.gz

tabix -p vcf ./vcf/hirundo_rustica.final.auto.vcf.gz 
tabix -p vcf ./vcf/hirundo_rustica.final.chrZ.vcf.gz
tabix -p vcf ./vcf/hirundo_rustica.final.chrW.vcf.gz 
```



## 11. Additional processing and filtering for imputation

### 11.1 Extract biallelic SNPs.
```
bcftools view --threads 16 -m2 -M2 -U -v snps -O z -o ./vcf/hirundo_rustica.final.auto.snps.vcf.gz ./vcf/hirundo_rustica.final.auto.vcf.gz 

bcftools view --threads 16 -m2 -M2 -U -v snps -O z -o ./vcf/hirundo_rustica.final.chrZ.snps.vcf.gz ./vcf/hirundo_rustica.final.chrZ.vcf.gz 

bcftools view --threads 16 -m2 -M2 -U -v snps -O z -o ./vcf/hirundo_rustica.final.chrW.snps.vcf.gz ./vcf/hirundo_rustica.final.chrW.vcf.gz 
```

# MIGHT BE TOO STRICT - SKIP
# 
# ### 11.2 Extract SNPs with < 20% missing data.
# ```
# bcftools view -i 'F_MISSING<0.2' -O z -o ./vcf/hirundo_rustica.final.auto.snps.miss02.vcf.gz ./vcf/hirundo_rustica.final.auto.snps.vcf.gz 
# 
# bcftools view -i 'F_MISSING<0.2' -O z -o ./vcf/hirundo_rustica.final.chrZ.snps.miss02.vcf.gz ./vcf/hirundo_rustica.final.chrZ.snps.vcf.gz
# 
# bcftools view -i 'F_MISSING<0.2' -O z -o ./vcf/hirundo_rustica.final.chrW.snps.miss02.vcf.gz ./vcf/hirundo_rustica.final.chrW.snps.vcf.gz
# ```


### 11.2 Impose minor allele frequency and minor allele count filters.
```
bcftools view --threads 16 -c 2 -q 0.05:minor -m2 -M2 -O z -o ./vcf/hirundo_rustica.final.auto.snps.maf05.vcf.gz ./vcf/hirundo_rustica.final.auto.snps.vcf.gz

bcftools view --threads 16 -c 2 -q 0.05:minor -m2 -M2 -O z -o ./vcf/hirundo_rustica.final.chrZ.snps.maf05.vcf.gz ./vcf/hirundo_rustica.final.chrZ.snps.vcf.gz

bcftools view --threads 16 -c 2 -q 0.05:minor -m2 -M2 -O z -o ./vcf/hirundo_rustica.final.chrW.snps.maf05.vcf.gz ./vcf/hirundo_rustica.final.chrW.snps.vcf.gz

tabix -p vcf ./vcf/hirundo_rustica.final.auto.snps.maf05.vcf.gz 
tabix -p vcf ./vcf/hirundo_rustica.final.chrZ.snps.maf05.vcf.gz 
tabix -p vcf ./vcf/hirundo_rustica.final.chrW.snps.maf05.vcf.gz  
```

NOTE: Search each file for candidate snps identifed in ventral color gwas from Schield et al. 2025 

### ventral color SNP list
* NC_053453.1	44285485
* NC_053453.1	46520417
* NC_053488.1	18954203
* NC_053488.1	19018140
* NC_053488.1	19332900
* NC_053488.1	33380365
* NC_053488.1	35903569
* NC_053488.1	44206228
* NC_053488.1	46065026
* NC_053488.1	46901341

### tail SNP list
* NC_053450.1	98076904
* NC_053450.1	98076904
* NC_053450.1	100834905

### File list
* hirundo_rustica.final.auto.snps.maf05.vcf.gz                    
* hirundo_rustica.final.auto.snps.vcf.gz
* hirundo_rustica.final.chrW.snps.maf05.vcf.gz                    
* hirundo_rustica.final.chrW.snps.vcf.gz
 
 
### 11.3 Filter sites with extreme read depths and on unassigned scaffolds

We want to remove sites that have extreme mean coverage across all samples as these will be the most likely to produce spurious imputation results.

* Use `vcftools` to calculate mean depth per site:
```
vcftools --gzvcf hirundo_rustica.final.auto.snps.maf05.vcf.gz  --site-mean-depth --out hirundo_rustica.final.auto.snps.maf05 

vcftools --gzvcf hirundo_rustica.final.chrZ.snps.maf05.vcf.gz  --site-mean-depth --out hirundo_rustica.final.chrZ.snps.maf05 

vcftools --gzvcf hirundo_rustica.final.chrW.snps.maf05.vcf.gz  --site-mean-depth --out hirundo_rustica.final.chrW.snps.maf05 
```

Note: there is likely a faster `bcftools` solution to this, using the `query` subcommand to output information from the FORMAT field.


### 11.4 Use AWK to calculate mean, SD, 2.5%, and 97.5% quantiles of read depth across all sites in the autosome, Z and W chromosomes

Note: column three is the 'MEAN_DEPTH' for each site 
```
awk 'NR > 1 { sum += $3; count++ } END { if (count > 0) print sum / count }' hirundo_rustica.final.auto.snps.maf05.ldepth.mean

awk 'NR > 1 {sum+=$3; sumsq+=$3*$3} END {print sqrt(sumsq/NR - (sum/NR)^2)}'  hirundo_rustica.final.auto.snps.maf05.ldepth.mean

sort -n -k3,3 hirundo_rustica.final.auto.snps.maf05.ldepth.mean | awk '{a[NR]=$3} END{print a[int(NR*0.025 + 0.5)]}'

sort -n -k3,3 hirundo_rustica.final.auto.snps.maf05.ldepth.mean | awk '{a[NR]=$3} END{print a[int(NR*0.975 + 0.5)]}'

awk 'NR > 1 { sum += $3; count++ } END { if (count > 0) print sum / count }' hirundo_rustica.final.chrZ.snps.maf05.ldepth.mean

awk 'NR > 1 {sum+=$3; sumsq+=$3*$3} END {print sqrt(sumsq/NR - (sum/NR)^2)}'  hirundo_rustica.final.chrZ.snps.maf05.ldepth.mean

sort -n -k3,3 hirundo_rustica.final.chrZ.snps.maf05.ldepth.mean | awk '{a[NR]=$3} END{print a[int(NR*0.025 + 0.5)]}'

sort -n -k3,3 hirundo_rustica.final.chrZ.snps.maf05.ldepth.mean | awk '{a[NR]=$3} END{print a[int(NR*0.975 + 0.5)]}'

awk 'NR > 1 { sum += $3; count++ } END { if (count > 0) print sum / count }' hirundo_rustica.final.chrW.snps.maf05.ldepth.mean

awk 'NR > 1 {sum+=$3; sumsq+=$3*$3} END {print sqrt(sumsq/NR - (sum/NR)^2)}'  hirundo_rustica.final.chrW.snps.maf05.ldepth.mean

sort -n -k3,3 hirundo_rustica.final.chrW.snps.maf05.ldepth.mean | awk '{a[NR]=$3} END{print a[int(NR*0.025 + 0.5)]}'

sort -n -k3,3 hirundo_rustica.final.chrW.snps.maf05.ldepth.mean | awk '{a[NR]=$3} END{print a[int(NR*0.975 + 0.5)]}
```

* Autosome mean depth = 2.38042, SD = 1.17345, 2.5th depth quantile = 1.74903, and 97.5th depth quantile = 2.65504. 

* Z-chromosome mean depth = 1.80115, SD = 0.198465, 2.5th depth quantile = 1.28973, and 97.5th depth quantile = 2.05523.

* W-chromosome mean depth = 1.2481, SD = 0.454868, 2.5th depth quantile = 0.423913, and 97.5th depth quantile = 2.32609. 


### 11.5 Recode sites with mean depth above 3.5 (which is 3 SD from autosomes) as missing data to avoid effects of potential paralogous mappings and remove sites on unassigned scaffolds using `bcftools`.

This command will output sites on scaffolds listed in `/genome_annotation/scaffold.assigned...bed` and recode sites with extremely high mean read depth:


```
bcftools view --threads 16 -R ./genome_annotation/scaffold.assigned-autosome.bed ./vcf/hirundo_rustica.final.auto.snps.maf05.vcf.gz | bcftools filter --threads 16 -e 'FMT/DP > 3.5' --set-GTs . -O z -o ./vcf/hirundo_rustica.final.auto.snps.maf05.depth.chrom.vcf.gz
tabix -p vcf ./vcf/hirundo_rustica.final.auto.snps.maf05.depth.chrom.vcf.gz

bcftools view --threads 16 -R ./genome_annotation/scaffold.assigned-chrZ.bed ./vcf/hirundo_rustica.final.chrZ.snps.maf05.vcf.gz | bcftools filter --threads 16 -e 'FMT/DP > 3.5' --set-GTs . -O z -o ./vcf/hirundo_rustica.final.chrZ.snps.maf05.depth.chrom.vcf.gz
tabix -p vcf ./vcf/hirundo_rustica.final.chrZ.snps.maf05.depth.chrom.vcf.gz

bcftools view --threads 16 -R ./genome_annotation/scaffold.assigned-chrW.bed ./vcf/hirundo_rustica.final.chrW.snps.maf05.vcf.gz | bcftools filter --threads 16 -e 'FMT/DP > 3.5' --set-GTs . -O z -o ./vcf/hirundo_rustica.final.chrW.snps.maf05.depth.chrom.vcf.gz
tabix -p vcf ./vcf/hirundo_rustica.final.chrW.snps.maf05.depth.chrom.vcf.gz
```

* Spot check a few examples where genotypes contain exrtremely high depth

```
bcftools view hirundo_rustica.final.auto.snps.maf05.vcf.gz NC_053453.1:44285485

bcftools view hirundo_rustica.final.auto.snps.maf05.depth.chrom.vcf.gz NC_053453.1:44285485
```


### 11.6 Copy read depth summaries to local drive for storage.
```
scp zchlaubach@Terminator:/data2/hirundo_plumage_MR/vcf/hirundo_rustica.final.auto.snps.maf05.ldepth.mean /Users/zacharylaubach/WD/Git/swallow/bs_MR/cov

scp zchlaubach@Terminator:/data2/hirundo_plumage_MR/vcf/hirundo_rustica.final.chrZ.snps.maf05.ldepth.mean /Users/zacharylaubach/WD/Git/swallow/bs_MR/cov

scp zchlaubach@Terminator:/data2/hirundo_plumage_MR/vcf/hirundo_rustica.final.chrW.snps.maf05.ldepth.mean /Users/zacharylaubach/WD/Git/swallow/bs_MR/cov
```


### 11.7 Cleaning up

These filtering steps produce many intermediate VCF files, especially in the initial parallel processing steps for raw and hard filtered VCFs. Delete these to save some disk space.
```
rm ./vcf/hirundo_rustica....*

```



## 12. Impute and call genotypes

### 12.1 Install Beagle - the software used for imputation
```
wget https://faculty.washington.edu/browning/beagle/beagle.28Jun21.220.jar
```


### 12.2 Format SNP data for imputation

* Concatenate autosomal and Z chromosome VCFs
```
bcftools concat -O z -o ./vcf/hirundo_rustica.final.auto+chrZ.snps.maf05.depth.chrom.vcf.gz ./vcf/hirundo_rustica.final.auto.snps.maf05.depth.chrom.vcf.gz ./vcf/hirundo_rustica.final.chrZ.snps.maf05.depth.chrom.vcf.gz
```

* Extract SNP data for males
```
bcftools view --threads 16 -S ./genome_lists/sample.ID.male.list -O v ./vcf/hirundo_rustica.final.auto+chrZ.snps.maf05.depth.chrom.vcf.gz | bcftools filter -e 'MAF <= 0.05' -O z -o ./vcf/hirundo_rustica.final.auto+chrZ.snps.maf05.depth.chrom.male.vcf.gz
```

* Extract SNP data for females
```
bcftools view --threads 16 -S ./genome_lists/sample.ID.female.list -O v ./vcf/hirundo_rustica.final.auto+chrZ.snps.maf05.depth.chrom.vcf.gz | bcftools filter -e 'MAF <= 0.05' -O z -o ./vcf/hirundo_rustica.final.auto+chrZ.snps.maf05.depth.chrom.female.vcf.gz

bcftools view --threads 16 -S ./genome_lists/sample.ID.female.list -O z -o ./vcf/hirundo_rustica.final.chrW.snps.maf05.depth.chrom.female.vcf.gz ./vcf/hirundo_rustica.final.chrW.snps.maf05.depth.chrom.vcf.gz

bcftools concat --threads 16 -O z -o ./vcf/hirundo_rustica.final.auto+chrZ+chrW.snps.maf05.depth.chrom.female.vcf.gz ./vcf/hirundo_rustica.final.auto+chrZ.snps.maf05.depth.chrom.female.vcf.gz ./vcf/hirundo_rustica.final.chrW.snps.maf05.depth.chrom.female.vcf.gz 

rm ./vcf/hirundo_rustica.final.chrW.snps.maf05.depth.chrom.female.vcf.gz 
rm ./vcf/hirundo_rustica.final.auto+chrZ.snps.maf05.depth.chrom.female.vcf.gz
```


### 12.3 Perform imputation using Beagle
```
java -Xmx96g -jar beagle.28Jun21.220.jar nthreads=16 gt=./vcf/hirundo_rustica.final.auto+chrZ.snps.maf05.depth.chrom.vcf.gz out=./vcf_imputed/hirundo_rustica.final.auto+chrZ.snps.maf05.depth.chrom.impute

java -Xmx96g -jar beagle.28Jun21.220.jar nthreads=16 gt=./vcf/hirundo_rustica.final.auto+chrZ.snps.maf05.depth.chrom.male.vcf.gz out=./vcf_imputed/hirundo_rustica.final.auto+chrZ.snps.maf05.depth.chrom.male.impute

java -Xmx96g -jar beagle.28Jun21.220.jar nthreads=16 gt=./vcf/hirundo_rustica.final.auto+chrZ+chrW.snps.maf05.depth.chrom.female.vcf.gz out=./vcf_imputed/hirundo_rustica.final.auto+chrZ+chrW.snps.maf05.depth.chrom.female.impute

tabix -p vcf ./vcf_imputed/hirundo_rustica.final.auto+chrZ.snps.maf05.depth.chrom.impute.vcf.gz

tabix -p vcf ./vcf_imputed/hirundo_rustica.final.auto+chrZ.snps.maf05.depth.chrom.male.impute.vcf.gz

tabix -p vcf ./vcf_imputed/hirundo_rustica.final.auto+chrZ+chrW.snps.maf05.depth.chrom.female.impute.vcf.gz
```


### 12.4 extract hard genotype calls for candidate loci using bcftools and copy to csv file
```
bcftools view hirundo_rustica.final.auto+chrZ.snps.maf05.depth.chrom.impute.vcf.gz NC_053453.1:44285485

bcftools view hirundo_rustica.final.auto+chrZ.snps.maf05.depth.chrom.impute.vcf.gz NC_053488.1:19332900

bcftools view hirundo_rustica.final.auto+chrZ.snps.maf05.depth.chrom.impute.vcf.gz NC_053488.1:4420622

bcftools view hirundo_rustica.final.auto+chrZ.snps.maf05.depth.chrom.impute.vcf.gz NC_053450.1:100834905
```


### 12.5 Copy imputed vcf file to local drive
```
scp zchlaubach@Terminator:/data2/hirundo_plumage_MR/vcf_imputed/hirundo_rustica.final.auto+chrZ.snps.maf05.depth.chrom.impute.vcf.gz /Users/zacharylaubach/WD/Git/swallow/bs_MR_plumage_traits_repro_success/raw_data

scp zchlaubach@Terminator:/data2/hirundo_plumage_MR/vcf_imputed/hirundo_rustica.final.auto+chrZ+chrW.snps.maf05.depth.chrom.female.impute.vcf.gz /Users/zacharylaubach/WD/Git/swallow/bs_MR_plumage_traits_repro_success/raw_data

scp zchlaubach@Terminator:/data2/hirundo_plumage_MR/vcf_imputed/hirundo_rustica.final.auto+chrZ.snps.maf05.depth.chrom.male.impute.vcf.gz /Users/zacharylaubach/WD/Git/swallow/bs_MR_plumage_traits_repro_success/raw_data

```





Maybe don't need to do this....

### 12.4 Extract VCFs and convert to Plink format
```
bcftools view --threads 16 -S ./phenotype_lists/list.hybrid.plum -O v ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.impute.vcf.gz | bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05' -O z -o ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.hybrid_all.impute.plum.vcf.gz


bcftools view --threads 16 -S ./phenotype_lists/list.hybrid.tail -O v ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.impute.vcf.gz | bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05' -O z -o ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.hybrid_all.impute.tail.vcf.gz
bcftools view --threads 16 -S ./phenotype_lists/list.full.plum -O v ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.impute.vcf.gz | bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05' -O z -o ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.impute.plum.vcf.gz
bcftools view --threads 16 -S ./phenotype_lists/list.full.tail -O v ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.impute.vcf.gz | bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05' -O z -o ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.impute.tail.vcf.gz
bcftools view --threads 16 -S ./phenotype_lists/list.male.plum -O v ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.male.impute.vcf.gz | bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05' -O z -o ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.male.impute.plum.vcf.gz
bcftools view --threads 16 -S ./phenotype_lists/list.male.tail -O v ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.male.impute.vcf.gz | bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05' -O z -o ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.male.impute.tail.vcf.gz
bcftools view --threads 16 -S ./phenotype_lists/list.female.plum -O v ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ+chrW.snps.miss02.maf05.female.impute.vcf.gz | bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05' -O z -o ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ+chrW.snps.miss02.maf05.female.impute.plum.vcf.gz
bcftools view --threads 16 -S ./phenotype_lists/list.female.tail -O v ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ+chrW.snps.miss02.maf05.female.impute.vcf.gz | bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05' -O z -o ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ+chrW.snps.miss02.maf05.female.impute.tail.vcf.gz
plink --vcf ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.hybrid_all.impute.plum.vcf.gz --make-bed --out ./input/gwas.hybrid.plum --allow-extra-chr
plink --vcf ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.hybrid_all.impute.tail.vcf.gz --make-bed --out ./input/gwas.hybrid.tail --allow-extra-chr
plink --vcf ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.impute.plum.vcf.gz --make-bed --out ./input/gwas.full.plum --allow-extra-chr
plink --vcf ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.impute.tail.vcf.gz --make-bed --out ./input/gwas.full.tail --allow-extra-chr
plink --vcf ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.male.impute.plum.vcf.gz --make-bed --out ./input/gwas.male.plum --allow-extra-chr
plink --vcf ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.male.impute.tail.vcf.gz --make-bed --out ./input/gwas.male.tail --allow-extra-chr
plink --vcf ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ+chrW.snps.miss02.maf05.female.impute.plum.vcf.gz --make-bed --out ./input/gwas.female.plum --allow-extra-chr
plink --vcf ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ+chrW.snps.miss02.maf05.female.impute.tail.vcf.gz --make-bed --out ./input/gwas.female.tail --allow-extra-chr






