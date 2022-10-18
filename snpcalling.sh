###snpcalling pipeline. You can contact yanglin(sheeplin0605@163.com) if you get any wrong.###

#software path, we here ensure that each software can run successfully in your terminal before runing yourself snpcalling pipeline, otherwise it will run failed.
fastp=/public/software/182_software/fastp-master/fastp
bwa=/usr/bin/bwa
samtools=/usr/bin/samtools
depth_coverage=/public/home/sheeplin/04.jajuar.vcf/new_cal_individual_depth.pl
picard=/public/software/182_software/picard.jar
gatk=/public/software/182_software/gatk-4.1.2.0/gatk
ParaFly=/usr/bin/ParaFly
vcftools=/public/software/182_software/vcftools/vcftools

#reference, here you need change to yourself genome reference path.
ref=/public/home/sheeplin/04.jajuar.vcf/cat.fa

#index genome reference
$bwa index $ref
$samtools faidx $ref
java -jar $picard CreateSequenceDictionary R=$ref O=$ref.dict

#get each individual's commands from rawdata(fq.gz) to g.vcf.gz.
for i in `ls *_1.fq.gz` 
do
i=${i/_1.fq.gz/}
cat <<EOF >$i.snpcalling.sh
$fastp -i ${i}_1.fq.gz -o ${i}_1.clean.fq.gz -I ${i}_2.fq.gz -O ${i}_2.clean.fq.gz
$bwa mem -t 16 -M -R '@RG\tID:$i\tSM:$i\tPL:Illumina' $ref  ${i}_1.clean.fq.gz ${i}_2.clean.fq.gz|samtools view -b -S - >$i.bam
rm ${i}_1.clean.fq.gz ${i}_2.clean.fq.gz 
$samtools sort -@ 5 $i.bam -o $i.sort.bam
rm $i.bam
$samtools flagstat -@ 5 $i.sort.bam >$i.flagstat.txt
$samtools depth $i.sort.bam >$i.depth.txt
$depth_coverage $i.depth.txt >$i.depth-coverage.txt #get each individual's sequencing statistic infomation
rm $i.depth.txt
java -jar $picard MarkDuplicates REMOVE_DUPLICATES=true I=$i.sort.bam O=$i.sort.dd.bam M=$i.sort.dd.metrics 2>$i.picard.log
rm $i.sort.bam
$samtools index -@ 5 $i.sort.dd.bam
$gatk HaplotypeCaller -R $ref -I $i.sort.dd.bam -ERC GVCF -O $i.g.vcf.gz 2>$i.HaplotypeCaller.log
EOF
done
#tips: if you want to speed the runnig of gatk HaplotypeCaller, you can run by splitting your ref genome into chromosomes with parameter "--intervals chr_ID". For example: gatk HaplotypeCaller -R cat.fa -I SRR11097154.sort.dd.bam -ERC GVCF --intervals CM001378.3 -O SRR11097154.CM001378.3.g.vcf.gz 2>SRR11097154.CM001378.3.HaplotypeCaller.log. After run all chr successfully, you need to merge all chr's g.vcf.gz like this: "gatk MergeVcfs -I SRR11097154.CM001378.3.g.vcf.gz -I SRR11097154.CM001379.3.g.vcf.gz -I SRR11097154.CM001380.3.g.vcf.gz -I ... -O SRR11097154.g.vcf.gz". 

#get run.snpcalling.list
for i in `ls *.snpcalling.sh`
do
echo "bash $i"
done >run.snpcalling.list

#parallel run run.snpcalling.list
$ParaFly -c run.snpcalling.list -CPU 10

#get gvcf merge command below, then run the gvcf merge command.
for i in `ls *.g.vcf.gz`; do echo "-V $i"; done |perl -e 'print"gatk CombineGVCFs -R \$ref ";while(<>){chomp;print"$_ ";}print"-O all_individual.merge.vcf.gz 2>gvcf.merge.log";'

#get gatk hard filter snp dataset
$gatk GenotypeGVCFs -R $ref -V all_individual.merge.vcf.gz --sample-ploidy 2 -O all_individual.raw.vcf.gz 2>GenotypeGVCFs.log
$gatk SelectVariants -select-type SNP -V all_individual.raw.vcf.gz -O all_individual.raw.snp.vcf.gz
$gatk SelectVariants -select-type INDEL -V all_individual.raw.vcf.gz -O all_individual.raw.indel.vcf.gz
$gatk VariantFiltration -V all_individual.raw.snp.vcf.gz  --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"  --filter-name "Filter"  -O all_individual.hardfilter.snp.vcf.gz 2>gatk.hardfilter.log

#vcftools filter and get clean snp dataset of population genetics 
$vcftools --gzvcf all_individual.hardfilter.snp.vcf.gz --max-missing 0.8 --maf 0.05 --minDP 4 --recode --remove-filtered-all --recode-INFO-all -c |gzip -c > all_individual.max-missing0.8maf0.05minDp4.vcf.gz
zcat all_individual.max-missing0.8maf0.05minDp4.vcf.gz |perl -e 'while(<>){print if(/^#/);if(/^chr.*\t\d+\t.\t\w\t\w\t/){print}}' |gzip -c >all_individual.final.vcf.gz

