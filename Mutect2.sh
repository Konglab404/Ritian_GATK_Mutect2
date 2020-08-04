#PBS -o ~/data/cmds/pbs_out/MuTect2.o
#PBS -e ~/data/cmds/pbs_out/MuTect2.e 
#PBS -N Mutect2
#PBS -q workq  
#PBS -l mem=150gb,walltime=25:00:00 
#PBS -l nodes=1:ppn=20 



for i in `seq -f '%04g' 0 19`
do
    gatk Mutect2 -R ~/data/BasicResource/Ref/Homo_sapiens_assembly38_gatk.fasta -I ~/data/Data/0_bwa_bam/${sample}.dedup.BQSR.sorted.bam -tumor ${sample} -L ~/data/BasicResource/GATK_bundle/hg38_scatter_interval_files/$i-scattered.interval_list --germline-resource ~/data/BasicResource/GATK_bundle/af-only-gnomad.hg38.vcf.gz --af-of-alleles-not-in-resource 0.00003125 --panel-of-normals ~/data/BasicResource/GATK_bundle/1000g_pon.hg38.vcf.gz -O ~/data/Data/0_bwa_bam/${sample}_somatic.$i.vcf.gz &
done
wait

gatk GatherVcfs -O ~/data/Data/1_vcfs/${sample}_somatic.vcf.gz \
    -I ~/data/Data/0_bwa_bam/${sample}_somatic.0000.vcf.gz \
    -I ~/data/Data/0_bwa_bam/${sample}_somatic.0001.vcf.gz \
    -I ~/data/Data/0_bwa_bam/${sample}_somatic.0002.vcf.gz \
    -I ~/data/Data/0_bwa_bam/${sample}_somatic.0003.vcf.gz \
    -I ~/data/Data/0_bwa_bam/${sample}_somatic.0004.vcf.gz \
    -I ~/data/Data/0_bwa_bam/${sample}_somatic.0005.vcf.gz \
    -I ~/data/Data/0_bwa_bam/${sample}_somatic.0006.vcf.gz \
    -I ~/data/Data/0_bwa_bam/${sample}_somatic.0007.vcf.gz \
    -I ~/data/Data/0_bwa_bam/${sample}_somatic.0008.vcf.gz \
    -I ~/data/Data/0_bwa_bam/${sample}_somatic.0009.vcf.gz \
    -I ~/data/Data/0_bwa_bam/${sample}_somatic.0010.vcf.gz \
    -I ~/data/Data/0_bwa_bam/${sample}_somatic.0011.vcf.gz \
    -I ~/data/Data/0_bwa_bam/${sample}_somatic.0012.vcf.gz \
    -I ~/data/Data/0_bwa_bam/${sample}_somatic.0013.vcf.gz \
    -I ~/data/Data/0_bwa_bam/${sample}_somatic.0014.vcf.gz \
    -I ~/data/Data/0_bwa_bam/${sample}_somatic.0015.vcf.gz \
    -I ~/data/Data/0_bwa_bam/${sample}_somatic.0016.vcf.gz \
    -I ~/data/Data/0_bwa_bam/${sample}_somatic.0017.vcf.gz \
    -I ~/data/Data/0_bwa_bam/${sample}_somatic.0018.vcf.gz \
    -I ~/data/Data/0_bwa_bam/${sample}_somatic.0019.vcf.gz \
    
gatk SortVcf -I ~/data/Data/1_vcfs/${sample}_somatic.vcf.gz \
    -O ~/data/Data/1_vcfs/${sample}_somatic.sorted.vcf.gz
    
rm ~/data/Data/0_bwa_bam/${sample}*00*



#GetPileupSummaries
gatk GetPileupSummaries \
    -L ./chromosome_intervals.intervals\
    -I ~/data/Data/0_bwa_bam/${sample}.dedup.BQSR.sorted.bam \
    -V ~/data/BasicResource/GATK_bundle/af-only-gnomad.hg38.SNP_biallelic.vcf.gz \
    -O ${sample}_somatic.pileups.table
    
#
gatk CalculateContamination \
    -I ${sample}_somatic.pileups.table \
    -O ${sample}_con.table &
    
gatk FilterMutectCalls \
    -R ~/data/BasicResource/Ref/Homo_sapiens_assembly38_gatk.fasta \
    -V ${sample}_somatic.vcf.gz \
    -contamination-table ${sample}_con.table \
    -O ${sample}_somatic.filtered.vcf.gz
    
bcftools filter -Oz -o ${sample}_somatic.filtered.PASS.vcf.gz -i 'FILTER=="PASS" & GERMQ>=30 & INFO/DP>=20' ${sample}_somatic.filtered.vcf.gz

