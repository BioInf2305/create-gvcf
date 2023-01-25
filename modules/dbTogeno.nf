#!/usr/bin/env nextflow

process variantCalling {
   publishDir(params.OUTPUT2, pattern: '*.raw.{vcf.gz,vcf.gz.tbi}')
   tag { "${chrm}_raw" }
   label "oneCpu"

   input:        
    tuple val(chrm),val(referencePrefix),file(reference)

   output:        
    tuple val(chrm), file ("${chrm}.raw.vcf.gz"), file("${chrm}.raw.vcf.gz.tbi"),val(referencePrefix),file(reference)

   script:     
    def gatkDatabase=params.databasePath
   """
   mkdir /home/maulik/data/Shared/Maulik/nextflowScripts/tmp/${chrm}_var

/data/headnode01/medugorac/Shared/tools/gatk-4.1.8.1/gatk --java-options "-Xmx${task.memory.toGiga()}G " GenotypeGVCFs -R ${referencePrefix}.fna -V gendb:///${gatkDatabase}/${chrm} -O ${chrm}.raw.vcf.gz --tmp-dir /home/maulik/data/Shared/Maulik/nextflowScripts/tmp/${chrm}_var

  
  rm -r /home/maulik/data/Shared/Maulik/nextflowScripts/tmp/${chrm}_var


   """
}


process extractSamples {
   publishDir(params.OUTPUT3, pattern: '*.raw.selectedSamples.{vcf.gz,vcf.gz.tbi}',mode:"copy")
   tag { "${chrm}_selectSample" }
   label "oneCpu"

   input:        
    tuple val(chrm),file(variantF),file(variantIndex),val(referencePrefix),file(reference)

   output:        
    tuple val(chrm), file ("${chrm}.raw.selectedSamples.vcf.gz"), file("${chrm}.raw.selectedSamples.vcf.gz.tbi"),val(referencePrefix),file(reference)

   script:     
    def sampleNameF=params.sampleNameF

   """
   mkdir /home/maulik/data/Shared/Maulik/nextflowScripts/tmp/${chrm}_sampleSelection


/data/headnode01/medugorac/Shared/tools/gatk-4.1.8.1/gatk --java-options "-Xmx${task.memory.toGiga()}G " SelectVariants -V ${variantF} --sample-name ${sampleNameF} -O ${chrm}.raw.selectedSamples.vcf.gz --tmp-dir /home/maulik/data/Shared/Maulik/nextflowScripts/tmp/${chrm}_sampleSelection
  
  rm -r /home/maulik/data/Shared/Maulik/nextflowScripts/tmp/${chrm}_sampleSelection


   """
}



process snpCallingFilter{
    publishDir(params.OUTPUT3, pattern: '*.snp.filt.{vcf.gz,vcf.gz.tbi}', mode:"copy")
    tag { "${chrm}_snp_filtered" }
    label "oneCpu"

    input:
        tuple val(chrm),file(rawSnpF),file(rawSnpIndex),val(referencePrefix),file(reference)
    output:
        tuple val(chrm),file("${chrm}.snp.filt.vcf.gz"),file("${chrm}.snp.filt.vcf.gz.tbi")

    script:
        def varQd        = params.qd
        def varQual      = params.qual
        def varSor       = params.sor
        def varFs        = params.fs
        def varMq        = params.mq
        def varMqRankSum = params.mqRankSum
        def varReadSum   = params.readSum

        """

        mkdir /home/maulik/data/Shared/Maulik/nextflowScripts/tmp/${chrm}_snp

/data/headnode01/medugorac/Shared/tools/gatk-4.1.8.1/gatk --java-options "-Xmx${task.memory.toGiga()}G " SelectVariants -V ${rawSnpF} -select-type SNP -O ${chrm}.snp.filt.1.vcf.gz --tmp-dir /home/maulik/data/Shared/Maulik/nextflowScripts/tmp/${chrm}_snp

        
/data/headnode01/medugorac/Shared/tools/gatk-4.1.8.1/gatk --java-options "-Xmx${task.memory.toGiga()}G " VariantFiltration -V ${chrm}.snp.filt.1.vcf.gz -filter "QD < ${varQd}.0" --filter-name "QD${varQd}" -filter "QUAL < ${varQual}.0" --filter-name "QUAL${varQual}" -filter "SOR > ${varSor}.0" --filter-name "SOR${varSor}" -filter "FS > ${varFs}.0"  --filter-name "FS${varFs}" -filter "MQ < ${varMq}.0" --filter-name "MQ${varMq}" -filter "MQRankSum < -${varMqRankSum}" --filter-name "MQRankSum-${varMqRankSum}" -filter "ReadPosRankSum < -${varReadSum}.0" --filter-name "ReadPosRankSum-${varReadSum}" -O ${chrm}.filt.2.vcf.gz --tmp-dir /home/maulik/data/Shared/Maulik/nextflowScripts/tmp/${chrm}_snp


/data/headnode01/medugorac/Shared/tools/gatk-4.1.8.1/gatk --java-options "-Xmx${task.memory.toGiga()}G " SelectVariants -R ${referencePrefix}.fna -V ${chrm}.filt.2.vcf.gz --exclude-filtered --restrict-alleles-to BIALLELIC -O ${chrm}.snp.filt.vcf.gz --tmp-dir /home/maulik/data/Shared/Maulik/nextflowScripts/tmp/${chrm}_snp

        rm -r /home/maulik/data/Shared/Maulik/nextflowScripts/tmp/${chrm}_snp


        """

}


workflow DBTOGENO {
    take:
        chrmRef
    main:
        variantFile=variantCalling(chrmRef)
        if (params.sampleNameF!="NA"){
            sampleVariantFile=extractSamples(variantFile)
          }
        finalSnpFile=snpCallingFilter( params.sampleNameF!="NA" ? sampleVariantFile: variantFile)
    emit:
        outF=finalSnpFile
}
