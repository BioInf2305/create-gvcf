#!/usr/bin/env nextflow



process runGvcfParallel {
   tag { "${chrm}_${sample}" }
   label "twoCpus"

   input:        
    tuple val(chrm),val(sample),file(bam),val(referencePrefix),file(reference)

   output:        
    tuple val(sample), file ("${chrm}_${sample}.g.vcf")

   script:     
    def tmp_dirct = params.tmp_dir
   """
      mkdir ${tmp_dirct}/${chrm}_${sample}

     /data/headnode01/medugorac/Shared/tools/gatk-4.1.8.1/gatk --java-options "-Xmx${task.memory.toGiga()}G " HaplotypeCaller -I ${sample}.bam -O ${chrm}_${sample}.g.vcf -R ${referencePrefix}.fna --emit-ref-confidence GVCF -pairHMM FASTEST_AVAILABLE --native-pair-hmm-threads ${task.cpus} -L ${chrm} --tmp-dir ${tmp_dirct}/${chrm}_${sample}/

     rm -r ${tmp_dirct}/${chrm}_${sample}/


   """
}

process mergeGvcf {
 publishDir(params.OUTPUT1, pattern:"*.g.{vcf,vcf.idx}",mode:"copy")
 tag { "${sample}" }
 label "oneCpu"

 input:
    tuple val(sample), file(gvcfFiles),val(referencePrefix),file(reference)

 output:
     path "${sample}_merged.g.vcf",     emit: vcf_file
     path "${sample}_merged.g.vcf.idx", emit: vcfIdx_file
 
 script:
    gvcfUpdated=gvcfFiles.sort { it.name.tokenize('_')[1].tokenize('.')[0].toInteger() }.join(" -I ")
    def tmp_dirct = params.tmp_dir
    """
    mkdir ${tmp_dirct}/${sample};

   /data/headnode01/medugorac/Shared/tools/gatk-4.1.8.1/gatk --java-options "-Xmx${task.memory.toGiga()}G " GatherVcfs -R ${referencePrefix}.fna -I ${gvcfUpdated} -O ${sample}_merged.g.vcf --TMP_DIR ${tmp_dirct}/${sample}/ ;

   rm -r ${tmp_dirct}/${sample}/

    """
}


process bamTogatkDatabase {

 tag { "${chrom}" }
 label "oneCpu"

 input:
    file(mergedGvcfs)
    file(vcfIdx)
    val(chrom)
 
 output:
    //tuple val(chrom),file("${chrom}.log"),emit:log_tuple
    val(chrom)

 script:
    mergedGvcfsUpdated=mergedGvcfs.join(" -V ")
    def gatkFolder=params.databasePath
    def addGvcf_v = params.addGvcf
    def tmp_dirct = params.tmp_dir
    def batchSize = params.batchSize
    if (addGvcf_v=="N")
    """

    mkdir ${tmp_dirct}/${chrom}

    /data/headnode01/medugorac/Shared/tools/gatk-4.1.8.1/gatk --java-options "-Xmx${task.memory.toGiga()}G -Xms${task.memory.toGiga()}G" GenomicsDBImport --genomicsdb-workspace-path ${gatkFolder}/${chrom} --batch-size ${batchSize} -L ${chrom} -V ${mergedGvcfsUpdated} --tmp-dir ${tmp_dirct}/${chrom}/

    rm -r ${tmp_dirct}/${chrom}/

    """
    else

    """

    mkdir ${tmp_dirct}/${chrom};

     /data/headnode01/medugorac/Shared/tools/gatk-4.1.8.1/gatk --java-options "-Xmx${task.memory.toGiga()}G -Xms${task.memory.toGiga()}G" GenomicsDBImport --genomicsdb-update-workspace-path ${gatkFolder}/${chrom} --batch-size ${batchSize} -L ${chrom} -V ${mergedGvcfsUpdated} --tmp-dir ${tmp_dirct}/${chrom}/ ;

     rm -r ${tmp_dirct}/${chrom}/

    """
}

process gvcfTogatkDatabase {

 tag { "${chrom}" }
 label "oneCpu"

 input:
    val(chrom)
 
 output:
    val(chrom)

 script:
    def gatkFolder=params.databasePath
    def addGvcf_v = params.addGvcf
    def tmp_dirct = params.tmp_dir
    def batchSize = params.batchSize
    def sampleMap = params.sampleMap
    if (addGvcf_v=="N")
    """

    mkdir ${tmp_dirct}/${chrom};

    /data/headnode01/medugorac/Shared/tools/gatk-4.1.8.1/gatk --java-options "-Xmx${task.memory.toGiga()}G -Xms${task.memory.toGiga()}G" GenomicsDBImport --genomicsdb-workspace-path ${gatkFolder}/${chrom} --batch-size ${batchSize} -L ${chrom} --sample-name-map ${sampleMap} --tmp-dir ${tmp_dirct}/${chrom}/ ;

    rm -r ${tmp_dirct}/${chrom}/

    """
    else

    """

    mkdir ${tmp_dirct}/${chrom};

     /data/headnode01/medugorac/Shared/tools/gatk-4.1.8.1/gatk --java-options "-Xmx${task.memory.toGiga()}G -Xms${task.memory.toGiga()}G" GenomicsDBImport --genomicsdb-update-workspace-path ${gatkFolder}/${chrom} --batch-size ${batchSize} -L ${chrom} --sample-name-map ${sampleMap} ;

     rm -r ${tmp_dirct}/${chrom}/

    """
}

workflow BAMTODB {
    take:
        chrmBamRef
        ref
        chrm
    main:
        if (params.sampleMap=="NA"){
        result = runGvcfParallel(chrmBamRef)
        resultGroupped=result.groupTuple()
        resultRefCombined=resultGroupped.combine(ref)
        mergeGvcf(resultRefCombined)
        chromComplete=bamTogatkDatabase(mergeGvcf.out.vcf_file.collect(),mergeGvcf.out.vcfIdx_file.collect(),chrm)
        }
        else{
        chromComplete=gvcfTogatkDatabase(chrm)
        }
        refChromComplete=chromComplete.combine(ref)

   emit:
        refChromOut=refChromComplete
}
