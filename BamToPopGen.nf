#!/usr/bin/env nextflow

nextflow.enable.dsl=2



sampleGvcfs_outDir = params.sampleGvcfs_outDir
bamFilesPath       = params.bamFilesPath
databasePath       = params.databasePath
vcfFilesPath       = params.vcfFilesPath
refFile            = params.refFile
chrmNameF          = params.chrmNameF
rawVcf_outDir      = params.rawVcf_outDir
filtVcf_outDir     = params.filtVcf_outDir
phased_outDir      = params.phased_outDir
tmp_dir            = params.tmp_dir
addGvcf            = params.addGvcf
sampleMap          = params.sampleMap
snpCall            = params.snpCall
indelCall          = params.indelCall
mixedCall          = params.mixedCall
phasedSnps         = params.phasedSnps
batchSize          = params.batch_size
sampleNameF        = params.sampleNameF
qd                 = params.qd
qual               = params.qual
sor                = params.sor
fs                 = params.fs
mq                 = params.mq
mqRankSum          = params.mqRankSum
readSum            = params.readSum



Channel
    .fromFilePairs( bamFilesPath ){ file -> file.name.replaceAll(/.bam|.bai$/,'') }
    .set { bamInfo }

channel 
    .fromFilePairs( vcfFilesPath )
    //.map { vcfFileKeys,files -> tuple(vcfFileKeys.replaceAll(/.filt.step3/,''),files)}
    .set { vcfFiles }

channel
    .fromPath( chrmNameF )
    .splitText( by:1 )
    .map{it -> it.trim()}
    .set { chrm }


channel
    .fromFilePairs( refFile, size:-1 )
    .set{ ref }


chrmBam=chrm.combine(bamInfo)
chrmBamRef=chrmBam.combine(ref)
chrmRef=chrm.combine(ref)            

include { BAMTODB } from "${baseDir}/modules/bamTodb" addParams(
        OUTPUT1: sampleGvcfs_outDir,
        databasePath: databasePath,
        tmp_dir: tmp_dir,
        addGvcf: addGvcf,
        batchSize: batchSize,
        sampleMap: sampleMap)

include { DBTOGENO } from "${baseDir}/modules/dbTogeno" addParams(
        OUTPUT2: rawVcf_outDir, 
        databasePath: databasePath,
        OUTPUT3: filtVcf_outDir,
        snpCall: snpCall,
        qd: qd,
        qual: qual,
        sor: sor,
        fs: fs,
        mq: mq,
        mqRankSum: mqRankSum,
        readSum: readSum
        )


include { GENOTOPHASE } from "${baseDir}/modules/genoTophase" addParams(OUTPUT4: phased_outDir)

workflow {
    if (addGvcf!="NA"){
    chrmComplete=BAMTODB(chrmBamRef,ref,chrm)
    }
    if (snpCall=="Y"){
    filteredVcf=DBTOGENO( addGvcf=="NA" ? chrmRef:chrmComplete.refChromOut)
    }
    if(phasedSnps=="Y"){
     trial=GENOTOPHASE( snpCall=="Y" ? filteredVcf.outF:vcfFiles)
    }
}
