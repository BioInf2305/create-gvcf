###import the necessary modules
"""
sys --> for exiting gracefully
argparse ---> add the options to the script
pysam --> reading vcf file
pyfasta --> reding fasta file --> only necessary if the guidefile for vcf convesion do not have the top bottom annotation
                              --> needed to run only once for each new assembly
Bio.Seq --> reverse complementing should be only running when using pyfasta; check if that is not the case

"""
import sys
import argparse
from pysam import VariantFile
from pyfasta import Fasta
from Bio.Seq import Seq
from collections import OrderedDict

class VcfToPlink:

    def __init__(self,vcfIn,customFormat,fastaFile,excludeChrm,guideFile,addMetaInfo,outPrefix):

        '''
        TO DO --> after the script is successfully validated change the variable names --> expTopF, checkTop

        '''
        self.vcfIn=vcfIn
        self.customFormat=customFormat
        self.expTopF=guideFile
        self.outPrefix=outPrefix
        self.fastaF=fastaFile
        self.checkTop=addMetaInfo
        self.excludeChrm=excludeChrm

        '''
        define all the genotypes and associated terms in the dictionary.

        '''
        self.topBottomDict={"AMBI":["AT","TA","GC","CG"],"HOM":["AA","TT","GG","CC"],\
        "TOP":["AG","GA","AC","CA"],"BOT":["TG","GT","CT","TC"]}

        '''
        define all the variables to be used in the script

        '''
        self.excludeChrmList=[]

        self.updatedGuideList=[]

        self.customFile=False

        self.customDict={}


        #for the ped file

        self.sampleGenoDict=OrderedDict()
        self.destPed=open(self.outPrefix+".ped","w")

        #for the map file

        self.destMap=open(self.outPrefix+".map","w")
        self.mapDict=OrderedDict()

        ##for the tab delimited illumina AB allele format
        self.abAlleleOut=open(self.outPrefix+".ABallele.txt","w")
        self.abAlleleDict=OrderedDict()


'''

This function reads guideFile;
note that the guideFile by default contains only the top alleles at 4th and 5th column

'''


    def readGuideFile(self):
        markerIdList=[]
        self.populateExcludeChrmList()
        if self.customFormat!="NotSet":
            self.customFile=True
            self.readCustomFile()
        else:
            self.collectSamplesVcf()
        with open(self.expTopF) as source:
            for line in source:
                line=line.rstrip().split()
                if line[1] not in self.excludeChrmList and line[1]!="0":   # --> here it check by default if the chromosome id is "0"
                    if line[0] not in markerIdList:
                        markerIdList.append(line[0])
                    else:
                        print("ERROR "+line[0]+" marker id is present at least twice in the file")
                        sys.exit(1)
                    self.snpId=line[0]
                    self.chroms=line[1]
                    self.pos=int(line[2])
                    self.expTopGeno=line[3]+line[4]
                    self.ambiGeno=False
                    if len(line)!=12:
                        if self.checkTop!=0 and self.fastaF!="NotSet":
                            if self.expTopGeno in self.topBottomDict["AMBI"]:
                                self.ambiGeno=True
                            self.illuTopSeq3Prime=Seq(line[5][-60:])
                            self.illuTopSeq5Prime=Seq(line[5][:60])
                            strandInfo,refAllele,maxSimilarity=self.addMetaInfoTag()
                            line+=[strandInfo,refAllele,maxSimilarity]
                            self.updatedGuideList.append(line[:])
                        else:
                            print("meta information from the fasta file is not set, create the guidefile by running addMetaInfo tag and \
                            refrence file first")
                            sys.exit(1)
                    else:
                        if self.customFile==True:
                            self.illuStrand=line[6]
                            self.sourceStrand=line[7]
                            self.refStrand=line[8]
                            if self.refStrand not in "FR":
                                print ("ERROR reference strand is not defined at "+self.snpId)
                                sys.exit(1)
                            for sample in self.customDict:
                                self.sample=sample
                                if self.snpId in self.customDict[sample]:
                                    self.obsGeno=Seq(self.customDict[sample][self.snpId])
                                    self.determineGenoCustom()
                        if self.customFile==False:
                            self.readVcfAssignStrand()
        if self.checkTop!=0:
            self.updatedGuideFile()
        else:
            self.writePedMap()




    def populateExcludeChrmList(self):
        if self.excludeChrm!="NotSet":
            with open(self.excludeChrm) as source:
                for line in source:
                    line=line.rstrip().split()
                    self.excludeChrmList.append(line[0])



    def readCustomFile(self):
        with open(self.customFormat) as customFile:
            for line in customFile:
                line=line.rstrip().split()
                if line[0] not in self.customDict:
                    self.customDict[line[0]]=OrderedDict()
                    self.sampleGenoDict[line[0]]=[line[0],line[0],"0","0","0","0"]
                 '''
                TO DO --> generate AB-allele file in the custom format
                 '''
                self.customDict[line[0]][line[1]]=line[2]+line[3]


    def collectSamplesVcf(self):
        self.vcfInF=VariantFile(self.vcfIn)
        self.totalSamples=list(self.vcfInF.header.samples)
        for sample in self.totalSamples:
            self.sampleGenoDict[sample]=[sample,sample,"0","0","0","0"]
            self.abAlleleDict[sample]=[]


    def addMetaInfoTag(self):
        fastaF=Fasta(self.fastaF)
        neighboringSeq=fastaF[self.chroms][self.pos:self.pos+60]
        obsRefBase=Seq(fastaF[chromos][position])
        possibleSeqList=[str(self.illuTopSeq3Prime).upper(),str(self.illuTopSeq5Prime).upper(),\
        str(self.illuTopSeq3Prime.reverse_complement()).upper(),\
        str(self.illuTopSeq5Prime.reverse_complement()).upper()]
        maxIdentity=calcMaxIdentity(neighboringSeq.upper(),possibleSeqList)
        strand="NA"
        if self.ambiGeno:
            strand=determineStrandAmbi(self.fastaF,self.chroms,self.pos)
        return strand,obsRefBase,maxIdentity


    def updatedGuideFile(self):
        dest=open(self.expTopF+".updated","w")
        for record in self.updatedGuideList:
            dest.write("\t".join(record)+"\n")
        dest.close()

    def writePedMap(self):
        for sample in self.sampleGenoDict:
            numSnpPed=(len(self.sampleGenoDict[sample])-6)/2
            if numSnpPed!=len(self.mapDict):
                print("ERROR the number of SNPs in the map file do not match with ",\
                "that of inside the ped file",len(self.mapDict),numSnpPed)
                sys.exit(1)
            else:
                self.destPed.write(" ".join(self.sampleGenoDict[sample])+"\n")
        for snps in self.mapDict:
            self.destMap.write("\t".join(self.mapDict[snps]))
            self.destMap.write("\n")
        for sample in self.abAlleleDict:
            for sampleGeno in self.abAlleleDict[sample]:
                self.abAlleleOut.write(sample+"\t"+"\t".join(sampleGeno))
                self.abAlleleOut.write("\n")
        self.destPed.close()
        self.abAlleleOut.close()
        self.destMap.close()


    def determineGenoCustom(self):
        if self.illuStrand=="TOP":
            if self.refStrand=="R":
                self.obsGeno=str(self.obsGeno.reverse_complement())
        elif self.illuStrand=="BOT":
            if self.refStrand=="F":
                self.obsGeno=str(self.obsGeno.reverse_complement())
        else:pass
        if self.obsGeno in self.topBottomDict["HOM"]:
            if self.obsGeno[0] not in self.expTopGeno:
                print("ERROR "+"top alleles at "+self.snpId+" are not according to what is expected")
                sys.exit(1)
            else:
                self.alleleA=self.alleleB=self.obsGeno[0]
        elif self.obsGeno in self.topBottomDict["AMBI"]:
            if "A" in self.obsGeno:
                self.alleleA="A"
                self.alleleB="T"
            elif "C" in self.obsGeno:
                self.alleleA="C"
                self.alleleB="G"
            else:
                print("ERROR "+"top alleles at "+self.snpId+" are not according to what is expected")
                sys.exit(1)
        elif self.obsGeno in self.topBottomDict["TOP"]:
            if "A" not in self.obsGeno:
                print("ERROR "+"top alleles at "+self.snpId+" are not according to what is expected")
                sys.exit(1)
            else:
                self.alleleA="A"
                self.alleleB=self.obsGeno[abs(1-str(self.obsGeno).index("A"))]
        else:
                print("ERROR "+"top alleles at "+self.snpId+" are not according to what is expected")
                sys.exit(1)

        self.sampleGenoDict[self.sample].append(self.alleleA)
        self.sampleGenoDict[self.sample].append(self.alleleB)

        self.abAlleleOut.write(self.sample+"\t"+self.snpId+\
        "\t"+self.alleleA+"\t"+self.alleleB+"\n")

    def readVcfAssignStrand(self):
        for rec in self.vcfInF.fetch(self.chroms,self.pos-1,self.pos):
            self.recCopy=rec
            self.obsGeno="NN"
            self.alleleA="N"
            self.alleleB="N"
            self.ref=rec.ref
            truePos=1
            if len(rec.ref)==1 and rec.alts!=None:
                if len(rec.alts)==1 and len(rec.alts[0])==1 and \
                len(rec.alts)==1 and rec.alts[0]!="*":
                    self.obsGeno=Seq(rec.ref+rec.alts[0])
                    revObsGeno=Seq(rec.alts[0]+rec.ref)
                else:pass
            elif len(rec.ref)==1 and (rec.alts==None):
                self.obsGeno=Seq(rec.ref+rec.ref)
                revObsGeno=Seq(rec.ref+rec.ref)
            else:pass
            if self.obsGeno!="NN":
                obsGenoList=[str(self.obsGeno),str(revObsGeno),str(self.obsGeno.reverse_complement()),\
                str(revObsGeno.reverse_complement())]
                maxCommon=0
                for geno in obsGenoList:
                    countCommon=0
                    for allele in geno:
                        if allele in self.expTopGeno:
                            countCommon+=1
                    if countCommon>=maxCommon:
                        maxCommon=countCommon
                if maxCommon==2:
                    if self.checkTop!=0:
                        truePos=self.checkNeighboringSeq()
                    else:pass
                    if truePos!=1:
                        print("Warning the neighboring sequence not matching the top seq at skipping",\
                        str(self.pos),"at chrom ",self.chroms)
                    else:
                        self.determineStrand()
                        self.genoFromStrand()
                        self.writeGenoToDict()


    def determineStrand(self):
        if self.obsGeno in self.topBottomDict["HOM"]:
            if self.obsGeno[0] not in self.expTopGeno:
                self.strand="BOT"
            elif self.expTopGeno in self.topBottomDict["AMBI"]:
                self.strand=determineStrandAmbi(self.fastaF,self.chroms,self.pos)
            else:
                self.strand="TOP"
        elif self.obsGeno in self.topBottomDict["AMBI"]:
            self.strand=determineStrandAmbi(self.fastaF,self.chroms,self.pos)
        elif self.obsGeno in self.topBottomDict["TOP"]:
            self.strand="TOP"
        elif self.obsGeno in self.topBottomDict["BOT"]:
            self.strand="BOT"
        else:
            print("ERROR unexpected genotype at ",self.chroms,self.pos)
            sys.exit(1)
        return self.strand

    def genoFromStrand(self):
        if self.strand=="NA":
            print("ERROR: strand not assigned at "+str(self.chroms)+"\t"+str(self.pos))
            sys.exit(0)
        else:
            if self.strand=="TOP":
                if self.obsGeno in self.topBottomDict["HOM"]:
                    self.alleleB=self.alleleA=self.obsGeno[0]
                elif self.obsGeno in self.topBottomDict["AMBI"]:
                    self.alleleA,self.alleleB=determineAlleleAB(self.obsGeno,1)
                elif "A" not in self.obsGeno:
                        print ("ERROR base A not present in top strand for non-ambiguous genotypes at "+\
                        str(self.chroms)+"\t"+str(self.pos))
                        sys.exit(0)
                else:self.alleleA,self.alleleB=determineAlleleAB(self.obsGeno,0)
                self.refAllele=self.obsGeno[0]
                self.altAllele=self.obsGeno[1]
            elif self.strand=="BOT":
                self.obsGeno=str(self.obsGeno.reverse_complement())
                if self.obsGeno in self.topBottomDict["HOM"]:
                    self.alleleB=self.alleleA=self.obsGeno[0]
                elif "A" not in self.obsGeno and "C" not in self.obsGeno:
                    print ("BOT ERROR base A or C not present in top strand for non-ambiguous genotypes at "+\
                        str(self.chroms)+"\t"+str(self.pos))
                    sys.exit(1)
                else:
                    self.alleleA,self.alleleB=determineAlleleAB(self.obsGeno,0)
                self.refAllele=str(self.obsGeno[1])
                self.altAllele=str(self.obsGeno[0])
            else:
                print("ERROR: unexpected strand "+self.strand+" assigned at "+\
                str(self.chroms)+"\t"+str(self.pos))
                sys.exit(0)

    def writeGenoToDict(self):
        if (self.refAllele=="N" or self.altAllele=="N") or (self.alleleA=="N" or self.alleleB=="N"):
            print("ERROR: allele not assigned at "+str(self.chroms)+"\t"+str(self.pos)+" "+self.obsGeno+" "+self.strand)
            sys.exit(0)
        else:
            for sample in self.totalSamples:
                if self.recCopy.samples[sample]["GT"]==(0,0):
                    self.sampleGenoDict[sample].append(self.refAllele)
                    self.alleleA=self.AlleleB=self.refAllele
                    self.sampleGenoDict[sample].append(self.refAllele)
                elif self.recCopy.samples[sample]["GT"]==(0,1):
                    self.sampleGenoDict[sample].append(self.refAllele)
                    self.sampleGenoDict[sample].append(self.altAllele)
                elif self.recCopy.samples[sample]["GT"]==(1,1):
                    self.sampleGenoDict[sample].append(self.altAllele)
                    self.sampleGenoDict[sample].append(self.altAllele)
                    self.alleleA=self.alleleB=self.altAllele
                else:
                    self.sampleGenoDict[sample].append("0")
                    self.sampleGenoDict[sample].append("0")
                    self.alleleA=self.alleleB="0"
                self.abAlleleDict[sample].append([self.snpId,self.alleleA,self.alleleB])
            if self.snpId not in self.mapDict:
                self.mapDict[self.snpId]=[self.chroms,self.snpId,"0",str(self.pos)]
            else:
                print("ERROR duplicate SNP id ",self.snpId)
                sys.exit(1)

def determineStrandAmbi(fastaIn,chromos,position):
    position=position-1
    strand="NA"
    fastaF=Fasta(fastaIn)
    strandAmbiDict={"AG":"TOP","GA":"BOT","AC":"TOP","CA":"BOT","TG":"TOP","GT":"BOT","TC":"TOP","CT":"BOT"}
    movBase=0
    obsRefBase=Seq(fastaF[chromos][position])
    while strand=="NA":
        movBase+=1
        prime5Base=fastaF[chromos][position-movBase].upper()
        prime3Base=fastaF[chromos][position+movBase].upper()
        newGeno=prime5Base+prime3Base
        if newGeno in strandAmbiDict:
            strand=strandAmbiDict[newGeno]
    return strand


def calcMaxIdentity(seqCompare,seqList):
    maxIdentity=0
    for seqP in seqList:
        maxIdentity1=0
        totalN=seqP.count("N")+seqCompare.count("N")
        for i in range(len(seqCompare)):
            if seqCompare[i]==seqP[i]:
                maxIdentity1+=1
        if maxIdentity1>maxIdentity:
            maxIdentity=maxIdentity1
    return maxIdentity


def determineAlleleAB(geno,isAmbi):
    ambiGeno=["AT","TA","GC","CG"]
    alleleA="N"
    alleleB="N"
    if isAmbi==1:
        if geno in ambiGeno:
            if "A" in geno:
                alleleA="A"
                alleleB="T"
            else:
                alleleA="C"
                alleleB="G"
    else:
        if "A" in geno:
            alleleA="A"
        elif "C" in geno:
            alleleA="C"
        else:
            geno="NN"
        alleleB=geno[abs(1-str(geno).index(alleleA))]
    return alleleA,alleleB

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="This python script will convert vcf or custom format to top alleles in custom and ped format",\
    epilog="author: Maulik Upadhyay (Upadhyaya.maulik@gmail.com)")
    parser.add_argument('-v',"--vcfFile",metavar="File",help="bam file",default="NotSet",required=False)
    parser.add_argument('-d',"--customF",metavar="File",help="custom input",default="NotSet",required=False)
    parser.add_argument('-f',"--fastaFile",metavar="File",help="reference fasta file",default="NotSet",required=False)
    parser.add_argument('-e'."--excludeChrm",metavar="File",help="file contaning chromosomes to exclude",default="NotSet",required=False)
    parser.add_argument('-g',"--guideFile",metavar="File",help="expected top alleles file",required=True)
    parser.add_argument('-a',"--addMetaInfo",metavar="Int",help="add meta information from guide file from ref fasta file",default=0,required=False)
    parser.add_argument('-o',"--outFile",metavar="Str",help="prefix of the output files",required=True)
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    elif args.vcfFile=="NotSet" and args.customF=="NotSet":
        print("either vcf or custom formatted file is required")
        parser.print_help(sys.stderr)
        sys.exit(1)
    elif args.addMetaInfo!=0 and args.fastaFile=="NotSet":
        print("when setting the metaInfo flag, it is necessary to provide the reference sequence")
        sys.exit(1)
    else:
        VcfToPlinkObj=VcfToPlink(args.vcfFile,args.customF,args.fastaFile,args.excludeChrm,args.guideFile,args.addMetaInfo,args.outFile)
        VcfToPlinkObj.readGuideFile()
