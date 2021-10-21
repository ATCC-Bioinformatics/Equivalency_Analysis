#!/usr/bin/env python

import sys
import operator
import os
import re
import subprocess
import glob
from collections import defaultdict
import time
import csv
import shutil

#Calculate coverage at VCF sites across BAM file
def calCov(inputVCF,bamFile):
    vcfFile=open(inputVCF,"r")
    passN=0
    coverageV=0
    totalBase=0
    CovList={}
    for inline in vcfFile:
        if not inline.startswith("#"):
            infor=inline.strip().split()
            if infor[6]=="PASS":
                passN+=1
                endpos=len(infor[4])+int(infor[1])-1
                totalBase+=len(infor[4])
                regionInfor=infor[0]+":"+infor[1]+"-"+str(endpos)
                cmd="samtools depth -a -r \""+regionInfor+"\" "+bamFile
                results=subprocess.check_output([cmd],shell=True).strip().decode("utf-8").split("\n")
                totalRegion=0
                regionCov=0
                if len(results)>0:
                    for eachItem in results:
                        numberTmp=eachItem.split("\t")
                        if len(numberTmp)>0:
                            covTmp=numberTmp[-1].strip()
                            if len(covTmp)==0:
                                covTmp=0.0
                        else:
                            covTmp=0.0
                        totalRegion+=float(covTmp)
                        coverageV+=float(covTmp)
                    regionCov=float(totalRegion)/float(len(infor[4]))
                CovList["\t".join(infor[:7])]=str(regionCov)
    vcfFile.close()
    return (coverageV,passN,totalBase,CovList)

#Reformat GFF file to vep-compatible gff
def GFFtovepGFF(inputFile,outputFolder):
    handle=open(inputFile,"rU")
    parnetID="None"
    outputFileName=outputFolder+"/ref.gff3"
    output=open(outputFileName,"w")
    geneCount=1
    
    for inline in handle:
        if not inline.startswith("#"):
            infor=inline.strip().split("\t")
            if len(infor)<3:
                continue
            if infor[2]=="gene":
                originalID="ID=gene"+str(geneCount)
                parnetID=originalID.replace("ID=","")
                AllAnnList=infor[-1].split(";")
                RestAnn=None
                for eachEle in AllAnnList:
                    if eachEle.startswith("ID="):
                        originalID=eachEle
                        parnetID=eachEle.replace("ID=","")
                        continue
                    else:
                        if RestAnn:
                            RestAnn+=";"+eachEle
                        else:
                            RestAnn=eachEle
                print (inline.strip(),file=output)
                infor[2]="transcript"
                if "ID=gene" in originalID: 
                    transcriptID=originalID.replace("gene","transcript")
                    exonID=originalID.replace("gene","exon")
                else:
                    transcriptID="ID=transcript-"+parnetID
                    exonID="ID=exon-"+parnetID
                print ("\t".join(infor[:-1])+"\t"+transcriptID+";"+RestAnn+";Parent="+parnetID+";biotype=protein_coding",file=output)
                infor[2]="exon"
                parnetID=transcriptID.replace("ID=","")
                print ("\t".join(infor[:-1])+"\t"+exonID+";"+RestAnn+";Parent="+parnetID,file=output)
            elif infor[2]=="CDS" or infor[2]=="AMR":
                AllAnnList=infor[-1].split(";")
                for eachindex,eachEle in enumerate(AllAnnList):
                    if eachEle.startswith("Parent="):
                        AllAnnList[eachindex]="Parent="+parnetID
                print ("\t".join(infor[:-1])+"\t"+";".join(AllAnnList),file=output)
            else:
                print (inline.strip(),file=output)
        else:
            print (inline.strip(),file=output)
                
    handle.close()
    output.close()
    return (outputFileName)

#Align reads to reference with BWA
def runbwa(refGenome,readsFileList,outputFolderPrefix,threadN,singleton):
    samList=[]
    samOutput=outputFolderPrefix+".sam"
    samOutputS=outputFolderPrefix+".singleton.sam"
    indexRef="bwa index "+refGenome
    bwaCmd="bwa mem "+refGenome+" "+(" ").join(readsFileList)+" -o "+samOutput+" -t "+str(threadN)+" 2> "+outputFolderPrefix+".bwa.log"
    subprocess.call([indexRef],shell=True)
    subprocess.call([bwaCmd],shell=True)
    samList.append(samOutput)
    if singleton:
        bwaCmd="bwa mem "+refGenome+" "+singleton+" -o "+samOutputS+" -t "+str(threadN)+" 2> "+outputFolderPrefix+".bwa.singleton.log"
        subprocess.call([bwaCmd],shell=True)
        samList.append(samOutputS)
    return (samList)

#Use GATK to generate VCF from BAM files according to ploidy, adding read groups and marking duplicates
def ProcessMapping(samFileList,refGenome,genomeType,vizFile):
    newBamList=[]
    for eachSam in samFileList:
        samTobam="samtools view -h -b -S "+eachSam+" > "+eachSam.replace(".sam",".bam")
        subprocess.call([samTobam],shell=True)
        newBamList.append(eachSam.replace(".sam",".bam"))
    samFile=samFileList[0]
    if len(newBamList)>1:
        samFile=samFileList[0].replace(".sam",".all.sam")
        combineBam="samtools merge "+samFile.replace(".sam",".bam")+" "+(" ").join(newBamList)
        subprocess.call([combineBam],shell=True)

    sortCmd="samtools sort "+samFile.replace(".sam",".bam")+" > "+samFile.replace(".sam",".sort.bam")
    indexBam="samtools index "+samFile.replace(".sam",".sort.bam")

    subprocess.call([sortCmd],shell=True)
    subprocess.call([indexBam],shell=True)
    vizFile.append(samFile.replace(".sam",".sort.bam"))
    vizFile.append(samFile.replace(".sam",".sort.bam.bai"))
    refGindex=".".join(refGenome.split(".")[:-1])+".dict"
    if not os.path.isfile(refGindex):
        gatkIndext="gatk CreateSequenceDictionary --QUIET true -R "+refGenome
        subprocess.call([gatkIndext],shell=True)
    if not os.path.isfile(refGenome+".fai"):
        samIndex="samtools faidx "+refGenome
        subprocess.call([samIndex],shell=True)

    bamFile=samFile.replace(".sam",".sort.bam")
    bamRGfile=bamFile.replace(".sort.bam",".rg.sort.bam")
    finalBam=bamRGfile.replace(".rg.sort.bam",".final.bam")
    rawVCF=finalBam.replace(".bam",".vcf")
    matrixDup=finalBam.replace(".final.bam",".duplicates.metrix")
    AddOrReplaceReadGroups="gatk AddOrReplaceReadGroups -I "+bamFile+" -O "+bamRGfile+" -LB library -PL illumina -PU unit1 -SM Sample1 --QUIET true"
    MarkDuplicates="gatk MarkDuplicates -I "+bamRGfile+" -O "+finalBam+" -M "+matrixDup+" --QUIET true"
    indexBam="samtools index "+finalBam
    if genomeType=="haploid":
        variantCall="gatk HaplotypeCaller --QUIET true -stand-call-conf 30 -ploidy 1 -I "+finalBam+" -O "+rawVCF+" -R "+refGenome
    elif genomeType=="diploid":
        variantCall="gatk HaplotypeCaller --QUIET true -stand-call-conf 30 -ploidy 2 -I "+finalBam+" -O "+rawVCF+" -R "+refGenome
    subprocess.call([AddOrReplaceReadGroups],shell=True)
    subprocess.call([MarkDuplicates],shell=True)
    subprocess.call([indexBam],shell=True)
    subprocess.call([variantCall],shell=True)
    return (rawVCF,vizFile,finalBam)

#Filter VCF file for assorted minimum criteria - different criteria apply to SNPs and InDels
def FilterVCF(rawVCF,referenceFa):
    snpVCF=rawVCF.replace(".vcf",".SNP.vcf")
    indelVCF=rawVCF.replace(".vcf",".indel.vcf")
    filterSNPvcf=rawVCF.replace(".vcf",".filter.SNP.vcf")
    filterIndelvcf=rawVCF.replace(".vcf",".filter.indel.vcf")
    ExtractSNP="gatk SelectVariants -V "+rawVCF+" -O "+snpVCF+" -R "+referenceFa+" -select-type SNP "
    ExtractIndel="gatk SelectVariants -V "+rawVCF+" -O "+indelVCF+" -R "+referenceFa+" -select-type INDEL "
    filterSNP="gatk VariantFiltration -V "+snpVCF+" -R "+referenceFa+" -O "+filterSNPvcf+" --QUIET true -filter \" QD < 2.0 || FS > 60.0 || MQ < 40.0 \" --filter-name \"SNPfilter-failed\""
    filterINDEL="gatk VariantFiltration -V "+indelVCF+" -R "+referenceFa+" -O "+filterIndelvcf+" --QUIET true  -filter \" QD < 2.0 || FS > 200.0  \" --filter-name \"INDELfilter-failed\""

    subprocess.call([ExtractSNP],shell=True)
    subprocess.call([ExtractIndel],shell=True)
    subprocess.call([filterSNP],shell=True)
    subprocess.call([filterINDEL],shell=True)
    return (filterSNPvcf,filterIndelvcf)

#Use VEP to generate *.variant.ann file with functional impact of genomic variation
def annVCF(vcf,gffFile,refGenome):
    gffGZ=gffFile+".gz.tbi"
    if not os.path.isfile(gffGZ): 
        cmd="grep -v \"#\" "+gffFile+" | sort -k1,1 -k4,4n -k5,5n -t$'\t'|sed '/^$/d'|bgzip -c >"+gffFile+".gz"
        subprocess.call([cmd],shell=True)
        tabixCMD="tabix -p gff "+gffFile+".gz"
        subprocess.call([tabixCMD],shell=True)
    vepCMD="vep -i "+vcf+" -gff "+gffFile+".gz -fasta "+refGenome+" -o "+vcf.replace(".vcf",".variant.ann")+" --format vcf"
    print("#### VEP command ####")
    print(vepCMD)
    subprocess.call([vepCMD],shell=True)
    return (vcf.replace(".vcf",".variant.ann"))

#Collate and process vcf and variant.ann file for functional impact of genomic variation
def VarAnnCsv(vcfFile,variantFile):
    handle=open(vcfFile,"rU")
    vcfDict={}
    vcfList=[]
    vcfAnnImport=defaultdict(list)
    vcfAnnSupp=defaultdict(list)
    annLine=None
    for inline in handle:
        if not inline.startswith("#"):
            infor=inline.strip().split("\t")
            vcfDict[infor[0]+":"+infor[1]]="\t".join(infor[:7])
            vcfList.append(infor[0]+":"+infor[1])
        elif inline.startswith("#CHROM"):
            infor=inline.strip().split("\t")
            annLine="\t".join(infor[:7])
    handle.close()
    if variantFile:
        handle=open(variantFile,"rU")
        maxLen=0
        for inline in handle:
            if inline.startswith("#Uploaded_variation"):
                infor=inline.strip().split("\t")
                annLine+="\t"+"\t".join(infor[3:])
            elif not inline.startswith("#"):
                infor=inline.strip().split("\t")
                AddInfor=infor[-1].split(";")
                
                if ("upstream" not in infor[6]) and ("downstream" not in infor[6]):
                    vcfAnnImport[infor[1]].append("\t".join(infor[3:]))
                else:
                    for eachEle in AddInfor:
                        if eachEle.startswith("DISTANCE="):
                            distance=int(eachEle.split("=")[-1])
                            if distance<=100:
                                vcfAnnSupp[infor[1]].append("\t".join(infor[3:]))
                                break
        handle.close()
    
    return (vcfDict,vcfList,vcfAnnImport,vcfAnnSupp,annLine)

#Move files to folder for visualization with broad institute
def vizFolder(vizFileList,outputFolder):
    allFiles=os.listdir(outputFolder)
    movingList=""
    for eachFile in allFiles:
        if (outputFolder+"/"+eachFile) not in vizFileList:
            movingList+=" "+outputFolder+"/"+eachFile
    tmpFolder=outputFolder+"/tmp/"
    igvFolder=outputFolder+"/igv/"
    if os.path.exists(tmpFolder):
        timestamp=str(int(time.time()))
        tmpFolder=outputFolder+"/tmp-"+timestamp+"/"
        igvFolder=outputFolder+"/igv-"+timestamp+"/"
    os.makedirs(tmpFolder)
    os.makedirs(igvFolder)
    cmd="mv "+movingList+" "+tmpFolder
    subprocess.call([cmd],shell=True)
    cmd="mv "+" ".join(vizFileList)+" "+igvFolder
    subprocess.call([cmd],shell=True)

#Combine VCFs and call consensus from Reference Genome
def vcfToConsensus(vcfList,refGenome,consensusSeq,outputFolder):
    allVCF=outputFolder+"/all.vcf"
    for eachVCF in vcfList:
        comSNP="bgzip -c "+eachVCF+" > "+eachVCF+".gz"
        subprocess.call([comSNP],shell=True)
        SNPindexCMD="bcftools index "+eachVCF+".gz"
        subprocess.call([SNPindexCMD],shell=True)
    gzVCF=list(map(lambda x: x+".gz",vcfList))
    mergeVCF="bcftools merge --merge all --force-samples "+" ".join(gzVCF)+" -O v | bcftools view --min-ac 1 > "+allVCF
    compAll="bgzip "+allVCF
    allindexCMD="bcftools index "+allVCF+".gz"
    consensus="cat "+refGenome+" | bcftools consensus -s Sample1 "+allVCF+".gz > "+consensusSeq
    subprocess.call([mergeVCF],shell=True)
    subprocess.call([compAll],shell=True)
    subprocess.call([allindexCMD],shell=True)
    subprocess.call([consensus],shell=True)


def main():
    readsList=[]
    refGenome=None
    refType="haploid"
    outputFolder=None
    gffFile=None
    singletonFile=None
    outputPrefix="result"
    threadN="10"
    usage="Usage: python3.5 ReferenceAssembly.py -pe-1 <fastq R1> -pe-2 <fastq R2 (optional)> -pe-s <fastq singleton (optional)> -ref <refenerence genome> -genomeType <haploid (default) or diploid> -out <outputFolder> -gff <reference genome gff (optional)> -prefix <output file prefix. Default: result> -threads <threads number. default: 10>"
    for idx in range(len(sys.argv)):
        if (sys.argv[idx] == "-ref") and (len(sys.argv) > idx + 1):
            refGenome = sys.argv[idx + 1]
        elif (sys.argv[idx] == "-h") or (sys.argv[idx] == "-help"):
            print (usage)
            sys.exit()
        elif (sys.argv[idx] == "-genomeType") and (len(sys.argv) > idx + 1):
            refType=sys.argv[idx + 1]
        elif ("-pe-" in sys.argv[idx]) and (len(sys.argv) > idx + 1):
            if sys.argv[idx]=="-pe-s":
                singletonFile=sys.argv[idx + 1]
            else:
                readsList.append(sys.argv[idx + 1])
        elif (sys.argv[idx] == "-out") and (len(sys.argv) > idx + 1):
            outputFolder=sys.argv[idx + 1]
            if not os.path.exists(outputFolder):
                os.makedirs(outputFolder)
        elif (sys.argv[idx] == "-prefix") and (len(sys.argv) > idx + 1):
            outputPrefix=sys.argv[idx + 1]
        elif (sys.argv[idx] == "-gff") and (len(sys.argv) > idx + 1):
            gffFile=sys.argv[idx + 1]
        elif (sys.argv[idx] == "-threads") and (len(sys.argv) > idx + 1):
            threadN=str(sys.argv[idx + 1])
    if refGenome and readsList and outputFolder:
        vizFileList=[]
        copyRef="cp "+refGenome+" "+outputFolder+"/ "
        subprocess.call([copyRef],shell=True)
        refGenome=outputFolder+"/"+os.path.basename(refGenome)
        sys.stdout.flush() 
        print ("Mapping the reads to the reference genome...",end=" ")
        samOutput=runbwa(refGenome,readsList,outputFolder+"/"+outputPrefix,threadN,singletonFile)
        sys.stdout.flush()
        print ("Done")
        print ("Processing mapping results...",end=" ")
        (rawVCFresult,vizFileList,finalBam)=ProcessMapping(samOutput,refGenome,refType,vizFileList)
        sys.stdout.flush()
        print ("Done")
        print ("Filtering SNPs and indels...",end=" ")
        (filterSNPvcf,filterIndelvcf)=FilterVCF(rawVCFresult,refGenome)
        sys.stdout.flush()
        print ("Done")
        vizFileList.append(filterSNPvcf)
        vizFileList.append(filterIndelvcf)
        vizFileList.append(refGenome)
        vizFileList.append(refGenome+".fai")
        newConsensus=outputFolder+"/consensus.fa"
        vcfToConsensus([filterSNPvcf,filterIndelvcf],refGenome,newConsensus,outputFolder)
        SNPvariant=None
        indelVariant=None
        if gffFile:
            gffFile=GFFtovepGFF(gffFile,outputFolder)
            sys.stdout.flush()
            print ("Variant annotation...",end=" ")
            SNPvariant=annVCF(filterSNPvcf,gffFile,refGenome)
            indelVariant=annVCF(filterIndelvcf,gffFile,refGenome)
            sys.stdout.flush()
            print ("Done")
            vizFileList.append(gffFile)
        else:
            print ("Missing gff file, the variant annotation skipped!")


        SNPvcfDict={}
        SNPvcfList=[]
        SNPannImport=defaultdict(list)
        SNPannSupp=defaultdict(list)
        annHeaderSNP=None
        indelvcfDict={}
        indelvcfList=[]
        indelannImport=defaultdict(list)
        indelannSupp=defaultdict(list)
        annHeaderIndel=None
        SNPcoverageV=0.0
        SNPpassN=0
        SNPtotalBase=0
        IndelcoverageV=0.0
        IndelpassN=0
        IndeltotalBase=0
        SNPCovDict={}
        indelCovDict={}
        if filterSNPvcf and SNPvariant \
                and (os.path.isfile(filterSNPvcf)) \
                and (os.path.isfile(SNPvariant)):
            (SNPvcfDict,SNPvcfList,SNPannImport,SNPannSupp,annHeaderSNP)=VarAnnCsv(filterSNPvcf,SNPvariant)
        else:
            (SNPvcfDict,SNPvcfList,SNPannImport,SNPannSupp,annHeaderSNP)=VarAnnCsv(filterSNPvcf,SNPvariant)
        if filterIndelvcf and indelVariant \
                and (os.path.isfile(filterIndelvcf)) \
                and (os.path.isfile(indelVariant)):
            (indelvcfDict,indelvcfList,indelannImport,indelannSupp,annHeaderIndel)=VarAnnCsv(filterIndelvcf,indelVariant)
        else:
            (indelvcfDict,indelvcfList,indelannImport,indelannSupp,annHeaderIndel)=VarAnnCsv(filterIndelvcf,indelVariant)
        if filterSNPvcf and os.path.isfile(filterSNPvcf):
            (SNPcoverageV,SNPpassN,SNPtotalBase,SNPCovDict)=calCov(filterSNPvcf,finalBam)
        if filterIndelvcf and os.path.isfile(filterIndelvcf):
            (IndelcoverageV,IndelpassN,IndeltotalBase,indelCovDict)=calCov(filterIndelvcf,finalBam)
            
        if SNPtotalBase+IndeltotalBase != 0:
            averageCov=(float(SNPcoverageV+IndelcoverageV))/float(SNPtotalBase+IndeltotalBase)
        else:
            averageCov = 0
            print("WARNING: Cannot calculate averageCov. Resetting it to zero")
            print("SNPtotalBase:", SNPtotalBase)
            print("IndeltotalBase:", IndeltotalBase)

        vizFolder(vizFileList,outputFolder)
        vizFileList.remove(refGenome)
        vizFileList.remove(refGenome+".fai")
        CleanvizFileList=list(filter(lambda x: ".bai" not in x,vizFileList))
        newvizFileList=list(map(lambda x:"[url]/"+os.path.basename(x),CleanvizFileList))
        igvLink="http://www.broadinstitute.org/igv/projects/current/igv.php?sessionURL="+",".join(newvizFileList)+"&maxHeapSize=2000m&genome=[url]/"+os.path.basename(refGenome)
        outputLink=open(outputFolder+'/igvlink.txt','a')
        print (igvLink,file=outputLink)
        outputLink.close()
        QCreport=open(outputFolder+'/QCreport.txt','w')
        print ("# of SNPs:"+str(SNPpassN)+"\n# of Indels:"+str(IndelpassN)+"\nAverage coverage:"+str(averageCov),file=QCreport)
        print ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tCOVERAGE",file=QCreport)
        for eachItem,value in SNPCovDict.items():
            print (eachItem+"\t"+value,file=QCreport)
        for eachItem,value in indelCovDict.items():
            print (eachItem+"\t"+value,file=QCreport)
        QCreport.close()
        shutil.move(outputFolder+"/tmp/consensus.fa", outputFolder+"/consensus.fa")

        with open(outputFolder+'/variant-ann.csv', 'a', newline='') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',',quoting=csv.QUOTE_MINIMAL)
            if annHeaderSNP:
                spamwriter.writerow(annHeaderSNP.split("\t"))
            elif annHeaderIndel:
                spamwriter.writerow(annHeaderIndel.split("\t"))
            for idx,eachEle in enumerate(SNPvcfList):
                flag=0
                if eachEle in SNPannImport:
                    for eachAnn in SNPannImport[eachEle]:
                        spamwriter.writerow(SNPvcfDict[eachEle].split("\t")+eachAnn.split("\t"))
                        flag=1
                if eachEle in SNPannSupp:
                    for eachAnn in SNPannSupp[eachEle]:
                        spamwriter.writerow(SNPvcfDict[eachEle].split("\t")+eachAnn.split("\t"))
                        flag=1
                if flag==0:
                    spamwriter.writerow(SNPvcfDict[eachEle].split("\t")+["-"]*11)
                
        with open(outputFolder+'/variant-ann.csv', 'a', newline='') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',',quoting=csv.QUOTE_MINIMAL)
            for idx,eachEle in enumerate(indelvcfList):
                positionInfor=int(eachEle.split(":")[-1])
                altPre=eachEle.split(":")[0]+":"+str((positionInfor-1))
                altAfter=eachEle.split(":")[0]+":"+str((positionInfor+1))
                flag=0
                if eachEle in indelannImport:
                    flag=1
                    for eachAnn in indelannImport[eachEle]:
                        spamwriter.writerow(indelvcfDict[eachEle].split("\t")+eachAnn.split("\t"))
                
                if eachEle in indelannSupp:
                    flag=1
                    for eachAnn in indelannSupp[eachEle]:
                        spamwriter.writerow(indelvcfDict[eachEle].split("\t")+eachAnn.split("\t")+[indelCovDict[eachEle]])
                if altPre in indelannImport:
                    flag=1
                    for eachAnn in indelannImport[altPre]:
                        spamwriter.writerow(indelvcfDict[eachEle].split("\t")+eachAnn.split("\t"))
                if altPre in indelannSupp:
                    flag=1
                    for eachAnn in indelannSupp[altPre]:
                        spamwriter.writerow(indelvcfDict[eachEle].split("\t")+eachAnn.split("\t"))
                if altAfter in indelannImport:
                    flag=1
                    for eachAnn in indelannImport[altAfter]:
                        spamwriter.writerow(indelvcfDict[eachEle].split("\t")+eachAnn.split("\t"))
                if altPre in indelannSupp:
                    flag=1
                    for eachAnn in indelannSupp[altPre]:
                        spamwriter.writerow(indelvcfDict[eachEle].split("\t")+eachAnn.split("\t"))
                if flag==0:
                    spamwriter.writerow(indelvcfDict[eachEle].split("\t")+["-"]*11)
    else:
        print (usage)
        sys.exit()
        


if __name__ == '__main__':
    main()
