#!/usr/bin/env python
import math
import sys
import pysam
import time
import os
#import multiprocessing
from optparse import OptionParser
'''

		*Program:bamToBw.py
		*Author:		Xukai Ma, Guo-Wei Li
		*Function:  Covert BAM files to bigWig files
		*Usage:
				python bamToBw.py -b bamFile -d
				-b: input BAM file
				-d: consider +/- strand
				-l: consider every reads' length
'''

def parseBamName(bamFile):
	if os.path.exists(bamFile):
		return bamFile
		
def creatChromeSize(bamFileName):
	preffixName, suffixName = os.path.splitext(bamFileName)
	tmpChromeSizeFilename = preffixName + ".chromesize"
	ftmp = open(tmpChromeSizeFilename, "w")

	for line in pysam.idxstats(bamFileName).strip().split('\n'):
		line = line.strip().split()
		if line[0] != "*":
			ftmp.write(line[0]+"\t"+line[1]+"\n")
	ftmp.close()
	return tmpChromeSizeFilename

def getReadLength(AlignSeg):
	#AlignSeg should be pysam.AlignedSegment object.
	length = AlignSeg.query_length
	if length > 0:
		return length
	length = AlignSeg.infer_query_length()

	return length

def divideBam(bamFileName):
	preffixName = os.path.splitext(bamFileName)[0]
	suffixName = os.path.splitext(bamFileName)[1]
	bamPlusFileName = preffixName+".PLUS"+suffixName
	bamMinusFileName = preffixName+".MINUS"+suffixName
	
	#Note: for 10x Genomics, library strand information is: [ + -> -F 0x10 ] [ - -> -f 0x10 ]
	creatPlusCMD = "samtools view -@ 4 -F 16 -h -b -o "+bamPlusFileName+" "+bamFileName ; os.system(creatPlusCMD)
	creatMinusCMD = "samtools view -@ 4 -f 16 -h -b -o "+bamMinusFileName+" "+bamFileName ; os.system(creatMinusCMD)
	
	return (bamPlusFileName, bamMinusFileName)

def convertName(bamFileNames):
	suffixbedgraph = ".bg"
	suffixbigwiggle = ".bw"
	suffixsortedbedgraph = "_sorted.bg" 
	
	bgFileNames = []
	bwFileNames = []
	sortedbgFileNames = []
	
	for bamFileName in bamFileNames:
		preffixFileName, suffixBam = os.path.splitext(bamFileName)
			
		bgFileName = preffixFileName+suffixbedgraph
		bgFileNames.append(bgFileName)
	
		bwFileName = preffixFileName+suffixbigwiggle
		bwFileNames.append(bwFileName)
			
		sortedbgFileName = preffixFileName+suffixsortedbedgraph 
		sortedbgFileNames.append(sortedbgFileName) 
	
	return (bgFileNames, bwFileNames, sortedbgFileNames) 


def converOneBam(bamFile, ratio, FLAG_STRAND, FLAG_LEHGTH):
	TOTAL_RATIO = 1
	
	if not os.path.exists(bamFile + '.bai'):
		pysam.index(bamFile)
	
	#creat chromesize tmp file
	chromeSizeFileName = creatChromeSize(bamFile)	

	if FLAG_STRAND == True:
		bamFileDivided = divideBam(bamFile)
	else:
		bamFileDivided = [bamFile]
	
	bgFileNames, bwFileNames, sortedbgFileNames = convertName(bamFileDivided)

	for i in range(len(bamFileDivided)):
		bamFileName = bamFileDivided[i]
		bgFileName = bgFileNames[i]
		bwFileName = bwFileNames[i]
		sortedbgFileName = sortedbgFileNames[i]

		outCMD = "genomeCoverageBed -bg -split -ibam "+bamFileName+\
						" -g "+chromeSizeFileName+" -scale "+str(TOTAL_RATIO)+\
						" > "+bgFileName
		os.system(outCMD)
		
		bgFileBackupName = bgFileName+".back"
		bgFile = open(bgFileName, "r")
		bgBackupFile = open(bgFileBackupName, "w")
		for line in bgFile:
			line = line.strip().split()
			line[3] = str(int(round(float(line[3]))))
			line = "\t".join(line)+"\n"
			bgBackupFile.write(line)
		bgBackupFile.close()
		bgFile.close()
		os.remove(bgFileName)
		os.rename(bgFileBackupName, bgFileName)
		
		outCMD = "sort -k1,1 -k2,2n "+bgFileName+" > "+sortedbgFileName
		os.system(outCMD) 
		os.system("rm "+bgFileName)
	
		outCMD = "bedGraphToBigWig "+sortedbgFileName+" "+chromeSizeFileName+" "+bwFileName
		os.system(outCMD)
		os.system("rm "+sortedbgFileName)
	os.remove(chromeSizeFileName)


if __name__ == '__main__':
	import pysam
	import os
	import multiprocessing

	from optparse import OptionParser
	parser = OptionParser()

	parser.add_option("-b", dest="bamFile",
                      help="bamFile")
	parser.add_option("-d", dest="FLAG_STRAND", action="store_true",
                      default=False,
                      help="To consider the strand")
	parser.add_option("-l", dest="FLAG_LEHGTH", action="store_true",
                      default=False, help="Consider lengths of all reads")
	(options, args) = parser.parse_args()
	bamFile = parseBamName(options.bamFile)

	converOneBam(bamFile, 1, options.FLAG_STRAND, options.FLAG_LEHGTH)
