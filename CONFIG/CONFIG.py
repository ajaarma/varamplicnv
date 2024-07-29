#!/usr/bin/python

import re,sys,os
import argparse

class CONFIG:

    def __init__(self,elements=[]):
        self.__elements={}
        for e in elements:
            self.__elements[e]=1

    def parseGCCommandArgs(self):

        cmdDict = {}

        parser = argparse.ArgumentParser(
                usage = '\npython getGCAmplicon.py '
                                '-a <amplicon-bed-file> -o <out-gc-file> '
                                '-b <two-bit-program> -f <two-bit-fasta>',
                version = '1.0',
                description='Program to get GC content corresponding to amplicon coordinates'
        )
        
        parser.add_argument('-a','--ampliconFile',help='Position Sorted\
                            Amplicion design file',action='store',
                            dest='bedFile',required=True)
        parser.add_argument('-o','--outFile',help='Output file name',
                            dest='outFile',action='store',required=True)
        parser.add_argument('-b','--twoBitToFa',help='Binaries for converting\
                            two bit to Fasta',action='store',dest='twoBitToFa',
                            required=True)
        parser.add_argument('-f','--twoBitFasta',help='Bit format fasta sequence',
                            action='store',dest='twoBitFasta',required=True)
        #parser.add_argument('-v','--version',help='show program\'s version and exit\n\n\n',
        #                    action='version',version='%(prog)s 1.0')

        self.results = parser.parse_args()
        self.bedFile = self.results.bedFile
        self.outFile = self.results.outFile
        self.twoBit = self.results.twoBitToFa
        self.twoBitFasta = self.results.twoBitFasta

        cmd_dict = {'ampliconFile':self.bedFile,'outFile':self.outFile,
                    'twoBit':self.twoBit,'twoBitFasta':self.twoBitFasta
                }

        print "\n\n"
        print "Entered Amplicon Design File is: ",self.bedFile
        print "Entered Output File for writing GC content correspond to Targeted Region is: ",self.outFile
        print "Entered binary for twoBitFasta : ",self.twoBit
        print "Entered bit format of Fasta sequence : ",self.twoBitFasta
        print "\n\n"

        return cmd_dict

    def parseRmDupAmpCommandArgs(self):

        cmdDict = {}

        parser = argparse.ArgumentParser(
                usage = '\npython rmDupAmp.py '
                                '-a <amplicon-bed-file> -o <out-gc-file> ',
                version = '1.0',
                description='Program to remove duplicate amplicon coordinates'
        )
        
        parser.add_argument('-a','--ampliconFile',help='Position Sorted\
                            Amplicion design file',action='store',
                            dest='bedFile',required=True)
        parser.add_argument('-s','--snpFile',help='Position Sorted\
                            SNP Amplicon design file', action='store',
                            dest='snpFile',required=True) 
        parser.add_argument('-r','--regionFile',help='Position Sorted\
                            ROI file',action='store',dest='regionFile',
                            required=True)
        parser.add_argument('-b','--bedTools',help='BedTools binaries',
                            action='store',dest='bedTools',required=True)
        parser.add_argument('-o','--outDir',help='Output directory',
                            dest='outDir',action='store',required=True)

        self.results = parser.parse_args()
        self.bedFile = self.results.bedFile
        self.snpFile = self.results.snpFile
        self.regFile = self.results.regionFile
        self.bedTools = self.results.bedTools
        self.outDir = self.results.outDir

        cmd_dict = {'ampliconFile':self.bedFile,'outDir':self.outDir,
                    'snpFile':self.snpFile,'regFile':self.regFile,
                    'bedTools':self.bedTools}

        print "\n\n"
        print "Entered Amplicon Design File is: ",self.bedFile
        print "Entered SNP Amplicon design file is: ",self.snpFile
        print "Entered Region of Interest (ROI) file is: ",self.regFile
        print "Entered Output directory is: ",self.outDir
        print "Entered bedtools binaries path is: ",self.bedTools
        print "\n\n"

        return cmd_dict

    def parseBAMCommandArgs(self):
        
        ''' Function to define command line arguments for processing BAM files '''
        cmdDict = {}
        parser = argparse.ArgumentParser(
                usage = '\npython processBAM.py '
                            '-i <input-BAM-directory> ' 
                            '-a <amplicon-bed-file> ' 
                            '-o <out-directory-path> '
                            '-b <analysis-batch-name> ',
                version = '1.0',
                description='Program to process BAM files and get read counts' 
                            'by mapping reads assigned uniquely to amplicons  '
        )
        parser.add_argument('-a','--ampliconFile',help='Position Sorted\
                            Amplicion design file',action='store',
                            dest='bedFile',required=True)
        parser.add_argument('-i','--inputBAM',help='Input directory of BAM\
                            files', action='store', dest='bamDir',required=True) 
        parser.add_argument('-o','--outDir',help='Output directory',
                            dest='outDir',action='store',required=True)
        parser.add_argument('-b','--batchNum',help='Analysis BATCH Name/Number',
                            dest='batchNum',action='store',required=True)
        self.results = parser.parse_args()
        self.bedFile = self.results.bedFile
        self.bamDir = os.path.abspath(self.results.bamDir)
        self.outDir = os.path.abspath(self.results.outDir)
        self.batch = self.results.batchNum

        cmd_dict = {'ampliconFile':self.bedFile,'outDir':self.outDir,
                    'bamDir':self.bamDir,'batch':self.batch}

        print "\n\n"
        print "Entered Amplicon Design File is: ",self.bedFile
        print "Entered BAM files directory is: ",self.bamDir
        print "Entered Output directory is: ",self.outDir
        print "Entered Analysis batch Name/Number is: ",self.batch
        print "\n\n"

        return cmd_dict
 
