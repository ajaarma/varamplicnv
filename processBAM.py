#!/usr/bin/python

import re,sys,os

global script_path
script_path = os.path.dirname(os.path.abspath(__file__))

sys.path.append(script_path+'/SCRIPTS/')
sys.path.append(script_path+'/CONFIG/')

from CONFIG import *
from BAM import *

if __name__=='__main__':

    objC = CONFIG()
    cmd_dict = objC.parseBAMCommandArgs()
    
    amp_file = cmd_dict['ampliconFile']
    bam_dir  = cmd_dict['bamDir']
    out_dir  = cmd_dict['outDir']
    batch_num = cmd_dict['batch']

    objB = BAM()
    # Create output sub-directories per batch
    dirDict = objB.createDir(out_dir,batch_num)

    # Read amplicon desing file by start coordinates
    ampInfoDict  = objB.readAmpliconBedFile(dirDict,amp_file)
    objB.parseBAMFile(dirDict,bam_dir,out_dir,ampInfoDict,batch_num) 

    # Combine all the read counts (RC) per amplicons of all the samples
    #for Matrix

    all_sample_mat = objB.combineReadCounts(dirDict)
    objB.createRCMatrix(dirDict,all_sample_mat,script_path)
     


    
