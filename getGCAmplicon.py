#!/usr/bin/python

import re,sys,os

global script_path
script_path = os.path.dirname(os.path.abspath(__file__))

sys.path.append(script_path+'/CONFIG/')

from CONFIG import *

if __name__=="__main__":


    objC = CONFIG()
    cmd_dict = objC.parseGCCommandArgs()
    bed_file = cmd_dict['ampliconFile']; out_file = cmd_dict['outFile']
    two_bit = cmd_dict['twoBit'];two_bit_fa = cmd_dict['twoBitFasta']
    
    fh = open(bed_file)
    out_dir = os.path.dirname(bed_file)
    
    #two_bit = "/opt/NGS/binaries/twoBitToFa/default/twoBitToFa"
    #two_bit_fa = "/opt/NGS/References/hg19/2bit/hg19.2bit"

    wh = open(out_file,"w")

    print >>wh,"CHR,CHR_ST,CHR_EN,AMP_ID,G_Count,C_Count,GC_Count"
    for lines in fh:
        lines = lines.strip()
        if not re.search("^browser|^track",lines):
            strs = re.split("\t",lines)
            chr_num = strs[0].strip()
            chr_st = int(strs[1].strip())+1
            chr_end = int(strs[2].strip())
            amp_id = strs[3].strip()
            command = two_bit+" "+two_bit_fa+":"+chr_num+":"+str(chr_st)+"-"+str(chr_end)+" -noMask "+out_dir+"/test.fa"
            os.system(command)

            fh2 = open(out_dir+"/test.fa")
            gc_count =0
            dna_list = [] 
            for lines in fh2:
                lines = lines.strip()
                if not re.search(">",lines):
                    dna_list.append(lines)
                    gc_count = gc_count+1
            count_g = "".join(dna_list).count("G")
            count_c = "".join(dna_list).count("C")
            gc_count = count_g+count_c
            fh2.close()
            print_list = [chr_num]+[chr_st-1]+[chr_end]+[amp_id]+[count_g]+[count_c]+[gc_count]
            print  >>wh,",".join(["%s" % x for x in print_list])
		
    wh.close()
