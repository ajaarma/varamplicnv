#!/usr/bin/python

import re,sys,os
from collections import OrderedDict

global script_path
script_path = os.path.dirname(os.path.abspath(__file__))

sys.path.append(script_path+'/CONFIG/')

from CONFIG import *

if __name__=="__main__":

    objC = CONFIG()
    cmd_dict = objC.parseRmDupAmpCommandArgs()
    
    amp_file = cmd_dict['ampliconFile']
    snp_file = cmd_dict['snpFile']
    reg_file = cmd_dict['regFile']
    out_dir = cmd_dict['outDir']
    bed_tools = cmd_dict['bedTools']

	#amp_file = sys.argv[1]
	#out_file = sys.argv[2]

    # Processing amplicon design file
    amp_hash = OrderedDict()
    fh = open(amp_file)

    for lines in fh:	
        lines = lines.strip()
        if not re.search("^browser|track",lines):
            strs = re.split("\t",lines); strs = [x.strip() for x in strs ]
            chr_num = strs[0]
            chr_st = strs[1]
            chr_en = strs[2]
            amp_id = strs[3]
            key_str = chr_num+"_"+str(chr_st)+"_"+str(chr_en)
            
            if amp_hash.has_key(key_str):
                
                tmp = amp_hash[key_str]
                tmp.append(amp_id)
                amp_hash[key_str] = tmp
            else:	
                amp_hash[key_str] = [amp_id]

    fh.close()
    
    #Processing snp-amplicon design file
    snp_hash = {}
    fh = open(snp_file)

    for lines in fh:
        lines = lines.strip()
        
        if not re.search('^browser|^track',lines):
            strs = re.split('\t',lines); strs = [x.strip() for x in strs]
            chr_num = strs[0]
            st_coord = strs[1]
            en_coord = strs[2]
            key_str = chr_num+"_"+st_coord+"_"+en_coord
            snp_hash[key_str] = 1

    # Output the amplicon coordinates without SNPs
    for keys in amp_hash:
        if snp_hash.has_key(keys):
            del amp_hash[keys]
                
    # Remove duplicate coordinates 
    amp2_hash = OrderedDict()

    for key_str in amp_hash:
        amp_id = amp_hash[key_str]

        if amp2_hash.has_key(key_str):
            tmp = amp2_hash[key_str]
            tmp.append(amp_id)
            amp2_hash[key_str] = tmp
        else:
            amp2_hash[key_str] = amp_id
    
    # Output the amplicon design file after filtering for SNPs
    tmp_file = os.path.abspath(out_dir)+'/tmp_AmpRmSNPDup.bed'
    out_file = os.path.abspath(out_dir)+'/AmpRmSNPRmDup.bed'
    wh = open(tmp_file,'w')

    for keys in amp2_hash:
        strs = re.split('\_',keys)
        print >>wh,'\t'.join(strs)+'\t'+amp2_hash[keys][0]
    wh.close()

    # Sort the amplicon design file obtained in previous step
    sort_cmd = 'sort -n -k1.4 -k2,2n '+tmp_file+' > '+out_file
    os.system(sort_cmd)
    os.system('rm '+tmp_file)

    # Add regions to ROI file
    tmp_reg_file = os.path.abspath(out_dir)+'/tmpGeneRegion.bed'
    sort_reg_gene_file = os.path.abspath(out_dir)+'/tmpSortROI.bed'
    merge_reg_gene_file = out_dir+'/tmpSortMergeROI.bed'
    sort_merge_gene_region_file = out_dir+'/sortMergeROIGene.bed'

    # Processing Input region file
    fh = open(reg_file)
    wh = open(tmp_reg_file,'w')

    for lines in fh:
        lines = lines.strip()
        if not re.search('^browser|^track',lines):
            strs = re.split("\t",lines);strs = [x.strip() for x in strs]
            print >>wh,'\t'.join([strs[0],strs[1],strs[2],strs[3]])
    wh.close()

    sort_cmd = 'sort -V -k1,1 -k2,2n '+tmp_reg_file+' > '+sort_reg_gene_file
    merge_cmd = bed_tools+' merge -i '+sort_reg_gene_file+' -c 4 -o distinct > '+merge_reg_gene_file

    os.system(sort_cmd)
    os.system(merge_cmd)
    
    fh = open(merge_reg_gene_file)

    gene_exons = OrderedDict()

    # Start processing the Gene-Region file and using hasing technique to store the coordinated per Gene
    for lines in fh:
        lines = lines.strip()
        strs = re.split('\t',lines); strs = [x.strip() for x in strs]
        gene_id = strs[3]
        if gene_exons.has_key(gene_id):
            tmp = gene_exons[gene_id]
            tmp.append(lines)
            gene_exons[gene_id] = tmp
        else:
            gene_exons[gene_id] = [lines]

    fh.close()
    
    # Iterate through the Gene-Exon dictionary to add sequentially the Regions to each of these coordinates of Region-Of-Interest
    wh = open(sort_merge_gene_region_file,'w')
    for keys in gene_exons:
        exon_length = len(gene_exons[keys])
        for i in range(0,exon_length):
            strs = re.split('\t',gene_exons[keys][i])
            print >>wh,strs[0].strip()+"\t"+str(strs[1].strip())+"\t"+str(strs[2].strip())+"\t"+strs[3].strip()+"|Region_"+str(i)
    wh.close()

    os.system('rm '+tmp_reg_file)
    os.system('rm '+sort_reg_gene_file)
    os.system('rm '+merge_reg_gene_file)

        
    
