#!/usr/bin/python

import re,sys,os,subprocess
import pysam,operator,csv

class BAM:

    def __init__(self,elements=[]):
        self.__elements={}
        for e in elements:
           self.__elements[e]=1

    
    def createDir(self,out_dir,batch_num):

        ''' Function to create subdirectories:NOT_MAPPED, RC, MAT_RC '''
        
        nm_dir = os.path.abspath(out_dir+'/'+batch_num+'/NOT_MAPPED/')
        rc_dir = os.path.abspath(out_dir+'/'+batch_num+'/RC/')
        mat_rc_dir = os.path.abspath(out_dir+'/'+batch_num+'/MAT_RC/')

        os.system('mkdir -p '+nm_dir)
        os.system('mkdir -p '+rc_dir)
        os.system('mkdir -p '+mat_rc_dir)

        dirDict = {'nmDir':nm_dir,'rcDir':rc_dir,'matDir':mat_rc_dir}

        return dirDict

    def readAmpliconBedFile(self,dirDict,amp_file):

        ''' Function to read amplicon design file by start/End position '''
        print ' -- Processing Amplicon design file '
        bed_hash ={}
        amplicon_by_start = {}
        amplicon_by_end = {}
        double_hash = {}
        amplicon_counts = {}

        fh1 = open(amp_file)

        for lines in fh1:
            lines = lines.strip()
            if not re.search('^browser|^track',lines):
                strs = re.split('\t',lines); strs = [x.strip() for x in strs]
                chr_str = strs[0]+"_"+str(strs[1])+"_"+str(strs[2])
                amp_id = strs[3]
                bed_hash[chr_str] = amp_id

                # Store by start position
                if not amplicon_by_start.has_key(str(strs[0])):
                    amplicon_by_start[str(strs[0])] = {}
                if not amplicon_by_start[str(strs[0])].has_key(int(strs[1])):
                    amplicon_by_start[str(strs[0])][int(strs[1])] = {}
                if amplicon_by_start[str(strs[0])][int(strs[1])].has_key(int(strs[2])):
                    doubles_hash[amp_id] = chr_str
                    continue
                amplicon_by_start[str(strs[0])][int(strs[1])][int(strs[2])] = amp_id

                # Store by end position
                if not amplicon_by_end.has_key(str(strs[0])):
                    amplicon_by_end[str(strs[0])] = {}
                if not amplicon_by_end[str(strs[0])].has_key(int(strs[2])):
                    amplicon_by_end[str(strs[0])][int(strs[2])] = {}
                amplicon_by_end[str(strs[0])][int(strs[2])][int(strs[1])] = amp_id

                # Counter
                amplicon_counts[amp_id] = 0

        fh1.close()
        
        ampInfoDict = {'amp':bed_hash, 'ampStart':amplicon_by_start,
                       'ampEnd': amplicon_by_end,'ampCount': amplicon_counts}
        
        return ampInfoDict

    def parseBAMFile(self,dirDict,bam_dir,out_dir,ampInfoDict,batch_num):

        ''' Function to process BAM files and getting read counts '''

        amplicon_counts = ampInfoDict['ampCount']
        amplicon_by_start = ampInfoDict['ampStart']
        amplicon_by_end = ampInfoDict['ampEnd']

        
        bam_list = os.listdir(bam_dir) 
        nm_dir = dirDict['nmDir']
        rc_dir = dirDict['rcDir']
        mat_rc_dir = dirDict['matDir']

        count_bam = 1

        wh1 = open(out_dir+'/'+batch_num+'/Alignment_Stats.txt','w')

        for ele in bam_list:
            
            if re.search('bam$',ele):
                sample_name = re.split('\.bam',ele)[0]
                sample_path = bam_dir+'/'+ele
                print ' -- Processing BAM file: ',ele,' Iteration count: ',count_bam
                bam_file = pysam.AlignmentFile(sample_path,'rb')
                #out_bam = pysam.AlignmentFile(nm_dir+'/'+sample_name+'_not_mapped.bam','wb',template=bam_file)
                
                read_hash = {}
                for amp_id in amplicon_counts:
                    amplicon_counts[amp_id] = 0
                #ampInfoDict = self.readAmpliconBedFile(dirDict,amp_file)
                #amplicon_counts = ampInfoDict['ampCount']

                # go through BAM using pysam
                bam_file = pysam.AlignmentFile(sample_path,'rb')
                exact_count = fuzzy_count = no_match = unmapped = align_problems = 0

                for read in bam_file.fetch():
                    #skip unmapped reads (or mate unmapped)
                    if read.is_unmapped or read.mate_is_unmapped:
                        unmapped += 1
                        continue

                    # some abnormal alignments : invalid chromosome names
                    if read.reference_id < 0 or read.next_reference_id < 0:
                        align_problems += 1
                        continue

                    # some abnormal alignments : different chromosomes
                    if read.reference_name != read.next_reference_name:
                        align_problems += 1
                        continue

                    # no probes on this chromosome
                    chr_txt = str(read.reference_name)
                    if not amplicon_by_start.has_key(chr_txt):
                        align_problems +=1
                        continue

                    start = read.reference_start
                    end   = read.reference_end

                    # This is the first read : cache
                    if not read_hash.has_key(read.query_name):
                        read_hash[read.query_name] = read
                        continue

                    # this is the mate : process
                    mate_start = read_hash[read.query_name].reference_start
                    mate_end   = read_hash[read.query_name].reference_end

                    # some abnormal alignments : same orientation
                    if read.is_reverse == read_hash[read.query_name].is_reverse:
                        align_problems += 1
                        read_hash.pop(read.query_name,None)
                        continue

                    min_start = min(start, mate_start)
                    max_end = max(end, mate_end)

                    if min_start in amplicon_by_start[chr_txt] and max_end\
                       in amplicon_by_start[chr_txt][min_start]:
                        amp_id = amplicon_by_start[chr_txt][min_start][max_end]
                        amplicon_counts[amp_id] += 1
                        exact_count += 1
                        read_hash.pop(read.query_name,None)
                    else:
                        #out_bam.write(read)
                        #out_bam.write(read_hash[read.query_name])
                        read_hash.pop(read.query_name,None)
                        no_match += 1

                bam_file.close()
                #out_bam.close()

                total_read_count = exact_count + fuzzy_count + no_match + unmapped + align_problems
                perc_exact = float(exact_count)/float(total_read_count)
               
                print >>wh1,ele,' The iternation count: ',count_bam
                print >>wh1,"Total Reads: %s " % (exact_count + fuzzy_count + 
                                                  no_match + unmapped + align_problems)

                print >>wh1,"Exact match : ",exact_count
                print >>wh1,"Fuzzy match : ",fuzzy_count
                print >>wh1,"No match    : ",no_match
                print >>wh1,"Not mapped  : ",unmapped
                print >>wh1,"Abnormal Alignments: ",align_problems
                print >>wh1,"Percentage Exact: ",perc_exact
                print >>wh1,"\n","###################################\n"
                
                wh2 = open(rc_dir+'/'+ele+'_counts.newflow.csv','w')
                for keys in amplicon_counts:
                    print >>wh2,keys+'\t'+str(amplicon_counts[keys])

                wh2.close()
                count_bam = count_bam+1

        wh1.close()

    def combineReadCounts(self,dirDict):

        ''' Function to combine read counts (Rc) per amplicons of all the samples '''
        
        print ' -- Combining read counts (RC) for all the samples '
        
        rc_dir = dirDict['rcDir']
        rc_dir_files = os.listdir(rc_dir)
        mat_rc_dir = dirDict['matDir']

        all_sample_mat = mat_rc_dir+'/ALLSampleMatrix.txt'
        wh = open(all_sample_mat,'w')

        for sample_name in rc_dir_files:

            sample_strs = re.split('\_counts',sample_name)
            sample_path = rc_dir+'/'+sample_name

            fh = open(sample_path)
            for lines in fh:
                lines = lines.strip()
                strs = re.split('\t',lines);strs = [x.strip() for x in strs]
                amp_id = strs[0]
                count = strs[1]
                print >>wh,amp_id+'\t'+str(count)+'\t'+sample_strs[0].strip()
            fh.close()

        wh.close()

        return all_sample_mat

       
    def createRCMatrix(self,dirDict,all_sample_mat,script_path):

        ''' Function to launch R script to create read count (RC) matrix '''

        print ' -- Converting read count (RC) to matrix '

        mat_dir = dirDict['matDir']

        cmd_dir = 'Rscript '+script_path+'/createMatrix.R  '+all_sample_mat+' '+mat_dir
        
        os.system(cmd_dir)

