#!/usr/bin/env /Local/genomics/virtualenvs/py2711/bin/python
"""
requiers HTseq (http://www-huber.embl.de/users/anders/HTSeq/)
parse sam/bam alignment and extract reads/hits by quality or gene they were aligned to.
counting is adopted from HTseq-count script and uses the "Union" mode only.
requiers SAM file to have and MD string (either inserted by aligner or added using samtools calmd)
and (http://code.google.com/p/pysam/).
Input is a bam/sam file sorted or a stream.
gff file is not required.
currently: count +1 if EITHER of the PE reads aligned. DOES NOT count BOTH. 
add different options for counting? add the other arg parse options
"""

import os, sys, string
import csv
import re
#sys.path.insert(0, '/srv01/technion/philosof/.local/lib/python2.7/site-packages/HTSeq-0.6.1-py2.7-linux-x86_64.egg/HTSeq/')
import HTSeq #from where, on atlas?
import itertools
import warnings, traceback, argparse

current_dir = os.getcwd()

def main():
    exe_parser = argparse.ArgumentParser()
    exe_parser.add_argument('infile', type=str, help='<input file> [(full path), -b/-s required]')
    exe_parser.add_argument("-s", "--samout", help = "output not aligned reads to [file path].", type=str)
    exe_parser.add_argument("-v", "--verbose", help="verbose. (default = TRUE).", action = "store_true")
    exe_parser.add_argument("-i", "--min_id", help = "minimal percent id of hit to consider. if 0 then ignored. (default = 0).", type=int)
    exe_parser.add_argument("-l", "--min_len", help = "minimal read_2_length. (default = 90).", type=int)
    exe_parser.add_argument("-c", "--max_clip", help = "proportion of bases clipped from read for alignment. (default = 0.3).", type=float)
    exe_parser.add_argument("-m", "--mode", help="(p)aired end or (s)ingle end (default=n)", type=str)
    exe_parser.add_argument("--order", help="order of bam/sam: (n)ame or p(osition). default=n", type=str)

    exe_parser.add_argument("-o", "--out", help="name of counts output file.", type=str)
    
    args = exe_parser.parse_args()

    paired_end = True
    pe_order = 'n'
    #min_read_len = 90.0  
    
    if args.verbose:
        verbose = True
    else:
        verbose = False
    

    #max_clip_ = float(0.1) 
    min_id = float(98.0)
    min_score = 0


    ####
    if args.samout:
        samoutfile = open(args.samout, "w")
    else:
        samoutfile = None
    if args.out:
        outfile = open(args.out, "w")
    else:
        outfile = None
    
    if args.min_len:
        min_read_len  = args.min_len
    else:
        min_read_len = 80.0
    if args.min_id:
        min_id = args.min_id
    else:
        min_id = float(0)

    if args.max_clip:
        max_clip_ = args.max_clip
    else:
        max_clip_ = float(0.1) 


    hit_dict = dict()
    notaligned = 0
    
    if args.infile:
        try: #try opening the file
            if args.infile == '-': #get sam on a stream
                seqfile = HTSeq.SAM_Reader(sys.stdin)
                if paired_end:
                    read_seq_iter = iter(seqfile)
                    first_read = read_seq_iter.next()
                    if pe_order == 'p':
                        reader = HTSeq.pair_SAM_alignments_with_buffer(seqfile)
                    elif pe_order == 'n':
                        reader = HTSeq.pair_SAM_alignments(seqfile)#(read_seq)
                else:
                    reader = seqfile
            elif args.infile != '-':
                    seqfile= HTSeq.SAM_Reader(args.infile)#SAM_Reader(args.infile)
                    if args.paired_end_mode:
                        read_seq_iter = iter(seqfile)
                        first_read = read_seq_iter.next()
                        read_seq = itertools.chain( [ first_read ], read_seq_iter)
                        reader = HTSeq.pair_SAM_alignments(read_seq)
                        if pe_order == 'p':
                            reader = HTSeq.pair_SAM_alignments_with_buffer(reader)
                        elif pe_order == 'n':
                            reader = HTSeq.pair_SAM_alignments(reader)
                    else:
                        reader = seqfile
                    #fread_seq_iter = iter(reader)
                    #first_read = iter(read_seq).next()
            elif args.infile == '':
                print "no input file type given. exiting..."
                sys.exit(1)
        except:
            print "failed processing SAM/BAM file"
            raise
    elif not args.infile:
        print "no input file given. exiting..."
        sys.exit(1)
    read_counter = 0
    for alignment in reader: 
        read_counter += 1
        read1_name = ''
        read_2_name = '' 
        if read_counter % 1000000 == 0 and verbose:
            if verbose:
                print read_counter, 'alignment pairs processed'
        if (alignment[0] is None) or not alignment[0].aligned:
            if (alignment[1] is None) or not alignment[1].aligned:
                notaligned += 1
        else:
            iv_seq = tuple()
            if (alignment[0] is not None) and (alignment[0].aligned):
                read_1_name = alignment[0].read.name
                read_1_length = len(alignment[0].read.seq)
                #print 'read_1_length', read_1_length, alignment[0].read.seq, alignment[0].original_sam_line
                opt_1_fields = alignment[0].optional_fields
                flag_1 = alignment[0].flag
                contig = alignment[0].iv.chrom
                id_check_1 = False
                if min_id > 0:
                    cigar_1_string = parse_cigar(alignment[0].original_sam_line.split('\t')[5]) #just the cigar string without the fancy HTseq additions
                    cigar_1_soft_clipped, cigar_1_M, cigar_1_insertions, cigar_1_deletions, cigar_1_insertions = parse_cigar_alignment(cigar_1_string)
                    score_1, md_1_matches, md_1_deletions, md_1_mismatches = parse_opt_fields(opt_1_fields) #get alignment data from md string
                    clipped_1 = (float(cigar_1_soft_clipped)/float(read_1_length))
                    percent_1_id = (100.0 * ((float(md_1_matches) / (float(read_1_length - cigar_1_soft_clipped + cigar_1_insertions + cigar_1_deletions)))))
                    if float(percent_1_id) >= min_id:
                        if (float(cigar_1_soft_clipped)/float(read_1_length)) <= float(max_clip_): #check the clipping percentag
                            iv_seq_good_1 = True
                            id_check_1 = True
                        else:
                            if args.samout:
                                write_to_samout(samoutfile, paired_end, alignment, 'too_many_bases_clipped_from_read=' + str(cigar_1_soft_clipped))
                            iv_seq_good_1 = False
                    else:
                        id_check_1 = False
                        if args.samout:
                            write_to_samout(samoutfile, paired_end, alignment, "percent_lower_than_threshold=" + str(percent_1_id))
                else:
                    id_check_1 = True
                    iv_seq_good_1 = True
                
                if (int(read_1_length) >= int(min_read_len)) or (read_1_length == 1): #LOOK INTO THIS CONDITION? HOW TO GET LEN OF * IN SAM
                    #if float(percent_1_id) >= min_id:#float(args.min_id):
                    if id_check_1 == True:
                        iv_seq = (cigar_operation.ref_iv for cigar_operation in alignment[0].cigar if cigar_operation.type == "M" and cigar_operation.size > 0)
                        iv_seq_good_1 = True
                    #else:
                    #    iv_seq_good_1 = False
                else:
                    iv_seq_good_1 = False
                    if args.samout:
                        write_to_samout(samoutfile, paired_end, alignment, 'read_too_short=' + str(read_1_length))
            else:
                iv_seq_good_1 = False
            #now the second read
            if (alignment[1] is not None) and (alignment[1].aligned):
                read_2_name = alignment[1].read.name
                #print read_2_name, read_1_name
                read_2_length = len(alignment[1].read.seq)
                #print 'read_2_length', read_2_length, alignment[1].read.seq
                opt_2_fields = alignment[1].optional_fields
                flag_2 = alignment[1].flag
                contig = alignment[1].iv.chrom
                id_check_2 = False
                if min_id > 0:
                    cigar_2_string = parse_cigar(alignment[1].original_sam_line.split('\t')[5]) #just the cigar string without the fancy HTseq additions
                    cigar_2_soft_clipped, cigar_2_M, cigar_2_insertions, cigar_2_deletions, cigar_2_insertions = parse_cigar_alignment(cigar_2_string)
                    score_2, md_2_matches, md_2_deletions, md_2_mismatches = parse_opt_fields(opt_2_fields) #get alignment data from md string
                    clipped_2 = (float(cigar_2_soft_clipped)/float(read_2_length))
                    percent_2_id = (100.0 * ((float(md_2_matches) / (float(read_2_length - cigar_2_soft_clipped + cigar_2_insertions + cigar_2_deletions)))))
                    if float(percent_2_id) >= min_id:
                        if (float(cigar_2_soft_clipped)/float(read_2_length)) <= float(max_clip_): #check the clipping percentag
                            iv_seq_good_2 = True
                            id_check_2 = True
                        else:
                            if args.samout:
                                write_to_samout(samoutfile, paired_end, alignment, 'too_many_bases_clipped_from_read=' + str(cigar_2_soft_clipped))
                            iv_seq_good_2 = False
                    else:
                        id_check_2 = False
                        if args.samout:
                            write_to_samout(samoutfile, paired_end, alignment, "percent_lower_than_threshold=" + str(percent_2_id))
                else:
                    id_check_2 = True
                    iv_seq_good_2 = True
                
                if (int(read_2_length) >= int(min_read_len)) or (read_2_length == 1):
                    #if float(percent_1_id) >= min_id:#float(args.min_id):
                    if id_check_2 == True:
                        iv_seq = itertools.chain(iv_seq,(cigar_operation.ref_iv for cigar_operation in alignment[1].cigar if cigar_operation.type == "M" and cigar_operation.size > 0)) #get the alignment
                        iv_seq_good_2 = True
                    else:
                        iv_seq_good_2 = False
                else:
                    iv_seq_good_2 = False
                    if args.samout:
                        write_to_samout(samoutfile, paired_end, alignment, 'read_too_short=' + str(read_2_length))
            else:
                iv_seq_good_2 = False

            temp_iv = []
            try: ###^*(^*^&^*^*^*^*^* ANY CONDITOINS TO ADD??
                for iv in iv_seq:
                    temp_iv.append(iv.chrom) #get name of chromosom/contig both PE reads align to, append to list.
                #print temp_iv, 'TRY', iv_seq_good_1, iv_seq_good_2
                temp_iv = set(temp_iv)
                #if temp_iv[0]  ==  temp_iv[1]: #if the read maps to the same contig/chromosom then +1 only once
                for i in temp_iv:
			if iv_seq_good_1:
                        	if iv_seq_good_2:
					if percent_1_id >= percent_2_id:
						outfile.write('\t'.join(map(str,[read_1_name, i, percent_1_id]))+'\n')
					else:					
						outfile.write('\t'.join(map(str,[read_2_name, i, percent_2_id]))+'\n')
				elif alignment[1] is None:
		    			outfile.write('\t'.join(map(str,[read_1_name, i, percent_1_id]))+'\n')
            except Exception, e:
                print e, '*****EXCEPT'
                #if read_1_name:
                #    if id_check_1==True and iv_seq_good_1 == True:

                        #print temp_iv, 'EXCEPT1_giid'
                        #if len(temp_iv) > 0:
                 #       for i in temp_iv:
                  #          if temp_iv[0] not in hit_dict:
                   #             hit_dict[temp_iv[0]] = 1
                    #        else:
                     #           hit_dict[temp_iv[0]] += 1
                 
                #elif read_2_name:
                 #   if id_check_2==True and iv_seq_good_2== True:     
                  #      #print temp_iv, 'EXCEPT2'
                   #     if len(temp_iv) > 0:
                    #        if temp_iv[0] not in hit_dict:
                     #           hit_dict[temp_iv[0]] = 1
                      #      else:
                       #         hit_dict[temp_iv[0]] += 1

def parse_cigar_alignment(cigar_string):
    #:wq""
    #parse cigar string to get alignment data: total number of insertions, clipped bases, deletions
    #""
    cigar_soft_clipped = 0
    cigar_M = 0
    cigar_insertions = 0
    cigar_deletions = 0
    cigar_insertions = 0
    for c in cigar_string: #('5S, 12M, 3D, 4M, 1I, 45M, 3S')
        if 'S' in c:
            c = c.replace('S','')
            cigar_soft_clipped += int(c)
        elif 'M' in c:
            c = c.replace('M','')
            cigar_M += int(c)
        elif 'D' in c:
            c = c.replace('D','')
            cigar_deletions += int(c)
        elif 'I' in c:
            c = c.replace('I','')
            cigar_insertions += int(c)
    return cigar_soft_clipped, cigar_M, cigar_insertions, cigar_deletions, cigar_insertions
    #print 'cigar_soft_clipped -->', cigar_soft_clipped, 'cigar_M -->', cigar_M, 'cigar_deletions -->', cigar_deletions, 'cigaer_insertion -->', cigar_insertions
    #print 'cigar total -->', cigar_soft_clipped + cigar_M + cigar_deletions + cigar_insertions

def parse_opt_fields(opt_fields):
    #""
    #parse optional fields coulmn in sam file.
    #get MD string alignment data.
    #get score --> AS
    #""
    score = int()
    md_matches = int()
    md_deletions = int()
    md_mismatches = int()
    for field in opt_fields:
        if 'MD' in field[0]:
            #print 'MD -->', field[1],
            md_list = re.split('(\D|\W)', field[1])
            for i in md_list:
                if i:
                    if i.isdigit():
                        md_matches += int(i)
                    elif "^" in i:
                        md_deletions += 1
                    else:
                        md_mismatches += 1
            total_md = md_mismatches + md_deletions + md_matches
        if 'AS' in field[0]:
            score = field[1]
    return score, md_matches, md_deletions, md_mismatches

def invert_strand(iv):
    iv2 = iv.copy()
    if iv2.strand == "+":
        iv2.strand = "-"
    elif iv2.strand == "-":
        iv2.strand = "+"
    else:
        raise ValueError, "Illegal strand"
    return iv2


def parse_cigar(s):
    s = s #"33S29M83S"
    l = list()
    ll = list()
    for i in s:
        if i.isdigit():
            ll.append(i)
        else:
            ll.append(i)
            l.append(''.join(ll))
            ll = list()
    return l

def write_to_samout(samoutfile, paired_end, alignment, assignment):
      if samoutfile is None:
         return
      if not paired_end:
         alignment = (alignment,)
      for read in alignment:
         if read is not None:
            samoutfile.write(read.original_sam_line.rstrip() + "\tXF:Z:" + assignment + "\n" )
            
            
main()
