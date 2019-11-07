#!/hpcf/apps/python/install/2.7.13/bin/python
import sys
import os
import pandas as pd
import datetime
import getpass
import uuid
import argparse

"""
module load bedtools
module load python/2.7.13
module load cas-offinder

"""

def my_args():

	mainParser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	mainParser.add_argument('-f','--bed_file',  help="input regions to look for gRNAs", required=True)
	mainParser.add_argument('-e','--extend',  help="extend search area to left and right", default=100)
	mainParser.add_argument('-l','--sgRNA_length',  help="sgRNA_length", default=20)
	mainParser.add_argument('-g','--genome_fa',  help="genome fasta", default="/home/yli11/Data/Human/hg19/fasta/hg19.fa")
	mainParser.add_argument('--PAM',  help="PAM sequence", default="NGG")
	mainParser.add_argument('-o','--output',  help="output bed file name", default="candidate_gRNA.bed")
	
	##------- add parameters above ---------------------
	args = mainParser.parse_args()	
	return args
def write_file(file_name,message):
	out = open(file_name,"wt")
	out.write(message)
	out.close()
def extend_bed(args):
	df = pd.read_csv(args.bed_file,sep="\t",header=None)
	df[1] = df[1].astype(int)-args.extend
	df[2] = df[2].astype(int)+args.extend
	outfile = str(uuid.uuid4()).split("-")[-1]
	df[[0,1,2]].to_csv(outfile,sep="\t",header=False,index=False)
	return outfile
	
def get_fasta(genome_fa,extended_file):
	out = extended_file+".fa"
	command = "bedtools getfasta -fi %s -bed %s -fo %s"%(genome_fa,extended_file,out)
	os.system(command)
	print (command)
	return out
	
def run_casOFFinder(genome_fasta,PAM,your_seq_list,nMisMatch=0):
	cas_input = str(uuid.uuid4()).split("-")[-1]
	cas_output = str(uuid.uuid4()).split("-")[-1]
	pattern = "N"*len(your_seq_list[0])+PAM
	config = [genome_fasta,pattern]
	for i in your_seq_list:
		config.append(i+PAM+" %s"%(nMisMatch))
	write_file(cas_input,"\n".join(config))
	command = "cas-offinder %s C %s;rm %s"%(cas_input,cas_output,cas_input)	
	os.system(command)
	print (command)
	return cas_output
def row_apply(x):
	if x[4] == "-":
		start = x[2]+3
	else:
		start = x[2]
	return start
	
def get_GC(x):
	GC_count = 0
	for i in x:
		if i.upper() in ['G','C']:
			GC_count+=1
	return GC_count/float(len(x))

def cas_to_bed(x,PAM,output):
	df = pd.read_csv(x,sep="\t",header=None)
	df['start'] = df.apply(row_apply,axis=1)
	df['end'] = df['start']+20
	df['seq'] = [x.replace(PAM,"") for x in df[0].tolist()]
	# print (df.head())
	## get GC%
	df[5] = df['seq'].apply(get_GC)
	df[[1,'start','end','seq',5,4]].to_csv(output,sep="\t",header=False,index=False)
	
def find_offtarget(locus_bed_file, gRNA_bed_file):
	outfile = str(uuid.uuid4()).split("-")[-1]
	command = "bedtools intersect -a %s -b %s -wa > %s"%(gRNA_bed_file,locus_bed_file,outfile)
	os.system(command)
	df = pd.read_csv(outfile,sep="\t",header=None)
	df2 = pd.read_csv(gRNA_bed_file,sep="\t",header=None)
	df2['name'] = df2[0]+":"+df2[1].astype(str)+"-"+df2[2].astype(str)
	df['name'] = df[0]+":"+df[1].astype(str)+"-"+df[2].astype(str)
	df.index = df[3]
	df3 = pd.DataFrame(df2.groupby(3)['name'].agg(', '.join))
	df['other matches'] = df3['name']
	def remove_self_match(r):
		myList = r['other matches'].split(", ")
		myList.remove(r['name'])
		return ",".join(myList)

	df['other matches'] = df.apply(remove_self_match,axis=1)

	def get_numOfftarget(r):
		if r['other matches'] == "":
			return 0
		myList = r['other matches'].split(",")
		return len(myList)
	df['numOffTargets'] = df.apply(get_numOfftarget,axis=1)
	df[[0,1,2,3,4,5,'numOffTargets','other matches']].to_csv(gRNA_bed_file+".off_targets.info.csv",index=False)
	df[[0,1,2,3,4,5]].to_csv(gRNA_bed_file,index=False,header=False,sep="\t")
	os.system("rm %s"%(outfile))

def main():

	args = my_args()
	extended_file = extend_bed(args)
	extended_fa = get_fasta(args.genome_fa,extended_file)
	# run 1 
	cas_output = run_casOFFinder(extended_fa,args.PAM,["N"*args.sgRNA_length])
	# get sequneces
	df = pd.read_csv(cas_output,sep="\t",header=None)
	df = df.dropna()
	candidate_gRNAs = [x[:-len(args.PAM)] for x in df[3].tolist()]
	print ("Total number of possible gRNA is: %s"%(len(candidate_gRNAs)))
	os.system("rm %s"%(cas_output))
	# run 2 to get genomic coordinates
	cas_output = run_casOFFinder(args.genome_fa,args.PAM,candidate_gRNAs)
	cas_to_bed(cas_output,args.PAM,args.output)
	## check if bedtools getfasta can get exact sequence - YES
	find_offtarget(extended_file, args.output)
	os.system("rm %s"%(extended_fa))
	os.system("rm %s"%(extended_file))
	os.system("rm %s"%(cas_output))
	
	
if __name__ == "__main__":
	main()


















































