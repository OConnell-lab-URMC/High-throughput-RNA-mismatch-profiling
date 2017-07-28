"""

Bindnseq pipeline
	1. Get deduplicated reads
			Read input fastq file and count UMI kmers
			Use collections.Counter()
			Dump all entries where count==1, keep remaining with a flag
				kmer_counts[kmer] = [counts, flag]
			Read input fastq, get UMI
				If UMI is not in kmer_counts.keys, yield
				If UMI is in kmer_counts.keys and flag==False, yield
	2. Count deduplicated sequences

"""
import argparse
import collections
import subprocess
import itertools
import os
import gzip

def run_all():
	#make output directory
	if not os.path.exists(args.output_dir):
		os.makedirs(args.output_dir)
	output_dir = args.output_dir
	
	#run stuff
	print('Deduplicating reads and counting sequences')
	umi_counts, read_counts = count_reads_deduped(args.reads)
	
	print('Writing read counts and mismatches as text')
	read_data = get_mismatches(read_counts)
	write_read_data(read_data, '%s/read_counts.txt.gz' % output_dir)
	
	print('Writing UMI counts as text')
	write_kmer_counts(umi_counts, '%s/umi_counts.txt.gz' % output_dir)
	#write_kmer_counts(read_counts, '%s/read_counts.txt' % output_dir)

def write_read_data(read_data, output_loc):
	with gzip.open(output_loc, 'w') as writer:
		writer.write(b'kmer\tcount\tref_id\tnum_mismatches\tmismatches\n')
		for (kmer, data) in read_data.items():
			(read_count, ref_id, mismatch_count, mismatches) = data
			mismatch_str = ','.join([str(i) for i in mismatches])
			line = '%s\t%i\t%s\t%i\t%s\n' % \
				(kmer, read_count, ref_id, mismatch_count, mismatch_str)
			writer.write(line.encode('utf-8'))
	
def write_kmer_counts(kmer_counts, output_loc):
	with gzip.open(output_loc, 'w') as writer:
		writer.write(b'kmer\tcount\n')
		for(kmer, count) in kmer_counts.items():
			writer.write(
				('%s\t%i\n' % (kmer, count)).encode('utf-8'))
	
def get_mismatches(read_counts):
	output_data = {}
	
	for (kmer, count) in read_counts.items():
		mismatches = {}
		for ref_id, ref_seq in REF_SEQS.items():
			mismatches[ref_id] = []
			target = ref_seq[READ_START : READ_END]
			for i in range(0, len(target)):
				if(target[i] != kmer[i]):
					mismatches[ref_id].append(i)
				
		#find the closest match to this kmer
		min_distance = 1000000 #some huge numbe here
		closest_match = None
		for ref_id, mismatch_lst in mismatches.items():
			if (len(mismatch_lst) < min_distance):
				min_distance = len(mismatch_lst)
				closest_match = ref_id
		output_data[kmer] = (
			count,
			closest_match, 
			len(mismatches[closest_match]),
			mismatches[closest_match])
		
		#print(kmer, len(mismatches['L3']), len(mismatches['L4']))
		#print(output_data[kmer])
		
	return output_data	

def count_reads_deduped(reads_f):
	umi_counts = collections.Counter()
	read_counts = collections.Counter()
	
	reads_iter = read_fastq(reads_f)
	while(True):
		try:
			read = next(reads_iter)
		except StopIteration:
			break		
		
		#UMI is stored in the read header
		umi_entry = read[0].split(' ')[0].split(':')[-1]
		umi = umi_entry[UMI_START : UMI_END]
		umi_counts.update([umi])
					
		#if this is the first observation of this umi
		if(umi_counts[umi] == 1):
			seq = read[1][READ_START : READ_END]
			#add to read counts
			read_counts.update([seq])
	return(umi_counts, read_counts)
				
def read_fastq(fq_file_lst):
	"""
	Args:
		fq_file (string):
			Filename for a fastq file
	Yields:
		read (list): Data from file, 4 lines at a time
		offset (int): Character offset for this read
	"""	
	fq_file_lst = fq_file_lst.split(',')
	zcat_cmd = ['zcat'] + fq_file_lst
	zcat = subprocess.Popen(
		zcat_cmd, stdout = subprocess.PIPE)
	
	i = 0
	for lines in itertools.zip_longest(*[zcat.stdout] *4):
		lines = [(l.rstrip()).decode('utf-8') for l in list(lines)]
		i += 1
		yield(lines)
		#if(i > 1000000):
		#	break
	
def get_args():
	parser = argparse.ArgumentParser(
		description = 'This script deduplicates bindnseq data')
		
	#wow this is verbose
	parser.add_argument(
		'--reads', 
		type=str, 
		help='Fastq.gz file, or comma separated list of fastq.gz files', 
		required=True)
	parser.add_argument(
		'--output_dir', 
		type=str, 
		help='Directory where outputs are written', 
		required=True)
	
	args = parser.parse_args()
	#args.barcodes = args.barcodes.split(',')
	return args

if __name__ == '__main__':
	args = get_args()
	"""
	#TO DO
	#sanitize args / check inputs
	
	
	#
	"""
	#constants go here (later in params.json)
	NUM_PROCESSES = 32
	
	#8bp- skip
	#4bp- outside guide but interesting (PFS etc)
	#12bp- guide
	#4bp- outside guide but interesting
	#rest- skip
	
	READ_START = 9		#guide starts at position 12 of the seq below
							#disregard bases more than 3bp away to reduce data
							
	READ_END = 35		#guide is 20 bp long + 3bp extra
	UMI_START = 0 
	UMI_END = 13
	REF_SEQS = {
		'L4': 'NNNNNNGCAGATATAGCCTGGTGGTTCAGGCGGCGCATGCTTAAGATCGGA',
		'L3': 'NNNNNNTGGCTGGTGAACTTCCGATAGTGCGGGTGTTGAATCCAGATCGGA'}
	run_all()





