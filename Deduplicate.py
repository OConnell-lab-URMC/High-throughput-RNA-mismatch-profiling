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

def run_all():
	#make output directory
	if not os.path.exists(args.output_dir):
		os.makedirs(args.output_dir)
	output_dir = args.output_dir
	
	#run stuff
	print('Deduplicating reads and counting sequences')
	umi_counts, read_counts = count_reads_deduped(args.reads)
	print('Writing UMI counts as text')
	write_kmer_counts(umi_counts, '%s/umi_counts.txt' % output_dir)
	#write_kmer_counts(read_counts, '%s/read_counts.txt' % output_dir)
	
	read_data = get_mismatches(read_counts)
	write_read_data(read_data, '%s/read_counts.txt' % output_dir)

def write_read_data(read_data, output_loc):
	with open(output_loc, 'w') as writer:
		writer.write('kmer\tcount\tnum_mismatches\tmismatches\n')
		for (kmer, data) in read_data.items():
			(read_count, mismatch_count, mismatches) = data
			mismatch_str = ','.join([str(i) for i in mismatches])
			line = '%s\t%i\t%i\t%s\n' % (kmer, read_count, mismatch_count, mismatch_str)
			writer.write(line)
	
	

def write_kmer_counts(kmer_counts, output_loc):
	with open(output_loc, 'w') as writer:
		writer.write('kmer\tcount\n')
		for(kmer, count) in kmer_counts.items():
			writer.write('%s\t%i\n' % (kmer, count))
			#print('%s\t%i' % (kmer, count))
	

def get_mismatches(read_counts):
	output_data = {}
	
	target = REF_SEQ[READ_START : READ_END]
	for (kmer, count) in read_counts.items():
		mismatches = []
		for i in range(0, len(target)):
			if(target[i] != kmer[i]):
				mismatches.append(i)
		output_data[kmer] = (count, len(mismatches), mismatches)
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
		umi = read[0].split(' ')[0].split(':')[-1]
		#last base of umi is part of adaptor
		umi = umi[UMI_START : UMI_END]
		umi_counts.update([umi])
				
		#if this is the first observation of this umi
		if(umi_counts[umi] == 1):
			seq = read[1][READ_START : READ_END]
			#add to read counts
			read_counts.update([seq])
	return(umi_counts, read_counts)
			
	
def read_fastq(fq_file):
	"""
	Args:
		fq_file (string):
			Filename for a fastq file
	Yields:
		read (list): Data from file, 4 lines at a time
		offset (int): Character offset for this read
	"""	
	zcat = subprocess.Popen(
		#['gzcat', fq_file],
		['zcat', fq_file],
		stdout = subprocess.PIPE
	)
	for lines in itertools.zip_longest(*[zcat.stdout] *4):
		lines = [(l.rstrip()).decode('utf-8') for l in list(lines)]
		yield(lines)
	
def get_args():
	parser = argparse.ArgumentParser(
		description = 'This script deduplicates bindnseq data')
		
	#wow this is verbose
	parser.add_argument(
		'--reads', 
		type=str, 
		help='Fastq.gz file', 
		required=True)
	parser.add_argument(
		'--output_dir', 
		type=str, 
		help='Directory where outputs are written', 
		required=True)
	parser.add_argument(
		'--target_seq',
		type=str,
		help='Target sequence for the expt',
		default='NNNNNNGAGTGGCAGATATAGCCTGGTGGTTCAGGCAGATCGGAAGAGCAC')
	
	
	
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
	READ_START = 12
	READ_END = 36
	UMI_START = 0 
	UMI_END = 13
	REF_SEQ = args.target_seq
	run_all()


	