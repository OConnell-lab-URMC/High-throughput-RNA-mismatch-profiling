"""
Batch script to run
DeduplicateMultipleTargets.py
"""

import sys
import os
import subprocess
import itertools

raw_data = sys.argv[1]
output_dir = sys.argv[2]
dedupe = sys.argv[3]

if not os.path.exists(output_dir):
	os.mkdir(output_dir)


#get files as list
fq_files = list(os.listdir(raw_data))
get_prefix = lambda txt: txt.split('_')[0]

samples = []
for f1, f2 in itertools.combinations(fq_files, 2):
	if not f1.endswith('.fastq.gz'):
		pass
	if not f2.endswith('.fastq.gz'):
		pass
	
	#else
	if(get_prefix(f1) == get_prefix(f2)):
		samples.append((f1, f2))


for (f1, f2) in samples:
	sample_dir = '%s/sample_%s' % (output_dir, get_prefix(f1))
	fnames = '%s%s,%s%s' % (raw_data, f1, raw_data, f2)
	
	dedupe_echo = subprocess.Popen(['echo',
		'python3',
		dedupe,
		'--reads',			fnames,
		'--output_dir',	sample_dir])
	dedupe_echo.communicate()[0]
	
	
	dedupe_cmd = subprocess.Popen([
		'python3',
		dedupe,
		'--reads',			fnames,
		'--output_dir',	sample_dir])
	dedupe_cmd.communicate()[0]
	
	