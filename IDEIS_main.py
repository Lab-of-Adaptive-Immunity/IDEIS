import os, sys, glob
import gzip
import pathlib
import argparse
import inspect
import shutil
import subprocess

########################################################################
#
# License:
#
# IDEIS (c) by Lab of Adaptive Immunity from Institute of Molecular Genetics of the Czech Academy of Sciences
#
# IDEIS and all its parts are licensed under a
# Creative Commons Attribution 4.0 International License.
#
# You should have received a copy of the license along with this
# work. If not, see <https://creativecommons.org/licenses/by/4.0/>.
#
########################################################################
#
# Acknowledgements:
#
# This project was supported by the National Institute of Virology and 
# Bacteriology (Programme EXCELES, LX22NPO5103 to Ondrej Stepanek) - 
# funded by the European Union - Next Generation EU.
#
########################################################################
#
# Description:
#
# IDEIS is a program for detection of CD45 isoforms in human and murine data.
# For user cases and instructions, please see python IDEIS_main.py -h 
#
########################################################################

########################################################################
# CLASSES FOR HANDLING EXONS, RULES and TRANSCRIPTS, PARSER OVERRIDE
#
# These classes are used to handle prep. steps for transcriptome creation.
# There is also redefinition of Parser class (library argparse) to override error handling.

# Used to store exons from fasta file
# - header: header of fasta starting with '>' and having two fields split by whitespace
class Exon:
  def __init__(self, header):
    header_split = header.strip().split()
    if header_split[0] == '>':
      header_split.pop(0)
    self.symbol = header_split[1]
    self.name = header_split[0].replace('>','')
    self.exon_seq = ''
    self.exon_length = 0
  
  # extends current exon
  # - extesion: sequence to extend exon by
  def extend_exon(self, extension):
    self.exon_seq = self.exon_seq + extension
    self.exon_length = len(self.exon_seq)

  # Repr. operator (print)
  def __repr__(self):
    return('Exon:\nExon symbol: %s\nExon name: %s\nExon seq.: %s\n---\n'%(self.symbol, self.name, self.exon_seq))

# Used to store rules to create transcripts from exons
# - line: line, as read from file
class Rule:
  def __init__(self, line):
    line_split = line.strip().split(':')
    self.main_exon = line_split[0]
    self.all_exons = line_split[1].split(',')
    self.transcript_name = line_split[2]
    self.gene_name = line_split[3]
  
  # Repr. operator (print)
  def __repr__(self):
    return('Rule:\nTranscript name: %s\nMain exon: %s\nAll exons: %s\n---\n'%(self.transcript_name, self.main_exon, self.all_exons))

# Used to store generated transcripts
# - rule: rule used to generate transcript, with appropriate transcript and gene name
# - exon_dixt: dictionary of exons that corresponds to given rule
# - read_length: length of reads used
# - cutoff_value: cutoff from length of reads for flanks
class Transcript:
  def __init__(self, rule, exon_dict, read_length, cutoff_value):
    self.rule = rule
    self.flanks_length = read_length - cutoff_value
    self.transcript_seq = ''

    # construct transcript sequence
    main_exon_pos = 0
    if(self.rule.main_exon):
      main_exon_pos = self.rule.all_exons.index(self.rule.main_exon)
    else:
      main_exon_pos = self.rule.all_exons.index('-')
      
    left_flank = ''
    left_flank_exons = main_exon_pos - 1
    while (len(left_flank) < self.flanks_length  and left_flank_exons >= 0):
      left_flank = exon_dict[self.rule.all_exons[left_flank_exons]].exon_seq + left_flank
      left_flank_exons -= 1
        
    if(len(left_flank) < self.flanks_length):
      print("Not enough exons for left flank. This is fine if left flank is absent (exon at or close to the start of complete transcript) but if unexpected it may negatively affect the performance of the software.")
      
    else:
      left_flank = left_flank[(len(left_flank) - self.flanks_length):]
      
    right_flank = ''
    right_flank_exons = main_exon_pos + 1
    while (len(right_flank) < self.flanks_length  and right_flank_exons < len(self.rule.all_exons)):
      right_flank += exon_dict[self.rule.all_exons[right_flank_exons]].exon_seq
      right_flank_exons += 1
       
    if(len(right_flank) < self.flanks_length):
      print("Not enough exons for right flank. This is fine if right flank is absent (exon at or close to the end of complete transcript) but if unexpected it may negatively affect the performance of the software.")
    
    else:
      right_flank = right_flank[:self.flanks_length]       

    if(self.rule.main_exon):
      self.transcript_seq = left_flank + exon_dict[self.rule.main_exon].exon_seq + right_flank
    else:
      self.transcript_seq = left_flank + right_flank

  # Comparison operator
  def __eq__(self, transcript2):
    return(self.transcript_seq == transcript2.transcript_seq)

  # Repr. operator (print)
  def __repr__(self):
    return('Transcript\nTranscript name: %s\nMain exon: %s\nAll exons: %s \nTranscript seuqence: %s \n---\n'%(self.rule.transcript_name, self.rule.main_exon, self.rule.all_exons, self.transcript_seq))

  # writes generated transcript to file open in output_fh
  # - output_fh: file handle where transcript will be written
  def write_transcript(self, output_fh):
    output_fh.write('>' + self.rule.transcript_name + ' ' + self.rule.gene_name + ' cDNA:protein_coding\n') # header first
    seq_formatted = [self.transcript_seq[i:i+60] for i in range(0, len(self.transcript_seq), 60)]
    for seq_line in seq_formatted:
      output_fh.write(seq_line + '\n')

# Overloads error behavior for parser to show help when error is thrown (mostly because not all required arguments were given)
# - argparse.ArgumentParser: - function to create parser from package argparse
class Parser_Error_Override(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('Error: %s\n'%message)
        self.print_help()
        sys.exit(2)
    
########################################################################
# CLASS END
  
########################################################################
# FUNCTIONS START
#


# Load configuration script as a dictionary of soft:path
# - config_file: path to configuration file. Details: 
#  -> 'PATHS:' test is required
#  -> write all paths bellow that line (all lines above will be ignored)
#  -> '#' to comment lines below 'PATHS:' line to ignore it
def load_config(config_file):
	config_fh = open(config_file, 'r')
	config_mem = config_fh.readlines()
	config_fh.close()
	config_paths = {}
	i = 0
	while i < len(config_mem):
		if config_mem[i].strip() == 'PATHS:':
			i += 1
			while i < len(config_mem):
				print(config_mem[i])
				if config_mem[i].strip() != '' and config_mem[i][0] != '#': 
					(soft, path) = config_mem[i].split('=')
					(soft, path) = (soft.strip(), path.strip())
					config_paths[soft] = path
				i += 1	
		
		i += 1	
	
	return config_paths


# Checks for specific soft necessary to run this script
# - soft: software name to check
# - config_paths: list of imported configuration paths from config file
def check_software(soft, config_paths):
	# first check config paths
	# practical in case you want to use local installation instead
	if soft in config_paths:
		print("Status of " + soft + ": FOUND in config file")
	elif shutil.which(soft) is not None:
		print("Status of " + soft + ": FOUND installed on computer")
		config_paths[soft] = soft
	else:
		sys.exit("The software " + soft + " was not found. Please install it or set a path to local installation.")


# Check installation of software
# - config_file: path to internal configuration file.
def check_all(config_paths):
	check_software('Rscript', config_paths)
	check_software('samtools', config_paths)
	check_software('salmon', config_paths)
	print("Software check successfully performed. \n\n")


# Builds the initial parser object to handle parameters and inputs
def create_parser():
	new_parser = Parser_Error_Override(description='Identifies CD45 isoforms of gene Ptprc/PTPRC.\n This program needs list of exons and read length to build transcritome that will be used to map reads to them. It will generate a table of cells containing CD45 isoforms.')
	new_parser.add_argument('reference', metavar='R',
							help='Path to Salmon reference/indices with prefix. The directory will have addad suffix made from used parameters to specify transcriptome; don\'t add this parameter to call. If it does not exist, it will be created in this directory as long as fasta file is supplied. For GRCm36, BALBcj and Human, fasta files are supplied.')
	new_parser.add_argument('bam_path', metavar='B',
							help='Path to bamfile containing the alignment of reads of dataset of interest. The index file .bam.bai is expected in the same directory with the same prefix.')
	new_parser.add_argument('output', metavar='O',
							help='Path to output directory where the results should be stored, including the name of the directory itself. ')		
	new_parser.add_argument('--species', '-g', dest = 'species', default='GRCm38',
							help='The species of analyzed organism. Three base values are \'GRCm38\' (default), \'BALBcj\', \'Homo-Sapiens\'. You can specify other value, but then you must provide .fasta with exons yourself. The value will be then used to name auto-generated transcriptome.')
	new_parser.add_argument('--gene-range', dest = 'gene_range',
							help='The range of Ptprc gene. For provided species the detection is automatic, otherwise it must be provided here in format \'chr:range_start-range_stop\'. Can be approximative as long as all key exons are contained. The format is chr:start-end.')
	new_parser.add_argument('--fasta', dest = 'fasta_path',
							help='Path to fasta with exons. Needed only if reference for given parameters does not exist and the species is not one of three base options.')		
	new_parser.add_argument('--rule-set', '-r', dest = 'rule_set_path',
							help='Path to files with rules to build trancriptome. Each supported species hats its own rule set. You need to define this only if you use other Species or want to use your custom ruleset. For other species you can still choose one of default rule sets by specifying path to respective directory.')	
	new_parser.add_argument('--transcriptome', '-t', dest = 'transcriptome_path',
							help='Path where to save transcriptome. By default they will be saved to the \'Transcriptomes\' directory located within the root of the software.')
	new_parser.add_argument('--read-length', '-l', dest = 'read_length', type = int,
							help='Length of reads used in mapping. In case of variable length maximum should be used. Will be auto-detected if not provided. Please keep in mind that in case of complicated designs with reads of different sizes it is better to split them into batches by read length and process them separately.')							
	new_parser.add_argument('--cutoff-value', '-c', dest = 'cutoff_value', type = int, default = 20,
							help='Cutoff from read length to create flanks around exons during creation of trancriptome fasta from exons. This forces the reads to map at least with part to given exon. Higher values mean bigger part of read has to map on sequence but discards more reads. Default is 20.')		
	new_parser.add_argument('--sequencing-type', '-s', dest = 'sequencing_type', default = '5-prime',
							help='Sequencing direction. It can be either 5-prime or 3-prime. Default is 5-prime.')
	new_parser.add_argument('--whitelist', '-w', dest = 'whitelist',
							help='Path to the whitelist of barcodes (barcode per row, no header). Mutually exclusive with expect-cells. Recommended to use if available.')
	new_parser.add_argument('--force-cells', '-f', dest= 'force_cells', default = '6000',
							help='Force number of cells. Ignored if whitelist is provided. Default 6000. If you have whitelist, then consider using it, as it is better option in general.')
	new_parser.add_argument('--data-set', '-d', dest = 'data_set',
							help='Path to .rds data set containing SeuratObject generated from same data as thos used to identify isoforms. If used, the whitelist will be used from these data, and new SeuratoOject will be created with new assay containing raw counts. If used, options --whitelist and --force-cells.')
	new_parser.add_argument('--protocol', dest = 'protocol', default = 'Chromium',
							help='Protocol used. Allowed options are \'Chromium\', \'ChromiumV3\' and \'DropSeq\', Default is \'Chromium\'.')
	new_parser.add_argument('--ncores', dest = 'ncores', type = int, default = 1,
							help='Number of cores to use by Alevin. Default is 1.')
	new_parser.add_argument('--version', action='version', version='2.0')
	return(new_parser)


# Converts the bam file to fastq -> uses bamtofastq provided by 10X
# bam is filtered only to relevant locus, for which its coordinates and species are needed
# - bam_path: path to .bam to convert
# - gene_range: range of genes if supported species is not analyzed
# - species: species to analyze (supported: GRCm38, BALBcj and Homo-Sapiens)
# - script_path: location of this script
# - config_paths: paths to required software
def convert2fastq(bam_path, gene_range, species, output_path, script_path, config_paths):
	bam_check = subprocess.Popen([config_paths['samtools'], 'view', bam_path], stdout=subprocess.PIPE)
	bam_head = subprocess.check_output(['head', '-n', '1'], stdin = bam_check.stdout, text=True)
	bam_chr = bam_head.split('\t')[2]
	
	# determine range
	if gene_range is not None:
	  print('Custom gene range detected: %s will be used for gene range'%gene_range)
	elif species == 'GRCm38':
		gene_range = bam_chr + ':138112191-138175126'
	elif species == 'BALBcj':
		gene_range = bam_chr + ':136956138-137020033'
	elif species == 'Homo-Sapiens':
		gene_range = bam_chr + ':198639226-198708261'
	else:
	  sys.exit('ERROR! Using custom reference without specified gene range. Please specify gene range in format chr:range_start-range_stop.')

	output_path = os.path.join(output_path, "")
	fastq_path = os.path.join(output_path, 'Fastqs')
	if os.path.exists(fastq_path):
		subprocess.run(['rmdir', fastq_path])
	
	print('Converting to fasta files ...')
	try:
	  process_bamtofastq = subprocess.run([os.path.join(script_path, 'bamtofastq-1.4.1'), bam_path, '--locus=' + gene_range, fastq_path],
	                                       stdout=subprocess.PIPE, text=True)
	except FileNotFoundError:
	  errmsg = 'bamtofastq-1.4.1 is missing. Please get it here: https://github.com/10XGenomics/bamtofastq/releases/tag/v1.4.1.'
	  sys.exit(errmsg)
  
	
	if 'bamtofastq error:' in process_bamtofastq.stdout:
	  if 'error creating output directory:' in process_bamtofastq.stdout:
	    errmsg = 'Error detected in bamtofastq: The target directory already exists. '
	    errmsg += 'Most likely this is because the pipeline was already run and completed beforehand. The software will continue but the previous results will be overwritten.'
	    print(errmsg)
	  else:
	    errmsg = 'Unspecified error detected in bamtofastq. This might be because of bad format of .bam. This software only supports .bam generated by Cell Ranger from 10X data. '
	    errmsg += 'If you downloaded these data from SRA and are sure the original pipeline used Cell Ranger, it might be possible that they were processed and consequently .bam '
	    errmsg += 'no longer contain required flags. In this case try to download original files, if they are available.'
	    sys.exit(errmsg)
	print(process_bamtofastq.stdout)

	return(fastq_path)


# Finds 50 first reads in file and detects read length based on them (average rounded down if inequal)
# - fastq_path: path to fastq
def detect_read_length(fastq_path):
  first_fastq_R2 = glob.glob(os.path.join(fastq_path, '*', '*R2*'))[0]
  readlens = []
  
  print(first_fastq_R2)
  with gzip.open(first_fastq_R2,'rt') as fastq_gzip:
    count = 0
    for line in fastq_gzip:
      count += 1
      if count%4 == 2: # in fastq second line is read itself
        readlens.append(len(line.strip()))
      if count == 198:
        break
  
  # return average length of read    
  return(int(sum(readlens)/len(readlens)))


# Builds up a fasta transcriptome file tailored to analyzed data
# - exon_fasta_path: fasta listing all exons, with header >A:B, where B matches symbol of exon
# - rule_set_path: transcripts to build; for details see Readme
# - script_path: path to this script
# - transcript_path: path where built fasta of transcript will be saved
# - read_length: length of read to be used; if not available will be estimated manually based on length of sequence (in construction)
# - cutoff_value: cutoff from length of read to delimit sequence around the read of interest
# - species: species to save transcriptome as
def build_transcriptome_fasta(exon_fasta_path, rule_set_path, script_path, transcript_path, species,
                              read_length, cutoff_value):

	if cutoff_value > read_length - 20:
	  cutoff_value = read_length - 20
	  print('Warning: Cutoff value longer than read length - 20. This is not advised as it will affect even reads fiably mapping directly to exons. Setting length of cutoff to read length - 20.')

	if exon_fasta_path is None:
	  fasta_dir_path = os.path.join(script_path, 'Fasta_exons')
	  if species == 'GRCm38':
		  exon_fasta_path = os.path.join(fasta_dir_path, 'Mus_Musculus.GRCm38.102.exons.fa')
	  elif species == 'BALBcj':
		  exon_fasta_path = os.path.join(fasta_dir_path, 'Mus_musculus_balbcj.BALB_cJ_v1.dna.PTPRC_A.fa')
	  elif species == 'Homo-Sapiens':
		  exon_fasta_path = os.path.join(fasta_dir_path, 'Homo_Sapiens.GRCh38.102.exons.fa')
	  else:
		  sys.exit('ERROR! Reference not found and fasta with exons not provided. Please provide references or exon fasta files to create reference.')
	else:
	  if not os.path.isfile(exon_fasta_path):
	    sys.exit('ERROR! Fasta file does not exist.')

	if rule_set_path is None:
	  rule_dir_path = os.path.join(script_path, 'Rule_sets')
	  if species == 'GRCm38':
		  rule_set_path = os.path.join(rule_dir_path, 'Mus_Musculus.GRCm38.ruleset.txt')
	  elif species == 'BALBcj':
		  rule_set_path = os.path.join(rule_dir_path, 'Mus_Musculus.GRCm38.ruleset.txt')
	  elif species == 'Homo-Sapiens':
		  rule_set_path = os.path.join(rule_dir_path, 'Homo_Sapiens.GRCh38.ruleset.txt')
	  else:
		  sys.exit('ERROR! Rule-set not found.')
	else:
	  if not os.path.isfile(rule_set_path):
	    sys.exit('ERROR! Rule-set does not exist.')    
	
	# first load all exons
	exon_dict = {} # list of Exons objects
	with open(exon_fasta_path, 'r') as exon_rh:
	  line = exon_rh.readline()
	  new_exon = ''
	  while(line):
	    if line[0] == '#':
	      continue
	    if '>' in line:       # header
	      if(new_exon != ''): # new sequence, dump previous to table
	        exon_dict[new_exon.symbol] = new_exon
	        new_exon = ''
          
	      new_exon = Exon(line)

	    elif new_exon  != '': # if there is exon, extend it
	      new_exon.extend_exon(line.strip())
    
	    line = exon_rh.readline()
	  
	  exon_dict[new_exon.symbol] = new_exon
  
  # load list of rules
	rule_list = []
	with open(rule_set_path, 'r') as rule_list_rh:
	  for line in rule_list_rh.readlines():
	    if line[0] == '#':
	      continue
	    rule = Rule(line)
	    rule_list.append(rule)
      
  # build transcriptomes
	with open(transcript_path, 'w') as transcript_fh:
	  transcripts = []
	  for rule in rule_list:
	    new_transcript = Transcript(rule, exon_dict, read_length, cutoff_value)
	    if(new_transcript not in transcripts): # not the same sequence
	      transcripts.append(new_transcript)
	      new_transcript.write_transcript(transcript_fh)
	    else:
	      print('Message: The transcript %s is similar, under given conditions, to other transcripts already in list and will be ignored.'%(new_transcript.rule.transcript_name))
	     
    
# Builds reference/index for salmon if does not exist already at designed path
# - reference: name of reference
# - transcriptome_path: path to .fasta of transcriptome
# - ncores: number of cores to use
# - config_paths: paths to required software
def create_reference(reference, transcriptome_path, ncores, config_paths):
	print('Creating salmon indices ...')

	# now create index
	if not os.path.exists(os.path.dirname(reference)):
		print('Creating directories with references ...')
		os.makedirs(os.path.dirname(reference))
	
	threads = ''
	if ncores != None:
		threads = ' -p '+str(ncores)
	
	try:
	  subprocess.run([config_paths['salmon'], 'index', '-t', transcriptome_path, '-i', reference, threads], check = True)
	
	except subprocess.CalledProcessError:
	  errmsg = '\nError during creation of reference index. This might be because of wrong format of exon files, wrong lengths of read or cutoff or wrong rule set.'
	  errmsg += 'Check the input files and try again.'
	  sys.exit(errmsg)
	  
  # finally, we need to create tgMap file	
	with open(os.path.join(reference, 'tgMap.tsv'), 'w') as tgmap_file:	
	  with open(transcriptome_path, 'r') as fasta_file:
	    for line in fasta_file:
	      if '>' in line[0]:	
	  	    transcript_name = line.split()[0].replace('>','')
	  	    gene_name = line.split()[1]
	  	    tgmap_file.write(transcript_name + '\t' + gene_name + '\n')

	print('Indexing Completed.')
	

# Creates all directories in output path:
# - output_path: where directory should be created
# - data_set: path to data_set. If not none, will create also directory wherenewly generated data_set will be stored
def create_out_hierarchy(output_path, data_set):
	if not os.path.exists(os.path.join(output_path, 'Results', 'Iso_Counts')):
		os.makedirs(os.path.join(output_path, 'Results', 'Iso_Counts'))
	if not os.path.exists(os.path.join(output_path, 'Results', 'Iso_Counts', 'raw_feature_matrix')):
		os.makedirs(os.path.join(output_path, 'Results', 'Iso_Counts', 'raw_feature_matrix'))	
	if not os.path.exists(os.path.join(output_path, 'Results', 'Datasets')) and data_set is not None:
		os.makedirs(os.path.join(output_path, 'Results', 'Datasets'))


# generates whitelist from rds file and returns path to it
# - data_set: .rds file to use for whitelisting of data
# - output_path: path to directory where generated whitelist will be stored 
# - script_path: path to scripts
# - config_paths: paths to required software
def get_whitelist_from_rds(data_set, output_path, script_path, config_paths):
  whitelist_path = os.path.join(output_path, 'Results', 'whitelist.csv')
  Rscript_path = os.path.join(script_path, 'Extract_whitelist_from_rds.R')
  try:
    subprocess.run([config_paths['Rscript'], '--no-restore', Rscript_path, data_set, whitelist_path], check = True)
  except subprocess.CalledProcessError:
    errmsg = 'Error during whitelist extraction from .rds. data set. This might be because the data set does not contain SeuratObject. '
    errmsg += 'Please check the format of the object or try again.'
    sys.exit(errmsg)
    
  return(whitelist_path)


# Performs counting using salmon allevin
# - fastq_path: path to directory with fastq files created by bamtofastq
# - output_path: path to directory where result hierarchy will be stored 
# - reference: path to salmon reference/indices
# - target_path: path directory where result will be stored
# - sequencing_type: direction of sequencing
# - script_path: path to scripts
# - data_set: .rds file to use for whitelisting of data
# - whitelist: whitelist to use; if not specified force_cells below is used instead
# - force_cells: number of expected cells
# - ncores: number of cores to use
# - config_paths: paths to required software
def run_allevin_analyses(fastq_path, output_path, reference, target_path, sequencing_type, script_path,
                         data_set, whitelist, force_cells, protocol, ncores, config_paths):
  
	print('Running cellranger analysis for ' + target_path + ' ...')
	
	# alevin requires list of fastq files. We list all R1 and R2 files
	fastqs_R1 = sorted(glob.glob(os.path.join(fastq_path, '*', '*R1*')))
	fastqs_R2 = sorted(glob.glob(os.path.join(fastq_path, '*', '*R2*')))
	
	library = 'ISF'
	if(sequencing_type == '3-prime'):
		library = 'ISR'

	cells_to_look = []
	if data_set is not None:
	  cells_to_look = ['--whitelist', get_whitelist_from_rds(data_set, output_path, script_path, config_paths)]
	elif whitelist is not None:
		cells_to_look = ['--whitelist', str(whitelist)]
	elif force_cells is not None:
		cells_to_look = ['--forceCells', str(force_cells)]
	
	protocol = '--chromium'
	if protocol == 'ChromiumV3':
		protocol = '--chromiumV3'
	elif protocol == 'Dropseq':
		protocol = '--dropseq'
		
	threads = '1'
	if ncores != None:
		threads = str(ncores)
    
	tgmap  = os.path.join(reference, 'tgMap.tsv')
	print([config_paths['salmon'], 'alevin', '-l', library, '-1'] + fastqs_R1 +['-2'] + fastqs_R2 + ['-i',
	                  reference, '-o', target_path, protocol, '--tgMap', tgmap] + cells_to_look + ['--dumpMtx', '-p', threads])
	try:
	  subprocess.run([config_paths['salmon'], 'alevin', '-l', library, '-1'] + fastqs_R1 +['-2'] + fastqs_R2 + ['-i',
	                  reference, '-o', target_path, protocol, '--tgMap', tgmap] + cells_to_look + ['--dumpMtx', '-p', threads],
	                  #'--bc-geometry', '1[1-16]', '--umi-geometry', '1[17-26]', '--read-geometry', '1[40-end]'],
	                  check = True)
    
	except subprocess.CalledProcessError:
	  errmsg = '\nThe application salmon alevin has thrown an error. The process was not finished. '
	  errmsg += 'Please check logs located in %s and try again. '%os.path.join(target_path, 'alevin')
	  errmsg += 'If you used --force-cells parameter, it might be too small or too high. Adjusting it may solve the issue.'
	  sys.exit(errmsg)

	return(os.path.join(output_path, 'Results', target_path))

########################################################################
# FUNCTIONS END HERE ^^
#
		
if __name__ == "__main__":

	# get path to directory containing this script, mostly to related files and software
	expath = inspect.getframeinfo(inspect.currentframe()).filename
	script_path = os.path.dirname(os.path.abspath(expath))
	script_path = os.path.join(script_path, "")
	
	# load configuration if necessary
	config_paths = load_config(script_path + 'config.txt')
	check_all(config_paths)
	print(config_paths)

	# read parameters
	new_parser = create_parser()
	main_parameters = new_parser.parse_args()
	
	# convert paths to absolute
	main_parameters.bam_path = os.path.abspath(main_parameters.bam_path)
	main_parameters.output = os.path.abspath(main_parameters.output)
	if(main_parameters.transcriptome_path is None):
	  main_parameters.transcriptome_path = os.path.join(script_path, 'Transcriptomes')
	  
	else:
	  main_parameters.transcriptome_path = os.path.abspath(main_parameters.transcriptome_path)
	
	if(main_parameters.fasta_path is not None):
	  main_parameters.fasta_path = os.path.abspath(main_parameters.fasta_path)

	if(main_parameters.rule_set_path is not None):
	  main_parameters.rule_set_path = os.path.abspath(main_parameters.rule_set_path)
	
	if(main_parameters.whitelist is not None):
	  main_parameters.whitelist = os.path.abspath(main_parameters.whitelist)

	if(main_parameters.data_set is not None):
	  main_parameters.data_set = os.path.abspath(main_parameters.data_set)
	  
	if main_parameters.ncores != None:
		if main_parameters.ncores < 0:
			sys.exit('ERROR! Invalid number of threads.')
					
	print(main_parameters, main_parameters.reference, main_parameters.bam_path)
	
	# 1.) Create hierarchy
	output_path = os.path.join(main_parameters.output, "")
	print('Creating output directory hierarchy.')
	create_out_hierarchy(output_path, main_parameters.data_set) # create all parent directories first

	# 2.) convert the above to fastq files for given CD45/PTPRC locus
	fastq_path = convert2fastq(main_parameters.bam_path, main_parameters.gene_range,
	                           main_parameters.species, output_path, script_path, config_paths)
	                           
	# 3.) detect length of read if not precised
	if main_parameters.read_length is None:
		main_parameters.read_length = detect_read_length(fastq_path)

	print('The read length that will be used: %d'%main_parameters.read_length)
	                           
	# 4.) reference creation (if does not exist for specified parameters yet)
	reference_suffix = '_%d_%d'%(main_parameters.read_length, main_parameters.cutoff_value)
	main_parameters.reference = os.path.abspath(main_parameters.reference + reference_suffix)
	
	# create transcriptome files and reference if reference is not created already (for given parameters)
	if not os.path.isdir(main_parameters.reference):
	  transcriptome_file_name =  '%s_%d_%d.fa'%(main_parameters.species, main_parameters.read_length, main_parameters.cutoff_value)
	  transcriptome_path = os.path.join(main_parameters.transcriptome_path, transcriptome_file_name)
	  build_transcriptome_fasta(main_parameters.fasta_path, main_parameters.rule_set_path, script_path, transcriptome_path, 
	                            main_parameters.species, main_parameters.read_length, main_parameters.cutoff_value)
	  create_reference(main_parameters.reference, transcriptome_path, main_parameters.ncores, config_paths)
  	
  # directories exist but are empty, remove them then continue with transcriptome and reference creation
	elif not os.listdir(main_parameters.reference):	
		subprocess.run(['rmdir', main_parameters.reference])
		transcriptome_file_name =  '%s_%d_%d.fa'%(main_parameters.species, main_parameters.read_length, main_parameters.cutoff_value)
		transcriptome_path = os.path.join(main_parameters.transcriptome_path, transcriptome_file_name)
		build_transcriptome_fasta(main_parameters.fasta_path, main_parameters.rule_set_path, script_path, transcriptome_path, 
		                          main_parameters.species, main_parameters.read_length, main_parameters.cutoff_value)
		create_reference(main_parameters.reference, transcriptome_path, main_parameters.ncores, config_paths)

	# 5.) run alevin to map all sequences to isoform transcripts
	fastq_path = run_allevin_analyses(fastq_path, output_path, main_parameters.reference, os.path.join(output_path, 'Results', 'alevin'), 
	                                  main_parameters.sequencing_type, script_path, main_parameters.data_set, main_parameters.whitelist, 
	                                  main_parameters.force_cells, main_parameters.protocol, main_parameters.ncores, config_paths)
	                                  
	# 6.) calls script to build final complete table using 10X whitelists, then export them as .mtx and .rds. If initial .rds with Seurat Object was provided, 
	#     it will also export data there
	print('Generating final matrix and aggregating all data ...')

	Rscript_path = os.path.join(script_path, 'Generate_final_mtx.R')	
	input_path = os.path.join(output_path, 'Results')
	count_output_path = os.path.join(output_path, 'Results', 'Iso_Counts')
	whitelist = str(main_parameters.whitelist)
	data_set = str(main_parameters.data_set)
	data_set_path = os.path.join(output_path, 'Results', 'Datasets', 'dataset_with_cd45.rds')
	
	try:
	  subprocess.run([config_paths['Rscript'], '--no-restore', Rscript_path, input_path, count_output_path, 
	                  main_parameters.sequencing_type, whitelist, data_set, data_set_path], check = True)
	
	except subprocess.CalledProcessError:
	  errmsg = 'Error during creation of final files. This might be problem with whitelists or R.'
	  errmsg += 'The most likely cause is the lack of packages for default installation of R or one precised in config file.'
	  errmsg += 'Please check your used installation of R, its packages and try again.'
	  sys.exit(errmsg)
	
	print('Done!')
