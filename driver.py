import os
import sys

#input gathering
if len(sys.argv) < 3:
	print("Usage: python driver.py [input FASTA] [output file]")
	sys.exit()

seqsfile = sys.argv[1]
outputfile = sys.argv[2]

# file titles
blast_output = 'blast_outputs/'+seqsfile+'_BLAST.xml'
#prosite_raw_xml = 'prosite_outputs/' + seqsfile+'_prosite.xml'
prosite_raw_xml = 'prosite_outputs/' + seqsfile+'_prosite.xml'
prosite_ids = 'prosite_outputs/' + seqsfile+'.vals'
inter_table_file = 'intermediate_tables/' + seqsfile + '.tab'

#Start running scripts
print("BLASTING")
os.system('py blaster.py '+seqsfile)
print("PARSING DATA...")
os.system('py init_table.py ' + blast_output + ' ' + inter_table_file)
print("Getting final values...")
os.system('py table_extender.py ' + inter_table_file + ' ' + prosite_ids+ ' '+outputfile)
