from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import sys
import argparse
import re
import os
import requests
import time
def Blastseqs(file1, file2):
	#file1 = sys.argv[1]
	#file2 = sys.argv[2]
	with open(file1, "r") as infile:
		for rec in SeqIO.parse(infile, "fasta"):
			seqs = rec.seq.split('>')
			for peptide in seqs:
				amino = ">" + str(rec.description) + "\n"  + str(peptide)
				result_handle = NCBIWWW.qblast("blastp", "swissprot", amino )
				blast_result = open(file2, "a")
				blast_result.write(result_handle.read())
				blast_result.close()
				result_handle.close()
				
def PrositeSearch(file1,file2):
#	base = "https://prosite.expasy.org/cgi-bin/prosite/PSScan.cgi?seq="
	base = "https://prosite.expasy.org/cgi-bin/prosite/PSScan.cgi"
	#o#utfilepath= "prosite_outputs/"
	#f = open(file2, 'w')
	count = 1
	#endofurl = "&output=xml"
	with open(file1, "r") as infile:
		for rec in SeqIO.parse(infile, "fasta"):
			seqs = rec.seq.split('>')
			for peptide in seqs:
				#url = base + peptide + endofurl
				p = str(peptide)
				data = {'seq':p, 'output':'xml'}
				r = requests.post(base, data)
				#print(p)
				outpt = r.content.decode('utf-8')
				#currfile = outfilepath + str(count) + '.txt' 
				prosite_result = open(file2, "a")
				prosite_result.write(outpt)
				prosite_result.close()
				#with open(currfile,'w') as f:
				#	print(outpt, file = f)
				#count+=1
				#@print(peptide)
	return file2
	
def processProsite(file1,file2):
	with open(file1) as f:
		rawdata = f.readlines()
	results = list()
	print(rawdata)
	#data = rawdata.split('</matchset>')[0:-1]
	for i in rawdata:
		if '<signature_ac>' in i:
			curr = i .split('<signature_ac>')[1].split('<')[0].strip()
			results.append(curr)
		else:
			results.append('.')
	#print(results)
	with open(file2,'w') as f:
		for i in results:
			print(i,file=f)
			
	
	
def pfamSearch(file1,file2):
	base = "https://pfam.xfam.org/search/sequence"
	with open(file1, "r") as infile:
		for rec in SeqIO.parse(infile, "fasta"):
			seqs = rec.seq.split('>')
			for peptide in seqs:
				#url = base + peptide + endofurl
				p = str(peptide)
				data = {'seq':p, 'output':'xml'}
				r = requests.post(base, data)
				#print(p)
				outpt = r.content.decode('utf-8')
				results_link =base+'/resultset/' + outpt.split('job_id="')[1].split('"')[0]
				time.sleep(15)
				#print(results_link)
				s = requests.post(results_link,{'output':'xml'})
				finaloutpt = s.content.decode('utf-8')
				pfam_result = open(file2, "a")
				pfam_result.write(outpt)
				pfam_result.close()

				#currfile = outfilepath + str(count) + '.txt' 
				
				#with open(currfile,'w') as f:
				#	print(finaloutpt, file = f)
				#count+=		1
				##return



seqsfile = sys.argv[1]	
blast_output = 'blast_outputs/_'+seqsfile+'_BLAST.xml'
prosite_raw_xml = 'prosite_outputs/' + seqsfile+'_prosite.xml'
prosite_ids = 'prosite_outputs/' + seqsfile+'.vals'



print("BLASTING...")
Blastseqs(seqsfile, blast_output)
print("QUERYING PROSITE...")
PrositeSearch(seqsfile, prosite_raw_xml)
print('Processing')
processProsite(prosite_raw_xml,prosite_ids)

