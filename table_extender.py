from seqs import * 
import requests
#from bioservices.kegg import KEGG
import time
from bioservices import QuickGO
import sys
k = KEGG()
HEADER = "Query\tDescription\tE-val\tPfam_Acc\tPfam\tProsite\tKEGGPathways\tGO\tComments"
#need to do it like this i think
#shell% curl -LH 'Expect:' -F seq='<a.seq' -F output=xml 'https://pfam.xfam.org/search/sequence' 
def searchPFAM(acc):
		#if curr.pfam != '.':
		#	return
		base = 'https://pfam.xfam.org/protein?output=xml&acc='
		url = base + acc
		response = requests.get(url)
		xmldata = response.content.decode('utf-8')
		if 'match accession' in xmldata:
			data = xmldata.split('<match accession=')[1].split("type")[0]
			data = data.split('"')
			return data[1],data[3]
		else:
			#print('wa')
			return '.','.'
#same idea
def findGO(accession):
	base = 'https://www.uniprot.org/uniprot/'
	url = base + accession+'.xml'
		#print(url)
	response = requests.get(url)
	#print(response.content)
	xmldata = response.content.decode('utf-8')
	splitstr = '<dbReference type="GO" id="'
	if splitstr in xmldata:
		chopped = xmldata.split(splitstr)[1:]
		vals = list()
		for i in chopped:
			if 'GO' in i:
				curr = i.split('"')[0]
				if 'GO' in curr:
					
					vals.append(curr)
		
		return vals
	else:
		return []


table_file = sys.argv[1]
prosite_file = sys.argv[2]
output_file = sys.argv[3]


c = controller(table_file,output_file)
#c.writeOut(True)

#pre computed prosite vals. Line # = query #
with open(prosite_file, "r") as infile:
	prositevals = infile.readlines()

c.prepareData()

#Calculate prosite vals,pfam and go data here
for i in c.seqs:
	#id is just query num atm , so multiple seqs may have same id since i save mult alignments for some queries
	
	
	# GOnna eventually uncomment this, just saving time since these values are filled atm
	#i.prosite = vals[int(i.id)-1].strip()	
	#if i.pfam == '.' and i.pfam != '':
	#if i.pfam == '.' or i.pfam_acc == '.':
	i.pfam_acc,i.pfam=searchPFAM(i.accession)
	i.go = findGO(i.accession)
	curr = k.get(i.kegg_id)
	kdata = k.parse(curr)
		#nums.append(i.id)
		#if has pathways add them to pathways
	if 'PATHWAY' in kdata:
	
		i.pathways = list(kdata['PATHWAY'].values())

	else:
		i.pathways = ('.')
bestseqs = c.getBestSeqs()

''' Gets slimGO values' '''
slimgo = {}
for i in bestseqs:
	for j in i.go:
		slimgo[j] = ''
with open('goslim_generic.obo.txt') as f:
	bigd = f.read().split('[Term]')
for i in bigd:
	id = i.split('id: ')
	if len(id) > 1: 
		id = id[1].split('\n')[0]
		if id in slimgo:
		#name = i.split('
			if 'name' in i:
			
			#slimgo[id] = i.split('name:')[1].split('\n')[0]
				temp = i.split('name:')
				if len(temp) > 1:
					slimgo[id] = temp[1].split('\n')[0]
	

s = QuickGO() 
prosite_val_index = 0
with open(output_file,'w') as f:
	print(HEADER, file =f)
	for i in bestseqs:
		#i.getKEGG()
		i.prosite = prositevals[prosite_val_index].strip()
		prosite_val_index +=1
		newgo =list()
		for j in i.go:
			if j in slimgo and slimgo[j] != '':
				newgo.append(j+'-'+slimgo[j])
			
			else:
				#print ("Manual search GO")
				goname = s.gosearch(j) 
				if goname['numberOfHits'] > 0:
					goname= goname['results'][0]['name']
					goname = j + '-' + goname
					newgo.append(goname)
		
		i.go = newgo
		print(i.id, i.description, i.eval, i.pfam_acc, i.pfam, i.prosite,','.join(i.pathways),','.join(i.go),sep='\t',file=f)
	
#c.writeOut()