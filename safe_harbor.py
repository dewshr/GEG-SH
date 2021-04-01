import pandas as pd
import argparse
import os
import sys
script_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(script_path+'/src')
from utils import *
from loguru import logger
import numpy as np


############### parameters for the program #################
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input',default = None, help='input file with mobile element id')
parser.add_argument('-t','--thresh', default = None, help= 'fdr threshold to filter out variants')
parser.add_argument('-eqtl','--eqtl_genes', default = False, help= 'eQTL genes')
parser.add_argument('-tad', '--tad_domain', default = None, help='custom tad domain file, else will use default file provided')
parser.add_argument('-rr', '--repressive_region', default ='./data/blood_repressive_marks.bed', help ='bed file containing regions with repressive mark')
parser.add_argument('-ar', '--active_region', default = None, help ='bed file containing regions with active transcription mark')
parser.add_argument('-gd', '--gene_density', default = None, help='gene density in the tad domains, mean gene density will be used as default')
parser.add_argument('-o', '--output', default='./results', help ='ouput folder name')
parser.add_argument('-af1', '--allele_freq1', default = None, help ='lower limit threshold of allele frequency')
parser.add_argument('-af2', '--allele_freq2', default = 0.9, help ='upper limit threshold of allele frequency')
parser.add_argument('-hic', '--hic_interaction', default = './data/blood_hic_interaction.bed', help='chromatin interaction bed file')
parser.add_argument('-l', '--nearby_cancer_genes', default = 50000, help='any variant with oncogenes or tumor repressor genes 50kb upstream or downstream will be removed')
parser.add_argument('-fname','--file_name', default='result.csv', help='output file name')

args = parser.parse_args()



if not os.path.exists(args.output):
	os.makedirs(os.path.join(args.output, 'temp_files'))

logger_path = os.path.abspath(args.output)
logger.add(logger_path+'/variant_filter{time}.log', rotation='10 MB')
logger.level("WARNING", color="<bold><red>")

dir_ = os.path.abspath(args.output) 


if args.input == None:
	logger.error('input file not provided. Exiting....')
	parser.print_help()
	sys.exit(1)

dist = int(args.nearby_cancer_genes)/1000


# tumor repressor and oncogenes list
tumor_repressor_gene = pd.read_csv('./data/oncogenes_and_tumor_repressor_genes.bed',sep='\t')
tumor_repressor_genes_list = tumor_repressor_gene.loc[3].tolist() #column 3 has name of genes


#list of dosage sensitive genes
dosage_sensitive_genes = list(filter(None, open('./data/dosage_sensitive_genes.txt','r').read().split('\n')))


# getting the list of TADS that consist of onco genes and tumor supressor genes
if args.tad_domain == None:
	cancer_tad = pd.read_csv('./data/cancer_genes_tad.bed',header=None, sep='\t')
	cancer_tad['tad_name'] = cancer_tad[5] + '-'+cancer_tad[6].map(str)+'-'+cancer_tad[7].map(str)
	
	#list of tad domains with tumor repressor or oncogenes
	cancer_tad_list = cancer_tad['tad_name'].tolist()

	gene_density = pd.read_csv('./data/gene_density_all_tad.csv')
	gene_density = gene_density.set_index('name')

	# tad domain information for genes
	genes_tad = pd.read_csv('./data/genes_tad.bed', sep='\t', header=None).iloc[:,[3,5,6,7]]
	genes_tad.columns = ['gene','chr','start','end']
	genes_tad['name'] = genes_tad.chr + '-'+ genes_tad['start'].map(str)+'-'+genes_tad['end'].map(str)
	genes_tad = genes_tad.groupby('gene').agg(lambda x: list(set(list(x))))
	

else:
	# TAD gene density calculation
	os.system('bedtools intersect -a ./data/sorted_gene_annotation.bed -b {} -wb > {}/genes_tad.bed'.format(args.tad_domain,dir_))
	os.system('bedtools intersect -a ./data/oncogenes_and_tumor_repressor_genes.bed -b {} -wb > {}/cancer_genes_tad.bed'.format(args.tad_domain,dir_))

	cancer_tad = pd.read_csv('./results/cancer_genes_tad.bed',header=None, sep='\t')
	cancer_tad['tad_name'] = cancer_tad[5] + '-'+cancer_tad[6].map(str)+'-'+cancer_tad[7].map(str)
	
	#list of tad domains with tumor repressor or oncogenes
	cancer_tad_list = cancer_tad['tad_name'].tolist()

	gene_density = pd.read_csv('./results/genes_tad.bed', sep='\t', header=None).iloc[:,[3,5,6,7]]
	gene_density.columns = ['gene','chr','start','end']
	gene_density['length'] = abs(gene_density['start'] - gene_density['end'])
	gene_density['name'] = gene_density.chr + '-'+ gene_density['start'].map(str)+'-'+gene_density['end'].map(str)
	
	# tad domain information for genes
	genes_tad = gene_density.copy()
	genes_tad = genes_tad.groupby('gene').agg(lambda x: list(set(list(x))))
	
	# tad domain with gene density information
	gene_density = gene_density.groupby('name').agg(lambda x: list(set(list(x))))
	gene_density['density'] = gene_density.apply(lambda x: (len(x.gene)/x.length[0])*1000000, axis=1)
	#gene_density.to_csv('/Users/dshresth/Downloads/safe_harbor/data/gene_density_all_tad.csv')



# calculating gene density if not provided by user
if args.gene_density == None:
	gd = round(np.mean(gene_density.density),2)
else:
	gd = float(args.gene_density)



#################################### Reading in input files and getting TAD domain information for variants ###########################

logger.info('reading in the input file: {}'.format(args.input))
input_data = pd.read_csv(args.input, sep='\t')	# reading in SNP eQTL data
columns = [x.lower() for x in input_data.columns.tolist()]
input_data.columns = columns
input_data['chr'], input_data['start'], input_data['stop'] = zip(*input_data.apply(get_bed_file, axis=1))

columns_names = input_data.columns


# extending the sequence to look for tumor repressor or oncogenes nearby
input_data['extended_start'] = input_data['start'].apply(lambda x: int(x)- int(args.nearby_cancer_genes) if int(x)>int(args.nearby_cancer_genes) else 0)
input_data['extended_stop'] = input_data['stop']+ int(args.nearby_cancer_genes)
input_data = input_data.sort_values(by=['chr','start'], ascending=[True, True])
#input_data['tissue'] ='.'

logger.info('creating bed file from the input data')
input_data.loc[:,['chr','start','stop','id']].to_csv(os.path.join(dir_,'sorted_variant_coordinates.bed'), header=False, sep='\t',index=False)		# writing bed file for snp coordinates
input_data.loc[:,['chr','extended_start','extended_stop','id']].to_csv(os.path.join(dir_,'sorted_extended_variant_coordinates.bed'), header=False, sep='\t',index=False)

logger.info('looking for nearby tumor repressor or oncogenes')
os.system('bedtools intersect -a {}/sorted_extended_variant_coordinates.bed -b ./data/oncogenes_and_tumor_repressor_genes.bed -wb > {}/oncogenic_tumor_repressor_variant.bed'.format(dir_,dir_))

variant_nearby_cancer = pd.read_csv('{}/oncogenic_tumor_repressor_variant.bed'.format(dir_), header=None, sep='\t')
variant_nearby_cancer.columns = ['chr','start','stop','id','c_chr','c_start','c_stop','genes','gene_id']
variant_nearby_cancer_list = variant_nearby_cancer['id'].tolist()


# variants overlap region with TAD domain
if args.tad_domain == None:
	os.system('bedtools intersect -a {}/sorted_variant_coordinates.bed -b ./data/merged_gm12878.bed  -wb > {}/variant_tad.bed'.format(dir_,dir_)) 
else:
	os.system('bedtools intersect -a {}/sorted_variant_coordinates.bed -b {} -wb > {}/variant_tad.bed'.format(dir_,args.tad_domain, dir_)) 




#logger.info('reading in the variants with TAD overlap information')
variant_tad = pd.read_csv(dir_+'/variant_tad.bed',header=None, sep='\t').iloc[:,[3,4,5,6]] # only taking the columns that represents snp, chr, start and end
variant_tad.columns = ['snp','chr','start','end']
variant_tad['tad_name'] = variant_tad['chr'] + '-'+variant_tad['start'].map(str)+'-'+variant_tad['end'].map(str)
variant_tad = variant_tad.groupby('snp').agg(lambda x: list(x))
variant_tad = pd.DataFrame(variant_tad.iloc[:,3])
#variant_tad['snp'] = variant_tad.index


logger.info('assigning tad domain information to variants')
input_data['tad_name'] = input_data['id'].apply(lambda x: get_tad_info(x, variant_tad))


logger.info('checking common TAD domain betweem variants and  tumor repressor/oncogenes')
input_data['same_cancer_tad'] = input_data['tad_name'].apply(lambda x: True if len(set(x)&set(cancer_tad_list))>0 else False) # checking if any of the tads overlap







####################################### Calculating gene density for variants associated tad region #####################################################

logger.info('getting gene density for variants associated TAD domain')

input_data['gene_density'] = input_data['tad_name'].apply(lambda x: calculate_gene_density(x, gene_density))







####################################### Getting hic-promoter interaction information for variants ##################################################


logger.info('running bedtools to get the information regarding overlap with hic-promoter interaction region')
os.system('bedtools intersect -a {}/sorted_variant_coordinates.bed -b {} -wao > {}/variant_promoter_interaction.bed'.format(dir_, args.hic_interaction, dir_))

variant_hic_promoter = pd.read_csv(dir_+'/variant_promoter_interaction.bed',sep='\t', header=None)
variant_hic_promoter.columns = ['chr','start','end','snp','hic_chr','hic_start','hic_end','hic_interacted_gene', 'overlap']
variant_hic_promoter = variant_hic_promoter.drop_duplicates()

logger.info('checking if the interacted gene and variants are in same TAD domain')
variant_hic_promoter['common_tad'] = variant_hic_promoter.apply(lambda x: check_tad(x.snp, x.hic_interacted_gene, variant_tad, genes_tad), axis=1)

#variant_hic_promoter = variant_hic_promoter[variant_hic_promoter['common_tad']==0] #removed same tad interaction

logger.info('checking if interacted gene falls in dosage_sensitive_genes, 1 or more value assigned depending on number of interaction, else if there is no interacted gene or does not fall in dosage_sensitive_genes then assignn 0')
variant_hic_promoter['label'] = variant_hic_promoter['hic_interacted_gene'].apply(lambda x: 1 if x in dosage_sensitive_genes else 0)


variant_hic_promoter_ = variant_hic_promoter.iloc[:,[3,-2,-1,7]].groupby('snp').agg(lambda x: list(x))
variant_hic_promoter_['dosage_sensitive_interaction'] = variant_hic_promoter_['label'].apply(lambda x: sum(x))
variant_hic_promoter_['common_tad_count'] = variant_hic_promoter_['common_tad'].apply(lambda x: sum(x))



all_data = input_data.merge(variant_hic_promoter_.loc[:,['hic_interacted_gene','common_tad_count','dosage_sensitive_interaction']], left_on='id', right_on=variant_hic_promoter_.index)

logger.info('checking if any of interacted genes are tumor repressor or oncogenes')
all_data['hic_interacted_gene_test'] = all_data['hic_interacted_gene'].apply(lambda x: filter_tumor_repressor_genes(x, tumor_repressor_genes_list))







######################################## Checking for heterochromatin region and nearby genes ########################


logger.info('running bedtools to get the information regarding overlap with heterochromatin region')
os.system('bedtools intersect -a {}/sorted_variant_coordinates.bed -b {} -wb > {}/variant_h3k27me3.bed'.format(dir_,args.repressive_region,dir_))

heterochromatin = pd.read_csv(dir_+'/variant_h3k27me3.bed',sep='\t', header=None)
heterochromatin_region = heterochromatin.iloc[:,3].tolist() # variants id overlapping with repressive region

all_data['repressive_region'] = all_data['id'].apply(lambda x: check_active_status(x, heterochromatin_region))



logger.info('checking for variants with nearby oncogenes or tumor repressor genes')

#print(filter7.head())
all_data['nearby_cancer_genes'] = all_data['id'].apply(lambda x: filter_nearby_cancer_genes(x, variant_nearby_cancer_list))
all_data['nearby_cancer_gene_names'] = all_data['id'].apply(lambda x: get_nearby_genes(x, variant_nearby_cancer))


	

######################################## Checking for active chromatin region and nearby genes ########################

if args.active_region != None:
	logger.info('running bedtools to get the information regarding overlap with active transcription region')
	os.system('bedtools intersect -a {}/sorted_variant_coordinates.bed -b {} -wb > {}/variant_active.bed'.format(dir_,args.active_region,dir_))
else:
	logger.info('running bedtools to get the information regarding overlap with active transcription region using default file "blood_active_transcription_marks.bed".')
	os.system('bedtools intersect -a {}/sorted_variant_coordinates.bed -b  ./data/blood_active_transcription_marks.bed -wb > {}/variant_active.bed'.format(dir_,dir_))

active_region = pd.read_csv(dir_+'/variant_active.bed',sep='\t', header=None)
active_region_list = active_region.iloc[:,3].tolist()

logger.info('tagging ME overlapping with active chromatin region as True')
#filter7['active_region'] = filter7['id'].apply(lambda x: check_active_status(x, active_region_list))
all_data['active_region'] = all_data['id'].apply(lambda x: check_active_status(x, active_region_list))





################################################ Filtering by FDR if provided by user ######################################################
if args.thresh != None:
	logger.info('filtering by FDR > {}'.format(args.thresh))

	all_data['fdr_test'] = all_data['fdr'].apply(lambda x: fdr_filter(str(x).split(','), args.thresh))
	#logger.info('filter6 shape (after removing ME with FDR) > {} : {}'.format(args.thresh, filter6.shape))
	filtered_data = all_data[all_data['fdr_test']== True]
else:
	if 'fdr' in columns:
		logger.warning('\n\n----------- !!! Warning !!! -------------------\n\n there seems to be fdr column in input file, but fdr parameter is not used, so filtration steps is carried out without using fdr column.\n\n-----------------------------------------------\n')
	filtered_data = all_data.copy()






########################################### Filtering by allele frequency and eQTL genes if provided###############

if args.allele_freq1 != None:
	logger.info('filtering by allele frequency > {}'.format(args.allele_freq))
	#filter5 = filter4[filter4['AF'] > float(args.allele_freq)]
	all_data['af'] = all_data['af'].apply(lambda x: get_round_value(x))
	filtered_data = filtered_data[(filtered_data['af']> float(args.allele_freq1))&(filtered_data['af']< float(args.allele_freq2))]
	#logger.info('filter5 shape (after removing ME with AF) > {} : {}'.format(args.allele_freq, filter5.shape))
else:
	if 'af' in columns:
		logger.warning('\n\n----------- !!! Warning !!! -------------------\n\n there seems to be allele frequency column in input file, but "af" parameter is not used, so filtration steps is carried out without using "af" column.\n\n-----------------------------------------------\n')



if args.eqtl_genes != False:
	all_data['eQTL_test'] = all_data['eqtl'].apply(lambda x: filter_tumor_repressor_genes(x.split(','), tumor_repressor_genes_list))

	eqtl_passed_variants = all_data[all_data['eQTL_test']==False]['id'].tolist()
	#print(eqtl_passed_variants)

	logger.info('keeping only those variants if EQTL gene is not tumor repressor or oncogenes')
	filtered_data = filtered_data[filtered_data['id'].isin(eqtl_passed_variants)]




######################################### Filteration steps: applying all remaining filters ############################


logger.info('removing variants in same TAD domain as tumor repressor or oncogenes')
logger.info('removing variants with TAD domain having gene density < {}'.format(gd))
logger.info('removing variants interacting with dosage sensitive genes or tumor repressor or oncogenes or any interaction with genes in same tad domain')
logger.info('removing variants with nearby tumor repressor or oncogenes')

filtered_data = filtered_data[(filtered_data['same_cancer_tad']==False)&(filtered_data['gene_density']<gd)&(filtered_data['dosage_sensitive_interaction']==0)&(filtered_data['repressive_region']==False)&(filtered_data['nearby_cancer_genes']==False)&(filtered_data['common_tad_count']==0)&(filtered_data['hic_interacted_gene_test']==False)]





########################################### Generating output files ##############################


logger.info('generating output files....')
filtered_data['ucsc_link'] = 'http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position='+filtered_data['position'].map(str)
filtered_data.loc[:,['id','position','ucsc_link']].to_csv(dir_+'/filtered_'+args.file_name, index=False)

all_data['passed_all_filter'] = all_data['id'].apply(lambda x: check_final_list(x, filtered_data['id'].tolist()))
all_data['tad_name'] = all_data['tad_name'].apply(lambda x: ','.join(x))
all_data['hic_interacted_gene'] = all_data['hic_interacted_gene'].apply(lambda x: ','.join(x))

all_data.rename(columns={'gene_density':'gene_density (genes per million tad)', 'hic_interacted_gene_test':'hic_interacted_genes (oncogenic or tumor repressor)','nearby_cancer_genes':'nearby_cancer_genes ({}kb)'.format(dist)}, inplace=True)
all_data.to_csv(dir_+'/'+args.file_name)

logger.success('\n\nThe total number of variant id after removing dosage sensitive genes:  {}\n'.format(filtered_data.shape[0]))


try:
	os.system('mv {}/*.bed {}/temp_files/'.format(dir_, dir_))
	os.system('mv {}/*.log {}/temp_files/'.format(dir_, dir_))
except IOError:
	logger.info('either bed files or log files are not found, please check the folder again')



