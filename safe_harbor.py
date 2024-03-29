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
parser.add_argument('-t','--thresh', default = None, type =float, help= 'fdr threshold to filter out variants')
#parser.add_argument('-eqtl','--eqtl_genes', default = False, help= 'eQTL genes')
parser.add_argument('-tad', '--tad_domain', default = None, help='custom tad domain file, else will use default file provided')
parser.add_argument('-rr', '--repressive_region', default = script_path+'/data/blood_repressive_marks_state.bed', help ='bed file containing regions with repressive mark')
parser.add_argument('-ar', '--active_region', default = script_path+'/data/blood_active_transcription_marks_state.bed', help ='bed file containing regions with active transcription mark')
parser.add_argument('-gd', '--gene_density', default = None, help='gene density in the tad domains, mean gene density will be used as default')
parser.add_argument('-o', '--output', default='./results', help ='ouput folder name')
parser.add_argument('-br','--blacklist_region', default =None, help='bed file with the coordinates that the user does not want to include in output')
parser.add_argument('-af', '--allele_freq', default = None, type =float, help ='allele frequency threshold for the variant')
parser.add_argument('-hic', '--hic_interaction', default = script_path+'/data/blood_hic_interaction.bed', help='chromatin interaction bed file')
parser.add_argument('-l', '--nearby_cancer_genes', default = 0, help='default = 0, takes number as input representing the distance user wants to check for oncogenes or tumor repressor genes in upstream or downstream of the pMEI')
parser.add_argument('-fname','--file_name', default='result.csv', help='output file name')

args = parser.parse_args()



if not os.path.exists(args.output):
	os.makedirs(os.path.join(args.output, 'temp_files'))


logger_path = os.path.abspath(args.output)
logger.add(logger_path+'/safe_harbor.log', rotation='10 MB',mode='w')
logger.level("WARNING", color="<bold><red>")

logger.info('COMMAND USED:\npython ' + ' '.join(sys.argv) +'\n')
#sys.exit()

dir_ = os.path.abspath(args.output) 


if args.input == None:
	logger.error('input file not provided. Exiting....')
	parser.print_help()
	sys.exit(1)


dist = int(args.nearby_cancer_genes)/1000


info = "\n\nParameter information:\n\n" + "FDR threshold:\t{}\n".format(args.thresh) + "Variant Allele frequency:\t{}\n".format(args.allele_freq)

if args.tad_domain == None:
	info = info + "tad domain file used:\t{}\n".format(script_path+'/data/merged_gm12878.bed')
else:
	info = info + "tad domain file used:\t{}\n".format(args.tad_domain)

info = info + "Repressive region:\t{}\n".format(args.repressive_region) + "Active region:\t{}\n".format(args.active_region) + "Chromatin interaction data:\t{}\n".format(args.hic_interaction)



# tumor repressor and oncogenes list
tumor_repressor_gene = pd.read_csv(script_path+'/data/oncogenes_and_tumor_suppressor_genes.bed',sep='\t')
tumor_repressor_genes_list = tumor_repressor_gene.loc[3].tolist() #column 3 has name of genes


#list of dosage sensitive genes
#dosage_sensitive_genes = list(filter(None, open('./data/dosage_sensitive_genes.txt','r').read().split('\n')))
dosage_sensitive_genes = pd.read_csv(script_path+'/data/dosage_sensitive_genes.bed', header=None, sep='\t').iloc[:,3].tolist()


# getting the list of TADS that consist of onco genes and tumor supressor genes
if args.tad_domain == None:
	cancer_tad = pd.read_csv(script_path+'/data/cancer_genes_tad.bed',header=None, sep='\t')
	cancer_tad['tad_name'] = cancer_tad[5] + '-'+cancer_tad[6].map(str)+'-'+cancer_tad[7].map(str)
	
	#list of tad domains with tumor repressor or oncogenes
	cancer_tad_list = cancer_tad['tad_name'].tolist()

	dosage_tad = pd.read_csv(script_path+'/data/dosage_genes_tad.bed',header=None, sep='\t')
	dosage_tad['tad_name'] = dosage_tad[5] + '-'+dosage_tad[6].map(str)+'-'+dosage_tad[7].map(str)
	dosage_tad_list = dosage_tad['tad_name'].tolist()


	gene_density = pd.read_csv(script_path+'/data/gene_density_all_tad.csv')
	gene_density = gene_density.set_index('name')

	# tad domain information for genes
	genes_tad = pd.read_csv(script_path+'/data/genes_tad.bed', sep='\t', header=None).iloc[:,[3,5,6,7]]
	genes_tad.columns = ['gene','chr','start','end']
	genes_tad['name'] = genes_tad.chr + '-'+ genes_tad['start'].map(str)+'-'+genes_tad['end'].map(str)
	genes_tad = genes_tad.groupby('gene').agg(lambda x: list(set(list(x))))
	

else:
	# TAD gene density calculation
	os.system('bedtools intersect -a {}/data/grch37_ensembl_genes_104.bed -b {} -wb > {}/genes_tad.bed'.format(script_path,args.tad_domain,dir_))
	os.system('bedtools intersect -a {}/data/oncogenes_and_tumor_suppressor_genes.bed -b {} -wb > {}/cancer_genes_tad.bed'.format(script_path,args.tad_domain,dir_))
	os.system('bedtools intersect -a {}/data/dosage_sensitive_genes.bed -b {} -wb > {}/dosage_genes_tad.bed'.format(script_path,args.tad_domain,dir_))

	cancer_tad = pd.read_csv("./{}/cancer_genes_tad.bed".format(args.output),header=None, sep='\t')
	cancer_tad['tad_name'] = cancer_tad[5] + '-'+cancer_tad[6].map(str)+'-'+cancer_tad[7].map(str)

	dosage_tad = pd.read_csv('./{}/dosage_genes_tad.bed'.format(args.output),header=None, sep='\t')
	dosage_tad['tad_name'] = dosage_tad[5] + '-'+dosage_tad[6].map(str)+'-'+dosage_tad[7].map(str)
	
	#list of tad domains with tumor repressor or oncogenes
	cancer_tad_list = cancer_tad['tad_name'].tolist()
	dosage_tad_list = dosage_tad['tad_name'].tolist()

	gene_density = pd.read_csv('./{}/genes_tad.bed'.format(args.output), sep='\t', header=None).iloc[:,[3,5,6,7]]
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


info = info + "Nearby cancer gene distance:\t{}kb\n".format(dist) + "Gene density:\t{}, this value represents mean value if user parameter is None\n".format(gd)
info = info + "Output folder:\t{}\n".format(args.output) + "Output filename:\t{}\n".format(args.file_name)

logger.info(info)



#################################### Reading in input files and getting TAD domain information for variants ###########################

logger.info('reading in the input file: {}'.format(args.input))
input_data = pd.read_csv(args.input, sep='\t')	# reading in SNP eQTL data
columns = [x.lower() for x in input_data.columns.tolist()]
input_data.columns = columns
input_data['chr'], input_data['start'], input_data['stop'] = zip(*input_data.apply(get_bed_file, axis=1))

columns_names = input_data.columns


# extending the sequence to look for tumor repressor or oncogenes nearby
input_data = input_data.sort_values(by=['chr','start'], ascending=[True, True])
#input_data['tissue'] ='.'

logger.info('creating bed file from the input data')
input_data.loc[:,['chr','start','stop','id']].to_csv(os.path.join(dir_,'sorted_variant_coordinates.bed'), header=False, sep='\t',index=False)		# writing bed file for snp coordinates

if int(args.nearby_cancer_genes) > 0:
	input_data['extended_start'] = input_data['start'].apply(lambda x: int(x)- int(args.nearby_cancer_genes) if int(x)>int(args.nearby_cancer_genes) else 0)
	input_data['extended_stop'] = input_data['stop']+ int(args.nearby_cancer_genes)
	input_data.loc[:,['chr','extended_start','extended_stop','id']].to_csv(os.path.join(dir_,'sorted_extended_variant_coordinates.bed'), header=False, sep='\t',index=False)

	logger.info('looking for nearby tumor repressor or oncogenes')
	os.system('bedtools intersect -a {}/sorted_extended_variant_coordinates.bed -b ./data/oncogenes_and_tumor_suppressor_genes.bed -wb > {}/oncogenic_tumor_repressor_variant.bed'.format(dir_,dir_))

	variant_nearby_cancer = pd.read_csv('{}/oncogenic_tumor_repressor_variant.bed'.format(dir_), header=None, sep='\t')
	variant_nearby_cancer.columns = ['chr','start','stop','id','c_chr','c_start','c_stop','genes','gene_id']
	variant_nearby_cancer_list = variant_nearby_cancer['id'].tolist()


# variants overlap region with TAD domain
if args.tad_domain == None:
	os.system('bedtools intersect -a {}/sorted_variant_coordinates.bed -b {}/data/merged_gm12878.bed  -wb > {}/variant_tad.bed'.format(dir_,script_path,dir_)) 
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
input_data['same_dosage_tad'] = input_data['tad_name'].apply(lambda x: True if len(set(x)&set(dosage_tad_list))>0 else False) # checking if any of the tads overlap







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
#variant_hic_promoter.to_csv('./results/temp_variant_tad_count.csv')

#variant_hic_promoter = variant_hic_promoter[variant_hic_promoter['common_tad']==0] #removed same tad interaction

logger.info('checking if interacted gene falls in dosage_sensitive_genes, 1 or more value assigned depending on number of interaction, else if there is no interacted gene or does not fall in dosage_sensitive_genes then assignn 0')
variant_hic_promoter['label'] = variant_hic_promoter['hic_interacted_gene'].apply(lambda x: 1 if x in dosage_sensitive_genes else 0)


variant_hic_promoter_ = variant_hic_promoter.iloc[:,[3,-2,-1,7]].groupby('snp').agg(lambda x: list(x))
variant_hic_promoter_['dosage_sensitive_interaction'] = variant_hic_promoter_['label'].apply(lambda x: sum(x))
variant_hic_promoter_['common_tad_count'] = variant_hic_promoter_['common_tad'].apply(lambda x: sum(x))



all_data = input_data.merge(variant_hic_promoter_.loc[:,['hic_interacted_gene','common_tad_count','dosage_sensitive_interaction']], left_on='id', right_on=variant_hic_promoter_.index)

logger.info('checking if any of interacted genes are tumor repressor or oncogenes')
all_data['hic_interacted_gene_test'] = all_data['hic_interacted_gene'].apply(lambda x: filter_tumor_repressor_genes(x, tumor_repressor_genes_list))






######################################## Checking for repressive region and nearby genes ########################


logger.info('running bedtools to get the information regarding overlap with repressive region')
os.system('bedtools intersect -a {}/sorted_variant_coordinates.bed -b {} -wo > {}/variant_repressive_marks.bed'.format(dir_,args.repressive_region,dir_))

#repressive_marks = pd.read_csv(dir_+'/variant_repressive_marks.bed',sep='\t', header=None)
#repressive_marks_list = repressive_marks.iloc[:,3].tolist() # variants id overlapping with repressive region

#all_data['repressive_region'] = all_data['id'].apply(lambda x: check_active_status(x, repressive_marks_list))

variant_repressive = pd.read_csv(dir_+'/variant_repressive_marks.bed', sep='\t', header=None).iloc[:,[3,7]]
variant_repressive.columns = ['id', 'state']
variant_repressive = variant_repressive.groupby('id').agg(lambda x : ','.join(list(set(list(x)))))

all_data[['repressive_region', 'repressive_region_info']] = all_data['id'].apply(lambda x: check_chromatin_status(x, variant_repressive))

logger.info('checking for variants with nearby oncogenes or tumor repressor genes')

#print(filter7.head())
if int(args.nearby_cancer_genes) > 0:
	all_data['nearby_cancer_genes'] = all_data['id'].apply(lambda x: filter_nearby_cancer_genes(x, variant_nearby_cancer_list))
	all_data['nearby_cancer_gene_names'] = all_data['id'].apply(lambda x: get_nearby_genes(x, variant_nearby_cancer))


	

######################################## Checking for active chromatin region ########################

logger.info('running bedtools to get the information regarding overlap with active transcription region')
os.system('bedtools intersect -a {}/sorted_variant_coordinates.bed -b {} -wo > {}/variant_active.bed'.format(dir_,args.active_region,dir_))


#active_region = pd.read_csv(dir_+'/variant_active.bed',sep='\t', header=None)
#active_region_list = active_region.iloc[:,3].tolist()

logger.info('tagging ME overlapping with active chromatin region as True')
#filter7['active_region'] = filter7['id'].apply(lambda x: check_active_status(x, active_region_list))
#all_data['active_region'] = all_data['id'].apply(lambda x: check_active_status(x, active_region_list))



active_region = pd.read_csv(dir_+'/variant_active.bed', sep='\t', header=None).iloc[:,[3,7]]
active_region.columns = ['id', 'state']
active_region = active_region.groupby('id').agg(lambda x : ','.join(sorted(list(set(list(x))))))

all_data[['active_region', 'active_region_info']] = all_data['id'].apply(lambda x: check_chromatin_status(x, active_region))



################################################ Filtering by FDR if provided by user ######################################################
if args.thresh != None:
	logger.info('filtering by FDR > {}'.format(args.thresh))

	all_data['fdr_test'] = all_data['fdr'].apply(lambda x: fdr_filter(str(x).replace(' ','').split(','), args.thresh))
	#logger.info('filter6 shape (after removing ME with FDR) > {} : {}'.format(args.thresh, filter6.shape))
	filtered_data = all_data[all_data['fdr_test']== True]
else:
	if 'fdr' in columns:
		logger.warning('\n\n----------- !!! Warning !!! -------------------\n\n there seems to be fdr column in input file, but fdr parameter is not used, so filtration steps is carried out without using fdr column.\n\n-----------------------------------------------\n')
	filtered_data = all_data.copy()






########################################### Filtering by allele frequency and eQTL genes if provided###############

if args.allele_freq != None:
	logger.info('filtering by allele frequency > {}'.format(args.allele_freq))
	#filter5 = filter4[filter4['AF'] > float(args.allele_freq)]
	all_data['af'] = all_data['af'].apply(lambda x: get_round_value(x))
	filtered_data = filtered_data[(filtered_data['af']> float(args.allele_freq))&(filtered_data['af']< (1-float(args.allele_freq)))]
	#logger.info('filter5 shape (after removing ME with AF) > {} : {}'.format(args.allele_freq, filter5.shape))
else:
	if 'af' in columns:
		logger.warning('\n\n----------- !!! Warning !!! -------------------\n\n there seems to be allele frequency column in input file, but "af" parameter is not used, so filtration steps is carried out without using "af" column.\n\n-----------------------------------------------\n')



#if args.eqtl_genes != False:
#	all_data['eQTL_test'] = all_data['eqtl'].apply(lambda x: filter_tumor_repressor_genes(x.replace(' ','').split(','), tumor_repressor_genes_list))

#	eqtl_passed_variants = all_data[all_data['eQTL_test']==False]['id'].tolist()
	#print(eqtl_passed_variants)

#	logger.info('keeping only those variants if EQTL gene is not tumor repressor or oncogenes')
#	filtered_data = filtered_data[filtered_data['id'].isin(eqtl_passed_variants)]


######################################### removing variants overlapping with the blacklist region provided by the user ###############################

if args.blacklist_region != None:
	os.system('bedtools intersect -a {}/sorted_variant_coordinates.bed -b {}  -wa > {}/variant_blacklist_overlap.bed'.format(dir_, args.blacklist_region, dir_))

	blacklist_region = pd.read_csv(dir_+ '/variant_blacklist_overlap.bed', sep='\t', header=None)
	blacklist_region_list = blacklist_region[3].tolist()
	
	logger.info('tagging variants if they overlap with blacklisted region')
	all_data['blacklist_region'] = all_data['id'].apply(lambda x: tag_blacklist_region(x, blacklist_region_list))

	non_blacklisted_ids = all_data[all_data['blacklist_region']==False]['id'].tolist()

	logger.info('removing any variants overlapping with blacklisted region')
	filtered_data = filtered_data[filtered_data['id'].isin(non_blacklisted_ids)]

	

######################################### Filteration steps: applying all remaining filters ############################


logger.info('removing variants in same TAD domain as tumor repressor or oncogenes')
logger.info('removing variants with TAD domain having gene density < {}'.format(gd))
logger.info('removing variants interacting with dosage sensitive genes or tumor repressor or oncogenes or any interaction with genes in same tad domain')
logger.info('removing variants with nearby tumor repressor or oncogenes')

if int(args.nearby_cancer_genes) >0:
	filtered_data = filtered_data[(filtered_data['same_cancer_tad']==False)&(filtered_data['same_dosage_tad']==False)&(filtered_data['gene_density']<gd)&(filtered_data['dosage_sensitive_interaction']==0)&(filtered_data['repressive_region']==False)&(filtered_data['nearby_cancer_genes']==False)&(filtered_data['common_tad_count']==0)&(filtered_data['hic_interacted_gene_test']==False)]
else:
	filtered_data = filtered_data[(filtered_data['same_cancer_tad']==False)&(filtered_data['same_dosage_tad']==False)&(filtered_data['gene_density']<gd)&(filtered_data['dosage_sensitive_interaction']==0)&(filtered_data['repressive_region']==False)&(filtered_data['common_tad_count']==0)&(filtered_data['hic_interacted_gene_test']==False)]




########################################### Generating output files ##############################


logger.info('generating output files....')
filtered_data['ucsc_link'] = 'http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position='+filtered_data['position'].map(str)
filtered_data.loc[:,['id','position','active_region','active_region_info','ucsc_link']].to_csv(dir_+'/filtered_'+args.file_name, index=False)

all_data['passed_all_filter'] = all_data['id'].apply(lambda x: check_final_list(x, filtered_data['id'].tolist()))
all_data['tad_name'] = all_data['tad_name'].apply(lambda x: ','.join(x))
all_data['hic_interacted_gene'] = all_data['hic_interacted_gene'].apply(lambda x: ','.join(x))

all_data.rename(columns={'gene_density':'gene_density ({}, genes per million tad)'.format(gd), 'hic_interacted_gene_test':'hic_interacted_genes (oncogenic or tumor repressor)','nearby_cancer_genes':'nearby_cancer_genes ({}kb)'.format(dist)}, inplace=True)
all_data.to_csv(dir_+'/'+args.file_name, index=False)

logger.success('\n\nThe total number of identified safe harbor candidates is :  {}\n'.format(filtered_data.shape[0]))


try:
	os.system('mv {}/*.bed {}/temp_files/'.format(dir_, dir_))
except IOError:
	logger.info('no bed files found to move to temp folder, please check the folder again')



