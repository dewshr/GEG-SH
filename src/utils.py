from loguru import logger
import pandas as pd


########################## function to calculate the gene density in a given tad domain ######
@logger.catch
def calculate_gene_density(tad, gene_density):
	x = [ i if 'chr' in i else None for i in tad ]
	x = list(filter(None, x))
	val =0.0
	if len(x)>0:
		for i in x:
			try:
				val = val + float(gene_density[gene_density.index == i]['density'])
			except:
				pass
		return round(val/len(x), 2)
	else:
		return 0



######################## function to filter given user defined FDR threshold ###############
#False means it didn't pass the Fdr threshold

'''
@logger.catch
def fdr_filter2(x, threshold):
	fdr_list =[]
	for fdr in x:
		if fdr.lower() == 'nan':
			fdr_list.append(True)
		elif fdr.lower() == 'not-tested':
			fdr_list.append(True)
		elif float(fdr) < float(threshold):
			fdr_list.append(False)
		else:
			fdr_list.append(True)

	if False in fdr_list:
		return False
	else:
		return True
'''

@logger.catch
def fdr_filter(x, threshold):
	fdr_list =[]
	for fdr in x:
		try:
			if float(fdr) < threshold:
				fdr_list.append(False)
			else:
				fdr_list.append(True)
		except:
			fdr_list.append(True)

	if False in fdr_list:
		return False
	else:
		return True



######################## removing the oncogenic eQTL genes ###############
@logger.catch
def filter_tumor_repressor_genes(x, tumor_repressor_genes_list):
	#print(x)
	length = len(set(x).intersection(set(tumor_repressor_genes_list)))
	#print(length)
	if length == 0:
		return False
	else:
		return True


################ function to check for nearby tumor repressor or oncogenes
@logger.catch
def filter_nearby_cancer_genes(x, mei_nearby_cancer_list):
	#print(length)
	#print(x)

	if x in  mei_nearby_cancer_list:
		return True
	else:
		return False



############# function to check if the ME and gene interacted are in same tad domain ######
@logger.catch
def check_tad(snp, gene, mei_tad, genes_tad):
	#print(x['snp'])
	#print(genes_tad.head())
	try:
		snp_tad = mei_tad[mei_tad.index == snp]['tad_name'].tolist()[0]
	except:
		return 0
		
	try:
		gene_tad = genes_tad[genes_tad.index== gene]['name'].tolist()[0]

	except:
		return 0

	if len(snp_tad)!=0 and len(gene_tad)!=0 and len(set(snp_tad).intersection(set(gene_tad))) >0:
		return 1
	else:
		return 0
	



##################### function to scan for repressive and active chromatin region ############
def check_chromatin_status(x, chromatin_region_df):
	try:
		val = chromatin_region_df[chromatin_region_df.index==x]['state'][0]
	except:
		val = 'None'

	if val != 'None':
		return pd.Series([True, val])
	else:
		return pd.Series([False, val])



################## functuion to assign the tad domain information to ME #########
@logger.catch
def get_tad_info(x, mei_tad):
	try:
		return (mei_tad[mei_tad.index == x])['tad_name'][0]
	except:
		return ['-']




########################## rounding off the decimal values ##################
@logger.catch
def get_round_value(x):
    try:
        return round(float(x),3)
    except:
        return x



########################### tagging the raw data file if it passed all the filters or not ####################
@logger.catch
def check_final_list(x, final_list):
	if x in final_list:
		return True
	else:
		return False


################################# generating bed files using position column #######################
@logger.catch
def get_bed_file(df):
	x = df.position
	x_list = x.replace(':','-').split('-')
	#if 'INS' in x:
	#	return x_list[0], int(x_list[1]), int(x_list[1])
	#else:
	return x_list[0], int(x_list[1]), int(x_list[2])



######################## extracting nearby genes for the given MEIs ######################
@logger.catch
def get_nearby_genes(x, mei_nearby_cancer):
	try:
		return ', '.join(set(mei_nearby_cancer[mei_nearby_cancer.id==x]['genes'].tolist()))
	except:
		return ''


######################## checking if the variant position overlaps with blacklisted position ######################
@logger.catch
def tag_blacklist_region(x, blacklist_region_list):
	if x in blacklist_region:
		return True
	else:
		return False




