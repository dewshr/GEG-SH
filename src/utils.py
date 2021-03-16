from loguru import logger


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
@logger.catch
def fdr_filter(x, threshold):
	fdr_list =[]
	for fdr in x:
		if fdr.lower() == 'non-significant':
			fdr_list.append(True)
		elif float(fdr) < float(threshold):
			fdr_list.append(False)
		else:
			fdr_list.append(True)

	if False in fdr_list:
		return False
	else:
		return True



######################## function to filter using user defined tissue ###############
@logger.catch
def tissue_filter(x):
	for tissue in x:
		if args.tissue.lower() in tissue.lower():
			return True

	return False



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
		#print(x['snp'])
		snp_tad = mei_tad[mei_tad.index == snp]['tad_name'].tolist()[0]
	except:
		return 0
		
	try:
		gene_tad = genes_tad[genes_tad['gene']== gene]['name'].tolist()[0]

	except:
		return 0

	if len(snp_tad)!=0 and len(gene_tad)!=0 and len(set(snp_tad).intersection(set(gene_tad))) >0:
		return 1
	else:
		return 0
	

###################### function to check if ME lies in active chromatin region ######
@logger.catch
def check_active_status(x, active_region_list):
	#print(x.name)
	if x in active_region_list:
		return True
	else:
		return False


################## functuion to assign the tad domain information to ME #########
@logger.catch
def get_tad_info(x, mei_tad):
	try:
		return (mei_tad[mei_tad.index == x])['tad_name'][0]
	except:
		return ['-']

@logger.catch
def get_round_value(x):
    try:
        return round(float(x),2)
    except:
        return x


@logger.catch
def check_final_list(x, final_list):
	if x in final_list:
		return True
	else:
		return False


@logger.catch
def get_bed_file(df):
	x = df.position
	x_list = x.replace(':','-').split('-')
	#if 'INS' in x:
	#	return x_list[0], int(x_list[1]), int(x_list[1])
	#else:
	return x_list[0], int(x_list[1]), int(x_list[2])

@logger.catch
def get_nearby_genes(x, mei_nearby_cancer):
	try:
		return ', '.join(set(mei_nearby_cancer[mei_nearby_cancer.id==x]['genes'].tolist()))
	except:
		return ''


@logger.catch
def check_ins(df):
	if 'chr3_37361602_INS_ME_ALU' in df.id.tolist():
		print('\n--------yes------\n')
	else:
		print('\n--------no------\n')