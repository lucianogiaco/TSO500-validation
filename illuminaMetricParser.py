#############################
# Script to parse the MetricsOutput.tsv
# provided by Illumina tso Local App
# Author: Luciano Giaco'
# Date: 2021.11.11
#############################

import os
import argparse


class metrics_fields:
	"""docstring for metrics_fields"""
	def __init__(self, meArray):
		super(metrics_fields, self).__init__()


		self.header = meArray[3:6]
		self.runQCmetrics = meArray[8:13]
		self.analysisStatus = meArray[14:18]
		self.DnaLibraryQC = meArray[20:24]
		self.DnaSVandTmbQC = meArray[25:30]
		self.DnaMsiQC = meArray[31:34]
		self.DnaCnvQC = meArray[35:39]
		self.DnaExpanded = meArray[40:56]
		self.RnaLibraryQC = meArray[57:62]
		self.RnaExpanded = meArray[63:69]
		
	def getHeader(self):
		return self.header

	def getRunQCmetrics(self):
		return self.runQCmetrics

	def getAnalysisStatus(self):
		return self.analysisStatus

	def getDnaLibraryQC(self):
		return self.DnaLibraryQC

	def getDnaSVandTmbQC(self):
		return self.DnaSVandTmbQC
		
	def getDnaMsiQC(self):
		return self.DnaMsiQC

	def getDnaCnvQC(self):
		return self.DnaCnvQC

	def getDnaExpanded(self):
		return self.DnaExpanded

	def getRnaLibraryQC(self):
		return self.RnaLibraryQC

	def getRnaExpanded(self):
		return self.RnaExpanded



def existingFile(inputFile):
	if os.path.isfile(inputFile) is True:
		return(True)
	else:
		#print("[WARNING] "+inputFile)
		#print("[WARNING] The file '"+inputFile+"' does not exist")
		#print("[INFO] Exit")
		os.sys.exit()

def get_contaminatio(DNA_metrics):
	# parserDnaLibraryQC contamination
	k = DNA_metrics[0]
	k = str.strip(k)
	k = k.strip('[')
	k = k.strip(']')
	k = k.replace(' ', '')
	k = k.replace('\t', '')
	#print('[INFO] Parsing metrics')
	#print(k)
	#print('\n')

		
	samples = DNA_metrics[1]
	samples = str.strip(samples)
	samples = samples.split('\t')
	samples = samples[3:]

	score = DNA_metrics[2]
	score = str.strip(score)
	score = score.split('\t')
	score = score[3:]

	pval = DNA_metrics[3]
	pval = str.strip(pval)
	pval = pval.split('\t')
	pval = pval[3:]

	i = 0
	l = len(samples)
	main_dict = dict()
	while i < l:
		tSmp = samples[i]
		tScr = score[i]
		tPvl = pval[i]
		#print(tSmp)

		if tScr == 'NA' and tPvl == 'NA':
			main_dict[tSmp] = ['NA', 'NA', 'NA']
			#print(tScr)
			#print(tPvl)
			#print('\n')
		else:
			tScr = int(tScr)
			tPvl = float(tPvl)
			#print(tScr)
			#print(tPvl)
			#print('\n')
			if tScr <= 3106:
				main_dict[tSmp] = [tScr, tPvl, 1]
				next
			elif tScr > 3106 and tPvl <= 0.049:
				main_dict[tSmp] = [tScr, tPvl, 1]
				next
			else:
				main_dict[tSmp] = [tScr, tPvl, 0]
		i = i + 1

	return(main_dict)

def	parserDnaSVandTmbQC(TMB_metrics, main_dict):
	k = TMB_metrics[0]
	k = str.strip(k)
	k = k.strip('[')
	k = k.strip(']')
	k = k.replace(' ', '')
	k = k.replace('\t', '')
	#print('[INFO] Parsing metrics')
	#print(k)
	#print('\n')

	samples = TMB_metrics[1]
	samples = str.strip(samples)
	samples = samples.split('\t')
	samples = samples[3:]
	#print(samples)
	median_size = TMB_metrics[2]
	median_size = str.strip(median_size)
	median_size = median_size.split('\t')
	median_size = median_size[3:]
	#print(median_size)
	median_ex_cover = TMB_metrics[3]
	median_ex_cover = str.strip(median_ex_cover)
	median_ex_cover = median_ex_cover.split('\t')
	median_ex_cover = median_ex_cover[3:]
	#print(median_ex_cover)
	pct_ex_50x = TMB_metrics[4]
	pct_ex_50x = str.strip(pct_ex_50x)
	pct_ex_50x = pct_ex_50x.split('\t')
	pct_ex_50x = pct_ex_50x[3:]
	#print(pct_ex_50x)


	i = 0
	l = len(samples)
	while i < l:
		tSmp = samples[i]
		tMsz = median_size[i]
		tMcv = median_ex_cover[i]
		tPce = pct_ex_50x[i]
		#print(tSmp)

		if tMsz == 'NA' and tMcv == 'NA' and tPce == 'NA':
			#print(tMsz)
			#print(tMcv)
			#print(tPce)
			#print('\n')
			v = main_dict[tSmp]
			v = v + ['NA','NA','NA', 'NA','NA','NA']
			main_dict[tSmp] = v
		else:
			tMsz = int(tMsz)
			tMcv = int(tMcv)
			tPce = float(tPce)
			#print(tMsz)
			#print(tMcv)
			#print(tPce)
			#print('\n')
			v = main_dict[tSmp]
			v = v + [tMsz, tMcv, tPce]
			main_dict[tSmp] = v
			if tMsz >= 70:
				v = main_dict[tSmp]
				v.append(1)
				main_dict[tSmp] = v
			else:
				v = main_dict[tSmp]
				v.append(0)
				main_dict[tSmp] = v
			if tMcv >= 150:
				v = main_dict[tSmp]
				v.append(1)
				main_dict[tSmp] = v
			else:
				v = main_dict[tSmp]
				v.append(0)
				main_dict[tSmp] = v
			if tPce >= 90.0:
				v = main_dict[tSmp]
				v.append(1)
				main_dict[tSmp] = v
			else:
				v = main_dict[tSmp]
				v.append(0)
				main_dict[tSmp] = v
		i = i + 1

	return(main_dict)

def get_msi(MSI_metrics, main_dict):
	

	k = MSI_metrics[0]
	k = str.strip(k)
	k = k.strip('[')
	k = k.strip(']')
	k = k.replace(' ', '')
	k = k.replace('\t', '')
	#print('[INFO] Parsing metrics')
	#print(k)
	#print('\n')

	samples = MSI_metrics[1]
	samples = str.strip(samples)
	samples = samples.split('\t')
	samples = samples[3:]
	#print(samples)
	unst_msi = MSI_metrics[2]
	unst_msi = str.strip(unst_msi)
	unst_msi = unst_msi.split('\t')
	unst_msi = unst_msi[3:]
	#print(unst_msi)

	i = 0
	l = len(samples)
	while i < l:
		tSmp = samples[i]
		tMsi = unst_msi[i]
		#print(tSmp)

		if tMsi == 'NA':
			#print(tMsi)
			#print('\n')
			v = main_dict[tSmp]
			v = v + ['NA', 'NA']
			main_dict[tSmp] = v

		else:
			tMsi = int(tMsi)
			#print(tMsi)
			#print('\n')
			v = main_dict[tSmp]
			v = v + [tMsi]
			main_dict[tSmp] = v
			if tMsi >= 40:
				v = main_dict[tSmp]
				v.append(1)
				main_dict[tSmp] = v
			else:
				v = main_dict[tSmp]
				v.append(0)
				main_dict[tSmp] = v

		i = i + 1
	return(main_dict)

def get_CNV(CNV_metrics, main_dict):
	k = CNV_metrics[0]
	k = str.strip(k)
	k = k.strip('[')
	k = k.strip(']')
	k = k.replace(' ', '')
	k = k.replace('\t', '')
	#print('[INFO] Parsing metrics')
	#print(k)
	#print('\n')

	samples = CNV_metrics[1]
	samples = str.strip(samples)
	samples = samples.split('\t')
	samples = samples[3:]

	cov_mad = CNV_metrics[2]
	cov_mad = str.strip(cov_mad)
	cov_mad = cov_mad.split('\t')
	cov_mad = cov_mad[3:]

	median_bin = CNV_metrics[3]
	median_bin = str.strip(median_bin)
	median_bin = median_bin.split('\t')
	median_bin = median_bin[3:]

	i = 0
	l = len(samples)
	while i < l:
		tSmp = samples[i]
		tCvm = cov_mad[i]
		tMdb = median_bin[i]
		#print(tSmp)

		if tCvm == 'NA' and tMdb == 'NA':
			#print(tCvm)
			#print(tMdb)
			#print('\n')
			v = main_dict[tSmp]
			v = v + ['NA', 'NA', 'NA', 'NA']
			main_dict[tSmp] = v

		else:
			tCvm = float(tCvm)
			tMdb = float(tMdb)
			#print(tCvm)
			#print(tMdb)
			#print('\n')
			v = main_dict[tSmp]
			v = v + [tCvm, tMdb]
			main_dict[tSmp] = v
			if tCvm <= 0.210:
				v = main_dict[tSmp]
				v.append(1)
				main_dict[tSmp] = v
			else:
				v = main_dict[tSmp]
				v.append(0)
				main_dict[tSmp] = v
			if tMdb >= 1.0:
				v = main_dict[tSmp]
				v.append(1)
				main_dict[tSmp] = v
			else:
				v = main_dict[tSmp]
				v.append(0)
				main_dict[tSmp] = v

		i = i + 1
	return(main_dict)

def get_rnaQCmetrics(RNA_mestrics, main_dict):
	k = RNA_mestrics[0]
	k = str.strip(k)
	k = k.strip('[')
	k = k.strip(']')
	k = k.replace(' ', '')
	k = k.replace('\t', '')
	#print('[INFO] Parsing metrics')
	#print(k)
	#print('\n')

	samples = RNA_mestrics[1]
	samples = str.strip(samples)
	samples = samples.split('\t')
	samples = samples[3:]
	#print(samples)
	median_cv_size = RNA_mestrics[2]
	median_cv_size = str.strip(median_cv_size)
	median_cv_size = median_cv_size.split('\t')
	median_cv_size = median_cv_size[3:]
	#print(median_cv_size)
	total_on_target_size = RNA_mestrics[3]
	total_on_target_size = str.strip(total_on_target_size)
	total_on_target_size = total_on_target_size.split('\t')
	total_on_target_size = total_on_target_size[3:]
	#print(total_on_target_size)
	median_ins_size = RNA_mestrics[4]
	median_ins_size = str.strip(median_ins_size)
	median_ins_size = median_ins_size.split('\t')
	median_ins_size = median_ins_size[3:]
	#print(median_ins_size)

	i = 0
	l = len(samples)
	while i < l:
		tSmp = samples[i]
		tMsz = median_cv_size[i]
		tMcv = total_on_target_size[i]
		tPce = median_ins_size[i]
		#print(tSmp)

		if tMsz == 'NA' and tMcv == 'NA' and tPce == 'NA':
			#print(tMsz)
			#print(tMcv)
			#print(tPce)
			#print('\n')
			v = main_dict[tSmp]
			v = v + ['NA','NA','NA','NA','NA','NA']
			main_dict[tSmp] = v
		else:
			tMsz = float(tMsz)
			tMcv = int(tMcv)
			tPce = int(tPce)
			#print(tMsz)
			#print(tMcv)
			#print(tPce)
			#print('\n')
			v = main_dict[tSmp]
			v = v + [tMsz, tMcv, tPce]
			main_dict[tSmp] = v
			if tMsz <= 93:
				v = main_dict[tSmp]
				v.append(1)
				main_dict[tSmp] = v
			else:
				v = main_dict[tSmp]
				v.append(0)
				main_dict[tSmp] = v
			if tMcv >= 9000000:
				v = main_dict[tSmp]
				v.append(1)
				main_dict[tSmp] = v
			else:
				v = main_dict[tSmp]
				v.append(0)
				main_dict[tSmp] = v
			if tPce >= 80:
				v = main_dict[tSmp]
				v.append(1)
				main_dict[tSmp] = v
			else:
				v = main_dict[tSmp]
				v.append(0)
				main_dict[tSmp] = v
		i = i + 1

	return(main_dict)

def get_dna_exp(DNA_metrics_exp, main_dict):
	k = DNA_metrics_exp[0]
	k = str.strip(k)
	k = k.strip('[')
	k = k.strip(']')
	k = k.replace(' ', '')
	k = k.replace('\t', '')
	#print('[INFO] Parsing metrics')
	#print(k)
	#print('\n')
	samples = DNA_metrics_exp[1]
	samples = str.strip(samples)
	samples = samples.split('\t')
	samples = samples[3:]
	#print(samples)
	dna_exp = DNA_metrics_exp[2:]
	# dna_exp = str.strip(dna_exp)
	# dna_exp = dna_exp.split('\t')
	# dna_exp = dna_exp[3:]
	#print(dna_exp)
	i = 0
	l = len(dna_exp)
	while i < l:
		
		tDexp = dna_exp[i].split('\t')
		tDexp[-1] = str.strip(tDexp[-1])

		##print(tDexp)
		arr_to_app = []
		c = 0
		for t in tDexp[3:]:
			tSmp = samples[c]
			v = main_dict[tSmp]
			v = v + [t]
			main_dict[tSmp] = v
			# #print(t)
			c = c + 1

		i = i + 1

	return(main_dict)

def write_out(main_dict):
	##print('\t'.join(header))
	out = open('table.tab', 'a')

	for k, v in main_dict.items():
		out.write(k)
		out.write('\t')
		out.write('\t'.join([str(x) for x in v]))
		out.write('\n')
	out.close()	

def main(inputFile):

	existingFile(inputFile)

	mefile = open(inputFile,'r')
	meArray = mefile.readlines()
	mefile.close()

	MetricsField = metrics_fields(meArray)

	



	DNA_metrics = MetricsField.getDnaLibraryQC()

	main_dict = get_contaminatio(DNA_metrics)
	#print(main_dict)
	# os.sys.exit()
	TMB_metrics = MetricsField.getDnaSVandTmbQC()
	main_dict = parserDnaSVandTmbQC(TMB_metrics, main_dict)
	#print(main_dict)

	MSI_metrics = MetricsField.getDnaMsiQC()
	main_dict = get_msi(MSI_metrics, main_dict)
	#print(main_dict)

	CNV_metrics = MetricsField.getDnaCnvQC()
	main_dict = get_CNV(CNV_metrics, main_dict)
	#print(main_dict)

	RNA_mestrics = MetricsField.getRnaLibraryQC()
	main_dict = get_rnaQCmetrics(RNA_mestrics, main_dict)
	#print(main_dict)

	DNA_metrics_exp = MetricsField.getDnaExpanded()
	main_dict = get_dna_exp(DNA_metrics_exp, main_dict)

	RNA_metrics_exp = MetricsField.getRnaExpanded()
	main_dict = get_dna_exp(RNA_metrics_exp, main_dict)

	# parserRunQCmetrics
	k = MetricsField.getRunQCmetrics()[0]
	k = str.strip(k)
	k = k.strip('[')
	k = k.strip(']')
	k = k.replace(' ', '')
	

	for f in MetricsField.getRunQCmetrics()[2:]:
		
		f = str.strip(f)
		#print(f)
		s = f.split('\t')
		#print(s)
		gl = round(float(s[1]), 1)
		val = round(float(s[3]), 1)
		if val >= gl:
			for k, v in main_dict.items():
				v = v + [val, 1]
				main_dict[k] = v
		else:
			for k, v in main_dict.items():
				v = v + [val, 0]
				main_dict[k] = v

	return(main_dict)

	
		


	


if __name__ == '__main__':
	# parser variable
	parser = argparse.ArgumentParser(description='Illumina Metrics Parser')

	# arguments
	parser.add_argument('-i', '--inputFile', required=False,
						help='Insert MetricsOutput.tsv path file')

	args = parser.parse_args()
	inputFile = args.inputFile

	header = ['Sample', 'contamination_score', 'contamination_p_value', 'Contamination', 'median_insert_size', 'median_exon_coverage_count',
				'pct_exon_50x_perc', 'median_insert_size_guidelines', 'median_exon_coverage_guidelines',
				'pct_exon_50x_guidelines', 'unstable_msi_sites', 'unstable_msi_sites_guidelines','coverage_mad_count', 
				'median_bin_count_cnv_target_count','coverage_mad_guidelines', 
				'median_bin_count_cnv_target_guidelines',
				'rna_median_cv__gene_500x_perc', 'rna_total_on_target_reads_count', 
				'rna_median_insert_size_count',
				'rna_median_cv__gene_500x_guidelines', 'rna_total_on_target_reads_guidelines', 
				'rna_median_insert_size_guidelines', 
				'total_pf_reads', 'mean_family_size', 'median_target_coverage', 
				'pct_chimeric_reads', 'pct_exon_100x', 'pct_read_enrichment', 
				'pct_usable_umi_reads', 'mean_target_coverage', 'pct_aligned_reads',
				'pct_contamination_est', 'pct_pf_uq_reads', 'pct_target_04x_mean',
				'pct_target_100x', 'pct_target_250x', 'rna_pct_chimeric_reads',
				'rna_pct_on_target_reads', 'rna_scaled_median_gene_coverage', 'rna_total_pf_reads',
				'PCT_PF_READS_perc', 'PCT_PF_READS_guidelines', 'PCT_Q30_R1_perc',
				'PCT_Q30_R1_guidelines', 'PCT_Q30_R2_perc', 'PCT_Q30_R2_guidelines', 'run']


	metrics = ['210715_A01423_0008_AH35CWDRXY.tsv',
				'210729_A01423_0009_AH33WGDRXY.tsv',
				'211022_A01423_0010_AHGYFYDRXY.tsv',
				'211111_A01423_0011_AHH2Y2DRXY.tsv',
				'211122_A01423_0012_AH2YWCDRXY.tsv']

	main_dict_final = dict()
	out = open('table.tab', 'w')
	out.write('\t'.join(header))
	out.write('\n')
	out.close()	
	for inputFile in metrics:
		root, ext = os.path.splitext(inputFile)
		main_dict = main(inputFile)
		#print(main_dict)
		for k, v in main_dict.items():
			v = v + [root]
			main_dict[k] = v
		# main_dict_final = {main_dict_final, main_dict}
		write_out(main_dict)
		print(len(main_dict.keys()))


	

