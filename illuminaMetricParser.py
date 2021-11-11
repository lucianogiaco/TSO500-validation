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
		self.DnaLibraryQC = meArray[20:23]
		self.DnaSVandTmbQC = meArray[25:29]
		self.DnaMsiQC = meArray[31:33]
		self.DnaCnvQC = meArray[35:38]
		self.DnaExpanded = meArray[40:55]
		self.RnaLibraryQC = meArray[57:61]
		self.RnaExpanded = meArray[63:68]
		
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
		print("[WARNING] "+inputFile)
		print("[WARNING] The file '"+inputFile+"' does not exist")
		print("[INFO] Exit")
		os.sys.exit()

def main(inputFile):

	existingFile(inputFile)

	mefile = open(inputFile,'r')
	meArray = mefile.readlines()
	mefile.close()

	MetricsField = metrics_fields(meArray)

	# parserRunQCmetrics
	k = MetricsField.getRunQCmetrics()[0]
	k = str.strip(k)
	k = k.strip('[')
	k = k.strip(']')
	k = k.replace(' ', '')
	

	for f in MetricsField.getRunQCmetrics()[2:]:
		
		f = str.strip(f)
		s = f.split('\t')
		gl = round(float(s[1]), 1)
		val = round(float(s[3]), 1)
		if val >= gl:
			print(1)
		else:
			print(0)



	
	






if __name__ == '__main__':
	# parser variable
	parser = argparse.ArgumentParser(description='Illumina Metrics Parser')

	# arguments
	parser.add_argument('-i', '--inputFile', required=True,
						help='Insert MetricsOutput.tsv path file')

	args = parser.parse_args()
	inputFile = args.inputFile


	main(inputFile)