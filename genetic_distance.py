#!/usr/bin/env python

import numpy as np
from numpy import array
import sys
import csv

#genetic pairwise diversit per sample
def genetPairDiver(vcfIO):
	lines = 0
	filename = input("Give file name to save the matrix as 'name.txt':")
	for line in vcfIO:
		if line[1] == '#':
			pass
		else:
			if line[1] == 'C':
				line2 = line.split('\t')
				line2 = line2[9::]
				indexes = [x for x in range(0,len(line2))]
				num_diff = np.zeros((len(indexes), len(indexes))) #create the table that I can add if different the two indiv
				num_comp = np.zeros((len(indexes), len(indexes))) #create the table that I can add number of comparisons
				#print indexes, line2
			else:
				lines = lines + 1 #count the number of lines in the variant file
				line2 = line.split('\t')
				line2 = line2[9::]
				for indv in indexes:
					i = indv + 1
					gt1 = line2[indv].split(':') #get the correct index for each indv
					gt1 = gt1[0]
					while i < len(indexes): #to compare with the next indv in the line, The comparison with the ones before has already be done, and no need to compare with itself
						gt2 = line2[i].split(':')
						gt2 = gt2[0]
						#print lines,indv, i
						if gt1 == './.' or gt2 =='./.': #skip if missing info
						#	print 1
							pass 
						else:
							if gt1 != gt2:
								num_diff[indv, i] = num_diff[indv, i] + 1
								num_comp[indv, i] = num_comp[indv, i] + 1
								num_diff[i, indv] = num_diff[i, indv] + 1 #add in both columns and rows, the matrix for heatmap3 must be symmetrical
								num_comp[i, indv] = num_comp[i, indv] + 1
								#print num_comp
						#		print 2
							elif gt1 == gt2:
								num_comp[indv, i] = num_comp[indv, i] + 1
								num_comp[i, indv] = num_comp[i, indv] + 1 #count it only in the number of comparisons.
						#		print 3
						i = i + 1
	diff_1 = num_diff / num_comp
	diff_1 = diff_1 * float(lines)
	diff_1 = diff_1 / float(161845050)  #normalise by multiply total number of sites in variant sites file and divide by total number of lines in nonvariant file 
	np.savetxt('{}'.format(filename), diff_1, delimiter='\t')			

file = sys.argv[1]
vcfIO = open(file, 'r')
genetPairDiver(vcfIO)
sys.exit("Live long and prosper!")
