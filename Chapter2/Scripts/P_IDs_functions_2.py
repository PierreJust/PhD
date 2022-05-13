#! /usr/bin/env python 2.7
# coding: utf-8


#Pierre Justeau

import os
from math import sqrt
#################################################
##################Genpofad#######################
#################################################
def perId_uniparentalMarkers(Net, haplotypesFiles):
	
	linkSample1 = {}
	linkSample2 = {}
	dic_merged_links = {}
	dic_en_commun = {}
	dic_length = {}
	result = {}
	number_mut = {}
	
	for x in Net:                              
# we parse network
		for y in haplotypesFiles:              
# we parse haplotype files
			
			if x[0] == y[0]:			       
# if the first IDs of the network (first column x[0]) is == to the IDs of the haplotype files y[0]
				tmp= ','.join(x[0:2])	       
# we keep both IDs x[0:2] of the network that creates a link (e.g sampl1, sample2)
				linkSample1[str(tmp)] = y[1]   
# here we add the haplotypes that correspond to the sample from the first column x[0]
				
			elif x[1] == y[0]:				   
# if the second IDs of the network (second column x[1]) is == to the IDs of the haplotype files y[0]
				tmp2= ','.join(x[0:2])		   
# we keep both IDs x[0:2] of the network that creates a link (e.g sampl1, sample2)
				linkSample2[str(tmp2)] = y[1]  
# here we add the haplotypes that correspond to the samples from the second column x[1]
	
	
	
	
	for k,v in linkSample1.items():								
# here we parse our first dic
		for k2,v2 in linkSample2.items():                                                
# here we parse our second dic
			
			if k == k2:                                                                  
# if there IDs(links) are ==
				shared = [x for x in v.split(',') if x in v2.split(',')]              
# we look for each mutations that are ==
				              								
				dic_en_commun[str(k)] = shared                                           
# we create a dic with IDs(links) plus the mutations that the samples share
				dic_length[k] = len(v.split(',')), len(v2.split(','))                  	 
# we create a dic with the length of the haplotype (number mutations) for each links like : Sample1, Sample2 : {28, 32}
	

	
	for k,v in dic_en_commun.items():                           
# here we parse our dic with IDs and the mutations shared between samples
		for k2,v2 in dic_length.items():						
# here we parse our dic with the lengths			
			if k == k2:                                         
# if their IDs are ==
				result[str(k2)] = float(len(v))/max(v2)	
# we create a new dic with the IDs(links) of samples linked plus their genetic distance
				number_mut[str(k2)] = len(v)                    
# we create a new dic with the IDs(links) and the number of mutations shared
						
	return (result, number_mut)


def perId_mtDNA_MSY(perId_MSY, perId_mtDNA):
	
	perID_mt_y = {}
	result_perID_chrY_mt = {}
	
	for k,v in perId_MSY.items():                                   
# here we parse our dic for the genetic distance of MSY SNPs shared between IDs
		for k2,v2 in perId_mtDNA.items():                           
# here we parse our dic for the genetic distance of mtDNA SNPs shared between IDs
			
			if k == k2:						
# if they have the same links IDs 	
				perID_mt_y[k] = v,v2				
# we create a new dic with IDs plus the value of genetic distance for mtDNA plus for MSY
				
	for k,v in perID_mt_y.items():
		result_perID_chrY_mt[k] = sqrt(float(v[0])**2 + float(v[1])**2) 
# here we compute the geometric distance for MSY and mtDNA to have a single value (if they share same mtDNA Ychr HGs)
		
		
	return (result_perID_chrY_mt)

def split_perId_mtDNA_MSY(perId_score,list_Id_shared,nature_links):

	result_final_perId = {}
	
	for ids_key,distance_val in perId_score.items():
		if ids_key in list_Id_shared:
			tmp = [str(distance_val),nature_links]
			result_final_perId[ids_key]=','.join(tmp)

	return (result_final_perId)
	
	
	



