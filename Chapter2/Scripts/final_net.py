#! /usr/bin/env python 2.7
# -*- coding: utf-8 -*-

#Pierre Justeau

import csv
import os
from math import sqrt

def extraction_csv(csv_file):
	
    file_csv = []
    
    with open(csv_file, 'rb') as f:
        reader = csv.reader(f, delimiter=",")
        for rows in reader:
            file_csv.append(rows)
            
    return (file_csv)


def network(variable):
	
	net = []
	list_samples = list(variable.keys())
	n=len(list_samples)
	
	for i,key_i in enumerate(list_samples) :
		list1 = variable[key_i]
		for j in range(i + 1, n) :                                       
# j start at i + 1   
			key_j = list_samples[j]
			share = [elt for elt in list1 if elt in variable[key_j]]	 
# to get the data shared between samples e.g: ['H1mt']										                                 
#                                             ['H2mt', 'I2ay']
			if share :	                                                 
# we exclude samples without links
				 net.append( [key_i, key_j, '_'.join(share)] )	         
# create a new list with all the IDs (pairwise links) plus the data shared 
				 
	return net
	

def share(network1, network2):
	
	nets_compare=[]
	temp_list_compare=[]
	
	for x in network1:
		for y in network2:
			if x[0:2] == y[0:2]:               
# if both network share the same IDs
				temp = list(x)                  
# keep the information from network1
				temp.append(y[2])               
# to add more information (e.g: archaeological culture for the samples that share also the same sub haplogroups)
				nets_compare.append(temp)       
# we copy in a new list				
				temp2 = x[0:2]                  
# we copy the IDs as well in a new temp list (to use them in the next function)
				temp_list_compare.append(temp2)
					
	return (nets_compare,temp_list_compare)



def merge_networks(share_network1, temp_list_to_compare, share_network2):
	
	len_ind_net1 = [len(x) for x in share_network1]  
# to get the length of the indexes in share_network1
	len_ind_net2 = [len(x) for x in share_network2]  
# to get the length of the indexes in share_network2

	diff_len = max(len_ind_net2) - max(len_ind_net1) 
# to get the difference in length between both networks

	for x in share_network1:                    	 
# here we parse network1
		if x[0:2] not in temp_list_to_compare:       
# if the IDs of network1 are not in the temp_list 		
			temp = list(x)                           
# we create a temporary list which contains the information from network 1
			temp.extend(diff_len*('',))              
# here we extend our list based on the different length
			share_network2.append(temp)              
# here we add this temporary list to our last cluster 

	return (share_network2)
	
	
	
	
	
	
	
	
	
