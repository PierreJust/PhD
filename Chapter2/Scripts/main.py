#! /usr/bin/env python 2.7
# -*- coding: utf-8 -*-

#Pierre Justeau

import final_net
import P_IDs_functions_2
import re
import networkx as nx

###########################################################################################												##########################################Network##########################################
###########################################################################################

aDNA_database = final_net.extraction_csv('database_filtered.csv')

print "for mitochondrial haplogroup clustering, select number of the mitochondrial column: (for us) type 9"
arg1=int(raw_input())
print "for Y chromosome haplogroup clustering, select number of the Y chromosome column: (for us) type 7"
arg2=int(raw_input())
print "for mitochondrial sub haplogroup clustering, select the number of the mitochondrial column: (for us) type 10"
arg3=int(raw_input())
print "for Y chromosome sub haplogroup clustering, select the number of the Y chromosome column: (for us) type 8"
arg4=int(raw_input())
print "for archaeological culture clustering, select the number of the archaeological culture column (for us) type 2"
arg5=int(raw_input())

net_HGs = {}
net_sub_HGs = {}
net_archae_cult = {}
	
for i in aDNA_database:
	
	keys = i[0]
	
	HGs_mt = i[arg1]      
	HGs_MSY = i[arg2]    
	
	sub_HGs_mt = i[arg3] 
	sub_HGs_MSY = i[arg4]
	
	archae_cult = i[arg5] 
		
	net_HGs[keys] = HGs_mt, HGs_MSY
	net_sub_HGs[keys] = sub_HGs_mt, sub_HGs_MSY
	net_archae_cult[keys] = [archae_cult]


Hgs_cluster = final_net.network(net_HGs)
sub_HGs_cluster = final_net.network(net_sub_HGs)
cult_cluster = final_net.network(net_archae_cult)

########################################
### we compare the clusters to know: ###
########################################

# 1) Which samples share the same haplogroups and sub haplogroups
# 2) For samples share the same haplogroups and sub haplogroups, do they share the same culture?

share_HGs_subHGs, list_share1 = final_net.share(Hgs_cluster, sub_HGs_cluster)
share_HGs_subHGs_culture, list_share2 = final_net.share(share_HGs_subHGs, cult_cluster)

############################################################
### here we merge our cluster to create a single network ###
############################################################

merge_HGs_subHGsculture = final_net.merge_networks(Hgs_cluster, list_share1, share_HGs_subHGs_culture)
# here we merge the samples that share HGs,subHGs plus the samples that share only HGs
merge_HGs_subHGs_subHGsculture = final_net.merge_networks(share_HGs_subHGs, list_share2, merge_HGs_subHGsculture)
# here we add into the previous cluster the samples that share only HGs and sub HGs

################################################################################
### to create directed graph based on the C14 dating to have the links like: ###
################################################################################

 # 1) sample1(oldest), sample2(youngest) 
samples_edges_col1 = []
samples_edges_col2 = []

for x in merge_HGs_subHGs_subHGsculture:
# we parse our new network (Haplogroups, sub haplogroups, archaeological culture)
	for y in aDNA_database:	
# we parse our database
		
		if x[0] == y[0]:
# if the first IDs of the network x[0] is == to one of the IDs in the database
			temp2 = x[0], x[1], y[13]
# we add our IDs(links) from the network plus the c14 dating associated to the samples x[0]
			samples_edges_col1.append(temp2)     
			
		elif x[1] == y[0]:
# if the second IDs of the network x[1] is == to one of the IDs in the database
			temp3 = x[0], x[1], y[13]
# we add our IDs(links) from the network plus the c14 dating associated to the samples y[2]
			samples_edges_col2.append(temp3)     
			
IDs_samples_no_need_reversed = []
IDs_samples_need_reversed = []



for i,v in enumerate (samples_edges_col1):
# we parse samples_edges_col1
	for i2,v2 in enumerate (samples_edges_col2):
# we parse samples_edges_col2		
		if i == i2 and v > v2:
# if their indexes are == and the dating of v is > v2
			temp4 = v[0:2]									  
			IDs_samples_no_need_reversed.append(list(temp4))
# we copy the IDs(the links) in a list 
			
		elif i == i2 and v < v2:                              
# if their indexes are == and the dating of v is < v2
			temp5 = v2[0:2]
			IDs_samples_need_reversed.append(list(temp5))     
# we copy the IDs(the links) for which we will need to reverse the IDs (To have the oldest sample, the youngest: to create a directed graph

		
for x in range(len(merge_HGs_subHGs_subHGsculture)):	
	if merge_HGs_subHGs_subHGsculture[x][0:2] in IDs_samples_need_reversed:           
# we reverse the linked IDs directly in the network list and we add 'true' (boolean value to create directed edges(links)
		merge_HGs_subHGs_subHGsculture[x] = ''.join(merge_HGs_subHGs_subHGsculture[x][1]+','+merge_HGs_subHGs_subHGsculture[x][0]+','+','.join(merge_HGs_subHGs_subHGsculture[x][2:])+','+'true').split(',')

	elif merge_HGs_subHGs_subHGsculture[x][0:2] in IDs_samples_no_need_reversed:      
		merge_HGs_subHGs_subHGsculture[x] += 'true'.split()                           
# for the links that we do not need to reverse the IDs (for the one that are already in the correct order), we add 'true' (boolean value to create directed edges(links)
		
	elif merge_HGs_subHGs_subHGsculture[x][0:2] not in IDs_samples_no_need_reversed or merge_HGs_subHGs_subHGsculture[x][0:2] not in IDs_samples_need_reversed:
		merge_HGs_subHGs_subHGsculture[x] += 'false'.split()                          
# if the samples have the same dating, we add false to not consider directed edges in that situation

################################################################################################												##########################################Genpofad Distance#####################################
################################################################################################


haplotypes_file_mtDNA = final_net.extraction_csv('mtDNA_haplotype.csv')
haplotypes_file_MSY = final_net.extraction_csv('Ychr_haplotypes.csv')
																		 
# for mtDNA and MSY we compute genpofad distance  + the number of mutations shared

P_ID_mtDNA, number_mut_mtDNA = P_IDs_functions_2.perId_uniparentalMarkers(merge_HGs_subHGs_subHGsculture, haplotypes_file_mtDNA) 
P_ID_MSY, number_mut_MSY = P_IDs_functions_2.perId_uniparentalMarkers(merge_HGs_subHGs_subHGsculture, haplotypes_file_MSY) 
                                                                        
per_ID_MSY_mtDNA = P_IDs_functions_2.perId_mtDNA_MSY(P_ID_mtDNA, P_ID_MSY) 
# for MSY and mtDNA we compute geometric distance based on genpofad score to have a single value

ids_share_mt_y = []
ids_share_mt = []
ids_share_y = []
for ids in merge_HGs_subHGs_subHGsculture:
	if len(ids[2].split('_')) == 4:
		mt_y = ids[0]+','+ids[1]
		ids_share_mt_y.append(mt_y)

	elif len(ids[2].split('_')) < 4 and ids[2].split('_')[1] == 'Y':
		Y = ids[0]+','+ids[1]
		ids_share_y.append(Y)

	elif len(ids[2].split('_')) < 4 and ids[2].split('_')[1] == 'mt':
		mt = ids[0]+','+ids[1]
		ids_share_mt.append(mt)


mt_y = 'maternal_paternal'
per_ID_MSY_mtDNA_final = P_IDs_functions_2.split_perId_mtDNA_MSY(per_ID_MSY_mtDNA, ids_share_mt_y,mt_y)
mater = 'maternal'
P_ID_mtDNA_final = P_IDs_functions_2.split_perId_mtDNA_MSY(P_ID_mtDNA, ids_share_mt,mater)
pater='paternal'
P_ID_MSY_final = P_IDs_functions_2.split_perId_mtDNA_MSY(P_ID_MSY, ids_share_y,pater)
list_dics = [per_ID_MSY_mtDNA_final,P_ID_MSY_final,P_ID_mtDNA_final,number_mut_MSY,number_mut_mtDNA]

perIDs_number_mutations = {}
for dics in list_dics:
	for k,v in dics.items():
		perIDs_number_mutations.setdefault(k, []).append(str(v)) 
# in a new dic we merge all the values (percentage ID, number of mutations) associated to the same key
		

for k,v in perIDs_number_mutations.items():
	for x in range(len(merge_HGs_subHGs_subHGsculture)):
		
		if k == ','.join(merge_HGs_subHGs_subHGsculture[x][0:2]): 
# we add all the information from the dic perIDs_number_mutations to our network
			merge_HGs_subHGs_subHGsculture[x] += v

#################################################################################################
###########################################output files##########################################
#################################################################################################

nodes = []

for x in aDNA_database:
	nodes.append(','.join(x))

outfile_gdf=open("network.gdf","w")
outfile_edgeList=open("network_edges.txt","w")                                                  
# we write the header for the nodes (Basically all the information from the database)
outfile_gdf.write("nodedef>name VARCHAR, Publication VARCHAR, GroupID VARCHAR, Locality VARCHAR, Country VARCHAR, Latitude DOUBLE, Longitude DOUBLE, chrY_HGs VARCHAR, ChrY_subHGs VARCHAR, mtDNA_HGs VARCHAR, mtDNA_subHGs VARCHAR, calBC_start INT, calBC_end INT, calBC_average DOUBLE\n")
for x in nodes:
	outfile_gdf.write(str(x)+"\n")
# we write all nodes informations from the database	
# we write the header for the links 	
	
outfile_gdf.write("edgedef>node1 VARCHAR,node2 VARCHAR,link Haplogroup VARCHAR,link sub_haplogroup VARCHAR,link culture VARCHAR,directed BOOLEAN,weight float,link_types VARCHAR,MSY_mut_shared INT,mtDNA_mut_shared INT\n")	
for x in merge_HGs_subHGs_subHGsculture:
	outfile_gdf.write(','.join(x)+"\n")           
# we write the network (links)
	outfile_edgeList.write(','.join(x[0:2])+"\n") 
# we generate as well an links files for try and error

outfile_gdf.close()
outfile_edgeList.close()




#########################################################################################
###########################firstMetrics of the network###################################
#########################################################################################

G = nx.Graph()

filename = "network_edges.txt"
G = nx.read_edgelist(filename, delimiter=",",data = (("id",str),))
nbnodes = len(list(G.nodes()))
nbedges = len(list(G.edges()))

CCs = nx.connected_component_subgraphs(G);
count = 0
for CC in CCs:
	count += 1
	nbnodes = len(list(CC.nodes()))
	nbedges = len(list(CC.edges()))
	clust = 2.0*nbedges/(nbnodes*(nbnodes-1))
	print ("component",count,"nodes",nbnodes,"edges",nbedges,"clust",clust)
#nbedges = len(G.edges)

#for x in net_cult_sub_HGs_HGs_merged:
	#G = nx.Graph()
	#G = nx.read_edgelist(x[0:2])
	#CCs = nx.connected_components(G)
	#nbCCs = len(list(CCs))
	#print ('nbCCs', nbCCs)
