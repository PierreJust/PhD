#! /usr/bin/env python 3.7
# -*- coding: utf-8 -*-

import re
import csv
import os.path

final_tree_newick = []

def extraction_csv(csv_file):
	
    file_csv = []
    
    with open(csv_file, 'r') as fi:
        reader = csv.reader(fi, delimiter=",")
        for rows in reader:
            file_csv.append(rows)
            
    return (file_csv)


def extraction_tsv(csv_file):
	
    file_csv = []
    
    with open(csv_file, 'r') as f:
        reader = csv.reader(f, delimiter="\t")
        for rows in reader:
            file_csv.append(rows)
            
    return (file_csv)

isogg2018 = extraction_tsv("snps_isogg2018_fixed.txt")
samples = extraction_csv("list_samples_I2a.txt")


def extraction_multiple_coma_file(path_to_file):
	folder_path_MSY = path_to_file
	MSY_haplotype_file = []
	tmp_haplo_file = []
	for path, dirs, files in os.walk(folder_path_MSY):
		for filename in files:
			if filename.endswith(".txt"):
				name = os.path.basename(filename).split('.')#supprimer '.' ici pour supprimer extension.txt
				tmp = open(os.path.join(path, filename ),"r")
				reader = csv.reader(tmp, delimiter = ',')
				for line in reader:
					tmp2 = line	
					tmp3 = name[0].split() + tmp2
				
					tmp_haplo_file.append(tmp3)
	return(tmp_haplo_file)


list_firstSubHGsFromLastDerivedSNPs_isogg_2018 = {}
for x in samples:
	for y in isogg2018:
		if re.search(''.join(x).split('_')[1].split('~')[0],y[1].split('~')[0]) and len(y[1].split('~')[0]) == len(''.join(x).split('_')[1].split('~')[0])+1:
			keys_IDs = ''.join(x)
			if keys_IDs in list_firstSubHGsFromLastDerivedSNPs_isogg_2018.keys():
				list_firstSubHGsFromLastDerivedSNPs_isogg_2018[keys_IDs].append(y[1])
			else:
				list_firstSubHGsFromLastDerivedSNPs_isogg_2018[keys_IDs]=y[1].split()

list_count_anc_pos = {}
list_count_anc_pos_test = []
sampleIDs_number_match_isogg2018 = []
for cf,vf in list_firstSubHGsFromLastDerivedSNPs_isogg_2018.items():
	vf = dict([(n, vf.count(n)) for n in set(vf)])
	for key,val in vf.items():
		tmp = [key,str(val)]
		tmp2 = [cf,key,val]
		Keys_IDs_dic = cf
		list_count_anc_pos_test.append(tmp2)
		sampleIDs_number_match_isogg2018.extend([cf])
		if Keys_IDs_dic in list_count_anc_pos.keys():
			list_count_anc_pos[Keys_IDs_dic].append('_'.join(tmp))
		else:
			list_count_anc_pos[Keys_IDs_dic] = '_'.join(tmp).split()


final_list_count_anc_pos = []
for k,v in list_count_anc_pos.items():
	for p in list_count_anc_pos_test:
		if k == p[0]:
			tmp3 = [p[0],p[1],p[2],len(v)]
			final_list_count_anc_pos.append(tmp3)


FICHIER_list_firstSubHGsFromLastDerivedSNPs_isogg_2018 = open("list_firstSubHGsFromLastDerivedSNPs_isogg_2018.csv","w")
for d in final_list_count_anc_pos:
	tmp4 = d[0]+','+d[1]+','+str(d[2])+','+str(d[3])
	FICHIER_list_firstSubHGsFromLastDerivedSNPs_isogg_2018.write(tmp4+"\n")
FICHIER_list_firstSubHGsFromLastDerivedSNPs_isogg_2018.close()

fichier_list_count_anc_pos = open("list_count_anc_pos_2.csv","w")
for k10,v10 in list_count_anc_pos.items():
	relou = k10+','+','.join(v10)
	fichier_list_count_anc_pos.write(relou+"\n")
fichier_list_count_anc_pos.close()

samples_anc_states = extraction_multiple_coma_file("ancestral_state/")
samples_list_isogg2018 = extraction_csv("list_firstSubHGsFromLastDerivedSNPs_isogg_2018.csv")
samples_anc_states_plus1 = {}
for a in samples_anc_states:
	for b in samples_list_isogg2018:
		if a[0] == b[0].split('_')[0] and a[3] == b[1]:
			IDs_samples_anc = b[0]

			if IDs_samples_anc in samples_anc_states_plus1.keys():
				samples_anc_states_plus1[IDs_samples_anc].append(a[3])
			else:	
				samples_anc_states_plus1[IDs_samples_anc]=a[3].split()
				



list_count_anc_pos_samples ={}
for kk,vv in samples_anc_states_plus1.items():
	vv = dict([(n2, vv.count(n2)) for n2 in set(vv)])
	for key2,val2 in vv.items():
		tmpIDs_keys = kk
		tmp2 = [key2,str(val2)]
		if tmpIDs_keys in list_count_anc_pos_samples.keys():
			list_count_anc_pos_samples[tmpIDs_keys].append('_'.join(tmp2))
		else:
			list_count_anc_pos_samples[tmpIDs_keys]='_'.join(tmp2).split()

final = []
for k2,v2 in list_count_anc_pos_samples.items():
	for k3,v3 in list_count_anc_pos.items():
		if k2 == k3:
			tmp5 = [k2,','.join(v2),len(v2),len(v3)]
			final.append(tmp5)


for r in final:

	if r[2] == r[3]:
		print ('SURE1',r)
	elif r[1].count('~') == r[2]:
		print ('unshure',r)

	elif '~'not in r[1] and r[2] == 1:
		print ('half sure',r)
	
	elif r[1].count('~') >= r[2]/2:
		print ('half sure2',r)


	elif r[1].count('~') <= r[2]/2:
		print ('SURE2',r)

	else:
		print ('need better classification',r)



