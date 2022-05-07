#! /usr/bin/env python 3.7
# -*- coding: utf-8 -*-


#############################################################
# This script is a draft which needs a lot of clarification##
#############################################################

#Pierre Justeau
import random
from random import randint, sample
import re
import csv
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, CircleFace, TextFace, PieChartFace, COLOR_SCHEMES, add_face_to_node
import os.path
import glob
import collections 
from collections import Counter
import itertools

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



list_samples_I2a = extraction_csv("list_samples_I2a_annotated.txt")
published = []
KD = []
test = []
for a in list_samples_I2a:
	if re.match("KD",'_'.join(a).split('_')[0]):
		KD.append(''.join(a))
	else:
		published.append(''.join(a))


t = Tree("pathPhynder_output_tree.nwk")
	
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

list_samples_annotated  = []
list_samples_annotated_IdsHGs_only = []
list_main_HGs_allSamples = []
for leaf in t.traverse():
	if len(leaf.name) > 0:
		tempo = leaf.name.split('_')[0]+","+leaf.name
		list_samples_annotated.append(tempo)
		list_samples_annotated_IdsHGs_only.append(leaf.name)
		
list_main_HGs_allSamples = []		
for po in list_samples_annotated_IdsHGs_only:
	if po.split('_')[1][0:3] not in list_main_HGs_allSamples:
		list_main_HGs_allSamples.append(po.split('_')[1][0:3])
		
try1 = extraction_csv("guide_tree_ISOGG_I2a.txt")
dic_annotate_HGs_nodeTree_to_add_SNPs = {}

for bn in try1:
	for bn2 in bn[1:]:
		for choco in list_samples_annotated_IdsHGs_only:
			if re.match(bn2,choco.split('_')[1][0:]):#marche pas
				nodes_def = bn[0]
				if nodes_def in dic_annotate_HGs_nodeTree_to_add_SNPs.keys():
					dic_annotate_HGs_nodeTree_to_add_SNPs[nodes_def].append(choco)
				else:
					dic_annotate_HGs_nodeTree_to_add_SNPs[nodes_def]=[choco]
			
			elif len(bn2) > 3 and re.match(bn2,choco.split('_')[1][0:]):#bn2 == choco.split('_')[1][0:] :
				nodes_def2 = bn[0]
				if nodes_def2 in dic_annotate_HGs_nodeTree_to_add_SNPs.keys():
					dic_annotate_HGs_nodeTree_to_add_SNPs[nodes_def2].append(choco)
				else:
					dic_annotate_HGs_nodeTree_to_add_SNPs[nodes_def2]=[choco]
			if len(bn2) < 3 and re.match(bn2,choco.split('_')[1][0:]):#bn2 == choco.split('_')[1][0:] :
				nodes_def3 = bn[0]
#				print (nodes_def3)
				if nodes_def3 in dic_annotate_HGs_nodeTree_to_add_SNPs.keys():
					dic_annotate_HGs_nodeTree_to_add_SNPs[nodes_def3].append(choco)
				else:
					dic_annotate_HGs_nodeTree_to_add_SNPs[nodes_def3]=[choco]
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

#Ici on a un dico ou les clés sont les haplogroupes pour annoter les noeuds de l'arbre et les valeurs sont les liste d'individus associé avec get common_ancestor
#le but sera de matcher les mutations de chacun de mes individus avec les cles des ce dico

##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

samples_SNPs = extraction_multiple_coma_file("SNPs_noDeaminations/")#20736 SNPs
samples_SNPs_uniq_mutations = []

for uniqSNPs in samples_SNPs:
	if uniqSNPs[1:4] not in samples_SNPs_uniq_mutations:
		samples_SNPs_uniq_mutations.append(uniqSNPs[1:4]) #2316 uniq SNPs

list_deaminations = []
list_uncertain_positions = [] #~
uniq_HGs = []
SNPs_shared_and_notShared = {}
Samples_SNPs_list_pos = []

for f in samples_SNPs:
	Samples_SNPs_list_pos.append(f[1])
	if f[3] not in  uniq_HGs :
		uniq_HGs.append(f[3])

	keys_pos_branches_HGs = '_'.join(f[1:4])
	if keys_pos_branches_HGs in SNPs_shared_and_notShared.keys():
		 SNPs_shared_and_notShared[keys_pos_branches_HGs].append(f[0])
	else:
		SNPs_shared_and_notShared[keys_pos_branches_HGs] = f[0].split()

SNPs_to_add_on_the_tree = []
list_SNPs_shared_and_notShared = []
list_SNPs_shared = []
list_SNPs_not_Shared = []
for k,v in SNPs_shared_and_notShared.items():
	tmp = k.split()+v
	list_SNPs_shared_and_notShared.append(tmp)
	if len(v) > 1 :
		#print (k)
		tmp_SNPs_shared = k.split()+v
		list_SNPs_shared.append(tmp_SNPs_shared)
	elif len(v) == 1 :
		#print (k)
		tmp_SNPs_not_shared = k.split()+v
		list_SNPs_not_Shared.append(tmp_SNPs_not_shared)

	for k2,v2 in dic_annotate_HGs_nodeTree_to_add_SNPs.items():
		if k.split('_')[2] == k2:
			tmp1 = k.split()+v
			if tmp1 not in SNPs_to_add_on_the_tree:
				SNPs_to_add_on_the_tree.append(tmp1)

private_SNPs = [x for x in list_SNPs_shared_and_notShared if x not in SNPs_to_add_on_the_tree]
private_SNPs_shared = [x2 for x2 in list_SNPs_shared if x2 not in SNPs_to_add_on_the_tree]
private_SNPs_not_shared = [x3 for x3 in list_SNPs_not_Shared if x3 not in SNPs_to_add_on_the_tree]

if len(private_SNPs_shared) + len(private_SNPs_not_shared) + len(SNPs_to_add_on_the_tree) == len(list_SNPs_shared_and_notShared):

	print ("private_SNPs_shared :",len(private_SNPs_shared),"+ private_SNPs_not_shared",len(private_SNPs_not_shared),"+ SNPs_to_add_on_the_tree",len(SNPs_to_add_on_the_tree), "= list_SNPs_shared_and_notShared",len(list_SNPs_shared_and_notShared))
else:
	print ("fuck!,not working")



final_private_SNPs_shared = []
for x4 in private_SNPs_shared:

	for y in x4[1:]:
		tmp_private_snps_shared = y,x4[0]
		final_private_SNPs_shared.append(tmp_private_snps_shared)

final_private_SNPs_not_shared = []
for x5 in private_SNPs_not_shared:
	tmp_private_snps_not_shared =x5[1],x5[0]
	final_private_SNPs_not_shared.append(tmp_private_snps_not_shared)

final_private_SNPs = final_private_SNPs_shared+final_private_SNPs_not_shared

if len(final_private_SNPs_shared) + len(private_SNPs_not_shared) == len(final_private_SNPs):
	print ("final_private_SNPs_shared + private_SNPs_not_shared = final_private_SNPs")


##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

samples_geo_1000Y = extraction_csv("1000G_Y_samples_info.csv")
Y_samples_geo = {}
Y1000_samples_IDs_IDsHGs = []
for leaf in t.traverse():
	for b in samples_geo_1000Y:
		#print (x[0])
		if leaf.name.split('_')[0] == b[0]:
			Y1000_samples_IDs_IDsHGs.append(leaf.name)
			IDS = leaf.name
			#print (x,leaf.name)
			if IDS in Y_samples_geo.keys():
				Y_samples_geo[IDS] = str(b[1])
			else:
				Y_samples_geo[IDS] = str(b[1])


ancient_samples_geo = extraction_csv("ancient_samples_info.csv")
ancient_samples_geo_dic = dict(ancient_samples_geo)
print (ancient_samples_geo_dic)


##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
	
samples_anc_states = extraction_multiple_coma_file("ancestral_state/")
samples_list_isogg2018 = extraction_csv("list_firstSubHGsFromLastDerivedSNPs_isogg_2018.csv")
samples_anc_states_plus1 = {}
for a1 in samples_anc_states:
	for b1 in samples_list_isogg2018:
		if a1[0] == b1[0].split('_')[0] and a1[3] == b1[1]:
			IDs_samples_anc = b1[0]

			if IDs_samples_anc in samples_anc_states_plus1.keys():
				samples_anc_states_plus1[IDs_samples_anc].append(a1[3])
			else:	
				samples_anc_states_plus1[IDs_samples_anc]=a1[3].split()
				



list_count_anc_pos_samples ={}
for kk,vv in samples_anc_states_plus1.items():
	vv = dict([(n2, vv.count(n2)) for n2 in set(vv)])
	for key2,val2 in vv.items():
		tmpIDs_keys = kk
		tmp22 = [key2,str(val2)]
		if tmpIDs_keys in list_count_anc_pos_samples.keys():
			list_count_anc_pos_samples[tmpIDs_keys].append('_'.join(tmp22))
		else:
			list_count_anc_pos_samples[tmpIDs_keys]='_'.join(tmp22).split()


file_list_count_anc_pos = extraction_csv("list_count_anc_pos.csv")

list_count_anc_pos = {v_test[0]:v_test[1:] for v_test in file_list_count_anc_pos}
final_anc = []
for k22,v22 in list_count_anc_pos_samples.items():
	for k3,v3 in list_count_anc_pos.items():
		if k22 == k3:
			tmp5 = [k22,','.join(v22),len(v22),len(v3)]
			final_anc.append(tmp5)
			
final_anc_annotated = []
for r in final_anc:
	if r[1].count('~') == r[2]:
		tp = [r[0]+','+r[1]+','+str(r[2])+','+str(r[3])+','+"≈"] 
		final_anc_annotated.append(tp)
		#final_anc_annotated.append("unshure")
	elif '~'not in r[1] and r[2] == 1:
		tp2 = [r[0]+','+r[1]+','+str(r[2])+','+str(r[3])+','+"half sure"]
		final_anc_annotated.append(tp2)
	elif r[1].count('~') >= r[2]/2:
		tp3 = [r[0]+','+r[1]+','+str(r[2])+','+str(r[3])+','+"half sure"]
		final_anc_annotated.append(tp3)
	elif r[2] == r[3]:
		tp4 = [r[0]+','+r[1]+','+str(r[2])+','+str(r[3])+','+"sure"]
		final_anc_annotated.append(tp4)

	elif r[1].count('~') <= r[2]/2:
		tp5 = [r[0]+','+r[1]+','+str(r[2])+','+str(r[3])+','+"sure"]
		final_anc_annotated.append(tp5)
	else:
		print ('need better classification',r)

out_csv_file = open("final_anc_annotated.csv","w")
for gh in final_anc_annotated:
	out_csv_file.write(','.join(gh)+"\n")
	print (','.join(gh)) 

##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

schema_names = COLOR_SCHEMES.keys()
def layout(node):
	if node.name == "KD059_I2a1a2a1a2[?]_I2a":
		name_face = TextFace(node.name,fgcolor="red",fsize=10)
		node.add_face(name_face, column=0, position='branch-right')	

	elif node.name in KD:
		name_face = TextFace(node.name,fgcolor="green",fsize=10)
		node.add_face(name_face, column=0, position='branch-right')
	
	if node.name in published:
		name_face = TextFace(node.name,fgcolor="blue",fsize=10)
		node.add_face(name_face, column=0, position='branch-right')

	if node.name in Y1000_samples_IDs_IDsHGs:
		name_face = TextFace(node.name,fgcolor="dark",fsize=10)
		node.add_face(name_face, column=0, position='branch-right')

		node.add_features(Y_samples_geo=Y_samples_geo.get(node.name,"none"))
		N = AttrFace("Y_samples_geo", fsize=10)
		if 'Y_samples_geo' in node.features:
			faces.add_face_to_node(N, node, column=0, position="branch-top")
		
	if node.name in ancient_samples_geo_dic.keys():
		node.add_features(ancient_samples_geo_dic=ancient_samples_geo_dic.get(node.name,"none"))
		N2 = AttrFace("ancient_samples_geo_dic", fsize=10)
		if 'ancient_samples_geo_dic' in node.features:
			faces.add_face_to_node(N2, node, column=0, position="branch-top")



	if node.name in list_count_anc_pos.keys():
		Ancestral_State = {}
		Ancestral_State[node.name] = "Ancestral State from the most derived SNPs"
		node.add_features(Ancestral_State=Ancestral_State.get(node.name,"none"))
		N4 = AttrFace("Ancestral_State",fsize=12)
		N4.inner_border.width= 1
		N4.inner_border.type = 1
		if 'Ancestral_State' in node.features:
			faces.add_face_to_node(N4, node, column=2, position="branch-right")
	

def tree_style():
	t3 = Tree("pathPhynder_output_tree.nwk")


	I2a1a2a = []
	I2a1b = []
	for leaf in t3.iter_leaves():
		if re.match("I2a1a2a",leaf.name.split('_')[1]):
			I2a1a2a.append(leaf.name.split('_')[0])

		elif re.match("I2a1b",leaf.name.split('_')[1]):
			I2a1b.append(leaf.name.split('_')[0])
	
	for k,v in dic_annotate_HGs_nodeTree_to_add_SNPs.items():
		Hgs_on_node_tree = TextFace(k)
		Hgs_on_node_tree.fsize = 20
		nHgs_on_node_tree = t3.get_common_ancestor(v)

		legend = TextFace("SNPs that define the branch"+':'+k,fsize=15)
		legend.margin_top = 10
		legend.margin_right = 10
		legend.margin_left = 10
		legend.margin_bottom = 10
		legend.border.width = 1
		nHgs_on_node_tree.add_face(legend, column=0, position = "branch-top")

		nHgs_on_node_tree.img_style["size"] = 50
		nHgs_on_node_tree.img_style["fgcolor"] = "Black"
		nHgs_on_node_tree.add_face(Hgs_on_node_tree, column=1, position = "branch-right")

		for x in SNPs_to_add_on_the_tree: 
			if x[0].split('_')[2] == k:
				SNPs_shared_and_notShared_mapped_on_tree = TextFace(','.join(x))
				nHgs_on_node_tree.add_face(SNPs_shared_and_notShared_mapped_on_tree, column=0, position = "branch-top")

	for x2 in final_private_SNPs:
		if x2[0] in I2a1a2a and re.match("I2a1a2a",x2[1].split('_')[2]):
			for leaf in t3.iter_leaves():		
				if x2[0].split('_')[0] == leaf.name.split('_')[0]:
					private_SNPs = TextFace(','.join(x2))
					private_SNPs.fgcolor = "Dark"
					leaf.add_face(private_SNPs, column=1, position = "branch-right")
		elif x2[0] in I2a1b and re.match("I2a1b",x2[1].split('_')[2]):
			for leaf in t3.iter_leaves():		
				if x2[0].split('_')[0] == leaf.name.split('_')[0]: 
					private_SNPs = TextFace(','.join(x2))
					private_SNPs.fgcolor = "Dark"
					leaf.add_face(private_SNPs, column=1, position = "branch-right")



	for ancState in final_anc_annotated:
		for leaf3 in t3.iter_leaves():
			if ','.join(ancState).split(',')[0] == leaf3.name:
				ancestralState = TextFace(','.join(ancState).split(',')[-1])
				leaf3.add_face(ancestralState, column=2, position = "branch-right")			
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

	style = NodeStyle()
	style["fgcolor"] = "black"
	style["size"] = 0
	style["vt_line_color"] = "black"
	style["hz_line_color"] = "black"
	style["vt_line_width"] = 7
	style["hz_line_width"] = 7
	style["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
	style["hz_line_type"] = 2
	for n in t3.traverse():	
		n.set_style(style)

	style2 = NodeStyle()
	style2["fgcolor"] = "#000000"
	style2["shape"] = "circle"
	style2["vt_line_color"] = "#0000aa"
	style2["hz_line_color"] = "#0000aa"
	style2["vt_line_width"] = 7
	style2["hz_line_width"] = 7
	style2["vt_line_type"] = 2 # 0 solid, 1 dashed, 2 dotted
	style2["hz_line_type"] = 2
	t3.children[0].img_style = style2
	for redBranch in t3.iter_leaves():

		if redBranch.name == "mid002_I2a1a2a1a2a[?]":
			redBranch.set_style(style2)

	ts = TreeStyle()
	ts.layout_fn = layout
	ts.show_leaf_name = False

	ts.root_opening_factor = 1
	return t3, ts

	
if __name__ == "__main__":
	t3, ts = tree_style()
	i2a1 = []
	for l in t3:
		if re.match('I2',l.name.split('_')[1]):
			i2a1.append(l.name)
	I2a1_ancestor = t3.get_common_ancestor(i2a1)
	I2a1_ancestor.show(tree_style=ts)

