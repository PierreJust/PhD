#! /usr/bin/env python 2.7
# -*- coding: utf-8 -*-

#Deborah Diquelou and Pierre Justeau

import final_net
import copy
import P_IDs_functions_2



####################################################################
#############################Network################################
####################################################################


aDNA_database = final_net.extraction_csv('database_filtered_GCC.csv')

print("for mitochondrial haplogroup clustering, select number of the mitochondrial column: (for us) type 9")
arg1=int(input())
print("for Y chromosome haplogroup clustering, select number of the Y chromosome column: (for us) type 7")
arg2=int(input())
print("for mitochondrial sub haplogroup clustering, select the number of the mitochondrial column: (for us) type 10")
arg3=int(input())
print("for Y chromosome sub haplogroup clustering, select the number of the Y chromosome column: (for us) type 8")
arg4=int(input())
print("for archaeological culture clustering, select the number of the archaeological culture column (for us) type 2")
arg5=int(input())
print("for country, select the number of the country column (for us) type 4")
arg6=int(input())

net_HGs = {}          
# Create dictionnaries for haplo, sub haplo and archeo 
net_sub_HGs = {}
net_archae_cult = {}
net_country = {}
    
for i in aDNA_database:           
# For each individual in the databasis, the key is the its ID for the final dictionnaries created above.  Rows are lists, we take HGs and archeo of the individual i.
                                                               # Finally, the dictionaries contains for each individuals (=keys), their HG and archeo.
    keys = i[0]
    
    HGs_mt = i[arg1]     # 9, ninth column of the csv file
    HGs_MSY = i[arg2]    # 11, eleventh column of the csv file 
    
    sub_HGs_mt = i[arg3] # 10, tenth column of the csv file 
    sub_HGs_MSY = i[arg4]# 12, twelfth column of the csv file 
    
    archae_cult = i[arg5]# 2, second column of the csv file
    
    country = i[arg6]
        
    net_HGs[keys] = HGs_mt, HGs_MSY                                                # We have 3 dictionnaries : HG mt/Y, subHG mt/Y, archeo fot each individuals.  + country
    net_sub_HGs[keys] = sub_HGs_mt, sub_HGs_MSY
    net_archae_cult[keys] = [archae_cult]
    net_country[keys] = [country]


Hgs_cluster = final_net.network(net_HGs)                        # ID1,ID2, shared HG etc (list)
sub_HGs_cluster = final_net.network(net_sub_HGs)
cult_cluster = final_net.network(net_archae_cult)
country_cluster = final_net.network(net_country)      # ID1, ID2, shared country

                                                                        ### we compare the clusters to know: ###
                                                                          # 1) Which samples share the same haplogroups and sub haplogroups
                                                                          # 2) For samples share the same haplogroups and sub haplogroups, do they share the same culture?

share_HGs_subHGs, list_share1 = final_net.share(Hgs_cluster, sub_HGs_cluster)       # ID1, ID2, shared HG, shared subHG
share_HGs_subHGs_culture, list_share2 = final_net.share(share_HGs_subHGs, cult_cluster)       # ID1, ID2, shared HG, shared subHG, shared culture (only those from share1 that have the same culture)
#share_HGs_subHGs_culture_country, list_share3 = final_net.share(share_HGs_subHGs_culture, country_cluster)       # ID1, ID2, HG, subHG, culture, country    list3 = IDs with same country

                                                                        ### here we merge our cluster to create a single network ###


merge_HGs_subHGsculture = final_net.merge_networks(Hgs_cluster, list_share1, share_HGs_subHGs_culture)             # here we merge the samples that share HGs,subHGs plus the samples that share only HGs
merge_HGs_subHGs_subHGsculture = final_net.merge_networks(share_HGs_subHGs, list_share2, merge_HGs_subHGsculture)  # here we add into the previous cluster the samples that share only HGs and sub HGs
# ID1, ID2, HG, subHG, culture

                                                                        ### to create directed graph based on the dating carbon 14 to have the links like: ###
                                                                          # 1) sample1(oldest), sample2(youngest) 
samples_edges_col1 = []
samples_edges_col2 = []

for x in merge_HGs_subHGs_subHGsculture :         # we parse our new network (Haplogroups, sub haplogroups, archaeological culture)
    for y in aDNA_database:                      # we parse our database
        
        if x[0] == y[0]:                         # if the first IDs of the network x[0] is == to one of the IDs in the database
            temp2 = x[0], x[1], y[11]            # we add our IDs(links) from the network plus the c14 dating associated to the samples x[0]
            samples_edges_col1.append(temp2)     
            
        elif x[1] == y[0]:                       # if the second IDs of the network x[1] is == to one of the IDs in the database
            temp3 = x[0], x[1], y[11]            # we add our IDs(links) from the network plus the c14 dating associated to the samples y[2]
            samples_edges_col2.append(temp3)     
            
IDs_samples_no_need_reversed = []
IDs_samples_need_reversed = []


for i,v in enumerate (samples_edges_col1):                    # we parse samples_edges_col1
    for i2,v2 in enumerate (samples_edges_col2):              # we parse samples_edges_col2
        
        if i == i2 and v > v2:                                # if their indexes are == and the dating of v is > v2
            temp4 = v[0:2]                                    
            IDs_samples_no_need_reversed.append(list(temp4))  # we copy the IDs(the links) in a list 
            
        elif i == i2 and v < v2:                              # if their indexes are == and the dating of v is < v2
            temp5 = v2[0:2]
            IDs_samples_need_reversed.append(list(temp5))     # we copy the IDs(the links) for which we will need to reverse the IDs (To have the oldest sample, the youngest: to create a directed graph

        
for x in range(len(merge_HGs_subHGs_subHGsculture)):
    
    if merge_HGs_subHGs_subHGsculture[x][0:2] in IDs_samples_need_reversed:           # we reverse the linked IDs directly in the network list and we add 'true' (boolean value to create directed edges(links)
        merge_HGs_subHGs_subHGsculture[x] = ''.join(merge_HGs_subHGs_subHGsculture[x][1]+','+merge_HGs_subHGs_subHGsculture[x][0]+','+','.join(merge_HGs_subHGs_subHGsculture[x][2:])+','+'true').split(',')

    elif merge_HGs_subHGs_subHGsculture[x][0:2] in IDs_samples_no_need_reversed:      
        merge_HGs_subHGs_subHGsculture[x] += 'true'.split()                           # for the links that we do not need to reverse the IDs (for the one that are already in the correct order), we add 'true' (boolean value to create directed edges(links)
        
    elif merge_HGs_subHGs_subHGsculture[x][0:2] not in IDs_samples_no_need_reversed or merge_HGs_subHGs_subHGsculture[x][0:2] not in IDs_samples_need_reversed:
        merge_HGs_subHGs_subHGsculture[x] += 'false'.split()                          # if the samples have the same dating, we add false to not consider directed edges in that situation



for i in merge_HGs_subHGs_subHGsculture :    # add countries to the information (id_country) + id to the culture
    id1 = i[0]
    id2 = i[1]
    for j in aDNA_database :
        if j[0] == id1 :
            country = id1, j[10]
            country2 = '-'.join(country)
            i.append(country2)
    for j in aDNA_database :
        if j[0] == id2 :
            country = id2, j[10]
            country2 = '-'.join(country)
            i.append(country2)
    n = len(i[4])           ### APPEND ID TO CULTURE
    if n>1 :
        culture = i[4]
        ind1 = i[0]
        ind2 = i[1]
        id_culture = ind1, ind2, culture
        id_culture_join = '-'.join(id_culture)
        i[4]= id_culture_join


countries = []

for i in aDNA_database :
    if (i[arg6] not in countries) :                          
 # create one key per unique country
        countries.append(i[arg6])        

##################################################################
############################GCC construction######################
##################################################################        


gcc = {}
gcc_id_for_nodes = {}  
# to store ids of gcc components
n = 0
already_classed_HGs_subHGs = []   
# dictionnary for individuals in nodes sharing both HGs and both subHGs
already_classed_HGs_oneSubHG = []   
# individuals sharing both HGs and one subHG
already_classed_HGs = []      
# individuals sharing both HGs but no subHGs

for i in merge_HGs_subHGs_subHGsculture :
    for country in countries :
        if country in i[6] and country in i[7]:       
# for individuals in the same country
            if '_' in i[2]:
# and (len(i[3]) == 0):#or '_' not in i[2]: # :          ############## if id share both mt and Y,  but no subHG shared
                n = i[2]
                t = country, n
                cle = '-'.join(t)           
# the key is "country-HGmt_HGY"
                if cle not in gcc.keys() :        
# create a dictionnary with HGs as keys and rows with the HGs as values
                    if i[0] not in already_classed_HGs_subHGs and i[1] not in already_classed_HGs_subHGs :
                        gcc[cle] = []
                        gcc[cle].append(i)
                        already_classed_HGs.append(i[0])
                        already_classed_HGs.append(i[1])
                        ind1 = i[0]
                        ind2 = i[1]
                        ind = ind1, ind2
                        ind_gcc = '-'.join(ind)
                        gcc_id_for_nodes[cle] = []
                        gcc_id_for_nodes[cle].append(ind_gcc)
                else :
                    if i[0] not in already_classed_HGs_subHGs and i[1] not in already_classed_HGs_subHGs :
                        gcc[cle].append(i)
                        already_classed_HGs.append(i[0])
                        already_classed_HGs.append(i[1])
                        ind1 = i[0]
                        ind2 = i[1]
                        ind = ind1, ind2
                        ind_gcc = '-'.join(ind)
                        gcc_id_for_nodes[cle].append(ind_gcc)


gcc_merge = {}

for cle in gcc.keys() :
    if cle not in gcc_merge :
        gcc_merge[cle] = []
    gcc_merge[cle] = gcc[cle]
    if len(gcc[cle]) > 1 :
        values_HG = gcc_merge[cle]
        first_list = values_HG[0]
        name1 = first_list[0]
        name2 = first_list[1]
        hg = first_list[2]
        subhg = first_list[3]
        culture = first_list[4]
        country1 = first_list[6]
        country2 = first_list[7]
        n= len(values_HG)
        for j in range(1, n) :
            next_list = values_HG[j]
            first_list[0] = name1 + '_' + next_list[0]
            name1 = first_list[0]
            first_list[1] = name2 + '_' + next_list[1]
            name2 = first_list[1]
            if next_list[2] not in hg :
                first_list[2] = hg + '_' +  next_list[2]
                hg = first_list[2]
            if next_list[3] not in subhg :
                first_list[3] = subhg + '_' + next_list[3]
                subhg = first_list[3]
            if (len(culture) == 0 and len(next_list[4]) == 0) or (len(culture) > 0 and len(next_list[4]) == 0) :
                first_list[4] = culture
                culture = first_list[4]
            if (len(culture) == 0 and len(next_list[4]) > 0) or (len(culture) > 0 and len(next_list[4]) > 0) :
                first_list[4] = culture + '~' + next_list[4]
                culture = first_list[4]
            first_list[6] = country1 + '-' + next_list[6]
            country1 = first_list[6]
            first_list[7] = country2 + '-' + next_list[7]
            country2 = first_list[7]
            gcc_merge[cle] = values_HG[0]
    else :
        a = gcc_merge[cle]
        gcc_merge[cle] = a[0]

for i in gcc.keys() :        
# delete the gcc from merge_hg_subhg_culture
    for j in gcc[i] :
        if i in merge_HGs_subHGs_subHGsculture :
            merge_HGs_subHGs_subHGsculture.remove(i)

### replace id in the merge_hg_subhg_culture by the name of their gcc

counter = 0
id_to_replace = {}
for cle in gcc_merge.keys() :
    info = gcc_merge[cle]
    if len(info) == 1 :
        counter = counter + 1
        ind = info[0][0] + '_' + info[0][1]
        id_to_replace[counter] = (cle, ind)
    else :
        counter = counter + 1
        ind =  info[0] + '_' + info[1]
        id_to_replace[counter] = (cle, ind)       
# -> dictionnary with tuples (HGmt_Y, individuals of this gcc)


for i in merge_HGs_subHGs_subHGsculture :     
# compare id in the merge with id in the uplets
    has_been_renamed = 'no'
    for cle in id_to_replace :
        info = id_to_replace[cle]
        if i[0] in info[1] :
            i[0] = info[0]
            has_been_renamed = 'yes'
        if has_been_renamed == 'yes' :
            break 
# if the node is renamed, we break the loop to go to the next node, to avoid a node being renamed several times if the individual is present in several gcc
for i in merge_HGs_subHGs_subHGsculture :     
# compare id in the merge with id in the uplets
    has_been_renamed = 'no'
    for cle in id_to_replace :
        info = id_to_replace[cle]
        if i[1] in info[1] :
            i[1] = info[0]
            has_been_renamed = 'yes'
        if has_been_renamed == 'yes' :
            break 
# if the node is renamed, we break the loop to go to the next node, to avoid a node being renamed several times if the individual is present in several gcc

#########################################################################
###############################output files##############################
#########################################################################

nodes = []

for x in aDNA_database:
    nodes.append(','.join(x))

nodes2 = []
for n in nodes :     
# add a base node size of ten, with a second list nodes as intermediate (because didn't work without???)
    n+=',10,'
    nodes2.append(n)

nodes = copy.deepcopy(nodes2)

# create a list of id to delete in nodes
list_ID = []
for cle in id_to_replace :
    info = id_to_replace[cle]
    IDs = info[1]    
# ID1_ID2...
    if '_' not in IDs :
        list_ID.append(IDs)
    if '_' in IDs :
        names = IDs.split('_')
        for x in names :
            if x not in list_ID :
                list_ID.append(x)

for name in list_ID :
    for elt in nodes :
        database = elt.split(',')
        if name == database[0] :
            nodes.remove(elt)


list_cle = []  
# used next loop
for i in merge_HGs_subHGs_subHGsculture :
    key1 = i[0]
    key2 = i[1]
    if key1 not in list_cle :
        list_cle.append(key1)
    if key2 not in list_cle :
        list_cle.append(key2)
        
identity = {}
# add a line for each gcc
for cle in gcc_merge :
    if cle in list_cle :
        # we keep only the GCCs whom nodes have been renamed after (cf lines 412-429). Otherwise, we would have nodes with no connections.
        
        line_to_add = gcc_merge[cle]
        if len(line_to_add) == 1 :
            line_to_add = line_to_add[0]
        name = cle
        culture = ''
        
        # extract mt DNA and Y DNA
        if '_' in line_to_add[2] :             
# if two HGs are shared
            hgDNA = line_to_add[2]
            hg_info = hgDNA.split('_')     
# 0 = mtHG and 1 = Yhg
            mtDNAhg = hg_info[0]
            YDNA = hg_info[1]
        else :                                    
# if one HG is shared
            if 'mt' in line_to_add[2] :
                mtDNAhg = line_to_add[2]
                YDNA = ''
            if 'Y' in line_to_add[2] :
                YDNA = line_to_add[2]
                mtDNAhg = ''
                
        
        # extract mt subDNA and Y subDNA
        if len(line_to_add[3]) > 1 :       
# if the column is not empty -> subHG shared
            hgsubDNA = line_to_add[3]
            if 'mt' in hgsubDNA and 'Y' in hgsubDNA :
                subhg_info = hgsubDNA.split('_')     
# 0 = mtHG and 1 = Yhg
                submt_alone = subhg_info[0]
                suby_alone = subhg_info[1]
            else :
                if 'mt' in hgsubDNA :
                    submt_alone = hgsubDNA
                    suby_alone = ''
                if 'Y' in hgsubDNA :
                    submt_alone = ''
                    suby_alone = hgsubDNA
        else :
            submt_alone = ''
            suby_alone = ''
        
        mtDNAsubHG = submt_alone
        YDNAsubHG = suby_alone
        publication = ''
        
       
        geographical_area = ''
        
        #extract country
        ID_country = line_to_add[6]
        country_sep = ID_country.split('-')    
# 0 = ID1       1 = country
        country = country_sep[1]

        bpend = ''
        
        for dble in id_to_replace.values() :        
# find the corresponding gcc in id_to_replace to calculate the number of individuals and hence the node size
            if cle == dble[0] :
                #print cle,dble
                size_large = dble[1]
                size_id = size_large.split('_')        
                sizeInt = len(size_id) * 10
                size = str(sizeInt)
          
        # extract list of individuals in the gcc
        for cle2 in gcc_id_for_nodes :
            if cle == cle2 :
                individuals = gcc_id_for_nodes[cle]
                forID = individuals
                if len(individuals) > 1 :    
# when more than 2, they are in a list of links [ind1-ind2, ind3-ind4...]
                    samples = []
                    for ident in individuals :
                        names = ident.split('-')
                        name1 = names[0]
                        name2 = names[1]
                        if name1 not in samples :
                            samples.append(name1)
                        if name2 not in samples :
                            samples.append(name2)
                    individuals = '-'.join(samples)
                else :      
# individuals = [ind1-ind2]
                    samples = []
                    names = individuals[0].split('-')
                    name1 = names[0]
                    name2 = names[1]
                    indiv = name1, name2
                    individuals = '-'.join(indiv)
                    if name1 not in samples:
                        samples.append(name1)
                    if name2 not in samples:
                        samples.append(name2)

        
        # list hg and sub hg in the gcc
        list_hg_subhg = []
        list_hg_mt = []
        list_hg_Y = []
        list_subhg_mt = []
        list_subhgY = []
        for a in aDNA_database :
            for b in samples :
                if a[0]==b :
                    list_hg_subhg.append(a[3:7])
                    list_hg_mt.append(a[3])
                    list_hg_Y.append(a[5])
                    list_subhg_mt.append(a[4])
                    list_subhgY.append(a[6])
        list_hg_sub_str = str(list_hg_subhg)
        str2 = list_hg_sub_str.replace('[', '')
        str3 = str2.replace(']', '')
        str4 = str3.replace(', ', '_')
        hg_subhg = str4.replace("'", "")
        
        list_hg_mt_str = str(list_hg_mt)
        str22 = list_hg_mt_str.replace('[', '')
        str32 = str22.replace(']', '')
        str42 = str32.replace(', ', '_')
        hg_mt = str42.replace("'", "")
        
        list_hg_Y_str = str(list_hg_Y)
        str23 = list_hg_Y_str.replace('[', '')
        str33 = str23.replace(']', '')
        str43 = str33.replace(', ', '_')
        hg_Y = str43.replace("'", "")
        
        list_subhg_mt_str = str(list_subhg_mt)
        str24 = list_subhg_mt_str.replace('[', '')
        str34 = str24.replace(']', '')
        str44 = str34.replace(', ', '_')
        subhg_mt = str44.replace("'", "")
        
        list_subhgY_str = str(list_subhgY)
        str25 = list_subhgY_str.replace('[', '')
        str35 = str25.replace(']', '')
        str45 = str35.replace(', ', '_')
        subhgY = str45.replace("'", "")

        # find gps coordonates
        latt = []   #column 11
        lngi = []    #column 12
        names = individuals.split('-')      
# gives a list of IDs
        for i in names :
            for a in aDNA_database :
                if a[0] == i :
                    latt.append(float(a[8]))
                    lngi.append(float(a[9]))
        lat_sum = sum(latt) 
        lat = str(lat_sum/len(latt))
        lng_sum = sum(lngi)   
        lng = str(lng_sum/len(lngi)) 
        
        # datation
        datations = []
        for i in names :
			for a in aDNA_database :
				if a[0] == i :
					age = a[11]
					datations.append(float(age))

        n = len(datations)           
        if n == 0 :
            datation = ''
        if n >0 :
            date = sum(datations)/n
            datation=str(date)
        
        #period
        indiv_list = individuals.split('-')
        truc = []
        for indi in indiv_list :
            for a in aDNA_database :
                if a[0] == indi :
                    period_a = a[1]
                    if period_a not in truc :
                        truc.append(period_a)
                        truc.append(indi)
                    else :
                        truc.append(indi)
        period = '-'.join(truc)
        #geographical_area bpend
        to_join = name, period,publication, mtDNAhg, mtDNAsubHG, YDNA, YDNAsubHG,culture, lat, lng,country , datation, size, individuals, hg_mt, hg_Y, subhg_mt, subhgY
        complete_gcc = ','.join(to_join)
        nodes.append(complete_gcc)


outfile_gdf=open("networkGCC.gdf","w")
outfile_edgeList=open("Net_3_median_calBC_2_4_removedBad_mt_LN_K1a_J1c_T2b_subHGs_GCC_2.txt","w")
     
                                                  # we write the header for the nodes (Basically all the information from the database)
outfile_gdf.write("nodedef>name VARCHAR, Period VARCHAR,publications VARCHAR, mtDNAHGs VARCHAR, mtDNA_subHGs VARCHAR, chrYHGs VARCHAR, chrY_subHGs VARCHAR,Culture VARCHAR, Latitude DOUBLE, Longitude DOUBLE,Country VARCHAR,datation average DOUBLE, Size INT, Individuals VARCHAR, list_HGs_mt VARCHAR, list_HGs_Y VARCHAR, list_subHGs_mt VARCHAR, list_subHGs_Y VARCHAR\n")


for x in nodes:
    #print x
    outfile_gdf.write(str(x)+"\n")                # we write all nodes informations (from the database)
                                                  # we write the header for the links
                                                  
outfile_gdf.write("edgedef>node1 VARCHAR,node2 VARCHAR,link Haplogroup VARCHAR,link sub_haplogroup VARCHAR,link culture VARCHAR,directed BOOLEAN,country1 VARCHAR,country2 VARCHAR,link type VARCHAR\n")



def link_type(network) :      
# determine the type of the link between individuals to get different edge colours in Gephi
    edge_type = ""
    if '_' in network[2] and '_' in network[3] :
        edge_type = 'Both HGs and subHGs'
    if '_' in network[2] and len(network[3]) == 0 :
        edge_type = 'Only HGs'
    if '_' in network[2] and 'Y' in i[3] and 'mt' not in network[3] :
        edge_type = 'Both HGs and Y subHG'
    if '_' in network[2] and 'mt' in network[3] and 'Y' not in network[3] :
        edge_type = 'Both HGs and mt subHG'
    if 'mt' in network[2] and '_' not in network[2] and 'mt' in network[3] and '_' not in network[3] :
        edge_type = 'Maternal link for both HG and subHG'
    if 'Y' in network[2] and '_' not in network[2] and 'Y' in network[3] and '_' not in network[3] :
        edge_type = 'Paternal link for both HG and subHG'
    if 'Y' in network[2] and len(network[3]) == 0 :
        edge_type = 'Paternal HG'
    if 'mt' in network[2] and len(network[3]) == 0 :
        edge_type = 'Maternal HG'
        
    return(edge_type)
       
for x in merge_HGs_subHGs_subHGsculture:
    if x[0] != x[1] :
        outfile_gdf.write(','.join(x)+"\n")           
# we write the network (links)
        outfile_edgeList.write(','.join(x[0:2])+"\n") 
# we generate as well a links files for try and error

outfile_gdf.close()
outfile_edgeList.close()

