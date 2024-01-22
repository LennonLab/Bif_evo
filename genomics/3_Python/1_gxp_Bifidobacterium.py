##get_multiplicity_matrix() written originally by Will Shoemaker, modified by Roy Moger-Reischer

from __future__ import division
import os, math, pickle
#import mult_tools as mt
import numpy as np
import pandas as pd

from decimal import Decimal

from pathlib import Path
homedir = str(Path.home())
workingdir = homedir + '\\Documents\\GitHub\\Bifidobacterium'
os.chdir(workingdir)

output_to_keep = ['INS', 'DEL', 'SNP', 'SUB']
strains = ['Bif']

#%%
# two arguments: minimal or wildtype

#May want to  actually hard code this, depending on any indications of possible early cross contamination
def get_sites_to_remove(strain):
    if strain == 'Bif':
        ancestor_path = workingdir + '\\Pipeline\\1_Carbonate\\breseq\\breseq_25\\output\\output.gd'
        #C:\Users\rmoge\Documents\GitHub\Bifidobacterium\Pipeline\1_Carbonate\breseq\breseq_25\output
    sites_to_remove = []
    for i, line in enumerate(open(ancestor_path, 'r')):
        line_split = line.strip().split('\t')
        if line_split[0] in output_to_keep:
            sites_to_remove.append( line_split[3] + '_' + str(line_split[4]))
# =============================================================================
#     Below, you will add sites that did no show up in the ancestor, but are obvious problematic sites, either due to cross contamination , or due to a mutation in a founding population for one of the treatments specifically
# =============================================================================
    sites_to_remove+=['CP001853_146344','CP001853_1305547','CP001853_1437318']
# =============================================================================
#     
# =============================================================================
    return sites_to_remove
#%%
def get_gxp_matrix(output_name,include_polymorphisms=True):

    gene_by_pop_dict = {}#initialize GxP matrix
    for strain in strains:
        sites_to_remove = get_sites_to_remove(strain)
        if strain == 'Bif':
            dirs = ['breseq_1','breseq_2','breseq_3','breseq_4','breseq_5','breseq_6','breseq_7','breseq_8','breseq_9','breseq_10','breseq_11','breseq_12','breseq_13','breseq_14','breseq_15','breseq_16','breseq_17','breseq_18','breseq_19','breseq_20','breseq_21','breseq_22','breseq_23','breseq_24']
            #dirs = ['breseq_1','breseq_17']
            ref_path = workingdir + '\\Pipeline\\1_Carbonate\\BB12_Jensen_2021.gb'
        

        #effective_gene_lengths, effective_gene_lengths_syn, Lsyn, Lnon, substitution_specific_synonymous_fraction = mt.calculate_synonymous_nonsynonymous_target_sizes(ref_path)
        for dir in dirs:
            print(dir)
            pop = dir.split('_')[1]
            gene_count_dict_pop = {}#initialize the dict of genes for THIS pop
            gene_by_pop_dict[pop] = {}#GxP matrix add "rows", i.e., sub-dicts
            for i, line in enumerate(open(workingdir +'\\Pipeline\\1_Carbonate\\breseq\\' + dir + '\\output\\evidence\\annotated.gd', 'r')):
                line_split = line.strip().split('\t')
                if line_split[0] not in output_to_keep:
                    continue
                if line_split[3] + '_' + line_split[4] in sites_to_remove:
                    continue
                frequency = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])
                if frequency != 1 and include_polymorphisms==False:
                    continue
                if line_split[0] == 'SNP':
                    if [s for s in line_split if 'snp_type=' in s][0].split('=')[1] == 'nonsynonymous' or [s for s in line_split if 'snp_type=' in s][0].split('=')[1] == 'nonsense':
                        locus_tag = [s for s in line_split if 'locus_tag=' in s][0].split('=')[1]
                        frequency = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])
                        if ';' in locus_tag:
                            for locus_tag_j in locus_tag.split(';'):
                                if locus_tag_j not in gene_count_dict_pop:
                                    gene_count_dict_pop[locus_tag_j] = 0
                                gene_count_dict_pop[locus_tag_j] += frequency
                        else:
                            if locus_tag not in gene_count_dict_pop:
                                gene_count_dict_pop[locus_tag] = 0
                            gene_count_dict_pop[locus_tag] += frequency
                    else:
                        continue                            
                ######################################
                              
                else:
                    if len([s for s in line_split if 'gene_position=coding' in s]) >= 1:
                        locus_tag = [s for s in line_split if 'locus_tag=' in s][0].split('=')[1]
                        frequency = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])
                        if frequency != 1 and include_polymorphisms == False:
                            continue
                        if ';' in locus_tag:
                            for locus_tag_j in locus_tag.split(';'):
                                if locus_tag_j not in gene_count_dict_pop:
                                    gene_count_dict_pop[locus_tag_j] = 0
                                gene_count_dict_pop[locus_tag_j] += frequency

                        else:
                            if locus_tag not in gene_count_dict_pop:
                                gene_count_dict_pop[locus_tag] = 0
                            gene_count_dict_pop[locus_tag] += frequency
                    else:        
                        if line_split[0] == 'DEL':
                            if len([s for s in line_split if 'gene_position' in s]) == 0:
                                
                                if len([s for s in line_split if 'locus_tags_inactivated=' in s]) > 0:
                                    locus_tags_comma_unseparated = [s for s in line_split if 'locus_tags_inactivated=' in s][0].split('=')[1]
                                    frequency = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])
                                    #print(locus_tags_comma_unseparated)
                                    if frequency != 1 and include_polymorphisms==False:
                                        continue
                                    #print(type(locus_tags_comma_unseparated))
                                    if locus_tags_comma_unseparated[-1]=='"':
                                        #print('it is 17')
                                        for locus_tag_j in locus_tags_comma_unseparated[:-1].split(','):
                                            if locus_tag_j not in gene_count_dict_pop:
                                                gene_count_dict_pop[locus_tag_j] = 0
                                            gene_count_dict_pop[locus_tag_j] += frequency
                                    else:
                                        for locus_tag_j in locus_tags_comma_unseparated.split(','):
                                            if locus_tag_j not in gene_count_dict_pop:
                                                gene_count_dict_pop[locus_tag_j] = 0
                                            gene_count_dict_pop[locus_tag_j] += frequency                                        
                                
                                elif len([s for s in line_split if 'locus_tags_overlapping=' in s]) > 0:        
                                    locus_tags_comma_unseparated = [s for s in line_split if 'locus_tags_overlapping=' in s][0].split('=')[1]
                                    frequency = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])
                                    if frequency != 1 and include_polymorphisms==False:
                                        continue
                                    #print(locus_tags_comma_unseparated)
                                    for locus_tag_j in locus_tags_comma_unseparated.split(','):
                                        if locus_tag_j not in gene_count_dict_pop:
                                            gene_count_dict_pop[locus_tag_j] = 0
                                        gene_count_dict_pop[locus_tag_j] += frequency
                            else:
                                continue
                        else:
                            continue
    #return gene_count_dict_pop
#            for locus_tag_i in gene_count_dict_pop.keys():
#                gene_by_pop_dict[pop][locus_tag_i] = gene_count_dict_pop[locus_tag_i]

            for locus_tag_i in gene_count_dict_pop.keys():
#                 mult_i = gene_parallelism_statistics[locus_tag_i]['multiplicity']
#                 if mult_i > 0:
                locus_tag_i_num = locus_tag_i.split('_')[1]
                gene_by_pop_dict[pop][locus_tag_i_num] = gene_count_dict_pop[locus_tag_i]
    
    #return gene_by_pop_dict
                
            
                        

                        
# =============================================================================
# 
#             gene_parallelism_statistics = {}
#             for gene_i, length_i in effective_gene_lengths.items():
#                 gene_parallelism_statistics[gene_i] = {}
#                 gene_parallelism_statistics[gene_i]['length'] = length_i
#                 gene_parallelism_statistics[gene_i]['observed'] = 0
#                 gene_parallelism_statistics[gene_i]['multiplicity'] = 0
# 
#             # save number of mutations for multiplicity
#             for locus_tag_i, n_muts_i in gene_count_dict_pop.items():
#                 gene_parallelism_statistics[locus_tag_i]['observed'] = n_muts_i
# 
#             # save number of mutations for multiplicity
#             L_mean = np.mean(list(effective_gene_lengths.values()))
#             L_tot = sum(list(effective_gene_lengths.values()))
#             n_tot = sum(gene_count_dict_pop.values())
#             # go back over and calculate multiplicity
#             for locus_tag_i in gene_parallelism_statistics.keys():
#                 # double check the measurements from this
#                 gene_parallelism_statistics[locus_tag_i]['multiplicity'] = gene_parallelism_statistics[locus_tag_i]['observed'] *1.0/ effective_gene_lengths[locus_tag_i] * L_mean
#                 gene_parallelism_statistics[locus_tag_i]['expected'] = n_tot*gene_parallelism_statistics[locus_tag_i]['length']/L_tot
# 
#             # split locus tags
#             for locus_tag_i in gene_parallelism_statistics.keys():
#                 mult_i = gene_parallelism_statistics[locus_tag_i]['multiplicity']
#                 if mult_i > 0:
#                     locus_tag_i_num = locus_tag_i.split('_')[1]
#                     gene_by_pop_dict[pop][locus_tag_i_num] = mult_i
# 
# =============================================================================
    gene_by_pop_df = pd.DataFrame(gene_by_pop_dict)#outputs the GxP matrix
    gene_by_pop_df = gene_by_pop_df.T
    gene_by_pop_df.fillna(0, inplace=True)

    gene_by_pop_df_out = workingdir + '\\Pipeline\\3_Python\\' + output_name + '.csv'
    #gene_by_pop_df.to_csv(gene_by_pop_df_out, sep = '\t', index = True)
    gene_by_pop_df.to_csv(gene_by_pop_df_out, sep = ',', index = True)
    
    
    return gene_by_pop_dict


#%%
test1=get_gxp_matrix('Bifidobacterium_raw_gxp',True)
test2=get_gxp_matrix('Bifidobacterium_raw_gxp_fixed.only',False)
#%%

def get_gxp_matrix_signif_loci_only(output_name,include_polymorphisms=True,signifList=[]):

    gene_by_pop_dict = {}#initialize GxP matrix
    for strain in strains:
        sites_to_remove = get_sites_to_remove(strain)
        if strain == 'Bif':
            dirs = ['breseq_1','breseq_2','breseq_3','breseq_4','breseq_5','breseq_6','breseq_7','breseq_8','breseq_9','breseq_10','breseq_11','breseq_12','breseq_13','breseq_14','breseq_15','breseq_16','breseq_17','breseq_18','breseq_19','breseq_20','breseq_21','breseq_22','breseq_23','breseq_24']
            #dirs = ['breseq_1','breseq_17']
            ref_path = workingdir + '\\Pipeline\\1_Carbonate\\BB12_Jensen_2021.gb'
        

        #effective_gene_lengths, effective_gene_lengths_syn, Lsyn, Lnon, substitution_specific_synonymous_fraction = mt.calculate_synonymous_nonsynonymous_target_sizes(ref_path)
        for dir in dirs:
            #print(dir)
            pop = dir.split('_')[1]
            gene_count_dict_pop = {}#initialize the dict of genes for THIS pop
            gene_by_pop_dict[pop] = {}#GxP matrix add "rows", i.e., sub-dicts
            for i, line in enumerate(open(workingdir +'\\Pipeline\\1_Carbonate\\breseq\\' + dir + '\\output\\evidence\\annotated.gd', 'r')):
                line_split = line.strip().split('\t')
                if line_split[0] not in output_to_keep:
                    continue
                if line_split[3] + '_' + line_split[4] in sites_to_remove:
                    continue
                frequency = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])
                if frequency != 1 and include_polymorphisms==False:
                    continue
                if line_split[0] == 'SNP':
                    if [s for s in line_split if 'snp_type=' in s][0].split('=')[1] == 'nonsynonymous' or [s for s in line_split if 'snp_type=' in s][0].split('=')[1] == 'nonsense':
                        locus_tag = [s for s in line_split if 'locus_tag=' in s][0].split('=')[1]
                        frequency = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])
                        if ';' in locus_tag:
                            for locus_tag_j in locus_tag.split(';'):
                                if locus_tag_j not in gene_count_dict_pop:
                                    gene_count_dict_pop[locus_tag_j] = 0
                                gene_count_dict_pop[locus_tag_j] += frequency
                        else:
                            if locus_tag not in gene_count_dict_pop:
                                gene_count_dict_pop[locus_tag] = 0
                            gene_count_dict_pop[locus_tag] += frequency
                    else:
                        continue                            
                ######################################
                              
                else:
                    if len([s for s in line_split if 'gene_position=coding' in s]) >= 1:
                        locus_tag = [s for s in line_split if 'locus_tag=' in s][0].split('=')[1]
                        frequency = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])
                        if frequency != 1 and include_polymorphisms == False:
                            continue
                        if ';' in locus_tag:
                            for locus_tag_j in locus_tag.split(';'):
                                if locus_tag_j not in gene_count_dict_pop:
                                    gene_count_dict_pop[locus_tag_j] = 0
                                gene_count_dict_pop[locus_tag_j] += frequency

                        else:
                            if locus_tag not in gene_count_dict_pop:
                                gene_count_dict_pop[locus_tag] = 0
                            gene_count_dict_pop[locus_tag] += frequency
                    else:        
                        if line_split[0] == 'DEL':
                            if len([s for s in line_split if 'gene_position' in s]) == 0:
                                
                                if len([s for s in line_split if 'locus_tags_inactivated=' in s]) > 0:
                                    locus_tags_comma_unseparated = [s for s in line_split if 'locus_tags_inactivated=' in s][0].split('=')[1]
                                    frequency = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])
                                    #print(locus_tags_comma_unseparated)
                                    if frequency != 1 and include_polymorphisms==False:
                                        continue
                                    #print(type(locus_tags_comma_unseparated))
                                    if locus_tags_comma_unseparated[-1]=='"':
                                        #print('it is 17')
                                        for locus_tag_j in locus_tags_comma_unseparated[:-1].split(','):
                                            if locus_tag_j not in gene_count_dict_pop:
                                                gene_count_dict_pop[locus_tag_j] = 0
                                            gene_count_dict_pop[locus_tag_j] += frequency
                                    else:
                                        for locus_tag_j in locus_tags_comma_unseparated.split(','):
                                            if locus_tag_j not in gene_count_dict_pop:
                                                gene_count_dict_pop[locus_tag_j] = 0
                                            gene_count_dict_pop[locus_tag_j] += frequency                                        
                                
                                elif len([s for s in line_split if 'locus_tags_overlapping=' in s]) > 0:        
                                    locus_tags_comma_unseparated = [s for s in line_split if 'locus_tags_overlapping=' in s][0].split('=')[1]
                                    frequency = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])
                                    if frequency != 1 and include_polymorphisms==False:
                                        continue
                                    #print(locus_tags_comma_unseparated)
                                    for locus_tag_j in locus_tags_comma_unseparated.split(','):
                                        if locus_tag_j not in gene_count_dict_pop:
                                            gene_count_dict_pop[locus_tag_j] = 0
                                        gene_count_dict_pop[locus_tag_j] += frequency
                            else:
                                continue
                        else:
                            continue
    #return gene_count_dict_pop
#            for locus_tag_i in gene_count_dict_pop.keys():
#                gene_by_pop_dict[pop][locus_tag_i] = gene_count_dict_pop[locus_tag_i]

            for locus_tag_i in gene_count_dict_pop.keys():
#                 mult_i = gene_parallelism_statistics[locus_tag_i]['multiplicity']
#                 if mult_i > 0:
                locus_tag_i_num = locus_tag_i.split('_')[1]
                print(locus_tag_i_num)
                if locus_tag_i_num in signifList:
                    gene_by_pop_dict[pop][int(locus_tag_i_num)] = gene_count_dict_pop[locus_tag_i]
                else:
                    pass
    

    gene_by_pop_df = pd.DataFrame(gene_by_pop_dict)#outputs the GxP matrix
    gene_by_pop_df = gene_by_pop_df.T
    gene_by_pop_df.fillna(0, inplace=True)

    gene_by_pop_df_out = workingdir + '\\Pipeline\\3_Python\\' + output_name + '.csv'
    #gene_by_pop_df.to_csv(gene_by_pop_df_out, sep = '\t', index = True)
    gene_by_pop_df.to_csv(gene_by_pop_df_out, sep = ',', index = True)
    
    
    return gene_by_pop_dict

#%%
signifListBifidoInts=[489,532,625,684,1197,1492,778,2305,490,651,1639,936,1789,1346,1467,892,1831,1327,1616,1060,944,208,1023,327,1193,776,1402,573,2090];signifListBifidoBadLengths=[str(l) for l in signifListBifidoInts]; signifListBifido=[]
for l in signifListBifidoBadLengths:
    if len(l)==1:
        signifListBifido.append("0000"+l)
    elif len(l)==2:
        signifListBifido.append("000"+l)
    elif len(l)==3:
        signifListBifido.append("00"+l)
    elif len(l)==4:
        signifListBifido.append("0"+l)
    elif len(l)==5:
        signifListBifido.append(l)
    else:
        pass
    
test1=get_gxp_matrix_signif_loci_only('Bifidobacterium_raw_gxp_signif.only',True,signifListBifido)
test1=get_gxp_matrix_signif_loci_only('Bifidobacterium_raw_gxp_fixed.only_signif.only',False,signifListBifido)
print(signifListBifido)