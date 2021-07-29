from __future__ import division
import os, math, pickle
#import mult_tools as mt
import numpy as np
import pandas as pd

from decimal import Decimal

output_to_keep = ['INS', 'DEL', 'SNP', 'SUB']
strains = ['wildtype', 'minimal']


# two arguments: minimal or wildtype

#May want to  actually hard code this, depending on any indications of possible early cross contamination
def get_sites_to_remove(strain):
    if strain == 'minimal':
        acnestor_path = mt.get_path() + '/data/syn3B_minimal/mmW_3B.ancestor/output.gd'
    elif strain == 'wildtype':
        acnestor_path = mt.get_path() + '/data/syn1.0_wildtype/mm8_syn1.0.ancestor/output.gd'
    sites_to_remove = []
    for i, line in enumerate(open(acnestor_path, 'r')):
        line_split = line.strip().split('\t')
        if line_split[0] in output_to_keep:
            sites_to_remove.append( line_split[3] + '_' + str(line_split[4]))
    return sites_to_remove

def get_multiplicity_matrix():

    gene_by_pop_dict = {}#initialize GxP matrix
    for strain in strains:
        sites_to_remove = get_sites_to_remove(strain)
        if strain == 'minimal':
            dirs = ['syn3B_minimal/mm13', 'syn3B_minimal/mm11', 'syn3B_minimal/mm10', 'syn3B_minimal/mm9']
            ref_path = mt.get_path() + '/data/syn3B_minimal/reference/Synthetic.bacterium_JCVI-Syn3A.gb'
        elif strain == 'wildtype':
            dirs = ['syn1.0_wildtype/mm6', 'syn1.0_wildtype/mm4', 'syn1.0_wildtype/mm3', 'syn1.0_wildtype/mm1']
            ref_path = mt.get_path() + '/data/syn1.0_wildtype/reference/Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.gb'

        effective_gene_lengths, effective_gene_lengths_syn, Lsyn, Lnon, substitution_specific_synonymous_fraction = mt.calculate_synonymous_nonsynonymous_target_sizes(ref_path)
        for dir in dirs:
            pop = dir.split('/')[1]
            gene_count_dict_pop = {}
            gene_by_pop_dict[pop] = {}#GxP matrix add "rows", i.e., sub-dicts
            for i, line in enumerate(open(mt.get_path()+'/data/'+dir+'/annotated.gd', 'r')):
                line_split = line.strip().split('\t')
                if line_split[0] not in output_to_keep:
                    continue
                if line_split[3] + '_' + line_split[4] in sites_to_remove:
                    continue
                frequency = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])

                if line_split[0] == 'SNP':
                    if [s for s in line_split if 'snp_type=' in s][0].split('=')[1] == 'nonsynonymous':
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
                        if ';' in locus_tag:
                            for locus_tag_j in locus_tag.split(';'):
                                if locus_tag_j not in gene_count_dict_pop:
                                    gene_count_dict_pop[locus_tag_j] = 0
                                gene_count_dict_pop[locus_tag_j] += frequency

                        else:
                            if locus_tag not in gene_count_dict_pop:
                                gene_count_dict_pop[locus_tag] = 0
                            gene_count_dict_pop[locus_tag] += frequency
                            
                    if line_split[0] == 'DEL':
                         #locus_tags_inactivated=BIF_00532,BIF_01857,BIF_00807,BIF_01111,BIF_01110,BIF_00609,BIF_00200,BIF_00733
                        if [s for s in line_split if 'gene_position=' in s][0].split('=')[1][0:6] == 'interg':
                                                    continue
                        elif [s for s in line_split if 'gene_position=' in s][0].split('=')[1][0:6] == 'pseudo':
                            continue
                        else:
                            locus_tags_comma_unseparated = [s for s in line_split if 'locus_tags_inactivated=' in s][0].split('=')[1]
                            frequency = float([s for s in line_split if 'frequency=' in s][0].split('=')[1])
                            if len(locus_tags_comma_unseparated.split(',')) <= 1:
                                continue
                            else:
                                for locus_tag_j in locus_tags_comma_unseparated.split(','):
                                    if locus_tag_j not in gene_count_dict_pop:
                                        gene_count_dict_pop[locus_tag_j] = 0
                                    gene_count_dict_pop[locus_tag_j] += frequency                            
                            


                        

                        

            gene_parallelism_statistics = {}
            for gene_i, length_i in effective_gene_lengths.items():
                gene_parallelism_statistics[gene_i] = {}
                gene_parallelism_statistics[gene_i]['length'] = length_i
                gene_parallelism_statistics[gene_i]['observed'] = 0
                gene_parallelism_statistics[gene_i]['multiplicity'] = 0

            # save number of mutations for multiplicity
            for locus_tag_i, n_muts_i in gene_count_dict_pop.items():
                gene_parallelism_statistics[locus_tag_i]['observed'] = n_muts_i

            # save number of mutations for multiplicity
            L_mean = np.mean(list(effective_gene_lengths.values()))
            L_tot = sum(list(effective_gene_lengths.values()))
            n_tot = sum(gene_count_dict_pop.values())
            # go back over and calculate multiplicity
            for locus_tag_i in gene_parallelism_statistics.keys():
                # double check the measurements from this
                gene_parallelism_statistics[locus_tag_i]['multiplicity'] = gene_parallelism_statistics[locus_tag_i]['observed'] *1.0/ effective_gene_lengths[locus_tag_i] * L_mean
                gene_parallelism_statistics[locus_tag_i]['expected'] = n_tot*gene_parallelism_statistics[locus_tag_i]['length']/L_tot

            # split locus tags
            for locus_tag_i in gene_parallelism_statistics.keys():
                mult_i = gene_parallelism_statistics[locus_tag_i]['multiplicity']
                if mult_i > 0:
                    locus_tag_i_num = locus_tag_i.split('_')[1]
                    gene_by_pop_dict[pop][locus_tag_i_num] = mult_i

    gene_by_pop_df = pd.DataFrame(gene_by_pop_dict)#outputs the GxP matrix
    gene_by_pop_df = gene_by_pop_df.T
    gene_by_pop_df.fillna(0, inplace=True)

    gene_by_pop_df_out = mt.get_path() + '/data/mult_by_pop.txt'
    gene_by_pop_df.to_csv(gene_by_pop_df_out, sep = '\t', index = True)




#get_multiplicity_matrix()
# equal rate of evolution

#get_multiplicity()
