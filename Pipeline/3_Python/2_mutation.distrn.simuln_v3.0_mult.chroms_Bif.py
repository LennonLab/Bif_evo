# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 16:35:59 2019

@author: rmoge
"""
#%%
from __future__ import division
import re, os,sys, math, operator,random, copy,collections,time; import numpy as np; import pandas as pd; import csv
from itertools import groupby; import pprint as pp
import matplotlib.pyplot as plt
import statsmodels.stats.multitest as smt

from pathlib import Path
homedir = str(Path.home())
workingdir = homedir + '\\Documents\\GitHub\\Bifidobacterium'
os.chdir(workingdir)
#%%
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

#%%
# get all sequence records for the specified genbank file
recs = [rec for rec in SeqIO.parse(workingdir + "\\Pipeline\\1_Carbonate\\BB12_Jensen_2021.gb", "genbank")]

# print the number of sequence records that were extracted
print(len(recs))

# print annotations for each sequence record
for rec in recs:
	print(rec.annotations)
#%%    
myrec=recs[0]

# print the CDS sequence feature summary information for each feature in each
# sequence record
for rec in recs:
    feats = [feat for feat in rec.features if feat.type == "CDS"]
    for feat in feats:
        print(feat)
        
print(repr(myrec.seq))
#%%
gb_file=workingdir + "\\Pipeline\\1_Carbonate\\BB12_Jensen_2021.gb"
gb_record = SeqIO.read(open(gb_file,"r"), "genbank")
print("Name %s, %i features" % (gb_record.name, len(gb_record.features)))
print(repr(gb_record.seq))


#%%
print(myrec.features[2].location.start)
myrec[(myrec.features[2].location.start)]#this is how you access the 0th nt of the gb file...
#for the parallel stuff, we just need lengths
#%%
def gene_to_tag(gb_record, feature_type, qualifier) :
    answer = dict()
    for (index, feature) in enumerate(gb_record.features) :
        if feature.type=='gene':
            if qualifier in feature.qualifiers :
                #There should only be one locus_tag per feature, but there
                #are usually several db_xref entries
                for value in feature.qualifiers[qualifier] :
                    if value in answer :
                        print("WARNING - Duplicate key")
                        answer[value]['index'].append(index)
                    else :
                        answer[value] = dict()
                        answer[value]['index'] = list()
                        answer[value]['index'].append(index)
            if qualifier in feature.qualifiers:
                for value in feature.qualifiers[qualifier] :
                    answer[value]['locus_tag']=feature.qualifiers['locus_tag']    
    return answer

#test = index_genbank_features(myrec, "CDS", "gene")
GD = gene_to_tag(myrec, "locus_tag", "gene")

#%%
# I think it may make more sense to build the dict all in one loop, rather than piping an updated dict from one function to the next, though. I attempt to do that in this chunk.
#We are gonna build a badass dict of dicts. Each locus tag will be a key. Its value will be a dict. In that dict, each key will be some kind of feature (e.g., length) and its value will be the corresponding value for THAT LOCUS.
def index_genbank_features_multiplechroms(gb_records, feature_type, qualifier, qualifier_translation, AT_rel_rate, GC_rel_rate):
    
    maxperchromlist = list()
    for rec in gb_records:
        max_on_this_chrom = max([len(i) for i in rec.features[1:]])
        maxperchromlist.append(max_on_this_chrom)
    max_gene_length = max(maxperchromlist)
    rel_rate_MAX=max(AT_rel_rate,GC_rel_rate)
    print(str(max_gene_length) + " is max gene length")
    
    answer = dict()

    for rec in gb_records:
        for (index, feature) in enumerate(rec.features) :
            if feature.type==feature_type :
                if qualifier in feature.qualifiers :
                    #There should only be one locus_tag per feature, but there
                    #are usually several db_xref entries
                    for value in feature.qualifiers[qualifier] :
                        if value in answer :
                            print("WARNING - Duplicate key %s for %s features %i and %i" \
                               % (value, feature_type, answer[value], index))
                        else :
                            answer[value] = {}
                            answer[value]["index"] = index
                try:
                    mygg=feature.qualifiers[qualifier_translation]#attempt to access the gene name for this locus tag
                except KeyError:
                    mygg=feature.qualifiers[qualifier]#if there is no gene name, just call it by the locus tag
                if qualifier in feature.qualifiers:
                    for value in feature.qualifiers[qualifier] :
                        answer[value]['gene'] = mygg[0]
                if qualifier in feature.qualifiers:
                    for value in feature.qualifiers[qualifier] :
                        answer[value]['gene_length']=len(feature)
                if qualifier in feature.qualifiers:
                    for value in feature.qualifiers[qualifier] :        
                        answer[value]["gene_len_rel"]= (len(feature) / max_gene_length)
                if qualifier in feature.qualifiers:
                    for value in feature.qualifiers[qualifier]:
                        GCcount=0
                        for (index2,pos) in enumerate(gb_record.features[index]):
                            if gb_record[pos]=="C" or gb_record[pos]=="G":
                                GCcount+=1
                            else:
                                pass
                            #print(pos)
                        answer[value]["gene_GC_prop"]= (GCcount/len(feature))
                if qualifier in feature.qualifiers:
                    for value in feature.qualifiers[qualifier]:
                        answer[value]['mut_rel_rate']=((answer[value]['gene_GC_prop']*GC_rel_rate) + ((1 - answer[value]['gene_GC_prop'])*AT_rel_rate)) / rel_rate_MAX
                    
            
    return answer

#test = index_genbank_features(myrec, "CDS", "gene")
DD = index_genbank_features_multiplechroms(recs, "gene", "locus_tag", 'gene', AT_rel_rate = 1, GC_rel_rate=1)#at_rel_rate is the relative rate of mutation of an A:T nucleotide; gc_rel_rate mutatis mutandis
#%%
# I think it may make more sense to build the dict all in one loop, rather than piping an updated dict from one function to the next, though. I attempt to do that in this chunk.
#We are gonna build a badass dict of dicts. Each locus tag will be a key. Its value will be a dict. In that dict, each key will be some kind of feature (e.g., length) and its value will be the corresponding value for THAT LOCUS.
def index_genbank_features(gb_record, feature_type, qualifier, qualifier_translation, AT_rel_rate, GC_rel_rate):
    max_gene_length = max([len(i) for i in gb_record.features[1:]])
    rel_rate_MAX=max(AT_rel_rate,GC_rel_rate)
    print(str(max_gene_length) + " is max gene length")
    answer = dict()
    for (index, feature) in enumerate(gb_record.features) :
        if feature.type==feature_type :
            if qualifier in feature.qualifiers :
                #There should only be one locus_tag per feature, but there
                #are usually several db_xref entries
                for value in feature.qualifiers[qualifier] :
                    if value in answer :
                        print("WARNING - Duplicate key %s for %s features %i and %i" \
                           % (value, feature_type, answer[value], index))
                    else :
                        answer[value] = {}
                        answer[value]["index"] = index
            try:
                mygg=feature.qualifiers[qualifier_translation]#attempt to access the gene name for this locus tag
            except KeyError:
                mygg=feature.qualifiers[qualifier]#if there is no gene name, just call it by the locus tag
            if qualifier in feature.qualifiers:
                for value in feature.qualifiers[qualifier] :
                    answer[value]['gene'] = mygg[0]
            if qualifier in feature.qualifiers:
                for value in feature.qualifiers[qualifier] :
                    answer[value]['gene_length']=len(feature)
            if qualifier in feature.qualifiers:
                for value in feature.qualifiers[qualifier] :        
                    answer[value]["gene_len_rel"]= (len(feature) / max_gene_length)
            if qualifier in feature.qualifiers:
                for value in feature.qualifiers[qualifier]:
                    GCcount=0
                    for (index2,pos) in enumerate(gb_record.features[index]):
                        if gb_record[pos]=="C" or gb_record[pos]=="G":
                            GCcount+=1
                        else:
                            pass
                        #print(pos)
                    answer[value]["gene_GC_prop"]= (GCcount/len(feature))
            if qualifier in feature.qualifiers:
                for value in feature.qualifiers[qualifier]:
                    answer[value]['mut_rel_rate']=((answer[value]['gene_GC_prop']*GC_rel_rate) + ((1 - answer[value]['gene_GC_prop'])*AT_rel_rate)) / rel_rate_MAX
                    
            
    return answer

#test = index_genbank_features(myrec, "CDS", "gene")
DD2 = index_genbank_features(myrec, "gene", "locus_tag", 'gene', AT_rel_rate = 1, GC_rel_rate=1)#at_rel_rate is the relative rate of mutation of an A:T nucleotide; gc_rel_rate mutatis mutandis
#%%
#Import mutation frequency data
def freq_list(infile=""):
    mydf = pd.read_csv(infile)
    return mydf          
            




F_freqlist_nonsyn = freq_list(infile=workingdir+"\\Pipeline\\3_Python\\mutation.frequencies_nonsynonymous_F.csv")
F_freqs_nonsyn = tuple(F_freqlist_nonsyn.loc[:, 'freq'].values)

P_freqlist_nonsyn = freq_list(infile=workingdir+"\\Pipeline\\3_Python\\mutation.frequencies_nonsynonymous_P.csv")
P_freqs_nonsyn = tuple(P_freqlist_nonsyn.loc[:, 'freq'].values)

X_freqlist_nonsyn = freq_list(infile=workingdir+"\\Pipeline\\3_Python\\mutation.frequencies_nonsynonymous_X.csv")
X_freqs_nonsyn = tuple(X_freqlist_nonsyn.loc[:, 'freq'].values)


F_freqlist_syn = freq_list(infile=workingdir+"\\Pipeline\\3_Python\\mutation.frequencies_synonymous_F.csv")
F_freqs_syn = tuple(F_freqlist_syn.loc[:, 'freq'].values)

P_freqlist_syn = freq_list(infile=workingdir+"\\Pipeline\\3_Python\\mutation.frequencies_synonymous_P.csv")
P_freqs_syn = tuple(P_freqlist_syn.loc[:, 'freq'].values)

X_freqlist_syn = freq_list(infile=workingdir+"\\Pipeline\\3_Python\\mutation.frequencies_synonymous_X.csv")
X_freqs_syn = tuple(X_freqlist_syn.loc[:, 'freq'].values)

#%%
roogeld = dict()
roogeld["BIF_00360"]= DD["BIF_00360"]; roogeld["BIF_01452"] = DD["BIF_01452"]; roogeld["BIF_01454"] = DD["BIF_01454"]

def simulate_mutations_GC(mutation_freq_tup, indict, reps=10000):
    counter = 0; sim_master=dict(); sim_dict=dict()
    for k in indict.keys():
        sim_master[k] = list()
    while counter < reps:
        sim_temp=dict()
        for index, m in enumerate(mutation_freq_tup):
            while True:
                w = random.choice(list(indict));r = random.random();rGC=random.random()
                #print("w = " + str(w) + ", w len rel = " + str(indict[w]['gene_len_rel']) + " and r = " + str(r))
                if indict[w]['gene_len_rel'] >= r and indict[w]['mut_rel_rate'] >= rGC:
                    #do something###########################################
                    ##################################
                    #################################################
                    try:
                        sim_temp[w] += m
                    except KeyError:
                        sim_temp[w] = 0
                        sim_temp[w] += m
                    #print("value added = " + str(m))
                    break
                else:
                    continue
        sim_dict[counter] = sim_temp
        for k,v in sim_temp.items():
            sim_master[k].append(v)
        for km in sim_master.keys():
            if len(sim_master[km]) < counter+1:
                sim_master[km].append(0)
        counter += 1
        #pp.pprint(sim_temp)
    
    return sim_master

def simulate_mutations_noGC(mutation_freq_tup, indict, reps=10000):
    counter = 0; sim_master=dict(); sim_dict=dict()
    for k in indict.keys():
        sim_master[k] = list()
    while counter < reps:
        sim_temp=dict()
        for index, m in enumerate(mutation_freq_tup):
            while True:
                w = random.choice(list(indict));r = random.random()
                #print("w = " + str(w) + ", w len rel = " + str(indict[w]['gene_len_rel']) + " and r = " + str(r))
                if indict[w]['gene_len_rel'] >= r:
                    #do something###########################################
                    ##################################
                    #################################################
                    try:
                        sim_temp[w] += m
                    except KeyError:
                        sim_temp[w] = 0
                        sim_temp[w] += m
                    #print("value added = " + str(m))
                    break
                else:
                    continue
        sim_dict[counter] = sim_temp
        for k,v in sim_temp.items():
            sim_master[k].append(v)
        for km in sim_master.keys():
            if len(sim_master[km]) < counter+1:
                sim_master[km].append(0)
        counter += 1
        #pp.pprint(sim_temp)
    
    return sim_master
#%%
start = time.time()
Fsim = simulate_mutations_noGC((F_freqs_nonsyn),DD,reps = 100000)
end = time.time()
print(end - start)


#1 000 reps took 9 seconds. So I could do 100 000 reps in 15 minutes.

Fsimpd=pd.DataFrame.from_dict(Fsim)
Fsimpd.to_csv(index=True,path_or_buf=workingdir+"\\Pipeline\\3_Python\\simout_F.nonsyn_100000.csv")
#%%
start = time.time()
Psim = simulate_mutations_noGC((P_freqs_nonsyn),DD,reps = 100000)
end = time.time()
print(end - start)

Psimpd=pd.DataFrame.from_dict(Psim)
Psimpd.to_csv(index=True,path_or_buf=workingdir+"\\Pipeline\\3_Python\\simout_P.nonsyn_100000.csv")
#%%
start = time.time()
Xsim = simulate_mutations_noGC((X_freqs_nonsyn),DD,reps = 100000)
end = time.time()
print(end - start)

Xsimpd=pd.DataFrame.from_dict(Xsim)
Xsimpd.to_csv(index=True,path_or_buf=workingdir+"\\Pipeline\\3_Python\\simout_X.nonsyn_100000.csv")

#%%
start = time.time()
Fsim_syn = simulate_mutations_noGC((F_freqs_syn),DD,reps = 100000)
end = time.time()
print(end - start)


#1 000 reps took 9 seconds. So I could do 100 000 reps in 15 minutes.

Fsimpd_syn=pd.DataFrame.from_dict(Fsim_syn)
Fsimpd_syn.to_csv(index=True,path_or_buf=workingdir+"\\Pipeline\\3_Python\\simout_F.syn_100000.csv")
#%%
start = time.time()
Psim_syn = simulate_mutations_noGC((P_freqs_syn),DD,reps = 100000)
end = time.time()
print(end - start)

Psimpd_syn=pd.DataFrame.from_dict(Psim_syn)
Psimpd_syn.to_csv(index=True,path_or_buf=workingdir+"\\Pipeline\\3_Python\\simout_P.syn_100000.csv")
#%%
start = time.time()
Xsim_syn = simulate_mutations_noGC((X_freqs_syn),DD,reps = 100000)
end = time.time()
print(end - start)

Xsimpd_syn=pd.DataFrame.from_dict(Xsim_syn)
Xsimpd_syn.to_csv(index=True,path_or_buf=workingdir+"\\Pipeline\\3_Python\\simout_X.syn_100000.csv")

#%%

#You can import simulation results from a file.
start = time.time();print(start)
recs = [rec for rec in SeqIO.parse(workingdir + "\\Pipeline\\1_Carbonate\\BB12_Jensen_2021.gb", "genbank")]
myrec=recs[0]

GD = gene_to_tag(myrec, "locus_tag", "gene")
DD = index_genbank_features_multiplechroms(recs, "gene", "locus_tag", 'gene', AT_rel_rate = 1, GC_rel_rate=1)
#import statsmodels.stats.multitest as smt

#smopd = pd.read_csv(r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\simulations_syn3B_100000.csv")
#smo1pd = pd.read_csv(r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\simulations_syn1.0_100000.csv")
start = time.time();print(start)
Fsimpd = pd.read_csv(workingdir+"\\Pipeline\\3_Python\\simout_F.nonsyn_100000.csv")
Fsim = pd.DataFrame.to_dict(Fsimpd)
end = time.time()
print(end - start)
#2405 seconds for 100 000 reps

start = time.time();print(start)
Psimpd = pd.read_csv(workingdir+"\\Pipeline\\3_Python\\simout_P.nonsyn_100000.csv")
Psim = pd.DataFrame.to_dict(Psimpd)
end = time.time()
print(end - start)

start = time.time();print(start)
Xsimpd = pd.read_csv(workingdir+"\\Pipeline\\3_Python\\simout_X.nonsyn_100000.csv")
Xsim = pd.DataFrame.to_dict(Xsimpd)
end = time.time()
print(end - start)

start = time.time();print(start)
Fsimpd_syn = pd.read_csv(workingdir+"\\Pipeline\\3_Python\\simout_F.syn_100000.csv")
Fsim_syn = pd.DataFrame.to_dict(Fsimpd_syn)
end = time.time()
print(end - start)
#2405 seconds for 100 000 reps

start = time.time();print(start)
Psimpd_syn = pd.read_csv(workingdir+"\\Pipeline\\3_Python\\simout_P.syn_100000.csv")
Psim_syn = pd.DataFrame.to_dict(Psimpd_syn)
end = time.time()
print(end - start)

start = time.time();print(start)
Xsimpd_syn = pd.read_csv(workingdir+"\\Pipeline\\3_Python\\simout_X.syn_100000.csv")
Xsim_syn = pd.DataFrame.to_dict(Xsimpd_syn)
end = time.time()
print(end - start)

#%%
#Now do the analyses on the imported simulation results


#%%
def adjust_BH(roypvals=dict()):
    pvals_list=list()
    pgenes_list=list()
    for k,v in roypvals.items():
        pvals_list.append(v)
        pgenes_list.append(k)
        
    fdr_arrays=smt.multipletests(pvals=pvals_list,alpha=0.05,method='fdr_bh',is_sorted=False,returnsorted=False)
    
    fdr=dict()
    for index,i in enumerate(fdr_arrays[1]):
        fdr[pgenes_list[index]]=i
        
    return fdr
#%%
F_pvals = dict();reps=100000


locustagstotest_F=['BIF_00489','BIF_00490','BIF_00532','BIF_00625','BIF_00651','BIF_00684','BIF_00778','BIF_00892','BIF_00936','BIF_01346','BIF_01467','BIF_01492','BIF_01639','BIF_01789','BIF_02305','BIF_01197']
locustagssums_F=[6.18,2.0,3.731,3.778,1.671,5.777,2.0,.462,1.657,1.140,1.44,3.923,2,2,1.474,3]

for index,j in enumerate(locustagstotest_F):
    #print(index)
    #print(j)
    F_pvals[j]=max(1/reps,np.count_nonzero([i>=locustagssums_F[index] for i in Fsim[j]])/reps)
    
print(F_pvals)
F_fdr=adjust_BH(roypvals=F_pvals)
print(F_fdr)

F_fdr_out=dict()
for k,v in F_fdr.items():
    F_fdr_out[k]=list()
    F_fdr_out[k].append(v)

F_fdrpd=pd.DataFrame.from_dict(F_fdr_out)
F_fdrpd.to_csv(index=True,path_or_buf=workingdir+"\\Pipeline\\3_Python\\fdrout_F.nonsyn_100000.csv")
#%%
P_pvals = dict();reps=100000


locustagstotest_P=['BIF_00489','BIF_00490','BIF_00651','BIF_00778','BIF_00936','BIF_00944','BIF_01327','BIF_01467','BIF_01492','BIF_01616','BIF_01639','BIF_01789','BIF_01831','BIF_01060']
locustagssums_P=[7.564,4.0,5.0,2.0,6.0,2.0,1.332,2.0,3.568,1.765,3.0,4.0,4.0,2.0]

for index,j in enumerate(locustagstotest_P):
    #print(index)
    #print(j)
    P_pvals[j]=max(1/reps,np.count_nonzero([i>=locustagssums_P[index] for i in Psim[j]])/reps)
    
print(P_pvals)
P_fdr=adjust_BH(roypvals=P_pvals)
print(P_fdr)

P_fdr_out=dict()
for k,v in P_fdr.items():
    P_fdr_out[k]=list()
    P_fdr_out[k].append(v)

P_fdrpd=pd.DataFrame.from_dict(P_fdr_out)
P_fdrpd.to_csv(index=True,path_or_buf=workingdir+"\\Pipeline\\3_Python\\fdrout_P.nonsyn_100000.csv")

#%%
X_pvals = dict();reps=100000


locustagstotest_X=['BIF_00208','BIF_00327','BIF_00489','BIF_00490','BIF_00532','BIF_00573','BIF_00776','BIF_00936','BIF_00944','BIF_01193','BIF_01492','BIF_01616','BIF_01789','BIF_02090','BIF_01402','BIF_01023']
locustagssums_X=[5,2.745,13.125,1.722,3.004,2.0,2.001,.895,2.575,2.643,1.056,4.010,.698,1.508,2.013,3.641,1.248]

for index,j in enumerate(locustagstotest_X):
    #print(index)
    #print(j)
    X_pvals[j]=max(1/reps,np.count_nonzero([i>=locustagssums_X[index] for i in Xsim[j]])/reps)
    
print(X_pvals)
X_fdr=adjust_BH(roypvals=X_pvals)
print(X_fdr)

X_fdr_out=dict()
for k,v in X_fdr.items():
    X_fdr_out[k]=list()
    X_fdr_out[k].append(v)

X_fdrpd=pd.DataFrame.from_dict(X_fdr_out)
X_fdrpd.to_csv(index=True,path_or_buf=workingdir+"\\Pipeline\\3_Python\\fdrout_X.nonsyn_100000.csv")


#%%
Fsyn_pvals = dict();reps=100000


locustagstotest_Fsyn=[]
locustagssums_Fsyn=[]

for index,j in enumerate(locustagstotest_Fsyn):
    #print(index)
    #print(j)
    Fsyn_pvals[j]=max(1/reps,np.count_nonzero([i>=locustagssums_Fsyn[index] for i in Fsim_syn[j]])/reps)
    
print(Fsyn_pvals)
Fsyn_fdr=adjust_BH(roypvals=Fsyn_pvals)
print(Fsyn_fdr)

Fsyn_fdr_out=dict()
for k,v in Fsyn_fdr.items():
    Fsyn_fdr_out[k]=list()
    Fsyn_fdr_out[k].append(v)

Fsyn_fdrpd=pd.DataFrame.from_dict(Fsyn_fdr_out)
Fsyn_fdrpd.to_csv(index=True,path_or_buf=workingdir+"\\Pipeline\\3_Python\\fdrout_F.syn_100000.csv")
#%%
Psyn_pvals = dict();reps=100000


locustagstotest_Psyn=['BIF_00469']
locustagssums_Psyn=[2.0]

for index,j in enumerate(locustagstotest_Psyn):
    #print(index)
    #print(j)
    Psyn_pvals[j]=max(1/reps,np.count_nonzero([i>=locustagssums_Psyn[index] for i in Psim_syn[j]])/reps)
    
print(Psyn_pvals)
Psyn_fdr=adjust_BH(roypvals=Psyn_pvals)
print(Psyn_fdr)

Psyn_fdr_out=dict()
for k,v in Psyn_fdr.items():
    Psyn_fdr_out[k]=list()
    Psyn_fdr_out[k].append(v)

Psyn_fdrpd=pd.DataFrame.from_dict(Psyn_fdr_out)
Psyn_fdrpd.to_csv(index=True,path_or_buf=workingdir+"\\Pipeline\\3_Python\\fdrout_P.syn_100000.csv")

#%%
Xsyn_pvals = dict();reps=100000


locustagstotest_Xsyn=['BIF_02970']
locustagssums_Xsyn=[1.047]

for index,j in enumerate(locustagstotest_Xsyn):
    #print(index)
    #print(j)
    Xsyn_pvals[j]=max(1/reps,np.count_nonzero([i>=locustagssums_Xsyn[index] for i in Xsim_syn[j]])/reps)
    
print(Xsyn_pvals)
Xsyn_fdr=adjust_BH(roypvals=Xsyn_pvals)
print(Xsyn_fdr)

Xsyn_fdr_out=dict()
for k,v in Xsyn_fdr.items():
    Xsyn_fdr_out[k]=list()
    Xsyn_fdr_out[k].append(v)

Xsyn_fdrpd=pd.DataFrame.from_dict(Xsyn_fdr_out)
Xsyn_fdrpd.to_csv(index=True,path_or_buf=workingdir+"\\Pipeline\\3_Python\\fdrout_X.syn_100000.csv")
#use similar code to output the nonsynonymous tables
# =============================================================================
# #%%
# def compare_recs(gb_record1,gb_record2):
#     answer = dict()
#     bads=list()
#     count=0
#     for (index, feature) in enumerate(gb_record1.features) :
#         if feature.type=='CDS':
#             if "locus_tag" in feature.qualifiers :
#                 lt1 = feature.qualifiers['locus_tag']
#                 num1 = lt1[0][-5:]
#                 #print(num1)
#                 try:
#                     psq1 = feature.qualifiers['translation']
#                 except KeyError:
#                     pass
#                 #try:
#                  #   print(psq)
#                 #except UnboundLocalError:
#                     #print("")
#                 for (index2,feature2) in enumerate(gb_record2.features):
#                     if feature2.type=='CDS':
#                         if "locus_tag" in feature2.qualifiers:
#                             lt2 = feature2.qualifiers['locus_tag']
#                             num2 = lt2[0][-5:]
#                             if num1 == num2:
#                                 #print('match!')
#                                 count+=1
#                                 try:
#                                     psq2 = feature2.qualifiers['translation']
#                                     if psq1 != psq2:
#                                         print('not the same')
#                                         bads+=lt1
#                                     else:
#                                         print('same')
#                                 except:
#                                     pass
#                                        
#                                        
#     return bads
# 
# #test = index_genbank_features(myrec, "CDS", "gene")
# #cd=compare_recs(myrec, syn1rec)
# #print(cd)
# 
# #%%
# #A next step could be to add the gene LENGTH to the dict. Later when youre worrying about dNdS, you would add the gene (DNA) SEQUENCE to the dict.OR MAYBE THE gff3 file would have an easier way to access teh DNA sdequence; worth looking at.
#             
# def looptest(dicter, reps=2):
#     counter = 0
#     while counter < reps:
#         while True:
#             w = random.choice(list(dicter)); r = random.random()
#             print("w = " + str(w) + ", w len rel = " + str(dicter[w]['gene_len_rel']) + " and r = " + str(r))
#             if dicter[w]['gene_len_rel'] >= r:
#                 print("sucksess")
#                 break
#             else:
#                 continue
#         counter += 1    
#     
#     
# #looptest(DD)
#     
#     
#     
# #%%
# def output_genes(gb_record, outfile=""):
#     tempdf=pd.DataFrame()
#     counter=0
#     for (index,feature) in enumerate(gb_record.features):
#         if feature.type == 'gene':
#             if 'gene' in feature.qualifiers:
#                 for v in feature.qualifiers['gene']:
#                     tempdf[str(feature.qualifiers['gene'][0])] = ""
#             else:#~147 genes do not have assigned gene names and are known ONLY by the locus tag.
#                 for v in feature.qualifiers['locus_tag']:
#                     tempdf[str(feature.qualifiers['locus_tag'][0])] = ""
#             counter+=1
#     print(counter)            
# 
#     tempdf.to_csv(index=True,path_or_buf=outfile)    
#     return tempdf
# 
# gxp = output_genes(myrec,r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\gxp_headers_syn3B.csv")
# #pp.pprint(gxp)
# ##
# =============================================================================
