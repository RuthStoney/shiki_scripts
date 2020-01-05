
import os
import os.path
from pathlib import Path
import pandas as pd
from collections import Counter
import csv
import seaborn as sns; 
sns.set(color_codes=True)
import scipy.spatial as sp
import scipy.cluster.hierarchy as hc
import re
# =============================================================================
# feed in the files
# =============================================================================

 # feed in the files 'SILICO', 'UNIMAN'
# test git

group='EPFL'
paths = [x[0] for x in os.walk('results-master/'+group+'/')]

Epfl = pd.DataFrame()
for path in paths:
#        find the pathway files
    if Path(path+'/PATHWAYS.tsv').is_file():
    
        dfl=pd.read_csv(path+"/PATHWAYS.tsv", sep='\t')
        dfl.insert(0, "file", path.split("/")[-1])
        Epfl=pd.concat([Epfl, dfl])


group='UNIMAN'
paths = [x[0] for x in os.walk('results-master/'+group+'/')]

Uniman = pd.DataFrame()
for path in paths:
#        find the pathway files
    if Path(path+'/PATHWAYS.tsv').is_file():
    
        dfl=pd.read_csv(path+"/PATHWAYS.tsv", sep='\t')
        dfl.insert(0, "file", path.split("/")[-1])
        Uniman=pd.concat([Uniman, dfl])
     
        
group='SILICO'
paths = [x[0] for x in os.walk('results-master/'+group+'/')]

Silico = pd.DataFrame()
for path in paths:
#        find the pathway files
    if Path(path+'/pathways.csv').is_file():
    
        dfl=pd.read_csv(path+"/pathways.csv", sep='\t')
        dfl.insert(0, "file", path.split("/")[-1])
        Silico=pd.concat([Silico, dfl])

# =============================================================================
# bar charts
# =============================================================================
import matplotlib.pyplot as plt
import numpy as np

# count how many times each pathway is present      
count_man=Counter(Uniman['file'])        
count_epfl=Counter(Epfl['file'])
count_silico=Counter(Silico['file'])

# uniman counts
y_pos = np.arange(len(count_man))
plt.bar(y_pos, count_man.values())
plt.xticks(y_pos, count_man.keys(), rotation='vertical')


# for pathways mapped by multiple locations map the numbers
Epfl_Man = set(count_man.keys()).intersection(set(count_epfl.keys()))
Silico_Man = set(count_man.keys()).intersection(set(count_silico.keys()))
Silico_Epfl = set(count_epfl.keys()).intersection(set(count_silico.keys()))
all_comp = Epfl_Man|Silico_Man #|Silico_Epfl

def plot_multibar(count_man, count_epfl, count_silico):
    

    
    # comparison bar chart
    n_groups = len(all_comp)
    vals_man = [count_man[k] for k in all_comp]
    vals_epfl = [count_epfl[k] for k in all_comp]
    vals_silico = [count_silico[k] for k in all_comp]
    
    # create plot
    fig, ax = plt.subplots()
    index = np.arange(n_groups)
    bar_width = 0.3
    opacity = 0.8
    
    rects1 = plt.bar(index, vals_man, bar_width,
    alpha=opacity,
    color='b',
    label='UNIMAN')
    
    rects2 = plt.bar(index + bar_width, vals_epfl, bar_width,
    alpha=opacity,
    color='g',
    label='EPFL')
    
    rects3 = plt.bar(index + bar_width*2, vals_silico, bar_width,
    alpha=opacity,
    color='r',
    label='SILICO')
    
    plt.xlabel('Compound')
    plt.ylabel('Number of pathways')
    plt.title('Common compound pathways UNMAN/EPFL')
    plt.xticks(index + bar_width, all_comp, rotation='vertical')
    plt.yscale('log')
    plt.ylabel('log')
    plt.legend()
    
    plt.tight_layout()
    plt.show()


plot_multibar(count_man, count_epfl, count_silico)
# =============================================================================
# look at pathway contents
# =============================================================================
data=pd.DataFrame(columns=["MAN", "SILICO_before", "SILICO_after", "EPFL_before", "EPFL_after"], index=sorted(all_comp))




def sorted_nicely( l ):
    """ Sorts the given iterable in the way that is expected.
 
    Required arguments:
    l -- The iterable to be sorted.
 
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key = alphanum_key)

def jac(list1, list2):
    return 1-(len(set(list1).intersection(set(list2)))/len(set(list1).union(set(list2)))) 

count_epfl2= count_epfl.copy()
for target in sorted(all_comp):
    
    data.loc[target, "MAN"] = count_man[target]
    if Epfl[Epfl['file']==target].shape[0]>0:
        data.loc[target, "EPFL_before"] = Epfl[Epfl['file']==target].shape[0]
    
    # process the target 
    if Epfl[Epfl['file']==target].shape[0]==0:
        print (target +" not in epfl")
    elif Epfl[Epfl['file']==target].shape[0]==1:
        print (target +" appears once in epfl")
    else:
        print('\n' + target)
        comp = Epfl.loc[Epfl.file==target].copy() 
        comp['P_PR_INTERMEDIATES']=[list(filter(None, x.split(';'))) for x in comp['P_PR_INTERMEDIATES']]

        print("epfl total number of pathways in "+str(len(comp)))
        unique_intermed = list(set(tuple(sorted(sub)) for sub in comp['P_PR_INTERMEDIATES']))
        print("length once all identical intermediate lists removed "+ str(len(unique_intermed)) )
        data.loc[target, "EPFL_after"] = len(unique_intermed)
        print("number in manchester  "+ str(count_man[target] ))
        
        count_epfl2[target]=len(unique_intermed)
        
        if len(unique_intermed)>2:
            # =============================================================================
            # clustering
            # =============================================================================
            # CALCULATE JACCARD DISTANCES
            names=[str(x) for x in list(range(1,len(unique_intermed)+1))]
            # make the dict of pathways
            unique_intermed_dict=dict(zip(names, unique_intermed))
            with open('paths_epfl_' + target + '.csv', 'w') as f:
                for key in unique_intermed_dict.keys():
                    f.write("%s,%s\n"%(key,unique_intermed_dict[key]))
            
            
            dist2 = pd.DataFrame(index=names, columns=names)
            for k1 in list(unique_intermed_dict):
                for k2 in list(unique_intermed_dict):
                    dist2.loc[k1,k2]=jac(unique_intermed_dict[k1],unique_intermed_dict[k2])
                        
            # heatmap
            DF_dism = dist2.astype(float)
            linkage = hc.linkage(sp.distance.squareform(DF_dism), method='average')
            g=sns.clustermap(DF_dism, row_linkage=linkage, col_linkage=linkage, yticklabels=1, figsize=(5,DF_dism.shape[0]/5))
            g.savefig("EPFL_" + target+"_heatmap long.png")
            
            g=sns.clustermap(DF_dism, row_linkage=linkage, col_linkage=linkage, yticklabels=1)
            g.savefig("EPFL_" + target+"_heatmap.png")

        print("\n")


count_silico2= count_silico.copy()

for target in all_comp:
    data.loc[target, "SILICO_before"] = Silico[Silico['file']==target].shape[0]
    # process the target 
    if Silico[Silico['file']==target].shape[0]==0:
        print (target +" not in Silico")
    elif Silico[Silico['file']==target].shape[0]==1:
        print (target +" appears once in Silico")
    else:
        print('\n' + target + ' in silico')
        comp = Silico.loc[Silico.file==target].copy() 
        # check for nans

        comp['P_PR_INTERMEDIATES']=[list(filter(None, str(x).split(' | '))) if type(x)==str else 'NA'  for x in comp['P_PR_INTERMEDIATES']]

        print("silico total number of pathways is "+str(len(comp)))
        unique_intermed = list(set(tuple(sorted(sub)) for sub in comp['P_PR_INTERMEDIATES']))
        print("length once all identical intermediate lists removed "+ str(len(unique_intermed)) )
        data.loc[target, "SILICO_after"] = len(unique_intermed)
        print("number in manchester  "+ str(count_man[target] ))
        
        count_silico2[target]=len(unique_intermed)
        
        if len(unique_intermed)>2:
    
            # =============================================================================
            # clustering
            # =============================================================================
            # CALCULATE JACCARD DISTANCES
            names=[str(x) for x in list(range(1,len(unique_intermed)+1))]
            unique_intermed_dict=dict(zip(names, unique_intermed))
            with open('paths_epfl_' + target + '.csv', 'w') as f:
                for key in unique_intermed_dict.keys():
                    f.write("%s,%s\n"%(key,unique_intermed_dict[key]))
            
            dist2 = pd.DataFrame(index=names, columns=names)
            for k1 in list(unique_intermed_dict):
                for k2 in list(unique_intermed_dict):
                    dist2.loc[k1,k2]=jac(unique_intermed_dict[k1],unique_intermed_dict[k2])
                        
            # heatmap
            DF_dism = dist2.astype(float)
            linkage = hc.linkage(sp.distance.squareform(DF_dism), method='average')
            g=sns.clustermap(DF_dism, row_linkage=linkage, col_linkage=linkage)
            g.savefig("SILICO_" + target+"_heatmap.png")

       
        print("\n")



plot_multibar(count_man, count_epfl2, count_silico2)

data.to_csv("paths_SILICO_EPFL.csv")

# =============================================================================
# # dendrogram
# # http://datanongrata.com/2019/04/27/67/
# from scipy.cluster.hierarchy import linkage, dendrogram
# import scipy.spatial.distance as ssd
# 
# distArray = ssd.squareform(dist)
# data_link = linkage(distArray,  method='complete')
# B=dendrogram(data_link, truncate_mode="lastp",get_leaves=True, count_sort='ascending', show_contracted=True)
# =============================================================================