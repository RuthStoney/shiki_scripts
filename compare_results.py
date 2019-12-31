
import os
import os.path
from pathlib import Path
import pandas as pd
from collections import Counter
import math

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

# =============================================================================
# look at pathway contents
# =============================================================================

for target in all_comp:
    # process the target 
    if Epfl[Epfl['file']==target].shape[0]==0:
        print (target +" not in epfl")
    elif Epfl[Epfl['file']==target].shape[0]==1:
        print (target +" appears once in epfl")
    else:
        print('\n' + target)
        comp = Epfl.loc[Epfl.file==target].copy() 
        comp['P_PR_INTERMEDIATES']=[list(filter(None, x.split(';'))) for x in comp['P_PR_INTERMEDIATES']]

        print("total number of pathways in "+str(len(comp)))
        
        unique_intermed = list(set(tuple(sorted(sub)) for sub in comp['P_PR_INTERMEDIATES']))
        print("length once all identical intermediate lists removed "+ str(len(unique_intermed)) )
        
        if len(unique_intermed)>2:
    
            # =============================================================================
            # clustering
            # =============================================================================
            # CALCULATE JACCARD DISTANCES
            def jac(list1, list2):
                return 1-(len(set(list1).intersection(set(list2)))/len(set(list1).union(set(list2))))    
            # make a list
            l=[ jac(i1, i2) for i1 in unique_intermed for i2 in unique_intermed ]
            # turn it into a array
            dist=np.asarray(l).reshape(len(unique_intermed), len(unique_intermed))
                    
            # heatmap
            import seaborn as sns; sns.set(color_codes=True)
            g = sns.clustermap(dist).fig.suptitle(target)
        
        print("\n")



for target in all_comp:
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

        print("total number of pathways is "+str(len(comp)))
        
        unique_intermed = list(set(tuple(sorted(sub)) for sub in comp['P_PR_INTERMEDIATES']))
        print("length once all identical intermediate lists removed "+ str(len(unique_intermed)) )
        
        if len(unique_intermed)>2:
    
            # =============================================================================
            # clustering
            # =============================================================================
            # CALCULATE JACCARD DISTANCES
            def jac(list1, list2):
                return 1-(len(set(list1).intersection(set(list2)))/len(set(list1).union(set(list2))))    
            # make a list
            l=[ jac(i1, i2) for i1 in unique_intermed for i2 in unique_intermed ]
            # turn it into a array
            dist=np.asarray(l).reshape(len(unique_intermed), len(unique_intermed))
                    
            # heatmap
            import seaborn as sns; sns.set(color_codes=True)
            g = sns.clustermap(dist).fig.suptitle(target)
        
        print("\n")







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