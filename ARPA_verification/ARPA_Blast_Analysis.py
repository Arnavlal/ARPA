#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###################################################################
#ARPA: Alignment-Free Residue Pangenome Analysis
#Arnav Lal, Ahmed M. Moustafa and Paul J. Planet

#University of Pennsylvania, Philadelphia, PA, 19104
#Perelman School of Medicine, University of Pennsylvania, Philadelphia, PA 19104, USA
#Children's Hospital of Philadelphia, Philadelphia, PA, 19104
#Sackler Institute for Comparative Genomics, American Museum of Natural History, New York, NY 10024, USA

#Contact1: arnavlal@sas.upenn.edu
#Contact2: arnavlal@gmail.com
###################################################################

import matplotlib.pyplot as plt
import numpy as np
import pickle


#set your CD to the location of the BLAST analysis files.

per_id=pickle.load(open('per_id.pickle','rb'))
outer_per_id  = pickle.load(open('outer_per_id.pickle','rb'))

cov=pickle.load(open('cov.pickle','rb'))
outer_cov  = pickle.load(open('outer_cov.pickle','rb'))


per_id_sqz = [item for sublist in per_id for item in sublist]
for i in np.arange(len(outer_per_id)):
    if type(outer_per_id[i])==np.float64:
        outer_per_id[i] = list(np.array([outer_per_id[i]]))
outer_per_id_sqz = [item for sublist in outer_per_id for item in sublist]



cov_sqz = [item for sublist in cov for item in sublist]
for i in np.arange(len(outer_cov)):
    if type(outer_cov[i])==np.float64:
        outer_cov[i] = list(np.array([outer_cov[i]]))
outer_cov_sqz = [item for sublist in outer_cov for item in sublist]

#%%FIGURE 1E; Percent Identity

plt.figure(0)
plt.plot(np.arange(len(outer_per_id_sqz))/len(outer_per_id_sqz), outer_per_id_sqz, 'o', color = 'goldenrod',markersize=2)
plt.plot(np.arange(len(per_id_sqz))/len(per_id_sqz), per_id_sqz, 'o', color ='cornflowerblue',markersize=4)
plt.ylim(0,105)

plt.figure(1)
plt.hist(per_id_sqz, bins = np.linspace(0,100,101), color = 'cornflowerblue', orientation='horizontal')
plt.hist(outer_per_id_sqz, bins = np.linspace(0,100,101), color = 'goldenrod', orientation='horizontal')


#%%FIGURE 1G; Percent Identity vs. Query Coverage
plt.figure(3)
plt.plot(outer_per_id_sqz, outer_cov_sqz, 'o', color ='goldenrod', markersize=2)
plt.plot(per_id_sqz, cov_sqz, 'o', color ='cornflowerblue',markersize=4)
plt.xlabel('Percent Identity')
plt.ylabel('Query Coverage')

#%%Developing Query Coverage times Percent Identity Metric

outer_id = np.array(outer_per_id_sqz)*np.array(outer_cov_sqz)/100
inner_id = np.array(per_id_sqz)*np.array(cov_sqz)/100

#%%FIGURE 1F
plt.figure(4)
plt.plot(np.arange(len(outer_cov_sqz))/len(outer_cov_sqz), outer_id, 'o', color = 'goldenrod',markersize=2)
plt.plot(np.arange(len(cov_sqz))/len(cov_sqz), inner_id, 'o', color ='cornflowerblue',markersize=4)
#plt.ylabel('Query_cov times Percent_ID')

plt.figure(5)
plt.hist(inner_id, bins = np.linspace(0,100,101), color = 'cornflowerblue', orientation='horizontal')
plt.hist(outer_id, bins = np.linspace(0,100,101), color = 'goldenrod', orientation='horizontal')

