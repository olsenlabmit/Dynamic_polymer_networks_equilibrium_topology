#!/usr/local/bin/env python
# -*- coding: utf-8 -*-

"""
#######################################
#                                     #
#-- Depercolation properties of Networks --#
#------  Author: Devosmita sen --------#
#                                     #
#######################################

 Overall Framework (Steps):
     1. Generate a Network using KMC algorithm
     
     2. Do depercolation of network by breaking bonds
     3. Calculate current, other properties
"""
import os.path
import sys

import time
import math
import random
import matplotlib
import numpy as np

from numpy import linalg as LA
from scipy.optimize import fsolve
from matplotlib import pyplot as plt
import shutil
import ioLAMMPS
import networkx as nx



random.seed(a=None, version=2)
##print('First random number of this seed: %d'%(random.randint(0, 10000))) 

c0=np.array([6.5846,8.1319,11.5651,17.2674,28.6624,67.5382])
c0_str=np.array(['6.5846','8.1319','11.5651','17.2674','28.6624','67.5382'])
cs_list=np.array(['1.1408'])#,'1.4089','2.0037','2.9916','4.9659','11.7012'])
cs_list_num=np.array([1.1408,1.4089,2.0037,2.9916,4.9659,11.7012])


for step in [1,2,3]:
    mean_dgel_second_largest=[]
    mean_dgel_avg_wo_largest=[]
    std_dgel_second_largest=[]
    std_dgel_avg_wo_largest=[]
    avg_cluster_size_wo_largest_all=[]
    for i in cs_list:
        list_dgel_second_largest=[]
        list_dgel_avg_wo_largest=[]
        count=int(np.where(cs_list==str(i))[0])
        
        print('C0',c0[count])
        for run in [1,2,3,4,5,6,7,8,10]:
            print('Run',run)
            NRA=np.loadtxt("./Run"+str(run)+'/'  +'NRA.txt')

            NRA=int(NRA)
            NRB=int(NRA/2)
            
            G=nx.MultiGraph()
            netgen_flag = 0 # read from file- KMC network
            swell = 0
            if(netgen_flag==0):
              
               print('--------------------------')   
               print('----Reading Network-------')   
               print('--------------------------')
               
               filename ="./Run"+str(run)+'/'+str(step)+'_network_KMC.txt'

               G=ioLAMMPS.readLAMMPS_into_graph_from_KMC(G,filename)
               ##print('len(G.nodes())',len(G.nodes),len(G.edges))
            else:
               print('Invalid network generation flag')


            G_orig=G.copy()

            potential_node_1_list=[]

            for i in G:
                  if(len(list(G.neighbors(i)))>0):
                     potential_node_1_list.append(i)


            density=[]
            max_connected=[]
            second_largest_connected=[]
            avg_cluster_size_without_largest=[]
            path_exists_array=[]
            num_paths_array=[]
            fraction_paths_connected=[]
            frac_cleaved_array=np.arange(0.01,0.99,0.01)##[0.4,0.5,0.55,0.56,0.57,0.58,0.59,0.6,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.7,0.8,0.9]
            for frac_cleaved in frac_cleaved_array: 
               G_temp=G.copy()
              
               init_num_bonds=G_temp.number_of_edges()

               num_broken_bonds=0
               while(num_broken_bonds<int(frac_cleaved*init_num_bonds)):
                  for node_1 in potential_node_1_list:
                     if(num_broken_bonds<int(frac_cleaved*init_num_bonds)):
                        a=list(nx.all_neighbors(G_temp, node_1))
                        if(a!=[]):
                           node_2_list=list(G_temp.neighbors(node_1))
                           node_2=random.choice(node_2_list)
                           G_temp.remove_edge(node_1,node_2)

                           num_broken_bonds=num_broken_bonds+1
                 
               density.append(nx.density(G_temp))
               
               max_connected.append(max(len(cc) for cc in nx.connected_components(G_temp)))
               cc=list(nx.connected_components(G_temp))
               
               cc.sort(key=len)
               second_largest_connected.append(len(cc[len(cc)-2]))
               all_cluster_sizes=[len(x) for x in cc[0:-1]] # without largest
               cluster_sizes, number_distribution=np.unique(all_cluster_sizes,return_counts=True)

               numerator=0
               denominator=0
               for i in range(len(cluster_sizes)):
                   numerator=numerator+(cluster_sizes[i]**2)*number_distribution[i]
                   denominator=denominator+(cluster_sizes[i])*number_distribution[i]
                   
               avg_cluster_size_without_largest.append(numerator/denominator)#np.mean([len(x) for x in cc[0:-1]]))
               sys_size=len(G_temp)

            
            i=np.argmax(second_largest_connected)
            list_dgel_second_largest.append(frac_cleaved_array[i])
            #print('De-gel point from second largest connected',frac_cleaved_array[i])
            
            i=np.argmax(avg_cluster_size_without_largest)
            list_dgel_avg_wo_largest.append(frac_cleaved_array[i])
            print('De-gel point from average cluster size without largest',frac_cleaved_array[i])

               
            plt.figure(step)
            plt.plot(frac_cleaved_array,avg_cluster_size_without_largest,'o-',label=str(run))
            avg_cluster_size_wo_largest_all.append(avg_cluster_size_without_largest)
        
        list_dgel_second_largest=np.array(list_dgel_second_largest)
        list_dgel_avg_wo_largest=np.array(list_dgel_avg_wo_largest)

        mean_dgel_second_largest.append(np.mean(list_dgel_second_largest))
        mean_dgel_avg_wo_largest.append(np.mean(list_dgel_avg_wo_largest))
        std_dgel_second_largest.append(np.std(list_dgel_second_largest))
        std_dgel_avg_wo_largest.append(np.std(list_dgel_avg_wo_largest))
        plt.figure(step)
        plt.legend()
        
        np.savetxt(str(step)+'_deperc_data.txt',np.transpose(np.array([frac_cleaved_array,np.mean(avg_cluster_size_wo_largest_all,axis=0),np.std(avg_cluster_size_wo_largest_all,axis=0)])))
        np.savetxt(str(step)+'_deperc_thres.txt',np.array([mean_dgel_avg_wo_largest,std_dgel_avg_wo_largest]))
        
'''
plt.figure()
plt.errorbar(cs_list_num,np.sqrt(mean_dgel_second_largest),yerr=np.multiply(std_dgel_second_largest,0.5),fmt = 'o',color = 'green', 
        ecolor = 'green', elinewidth = 2,capsize=5)

plt.xlabel('Conc dimensionless')
plt.ylabel('gel pt from Size of SECOND largest connected component')
plt.savefig('1_gel pt. - size_second_largest_connected_component.png')

plt.figure()
plt.errorbar(1.0/cs_list_num,np.sqrt(mean_dgel_second_largest),yerr=np.multiply(std_dgel_second_largest,0.5),fmt = 'o',color = 'green', 
        ecolor = 'green', elinewidth = 2,capsize=5)
plt.xlabel('Conc dimensionless')
plt.ylabel('gel pt from Size of SECOND largest connected component')
plt.savefig('1_gel pt. - size_second_largest_connected_component_inv_dimC.png')


plt.figure()
plt.errorbar(cs_list_num,np.sqrt(mean_dgel_avg_wo_largest),yerr=np.multiply(std_dgel_avg_wo_largest,0.5),fmt = 'o',color = 'blue', 
        ecolor = 'blue', elinewidth = 2,capsize=5)
plt.xlabel('Conc dimensionless')
plt.ylabel('gel pt from weighted Average cluster size without largest')
plt.savefig('1_gel pt. - avg_cluster_size_withiut_largest.png')

plt.figure()
plt.errorbar(1.0/(cs_list_num),np.sqrt(mean_dgel_avg_wo_largest),yerr=np.multiply(std_dgel_avg_wo_largest,0.5),fmt = 'o',color = 'blue', 
        ecolor = 'blue', elinewidth = 2,capsize=5)
plt.xlabel('1/cR3')
plt.ylabel('gel pt from weighted Average cluster size without largest')
plt.savefig('1_gel pt. - avg_cluster_size_withiut_largest_inv_dimC.png')
'''
##np.savetxt('1_all_data_percolation.txt',np.transpose(np.array([cs_list_num,(np.sqrt(mean_dgel_avg_wo_largest)),(np.sqrt(mean_dgel_second_largest)),(np.multiply(std_dgel_avg_wo_largest,0.5))])))
plt.savefig('all_deperc_runs')
plt.show()
