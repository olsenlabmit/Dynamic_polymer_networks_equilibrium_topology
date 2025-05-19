#!/use/local/bin/env python
# -*- coding: utf-8 -*-

# Python script to write configuration for connectivity

import math
import numpy as np
import os.path

def readLAMMPS_into_graph_from_KMC(G,filename):#, vflag,frac_weak): # reads the file and also makes a networkX graph out of it
   chains= np.genfromtxt(filename,delimiter=',',skip_header=1)
   print('len(chains)',len(chains))
   secondary_loop=0
   for i in range(0,len(chains)):
      link_1=chains[i,0]
      link_2=chains[i,1]
      if(G.has_edge(link_1,link_2)):
         secondary_loop=secondary_loop+1
      G.add_edge(link_1,link_2) 
      
   return G



