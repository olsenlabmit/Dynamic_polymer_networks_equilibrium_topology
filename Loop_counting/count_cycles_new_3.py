#!/usr/local/bin/env python
# -*- coding: utf-8 -*-

"""
#######################################
#                                     #
#-- Counting loops using cocept of Nets --#
#------  Author: Devosmita Sen  --------#
#                                     #
#######################################

 Overall Framework (Steps):
     1. 
"""
import os.path
import sys
import time
import math
import random
import matplotlib
matplotlib.use('Agg')
import numpy as np
from timeit import default_timer as timer
import os
from matplotlib.pyplot import cm
import shortest_path_DS
import multiprocessing as mp

import networkx as nx
from itertools import combinations
##from relax import Optimizer
from numpy import linalg as LA
from scipy.optimize import fsolve
from matplotlib import pyplot as plt
##import param as p
import shutil

class node_part: 
    def __init__(self, M, y): 
        self.M = M 
        self.y = y


def generate_graph(netgen_flag,test, graph_algo, lattice_type,ite): # specify lattice type only if usig lattice, else: None
    directory = './original_files/'
    orig_dir = os.path.dirname(directory)
    files=os.listdir(orig_dir)

    directory = './'#+'network'+'/'
    if not os.path.isdir(directory):
      os.mkdir(directory)

    for fname in files:
     
    # copying the files to the
    # destination directory
       shutil.copy2(os.path.join(orig_dir,fname), directory)

    # now add path to frac_weak directory
    file_dir = os.path.dirname(directory)
    sys.path.append(file_dir)
    import ioLAMMPS
    import netgen


    G=nx.Graph()
    Gmult=nx.MultiGraph()
    if(test==False):
        if(graph_algo=='Gusev'):
            
            if(netgen_flag==0):
              print(graph_algo)
              vflag = 0
            ##   N = 12   
              print('--------------------------')   
              print('----Reading Network-------')   
              print('--------------------------')
              
              filename = "network_"+str(ite)+".txt"
              file_path = os.path.join(directory, filename)
              if not os.path.isdir(directory):
                 os.mkdir(directory)  
              [xlo, xhi, ylo, yhi, zlo, zhi, n_atoms, n_bonds, atoms, bonds, 
                          atom_types, bond_types, mass, loop_atoms,num_secondary_loops] = ioLAMMPS.readLAMMPS_into_graph(G,Gmult,file_path, vflag, frac_weak)
              n_chains=n_bonds
              print('len(loop_atoms)',len(loop_atoms))
              print('loop_atoms',loop_atoms)
            elif(netgen_flag==1): # Gusev network generation

              func=p.func
              l0   = 1
              prob = 1.0
              n_chains  = p.n_chains
              n_links   = int(2*n_chains/func)
              L=p.L
            
              num_secondary_loops=netgen.generate_network(G,prob, func, parameters,L, l0, n_chains, n_links, frac_weak)
              directory = './'+str(int(100*frac_weak))+'/' # while generating network- separate folder created so that existing network data doesn't get deleted
              filename = 'network.txt'
              file_path = os.path.join(directory, filename)
              if not os.path.isdir(directory):
                 os.mkdir(directory)  
              
              [xlo, xhi, ylo, yhi, zlo, zhi, n_atoms, n_bonds, atoms, bonds, 
                      atom_types, bond_types, mass, loop_atoms] = ioLAMMPS.readLAMMPS(file_path,0,frac_weak)
              n_chains=n_bonds
            else:
                print('Invalid network generation flag for graph algorithm: ',graph_algo)

        elif(graph_algo=='KMC'):
            
            if(netgen_flag==0):
              print(graph_algo)

              vflag = 0
            ##   N = 12   
              print('--------------------------')   
              print('----Reading Network-------')   
              print('--------------------------')

              if(ite==-1):
                  filename = "3_network_KMC.txt"
              else:            
                  filename = "3_network_KMC_"+str(ite)+".txt"
              file_path = os.path.join(directory, filename)
              if not os.path.isdir(directory):
                 os.mkdir(directory)  
            ##  [xlo, xhi, ylo, yhi, zlo, zhi, n_atoms, n_bonds, atoms, bonds, 
            ##              atom_types, bond_types, mass, loop_atoms,num_secondary_loops] = ioLAMMPS.readLAMMPS_into_graph(G,file_path, vflag, frac_weak)
              [loop_atoms, n_chains]=ioLAMMPS.readLAMMPS_into_graph_from_bond_only(G,Gmult,file_path)#, vflag, frac_weak)
##              print(935 in G)
##              stop
              print('len(G)',len(G))
##              stop
        elif(graph_algo=='lattice'):
            
            if(netgen_flag==0):
              print(graph_algo)

              vflag = 0
            ##   N = 12   
              print('--------------------------')   
              print('----Reading Network-------')   
              print('--------------------------')
              
              filename = lattice_type+"_connectivity.txt"
              file_path = os.path.join(directory, filename)
              if not os.path.isdir(directory):
                 os.mkdir(directory)  
            ##  [xlo, xhi, ylo, yhi, zlo, zhi, n_atoms, n_bonds, atoms, bonds, 
            ##              atom_types, bond_types, mass, loop_atoms,num_secondary_loops] = ioLAMMPS.readLAMMPS_into_graph(G,file_path, vflag, frac_weak)
              n_chains=ioLAMMPS.readLAMMPS_into_graph_lattice(G,file_path)#, vflag, frac_weak)    

            else:
              print('Invalid network generation flag for graph algorithm: ',graph_algo)
        else:
              print('Invalid network graph algorithm: ')


    elif(test==True):
####        TEST GRAPHS
        G = nx.Graph()
        G.clear()



        # test graphs containing primary and secondary loops
##        G.add_edge(1,2)
##        G.add_edge(2,1)
##        G.add_edge(1,1)
##        G.add_edge(2,2)


        if(netgen_flag==0):
              print(graph_algo)

              vflag = 0
            ##   N = 12   
              print('--------------------------')   
              print('----Reading Network-------')   
              print('--------------------------')
              
              filename = "network_test.txt"
              file_path = os.path.join(directory, filename)
              if not os.path.isdir(directory):
                 os.mkdir(directory)  
            ##  [xlo, xhi, ylo, yhi, zlo, zhi, n_atoms, n_bonds, atoms, bonds, 
            ##              atom_types, bond_types, mass, loop_atoms,num_secondary_loops] = ioLAMMPS.readLAMMPS_into_graph(G,file_path, vflag, frac_weak)
              [loop_atoms, n_chains]=ioLAMMPS.readLAMMPS_into_graph_from_bond_only(G,Gmult,file_path)#, vflag, frac_weak)
##              print(935 in G)
##              stop
              print('len(G)',len(G))


              
        '''
        nx.add_cycle(G,[1,2,5,4])
        nx.add_cycle(G,[1,4,3,2])
        nx.add_cycle(G,[1,4,6,2])
        '''

        
##        #square grid
##        nx.add_cycle(G,[0,1,4,5])
##        nx.add_cycle(G,[1,4,3,2])
##        nx.add_cycle(G,[4,3,8,7])
##        nx.add_cycle(G,[4,7,6,5])


##        nx.add_cycle(G,[1,2,5,4])
##        nx.add_cycle(G,[1,4,3,2])
##        nx.add_cycle(G,[1,4,6,2])
##
##        G.add_edge(1,2)
##        G.add_edge(1,2)
##        print(list(G.neighbors(1)))
##
##        # 4-7 ring- text example
##        nx.add_cycle(G,[1,2,3,4,5,6,7])
##        nx.add_cycle(G,[7,6,9,8])
##        nx.add_cycle(G,[6,5,10,9])
##        nx.add_cycle(G,[14,13,12,11,5,6,7])
##        G.add_edge(5,11)
##        G.add_edge(11,9)
##
##        nx.add_cycle(G,[6,5,10,9,8,7])
##
##        ##   Testing example 1
##        nx.add_cycle(G,[1,2,3,4])
##        nx.add_cycle(G,[1,5,3,6])
##        G.add_edge(1,3)
##
##        ##   Testing example 2
##        nx.add_cycle(G,[1,2,4,3])#[1,2,3,4]or [1,2,3,4]
##        nx.add_cycle(G,[3,5,2,4])
##        G.add_edge(1,3)
##
##        new test case:
##        nx.add_cycle(G,[1,2,3,4,5,6,7,8])
##        nx.add_cycle(G,[6,5,9,7])
##        G.add_edge(1,10)
##        G.add_edge(1,13)
##        G.add_edge(3,11)
##        G.add_edge(11,12)

####
####        loop_atoms=[]
####        for i in G.edges:
####            link_1=i[0]
####            link_2=i[1]
####            if(link_1==link_2):
####             loop_atoms.append(link_1)
##        DRAW NETWORK- commented out for larger systems since python cannot draw very big networks
        nx.draw(G, with_labels=True,node_color="tab:pink")
        plt.savefig('network.png',bbox_inches='tight')



##        plt.show()
    
    return G,Gmult, loop_atoms, n_chains



def count_cycles(node_list,G,max_loop_order,node_number,return_all_paths_list,return_main_node_list,return_tot_node_number,proc_count):
    print('inside function count_cycles')
    main_node_list=[]
    all_paths_list=[]
    node_number=[]
    tot_node_number=0 # gives the node given the index- node_number[3]=46 for example
    G_original=G
    count_outermost=0
    for i in node_list: # outer loop- loop 1
        ## CYCLES CALCULATION AROUND EACH OF THESE NODES CAN BE PARALLELIZED WRT COMPUTATION (NOT MEMORY)
      count_outermost=count_outermost+1
      if(count_outermost%100==0):
          print('node sr. num',count_outermost)
    ##      i=6
    ##      visited=np.zeros(len(G)+1) #IN A REAL NETWORK- MANY INDICES ARE SKIPPED. SO IF I JUST TAKE THE len(), THEN IT WILL BE WRONG- I WILL HAVE TO TAKE THE MAX
      visited=np.zeros(max(G)+1)
      printed=0

      if(len(list(G.neighbors(i)))>1): #if >1- then primary loops and molecules without any loops are excluded, else if !=0 or >0-#includes primary loop, dangling ends are also excluded
         internal_count=0
         printed=0
         main_node=i
         main_node_list.append(i)
         tot_node_number=tot_node_number+1
    ##     num_stars[str(main_node)]=0 #NOT REQUIRED HERE
         A=[]
####         f.write(str(main_node))
####         f.write(' ')
         ####print('-----------------------------')
    ##         print('main_node',main_node)
         node_number.append(main_node)
         temp_node_names=[]

         neigh=[n for n in G.neighbors(main_node)]
         new_path_required=0 #needed later
         inserted=0
    #CAN'T COUNT PRIMARY LOOPS HERE BECAUSE THE NETWORK HAS ALREADY BEEN PRUNED WRT PRIMARY AND SECONDARY LOOPS         
    ##         for i in neigh:
    ##             if(i==main_node):
    ##                 loop_dist[1-1]=loop_dist[1-1]+1; # primary loop detected
      ##   for n in neigh:
      ##      if(nx.shortest_path()==1 and nx.shortest_path()==1):
      ####IMPORTANT: I cannot track secondary loops or any structures which involve connections b/n only 2 nodes- will have to track them differently 
         res = list(combinations(neigh, 2))
         if(len(res)==0):
            print('primary loop- detected before for loop starting')
         for i_node_pair in range(0,len(res)):
    ##            print(res)
            Gcopy=G.copy()
            node_pair=res[i_node_pair]
    ########            print('node_pair',node_pair)
      ##      if(node_pair[0]==main_node or node_pair[1]==main_node):
      ##         print('primary loop detected')
      ##         continue
            
##            Gcopy.remove_node(main_node) #THIS IS REMOVING ALL THE EDGES OF MAIN_NODE!! BUT WE WANT TO ONLY REMOVE THE EDGES CONNECTING THE MAIN NODE TO THE TWO NODES IN CONSIDERATION
            Gcopy.remove_edge(node_pair[0],main_node)
            Gcopy.remove_edge(node_pair[1],main_node)
##            print('node_pair[0]',node_pair[0])
##            print('node_pair[1]',node_pair[1])
##            print('main_node',main_node)
    ##            for n in neigh:
    ##                if(n not in [node_pair[0],node_pair[1]]):
    ##                    Gcopy.remove_node(n) # this ensures that a path having shorted connection to the main node is eliminated
      ##      if(node_pair[0] in G or node_pair[1] in G):
      ##         print('Some error, primary loop NOT detected properly earlier')
      ##         stop
            try:
                min_path=nx.all_shortest_paths(Gcopy, source=node_pair[0], target=node_pair[1])
    ##                min_path_temp=nx.shortest_path(Gcopy, source=node_pair[0], target=node_pair[1])
                a=list(min_path)
    ##                print(a)
                a=[x for x in a if len(x)<max_loop_order]#!=path[index]]
    ##                print(a)
    ##                stop
                A.append(a)
            except nx.NetworkXNoPath:
                continue
         long_loops_removed_correctly=False   
         if(len(A)==0):
    ##             Anew=A # all these cases are written separately because each is different,. and some wont compute in other cases
             
             Anew=sorted(A,key=lambda l:len(l[0]))
    ##             stop
    ##         elif(len(A[0])==0):
    ##             Anew=sorted(A,key=lambda l:len(l[0]))
    ##             stop
         elif([] in A):
             Anew=[x for x in A if x!=[]]
             Anew=sorted(Anew,key=lambda l:len(l[0]))
    ##             stop
    ##         elif([] in A[0]):
    ##             Anew=sorted(A,key=lambda l:len(l[0]))
    ##             stop
    ##         if(False):
    ##             stop
         else: # eg. if A=[[]]- then there will be problem with sorting
             Anew=sorted(A,key=lambda l:len(l[0]))
             long_loops_removed_correctly=True

                  
    #sort A wrt path length and then do this analysis        
         Acopy=Anew.copy();
         path_count=-1
         
         count_out=0
         for path in Acopy:
             if(len(path)<max_loop_order):
                count_out=count_out+1
    ##                print('count_out',count_out)
    ##                path_count=path_count+1
    ##                print('path count', path_count)
    ##          else:
    ########                print('path',path)
                count_t=0
                for index in range(0,len(path)): # if there are degenerate/equivalent ones in the same path
                    count_t=count_t+1

    ##              index=0 # because always the first term in the path is onyl being tested- without this consition- a list with multiple paths which are non fundamental- onyl one path is removed,a nd not the other ones
                    #CHECKING IF ALL NODES ALONG THAT PATH HAVE BEEN VISITED
                    #IF YES- THEN THAT PATH IS REOMOVED AND WE WILL BE LOOKING FOR NEW PATH
                    if(index<len(path) and all(visited[test_node]==1 for test_node in path[index]) and (path[index][0] not in G.neighbors(path[index][-1]))): # means that all nodes of that path have already been visited
                           # the last condition of the loop is neccesary because the only case when 'checking whether all nodes in that path have been visited' - for removing that path is when the starting nodes themselves are connected
                           # in that case- that particular connection is not yet visited even though the nodes have been visited!
                                           
                       removed_path=path[index]
                       new_path_start=path[index][0]
                       new_path_end=path[index][-1]
                       len0=len(path)
    ##                            to_remove=path.copy()
    ##                       Acopy=Anew.copy()
                       for temp in range(0,len(Anew)):
                          if(path[index] in Anew[temp]):
                             
    ##                                    if(to_remove==[5,10,9,8,7]):
    ##                                       stop
    ##                                    print('path just before removing',path)
                             Anew[temp]=[x for x in Anew[temp] if x!=path[index]] # remove that path
    ##                                    Anew[temp].remove(to_remove[index])
                             #IMPORTANT: USING remove() MANIPULATES EVERYTHING- CREATES PROBLEMS!!!
                             # USE A MODIFIED LIST INSTEAD OF remove()!!!
    ##                                    print('path',path)
                             if(len(path)!=len0):
                                stop
                             break
                       

                       new_path_required=0
                       if([] in Anew):
    ##                          stop
                          new_path_required=1 #this is done so that we search for new paths only if we have
                          # exhausted every other path and all of them are redundant
                          Anew.remove([])
    ##                            elif(index<=len(path)-1):
    ##                                continue
                                        
    ##                       print(Anew)
    ##                            stop
    ##                            stop
    ##                            if(index==len(path)
                       if(new_path_required==1):
                          new_paths=shortest_path_DS.shortest_simple_paths_DS(Gcopy, source=new_path_start, target=new_path_end,max_loop_order=max_loop_order)
##                          print('before len(new_paths)')
##                          print('len(new_paths)',len(list(new_paths)))
                          new_paths_shortest=[]
##                          shortest_path_length=len(next(new_paths))
##                          print('shortest_path_length',shortest_path_length)
                          count_tmp=0
##                          for i in new_paths:
##                              print(i)
####                              break
##                          for i in new_paths:
##                              print('i2',i)
##                              break
                          for i in new_paths:
                              if(count_tmp==0):
##                                  print('i_count_tmp_0',i)
                                  new_paths_shortest.append(i)
                                  shortest_path_length=len(i)
##                                  print('shortest_path_length',shortest_path_length)
                                  count_tmp=count_tmp+1
                              else:
##                                  print('i_count_tmp_not_0',i)
                                  if(len(i)==shortest_path_length):
                                      new_paths_shortest.append(i)
                                  elif(len(i)<shortest_path_length):
                                      print('shortest_path_length',shortest_path_length)
                                      print('len(i)',len(i))
                                      print('i',i)
                                      print('new_paths_shortest',new_paths_shortest)
                                      print('count_tmp=',count_tmp)
                                      print('new_path_start=',new_path_start)
                                      print('new_path_end=',new_path_end)
                                      print('PROBLEM!!')
                                      for j in new_paths:
                                          print(j)
                                      stop
                                  else:
                                      break # can do this because the paths are arranged in increasing order of lengths. so if one path has higher length, all subsequent paths will have higher length!!
                                  count_tmp=count_tmp+1
##                          stop
##                          for tmp_path in new_paths:
##                              if(len(tmp_path)>10):
##                                  print('tmp_path: ',tmp_path)
                          
##                                  stop
                          loop_cnt=0
    ##                          length=len(new_paths)
    ##                          print('length=',length) ##!!!IMPORTANT: THERE IS SOME PROBLEM WITH PRINTING len(new_paths)
                          #MORE THAN ONCE- IT PRINTS IT OUT AS 0 LIST [] BUT IT IS ACTUALLY NOT SO!!
                          # THAT IS WHY IT IS GOING INTO THE LOOP BUT I AM NOT ABLE TO DEBUG!!

    ##                          new_paths=sorted(new_paths,key=lambda l:len(l[0]))
    ##                          count_t=0
                          T=0
##                          print('len(new_paths_shortest)=',len(list(new_paths_shortest)))
                          
##                          stop
                          for temp_path in (new_paths_shortest): # temp_path can contain more that one equivalent path
##                             print('len(tmp_path)=',len(tmp_path))
    ##                             print('(temp_paths)',(temp_path))
                             if(len(temp_path)>10):
                                 print(len(temp_path))
                                 print(temp_path)
##                                 stop
    ##                           if(all(len(x)<max_loop_orderlen(temp_path)<max_loop_order):
    ##                             if(count_outermost==55):
    ##                             temp_path=[x for x in temp_path if len(x)<max_loop_order]
                             T=T+1 

                             loop_cnt=loop_cnt+1
                         
                             for test_node in temp_path:# equivalent paths 
    ##                                 for test_node in eq_path: 
                                     if(all(visited[test_node]==1 for test_node in temp_path)): # means that all nodes of that path have already been visited
            ##                                    i=i+1 #this path is not valid
                                         #THIS HAS TO BE CHANGED
                                        continue
                                     else:
                                        inserted=0
                                        if(len(Anew[-1][0])<len(temp_path)): # new path has size greater than the paths earlier considered
                                            # here- need to be careful about finding a new loop of higher isze which is basically just an extension of the older loop which is removed
                                            # can be visualized very clearly in the square grid problem
                                                                            
                                            if(len(set(removed_path) & set(temp_path))<=2): # always will be 2 alt leasdt
                                               Anew.append([temp_path]) # path is valid, add to min_path list
        ########                                       print('inserted path of higher size')
                                               inserted=1
        ########                                       print('cond1,inserted=1',temp_path)
                                        else:
                                            for temp_idx in range(0,len(Anew)):
                                               if(len(Anew[temp_idx][0])==len(temp_path) and temp_path not in Anew[temp_idx]):
                                                  Anew[temp_idx].append(temp_path)
        ########                                          print('inserted another path of same existent size')
                                                  inserted=1
        ########                                          print('cond2,inserted=1',temp_path)


                                # the other case is when the new path has size which is less than the max size, but alsonot equal to any of the sizes present-
                                # ie. the size is in between- THIS CASE WILL NEVER ARISE- BECAUSE IF THE NEW PATH HAS LOWER SIZE THAN THE PATH BEING REMOVED,
                                #THEN IT WOULD ALREADY HAVE BEEN CONSIDERED EARLIER, AND THIS CASE WOULD NEVER HAVE COME UP IN THE FIRST PLACE

    ##                                if(new_path_required==1 and inserted==0):
    ##                                    num_stars[str(main_node)]=num_stars[str(main_node)]+1
                             for node in temp_path:
                                    visited[node]=1
                             printed=1

                                           # look through shortest_simple_path
                    elif(index<len(path)):            
                        for node in path[index]:
                            visited[node]=1



         if(printed==1):

                all_paths_list.append(Anew)


        # CHECKING AND COUNTING NUMBER OF * IN VERTEX SYMBOL

####                f.write(str(Anew))
    ##          for min_path in Anew:
    ##              node_part_My=node_part(min_path_length,len(a))# M and y
    ##          f.write('\n')

                                        
                                    # add to list of valid loops
                                    # path is valid
                  
                    

         if(printed==0 and len(list(G_original.neighbors(i)))>1):
    ########            print('written for node',i)
####            f.write(str(Anew))
            all_paths_list.append(Anew)

####         f.write('\n')
         num_neigh=len(res) # number of neighbors
             
####    f.close()
##    return tot_node_number
##    all_paths_list,main_node_list,tot_node_number,f,,i,return_dict
    return_all_paths_list[proc_count]=all_paths_list
    return_main_node_list[proc_count]=main_node_list
    return_tot_node_number[proc_count]=tot_node_number
    print('updated return dict()s for proc_count=',proc_count)
    


def symbol_str(a,main_node):
    # a is the vertex symbol of a node
    string=''
    for i in range(0, len(a)-1):
        string=string+str(a[i].M)+str(a[i].y).translate(subscript)+"."
    i=len(a)-1
    if(i>=0):
        string=string+str(a[i].M)+str(a[i].y).translate(subscript)
        for j in range(0,num_stars[str(main_node)]):
            string=string+'.*'
    return string

def symbol_str_no_subscript(a,main_node):
    # a is the vertex symbol of a node
    string=''
    for i in range(0, len(a)-1):
        string=string+str(a[i].M)+str(a[i].y).translate(subscript)+"."
    i=len(a)-1
    if(i>=0):
        string=string+str(a[i].M)+'_'+str(a[i].y)
        for j in range(0,num_stars[str(main_node)]):
            string=string+'.*'
    return string


if __name__=='__main__':
    
    ###############GENERATE NETWORK#######################
    start1=timer()
    random.seed(a=None, version=2)
    print('First random number of this seed: %d'%(random.randint(0, 10000))) 
    # This is just to check whether different jobs have different seeds
    ##global parameters
    ###parameters=np.zeros([2,6]) # N, b, K, fit_param, E_b,U0
    ###parameters[0,:]=np.array([p.N_low,p.b_low,p.K_low,p.fit_param_low,p.E_b_low,p.U0_low])
    ###parameters[1,:]=np.array([p.N_high,p.b_high,p.K_high,p.fit_param_high,p.E_b_high,p.U0_high])

    frac_weak=0.0


    netgen_flag = 0
    test=False
    graph_algo='KMC'
    lattice_type='None'
    func=4 # functionality of network
    num_nei=int(func*(func-1)/2)
    max_loop_order=60 # threshold loop order - will not count beyond this!
    fraction_to_count=1.0 # fraction of nodes on which the cycle counting should be done
    max_loop_order_to_consider=func*(func-1)/2.0 #ONLY VALID FOR LATTICES      
    #######################################
    NRA=np.genfromtxt('NRA.txt')
    ite_arr=np.arange(0, 2*NRA,NRA/10,dtype='int')
    ite_arr=np.append(ite_arr, -1)
    for ite in ite_arr:
      try:
        print('ite',ite)
        G,Gmult,loop_atoms,n_chains=generate_graph(netgen_flag,test, graph_algo,lattice_type,ite)# specify lattice type only if usig lattice, else: None
    
    #######################################
##    print(935 in G)
##    stop
        print('outside len(G)',len(G))
        end1 = timer()
        G_orig=G.copy()
    ##    print('outside: G.neighbor(4)',list(G_orig.neighbors(4)))

    ##    stop
    ##    print(len(G))
    ##    stop
        ######################STARTING LOOP COUNTING######################
        print('Starting LOOP COUNTING...')
        start=timer()
        ##   nx.draw(G, with_labels=True, font_weight='bold',node_color="tab:pink")
        '''
        ax = plt.gca()
        ax.set_title('Gussev network- all primary loops EXCLUDED, secondary loops shown as single connection')
        nx.draw(G, with_labels=True,font_size=7,node_color="tab:pink", node_size=85,ax=ax)
        plt.savefig('network_gussev.png')
        '''
        ##node_names=[]#np.zeros([G.size(),7]) # [0]- main_node, [1:6]- all the loop sizes for each pair of edges selected- this is not exactly in vertex symbol notation, but we can sort this later and make it vertex symbol
        

        
        ##main_node_list=[]
        ##all_paths_list=[]
        ##tot_node_number=0 # gives the node given the index- node_number[3]=46 for example
        ##num_stars={} # NUMBER OF STARS (*) IN VERTEX SYMBOL OF EACH NODE- python dictionary
        print('max_loop_order',max_loop_order)

        ##Removing nodes with degree less than 2 or coreness<2- these nodes will not be considered in the analysis
        ####just to reduce the number of computations
        removeNodes =set()

    ####    G_no_primary_loop=G.copy()
    ##    print('935 in G',935 in G)
    ##    print('935 in G_no_primary_loop',935 in G_no_primary_loop)
    ##    stop
    ####    for i in loop_atoms:
    ####        if(i in G_no_primary_loop): # if there are 2 primary loops on the same node, then this condition will not be satisfied, but that's okay because we anyways have to remove the node for coreness calculation
    ####            G_no_primary_loop.remove_node(i)
    ##        print('removed i=',i)
    ##    Coreness = nx.core_number(G_no_primary_loop)
        Coreness = nx.core_number(G)

    ####    for i in G_no_primary_loop:  #
    ####        if  G.degree(i)<=1 or Coreness[i] <= 1: # DONT USE CORENESS FOR CAYLEY TREES!!!
    ####            removeNodes.add(i)
    ##        if G.degree(i)<4 : ##DON'T DO THIS!!! IT BREAKS BONDS AND THUS CREATES MORE NODE WITH DEGREE<4
    ##            removeNodes.append(i)
        ##stop


        for i in G:  #
            if  G.degree(i)<=1 or Coreness[i] <= 1: # DONT USE CORENESS FOR CAYLEY TREES!!!
                removeNodes.add(i)
        G.remove_nodes_from(removeNodes) # HERE I AM ONLY REMOVING NODES WHICH ARE DEGREE=1 OR CORENESS=1, SO WON'T AFFECT THE CYCLE COUNTING EVEN IF IT MAKES THE NODE PARTIALL REACTED- WONT CHANGE CYCLE SIZE IN THIS CASE
        #BUT IN GENERAL: DONT REMOVE NODES!!! (ESPECIALLY FOR PRIMARY LOOP PRUNING!!) BECAUSE THAT CAUSES MORE NODES TO BECOME UNREACTED!!
        num_primary=0
        num_secondary=0
        for i in Gmult.edges:
            u=i[0]
            v=i[1]
            if(u==v):
                num_primary=num_primary+1
            
        print('num_primary',num_primary)
        
        edge_chain=list(Gmult.edges)
        for i in range(0,len(Gmult.edges)):
           for j in range(0,i):
               
               if ((edge_chain[i][0]==edge_chain[j][0]) and (edge_chain[i][1]==edge_chain[j][1])) or ((edge_chain[i][1]==edge_chain[j][0])and (edge_chain[i][0]==edge_chain[j][1])):
                   num_secondary=num_secondary+1
        num_secondary=num_secondary/2.0
        print('num_secondary',num_secondary)

        num_unreacted=0
        for i in Gmult:
            num_unreacted=num_unreacted+(4-Gmult.degree(i))
        tot_junctions=4*len(Gmult)
        frac_unreacted=num_unreacted/tot_junctions
        print('frac_unreacted',frac_unreacted)
        num_partially_react_nodes=0
        for i in Gmult:
            if(Gmult.degree(i)<4):
                num_partially_react_nodes=num_partially_react_nodes+1
        print('number of partially reacted nodes in Gmult=',num_partially_react_nodes)

        num_partially_react_nodes=0
        for i in G:
            if(G.degree(i)<4):
                num_partially_react_nodes=num_partially_react_nodes+1
        print('number of partially reacted nodes in G=',num_partially_react_nodes)

        num_partially_react_nodes=0
        for i in G_orig:
            if(G_orig.degree(i)<4):
                num_partially_react_nodes=num_partially_react_nodes+1
        print('number of partially reacted nodes in G_orig=',num_partially_react_nodes)
        
    ##    stop


        node_number=[]
        all_node_list=[]
        for i in G: # outer loop- loop 1
    ##      print('in loop 1')
    ##      print('degree=',G.degree(i))
            ## CYCLES CALCULATION AROUND EACH OF THESE NODES CAN BE PARALLELIZED WRT COMPUTATION (NOT MEMORY)
            if(graph_algo=='lattice'):
              if(G.degree(i)==func): #and all(G.degree(G.neighbors(i))==f)): ## CONSIDERING ONLY THOSE NODES WHICH HAVE FULL FUNCTIONALITY AND THOSE NODES WHICH ARE NOT ON BOUNDARY
                  #IF THE NODE IS ON THE BOUNDARY, THEN EITHER ITS DEGREE WILL BE <f OR IT WILL BE CONNECTED TO SOME NODE WHOSE DEGREE IS <f
        ##          print('in loop 2')
        ##          stop
                  count_node_i=True
                  for j in nx.descendants_at_distance(G, i, max_loop_order_to_consider-2):#node_connected_component(G, i):#G.neighbors(i):
                      # not only nearest neighbors, but also the nodes connected at distance k=(max loop size expected in lattice-2)=4 here
                      if(G.degree(j)<func):
                          count_node_i=False
        ##                  print('breaking')
                          break
                  if(count_node_i==True):
                      all_node_list.append(i)#list(G.nodes)

                      if(len(list(G.neighbors(i)))>1): #if >1- then primary loops and molecules without any loops are excluded, else if !=0 or >0-#includes primary loop, dangling ends are also excluded
                         node_number.append(i)
    ##        elif(graph_algo=='KMC'):
    ##          if(G.degree(i)==func): #and all(G.degree(G.neighbors(i))==f)): ## CONSIDERING ONLY THOSE NODES WHICH HAVE FULL FUNCTIONALITY AND THOSE NODES WHICH ARE NOT ON BOUNDARY
    ##              #IF THE NODE IS ON THE BOUNDARY, THEN EITHER ITS DEGREE WILL BE <f OR IT WILL BE CONNECTED TO SOME NODE WHOSE DEGREE IS <f
    ##    ##          print('in loop 2')
    ##    ##          stop
    ##              count_node_i=True
    ##              for j in nx.descendants_at_distance(G, i, max_loop_order_to_consider-2):#node_connected_component(G, i):#G.neighbors(i):
    ##                  # not only nearest neighbors, but also the nodes connected at distance k=(max loop size expected in lattice-2)=4 here
    ##                  if(G.degree(j)<func):
    ##                      count_node_i=False
    ##    ##                  print('breaking')
    ##                      break
    ##              if(count_node_i==True):
    ##                  all_node_list.append(i)#list(G.nodes)
    ##                  rand_num=random.random()
    ##                  if(len(list(G.neighbors(i)))>1 and rand_num<=fraction_to_count): #if >1- then primary loops and molecules without any loops are excluded, else if !=0 or >0-#includes primary loop, dangling ends are also excluded
    ##                     node_number.append(i)
            else:
                rand_num=random.random()
                if(len(list(G.neighbors(i)))>1 and rand_num<=fraction_to_count): #if >1- then primary loops and molecules without any loops are excluded, else if !=0 or >0-#includes primary loop, dangling ends are also excluded
                     node_number.append(i)
                all_node_list.append(i)

        start_mp=timer()
        ######### BEGIN MULTIPROCESSING #########################
        print('begin MULTIPROCESSING')
    ##    all_node_list=list(G.nodes) # equivalent to list(G)


        
        cpu = 20#os.cpu_count()-1  # No. of parallel processors
        N = len(list(G.nodes))     # Total no. of nodes in graph
        n = int(N/cpu)      # No. of nodes to be analyzed (for cycles) in each subprocess
        n_rem = N%cpu #remianing number of nodes to be given to another process
        extra_cpu=np.heaviside(n_rem,0) # will give value 0 or 1 - implies: whether the extra cpu is required or not
        cpu=cpu+int(extra_cpu)
        node_list = np.zeros((cpu),dtype=object)  # Total no. of node_list are sub-divided into multi-segment for mp
        for i in range(cpu-1):
                node_list[i] = all_node_list[i*n:(i+1)*n]
                i_final=i
        node_list[cpu-1]=all_node_list[(i_final+1)*n:len(all_node_list)]
    ##    print('len(G.nodes)',len(G.nodes))
    ##    print('len(node_list[0])',len(node_list[0]))
        manager = mp.Manager()
    ##    print('after manager')
        return_all_paths_list = manager.dict()
    ##    print('after manager.dict() 1')
        return_main_node_list = manager.dict()
    ##    print('after manager.dict() 2')
        return_tot_node_number = manager.dict()
    ##    print('after manager.dict() 3')


        #######################################
        processes = [mp.Process(target=count_cycles, args=(node_list[i],G,max_loop_order,node_number,return_all_paths_list,return_main_node_list,return_tot_node_number,i)) for i in range(cpu)]
        #######################################
        print('len(processes)=',len(processes))
        
        for process in processes:
            process.start()
            print('started process')
            # wait for all processes to finish
        for process in processes:
            process.join()
            print('joined process')
        print('completed processes')
        FINAL_all_paths_list=[]
        FINAL_main_node_list=[]
        FINAL_tot_node_number=0
        for i in range(cpu):
            FINAL_all_paths_list=FINAL_all_paths_list+return_all_paths_list[i]
            FINAL_main_node_list=FINAL_main_node_list+return_main_node_list[i]
            FINAL_tot_node_number=FINAL_tot_node_number+return_tot_node_number[i]
        ##    print('result[i]',result[i])

        print('end MULTIPROCESSING')    
        ##count_cycles(node_list,G,max_loop_order,return_all_paths_list,return_main_node_list,return_tot_node_number,i):
        #############END MULTIPROCESSING#######################

        end_mp=timer()
        print('Multiprocessing time=',end_mp-start_mp)



        num_stars={} # NUMBER OF STARS (*) IN VERTEX SYMBOL OF EACH NODE- python dictionary

        all_vertex_symbols=[]
        for Anew_node in FINAL_all_paths_list:
            temp_node_names=[]
            for minpath in Anew_node:

                node_part_My=node_part(len(list(minpath)[0])+1,len(list(minpath)))# M and y
                temp_node_names.append(node_part_My)
            all_vertex_symbols.append(temp_node_names)

        temp_count2=0
        all_vertex_symbols_temp=all_vertex_symbols.copy()
        for main_node2 in FINAL_main_node_list:
                neigh=[n for n in G.neighbors(main_node2)]
                num_stars[str(main_node2)]=len(list(combinations(neigh, 2)))-len(all_vertex_symbols_temp[temp_count2]) # WHY OS THIS LINE CREATING PROBLEMS!?????- HAVE TO CHECK!!
                temp_count2=temp_count2+1
        #print(len(all_vertex_symbols[0]))
        f=open('symbols.txt','w')
        count_node=0
        loop_order_count=np.zeros(max_loop_order+1)
        loop_order_count_mean_field=np.zeros(max_loop_order+1)

        subscript = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

        temp_count=0-1
        for symbol in all_vertex_symbols:
            f.write(str(FINAL_main_node_list[count_node]))
            f.write(' ')
            temp_count=temp_count+1
            ####print('FINAL_main_node_list[count_node]',FINAL_main_node_list[count_node])
            num_cycles_in_symbol=0 # number of actual cycles in the vertex symbol (without the *)
            
            neigh=[n for n in G.neighbors(FINAL_main_node_list[count_node])]
            
            junc_func=len(neigh) # number of connections of the node
            ####print('junc_func',junc_func)
            
            num_actual_cycles=len(all_vertex_symbols[count_node])
            ####print('num_actual_cycles',num_actual_cycles)
            
            ####print('num_stars[str(FINAL_main_node_list[count_node])]',num_stars[str(FINAL_main_node_list[count_node])])
            ysum=0
            for i in range(len(all_vertex_symbols[count_node])):           
                #M=all_vertex_symbols[count_node][i].M
                y=all_vertex_symbols[count_node][i].y
                ysum=ysum+y
                
            for i in range(len(all_vertex_symbols[count_node])):
        ##        for j in len(all_vertex_symbols[count]:
                f.write('(')
                f.write(str(all_vertex_symbols[count_node][i].M))
                f.write(',')
                f.write(str(all_vertex_symbols[count_node][i].y))
                f.write(') ')    
                
                M=all_vertex_symbols[count_node][i].M
                y=all_vertex_symbols[count_node][i].y
                #loop_order_count[M]=loop_order_count[M]+y
                if(y>0):
                    #loop_order_count[M]=loop_order_count[M]+((junc_func-2)/2)*1/num_actual_cycles
                    loop_order_count[M]=loop_order_count[M]+((junc_func-2)/2)*y/ysum
                else:
                    stop
                    
            '''
            for i in range(len(all_vertex_symbols[count_node])):
        ##        for j in len(all_vertex_symbols[count]:
                f.write('(')
                f.write(str(all_vertex_symbols[count_node][i].M))
                f.write(',')
                f.write(str(all_vertex_symbols[count_node][i].y))
                f.write(') ')    
                
                M=all_vertex_symbols[count_node][i].M
                y=all_vertex_symbols[count_node][i].y
                #loop_order_count[M]=loop_order_count[M]+y
                if(y>0):
                    loop_order_count[M]=loop_order_count[M]+1
            '''
            
            for j in range(0,num_stars[str(FINAL_main_node_list[count_node])]):
                f.write('.*')
            f.write('\n')              
            count_node=count_node+1



        #  mean field loop order 
        for i in range(1,len(loop_order_count)):
            loop_order_count_mean_field[i]=loop_order_count[i]#/i
            #loop_order_count_mean_field[i]=math.ceil(loop_order_count_mean_field[i])
            if (int(loop_order_count_mean_field[i] - 0.5) == loop_order_count_mean_field[i] - 0.5):  # this means that the number is something.5 (eg. 2.5, 3.5, etc), in this case, always round up
                loop_order_count_mean_field[i] = math.ceil(loop_order_count_mean_field[i])  # ceil
            else:
                loop_order_count_mean_field[i] = round(loop_order_count_mean_field[i])  # ceil


        loop_order_count_mean_field[1]=num_primary
        loop_order_count_mean_field[2]=num_secondary
        ##loop_order_count_mean_field=np.sqrt(loop_order_count)
        plt.figure()#figsize=(3, 3))
        ##loop_order_count_mean_field[1]=len(loop_atoms)
        ##loop_order_count_mean_field[2]=num_secondary_loops
        if(graph_algo=='KMC'):
          n_chains=np.genfromtxt('NRA.txt')
        else:
          n_chains=p.n_chains # this is the total number of chains taken in the system ( may include dangling ends)
        ####print('n_chains',n_chains)
        loop_fraction_mean_field=loop_order_count_mean_field/n_chains
    ##    plt.plot(np.arange(1,max_loop_order+1), (loop_order_count_mean_field[1:]),'bo-') # check this ceil() function!!- i am assuming this is the correct implementation, but not sure
        plt.plot(np.arange(1,max_loop_order+1), (loop_fraction_mean_field[1:]),'bo-') # check this ceil() function!!- i am assuming this is the correct implementation, but not sure
        # basically what is happwning is that for a specified loop order, it is not necessary that all the nodes along that path actually have that loop order in their own vertex symbol
        # this happens mostly for bigger loop orders
        # and hence, the experssion for calvulating the mean field loop order might not be valid in that case
        ##print((loop_order_count_mean_field))
        ####print('max frequency mean field:',np.argmax(loop_order_count_mean_field))
        f.close()
        plt.xlabel('Loop order')
        plt.ylabel('Frequency')
        plt.savefig('3_mean_field_loop_order_dist_new_'+str(ite),bbox_inches='tight')
        np.savetxt('3_loop_order_mean_field_new'+str(ite)+'.txt',np.array(np.transpose([np.arange(1,max_loop_order+1), (loop_order_count_mean_field[1:]),(loop_fraction_mean_field[1:])])),header='loop_order,count,fraction')

        # probability distribution of vertex symbols
        e = G.number_of_edges()
        v = G.number_of_nodes()
        c=nx.number_connected_components(G)
        print('e',e)
        print('v',v)
        print('c', c)
        cycle_rank = e - v + c
        num_cycles=sum(loop_order_count_mean_field[3:])
        print('number of cycles from  my algo (including primary and secondary)=',num_cycles+num_primary+num_secondary)
        print('number of cycles from  my algo (no primary and secondary)=', num_cycles)
        print('cycle rank +num_primary+num_secondary=',cycle_rank+num_primary+num_secondary)
        print('cycle rank=',cycle_rank)
        print('num_primary',num_primary)
        print('num_secondary',num_secondary)
      except Exception as e:
          print(e)

