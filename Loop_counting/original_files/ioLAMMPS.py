#!/use/local/bin/env python
# -*- coding: utf-8 -*-

# Python script to write configuration for LAMMPS

import math
import numpy as np
import os.path

def readLAMMPS(filename, vflag,frac_weak):

   f1=open(filename,"r")

   line1 = f1.readline()
   line2 = f1.readline()

   line3 = f1.readline()
   line3 = line3.strip()
   n_links = int(line3.split(" ")[0])
 
   line4 = f1.readline()
   line4 = line4.strip()
   atom_types = int(line4.split(" ")[0])

   line5 = f1.readline()
   line5 = line5.strip()
   n_chains = int(line5.split(" ")[0])

   line6 = f1.readline()\
           
   line6 = line6.strip()
   bond_types = int(line6.split(" ")[0])

   links_unsort  = np.zeros((n_links,4))
   links   = np.zeros((n_links,3), dtype = float)
   chains  = np.full((n_chains,4), -1, dtype = int)
   mass    = np.zeros(atom_types, dtype = float)

   line7 = f1.readline()
   line8 = f1.readline()
   line8 = line8.strip()
   xlo = float(line8.split(" ")[0])
   xhi = float(line8.split(" ")[1])

   line9 = f1.readline()
   line9 = line9.strip()
   ylo = float(line9.split(" ")[0])
   yhi = float(line9.split(" ")[1])

   line10 = f1.readline()
   line10 = line10.strip()
   zlo = float(line10.split(" ")[0])
   zhi = float(line10.split(" ")[1])


   for i in range (0, 3):
       f1.readline()
   
   for i in range(0, atom_types):
       line = f1.readline()
       line = line.strip()
       mass[i] = float(line.split(" ")[1])

   f1.close()


   links_unsort = np.genfromtxt(filename, usecols=(0,3,4,5), skip_header=18, max_rows=n_links)

   for i in range(0, n_links):
       index = int(links_unsort[i,0])
       links[index-1,:] = links_unsort[i,1:4]


##   chains[:,0] = N
#cnt,ctype,1,conn1,conn2
   if(vflag==0):
      chains[:,0:4] = np.genfromtxt(filename,usecols=(1,2,3,4), skip_header=17+n_links+3, max_rows=n_chains)
   elif(vflag==1):
      chains[:,0:4] = np.genfromtxt(filename,usecols=(1,2,3,4), skip_header=17+2*n_links+2*3, max_rows=n_chains)
   else:
      print("Invalid Velocity Flag")
##   print(chains)
   directory = './'+str(int(100*frac_weak))+'/'
   filename = 'primary_loops'
   file_path = os.path.join(directory, filename)
   if not os.path.isdir(directory):
      os.mkdir(directory)  
   loop_atoms = np.genfromtxt(file_path, usecols=(1), skip_header=0)
   loop_atoms.tolist() 

   return xlo, xhi, ylo, yhi, zlo, zhi, n_links, n_chains, links, chains, atom_types, bond_types, mass, loop_atoms



def readLAMMPS_into_graph_lattice(G,filename): # reads the file and also makes a networkX graph out of it

   chains= np.genfromtxt(filename,dtype='int',delimiter='\t')
##   print(chains)
##   print('len(chains)',len(chains))
   print(np.shape(chains))
##   stop
   secondary_loop=0
   for i in range(0,len(chains)):
      link_1=chains[i,0]
      link_2=chains[i,1]
      
##      print('link1',link_1,'link2',link_2)
      if(link_1!=link_2):  # pruning primary loops here
##         if(link_1==4):
##            print('link_1=4')
##         if(link_2==4):
##            print('link_2=4')
         G.add_edge(link_1,link_2)
##      else:
##         print('link_1=link_2=',link_1)
##      if(G.has_edge(link_1,link_2)):
##         secondary_loop=secondary_loop+1
##      G.add_edge(link_1,link_2) # no pruning of primary loops required here
##   print('G.neighbor(4)',list(G.neighbors(4)))
   return len(chains)



def readLAMMPS_into_graph(G,Gmult,filename, vflag,frac_weak): # reads the file and also makes a networkX graph out of it

   f1=open(filename,"r")

   line1 = f1.readline()
   line2 = f1.readline()

   line3 = f1.readline()
   line3 = line3.strip()
   n_links = int(line3.split(" ")[0])
 
   line4 = f1.readline()
   line4 = line4.strip()
   atom_types = int(line4.split(" ")[0])

   line5 = f1.readline()
   line5 = line5.strip()
   n_chains = int(line5.split(" ")[0])

   line6 = f1.readline()
           
   line6 = line6.strip()
   bond_types = int(line6.split(" ")[0])

   links_unsort  = np.zeros((n_links,4))
   links   = np.zeros((n_links,3), dtype = float)
   chains  = np.full((n_chains,4), -1, dtype = int)
   mass    = np.zeros(atom_types, dtype = float)

   line7 = f1.readline()
   line8 = f1.readline()
   line8 = line8.strip()
   xlo = float(line8.split(" ")[0])
   xhi = float(line8.split(" ")[1])

   line9 = f1.readline()
   line9 = line9.strip()
   ylo = float(line9.split(" ")[0])
   yhi = float(line9.split(" ")[1])

   line10 = f1.readline()
   line10 = line10.strip()
   zlo = float(line10.split(" ")[0])
   zhi = float(line10.split(" ")[1])


   for i in range (0, 3):
       f1.readline()
   
   for i in range(0, atom_types):
       line = f1.readline()
       line = line.strip()
       mass[i] = float(line.split(" ")[1])

   f1.close()


   links_unsort = np.genfromtxt(filename, usecols=(0,3,4,5), skip_header=18, max_rows=n_links)

   for i in range(0, n_links):
       index = int(links_unsort[i,0])
       links[index-1,:] = links_unsort[i,1:4]


##   chains[:,0] = N
#cnt,ctype,1,conn1,conn2
   if(vflag==0):
      chains[:,0:4] = np.genfromtxt(filename,usecols=(1,2,3,4), skip_header=17+n_links+3, max_rows=n_chains)
   elif(vflag==1):
      chains[:,0:4] = np.genfromtxt(filename,usecols=(1,2,3,4), skip_header=17+2*n_links+2*3, max_rows=n_chains)
   else:
      print("Invalid Velocity Flag")
   print('len(chains)',len(chains))
   for i in range(0,len(chains)):
      link_1=chains[i,2]
      link_2=chains[i,3]
      Gmult.add_edge(link_1,link_2)
      
      if(link_1==link_2):
         loop_atoms.append(link_1)
         G.add_node(link_1) # only add the node, don't add the bond in case of primary loop
      else:
         G.add_edge(link_1,link_2)
##      print('link1',link_1,'link2',link_2)
##      if(link_1!=link_2):  # pruning primary loops here
##         G.add_edge(link_1,link_2)
##      else:
##         print('primary loops detected')
##         stop
##         print('ADDED:','link1',link_1,'link2',link_2)
##      stop

   #counting secondary loops:
   num_secondary_loops=0
   for i in range(0,n_chains):
       for j in range(0,i):
           if (((chains[i,2]==chains[j,2]) and (chains[i,3]==chains[j,3])) or ((chains[i,2]==chains[j,3])and (chains[i,3]==chains[j,2]))):
               num_secondary_loops=num_secondary_loops+1

   
   directory = './'+'network'+'/'
   filename = 'primary_loops'
   file_path = os.path.join(directory, filename)
   if not os.path.isdir(directory):
      os.mkdir(directory)  
   loop_atoms = np.genfromtxt(file_path, usecols=(1), skip_header=0)
   loop_atoms.tolist() 
   for i in loop_atoms: # nodes involved in primary loop
      Gmult.add_edge(i,i)
   return xlo, xhi, ylo, yhi, zlo, zhi, n_links, n_chains, links, chains, atom_types, bond_types, mass, loop_atoms, num_secondary_loops





def readLAMMPS_into_graph_from_bond_only(G,Gmult,filename):#, vflag,frac_weak): # reads the file and also makes a networkX graph out of it


# here file is not in LAMMPS format- so no need to do all these formatting!
   loop_atoms=[]
   chains= np.genfromtxt(filename,dtype='int',delimiter=',',skip_header=1)
   print(chains)
   print('len(chains)',len(chains))
   print(np.shape(chains))
##   stop
   secondary_loop=0
   for i in range(0,len(chains)):
      link_1=chains[i,0]
      link_2=chains[i,1]
      Gmult.add_edge(link_1,link_2)
      
      if(link_1==link_2):
         loop_atoms.append(link_1)
         G.add_node(link_1) # only add the node, don't add the bond in case of primary loop
      else:
         G.add_edge(link_1,link_2)
##         print('primary loop at node ',link_1)
##      print('link1',link_1,'link2',link_2)
##      if(link_1!=link_2):  # pruning primary loops here
##         G.add_edge(link_1,link_2)
##      if(G.has_edge(link_1,link_2)):
##         secondary_loop=secondary_loop+1
##      G.add_edge(link_1,link_2) # no pruning of primary loops required here
      
   return loop_atoms, len(chains)


def writeLAMMPS(filename, xlo, xhi, ylo, yhi, zlo, zhi, links, chains, atom_types, bond_types, mass, loop_atoms):
   
   n_atoms  = len(links[:,0])
   n_bonds  = len(chains[:,0])
   n_loops  = len(loop_atoms)

   # LAMMPS input file begins here
   file1=open(filename,"w")
   file1.write("# Gels-Network -- Initial Configuration\n")
   file1.write("\n")
   file1.write("%i %s\n" % (n_atoms,"atoms"))
   file1.write("%i %s\n" % (atom_types,"atom types"))
   file1.write("%i %s\n" % (n_bonds,"bonds"))
   file1.write("%i %s\n" % (bond_types,"bond types"))
   file1.write("\n")
   file1.write("%7.4f %7.4f %s %s\n" % (xlo,xhi,"xlo","xhi"))
   file1.write("%7.4f %7.4f %s %s\n" % (ylo,yhi,"ylo","yhi"))
   file1.write("%7.4f %7.4f %s %s\n" % (zlo,zhi,"zlo","zhi"))
   
   file1.write("\n")
   file1.write("Masses\n")
   file1.write("\n")
   for i in range(0, atom_types):
       file1.write("%i %6.4f\n" % (i+1,mass[i]))
   
   file1.write("\n")
   file1.write("Atoms\n")
   file1.write("\n")
   cnt = 0
   cnt_loops = 0
   for i in range(0, n_atoms):
      # cnt = cnt + 1
      # file1.write("{:2} {:2} {:2} {:7} {:7} {:7}\n".format(cnt,1,1,links[i,0],links[i,1],links[i,2]))
       if(cnt_loops < n_loops and cnt == loop_atoms[cnt_loops]):
          cnt = cnt + 1
          cnt_loops = cnt_loops + 1
          file1.write("{:2} {:2} {:2} {:7} {:7} {:7}\n".format(cnt,1,2,links[i,0],links[i,1],links[i,2]))
       else:
          cnt = cnt + 1
          file1.write("{:2} {:2} {:2} {:7} {:7} {:7}\n".format(cnt,1,1,links[i,0],links[i,1],links[i,2]))
   
   
   # Bond writing starts here
   file1.write("\n")
   file1.write("Bonds\n")
   file1.write("\n")
   cnt = 0
   for i in range(0, n_bonds):
       cnt = cnt + 1
       #need to change the following lines according to new format!!
##       stop#- ifrst change the following line according to format
       file1.write("{:4} {:4} {:4} {:4}\n".format(cnt,chains[i,0]+1,chains[i,2], chains[i,3]))
   
   
   file1.close()



#def writeLAMMPSafternetgen(filename, xlo, xhi, ylo, yhi, zlo, zhi, n_links, n_chains, n_loops, n_dang, n_free, n_break, links, chains, atom_types, bond_types, mass):
def writeLAMMPSafternetgen(filename, xlo, xhi, ylo, yhi, zlo, zhi, links, chains, atom_types, bond_types, mass, loop_atoms):

#   n_atoms  = n_links
#   n_bonds  = n_chains - n_loops - n_dang - n_free -n_break
 
   n_atoms = len(links[:,0])  
   n_bonds = len(chains[:,0])  
   n_loops = len(loop_atoms)

   print(n_atoms, n_bonds, n_loops)
 
   # LAMMPS input file begins here
   file1=open(filename,"w")
   file1.write("# Gels-Network -- Initial Configuration\n")
   file1.write("\n")
   file1.write("%i %s\n" % (n_atoms,"atoms"))
   file1.write("%i %s\n" % (atom_types,"atom types"))
   file1.write("%i %s\n" % (n_bonds,"bonds"))# equivalent to n_chains
   file1.write("%i %s\n" % (bond_types,"bond types"))
   file1.write("\n")
   file1.write("%7.4f %7.4f %s %s\n" % (xlo,xhi,"xlo","xhi"))
   file1.write("%7.4f %7.4f %s %s\n" % (ylo,yhi,"ylo","yhi"))
   file1.write("%7.4f %7.4f %s %s\n" % (zlo,zhi,"zlo","zhi"))
   
   file1.write("\n")
   file1.write("Masses\n")
   file1.write("\n")
   for i in range(0, atom_types):
       file1.write("%i %6.4f\n" % (i+1,mass[i]))
   
   file1.write("\n")
   file1.write("Atoms\n")
   file1.write("\n")
   cnt = 0
   cnt_loops = 0
   for i in range(0, n_atoms):
       if(cnt_loops < n_loops and cnt == loop_atoms[cnt_loops]): 
          cnt = cnt + 1
          cnt_loops = cnt_loops + 1
          file1.write("{:2} {:2} {:2} {:7} {:7} {:7}\n".format(cnt,1,2,links[i,0],links[i,1],links[i,2]))
       else:
          cnt = cnt + 1
          file1.write("{:2} {:2} {:2} {:7} {:7} {:7}\n".format(cnt,1,1,links[i,0],links[i,1],links[i,2]))
   
  
   print(cnt, cnt_loops) 
   # Bond writing starts here
   file1.write("\n")
   file1.write("Bonds\n")
   file1.write("\n")
   cnt = 0
#   for i in range(0, n_chains):
#         if(chains[i,2]!=-1):
#            if(chains[i,1] != chains[i,2]):
#               cnt = cnt + 1
#               file1.write("{:4} {:4} {:4} {:4}\n".format(cnt,1,chains[i,1], chains[i,2]))
   
   for i in range(0, n_bonds):
       cnt = cnt + 1
       file1.write("{:4} {:4} {:4} {:4} {:4}\n".format(cnt,chains[i,0],1,chains[i,1], chains[i,2]))## cnt,ctype,1,con1,conn2
   
   file1.close()
