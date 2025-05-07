#!/use/local/bin/env python
# -*- coding: utf-8 -*-
##
##-------------------------------------------------
## Fast Inertial Relaxation Engine (FIRE) Optimizer
## Ref: Bitzek et al, PRL, 97, 170201 (2006)
##
## Author: Akash Arora
## Implementation is inspired from LAMMPS and ASE Master Code
##-------------------------------------------------

import math
import time
import random
import numpy as np
import scipy.optimize as opt
from numpy import linalg as LA
from scipy.optimize import fsolve
import os.path


class Optimizer(object):

    def __init__(self, atoms, bonds, xlo, xhi, ylo, yhi, zlo, zhi, r0, parameters, ftype):

        self.atoms = atoms
        self.bonds = bonds
        self.xlo = xlo
        self.xhi = xhi
        self.ylo = ylo
        self.yhi = yhi
        self.zlo = zlo
        self.zhi = zhi
        self.parameters = parameters
        self.r0 = r0          
##        self.N = N          
        self.ftype = ftype

    def bondlengths(self):
     
        atoms = self.atoms
        bonds = self.bonds
        Lx = self.xhi - self.xlo
        Ly = self.yhi - self.ylo
        Lz = self.zhi - self.zlo
        n_atoms = len(self.atoms[:,0])
        n_bonds = len(self.bonds[:,0])

        dist = np.zeros((n_bonds,4), dtype=float)

        for i in range (0, n_bonds):
          
              lnk_1 = bonds[i,2]-1
              lnk_2 = bonds[i,3]-1
              delr = atoms[lnk_1,:] - atoms[lnk_2,:]
              
              delr[0] = delr[0] - int(round(delr[0]/Lx))*Lx
              delr[1] = delr[1] - int(round(delr[1]/Ly))*Ly
              delr[2] = delr[2] - int(round(delr[2]/Lz))*Lz
                   
              dist[i,0:3] = delr
              dist[i,3] = LA.norm(delr)
    
        return dist

    
    def invlangevin(self, x):
        return x*(2.99942 - 2.57332*x + 0.654805*x**2)/(1-0.894936*x - 0.105064*x**2)

    def kuhn_stretch(self, lam, E_b):
       
        def func(x, lam, E_b):
            y = lam/x
            beta = self.invlangevin(y)
            return E_b*np.log(x) - lam*beta/x
   
        if lam == 0:
           return 1
        else:
           lam_b = opt.root_scalar(func,args=(lam, E_b),bracket=[lam,lam+1],x0=lam+0.05)
           return lam_b.root

    def get_bondforce(self, r,i):
        bonds=self.bonds
        ctype=bonds[i,0]
        parameters=self.parameters
        N=parameters[ctype,0]#bonds[i,0] gives the ctype
        b=parameters[ctype, 1]
        K=parameters[ctype, 2]
        
##        print(parameters)
##        stop
        fit_param=parameters[ctype, 3]
        E_b=parameters[ctype, 4]
##        K  = self.K
        r0 = self.r0
##        Nb = self.N # b = 1 (lenght scale of the system)
        
##        E_b = 1200
 
        x = (r-r0)/(N*b)
##        print('x ',x) 
        if(x<0.90):
##           print('get_bondforce, case 1')
           lam_b = 1.0
           fbkT  = self.invlangevin(x)
           fbond = -K*fbkT/r
        elif(x<1.4):
##           print('get_bondforce, case 2')
           lam_b = self.kuhn_stretch(x, E_b)
           fbkT  = self.invlangevin(x/lam_b)/lam_b
           fbond = -K*fbkT/r
        else:
##           print('get_bondforce, case 3')
           lam_b = x + 0.05
           fbkT  = 325 + 400*(x-1.4)            
           fbond = -K*fbkT/r
 
        return fbond, lam_b  
          

    def get_force(self):
       
##        N = self.N
##        E_b = 1200
        atoms = self.atoms
        bonds = self.bonds
        ftype = self.ftype
        Lx = self.xhi - self.xlo
        Ly = self.yhi - self.ylo
        Lz = self.zhi - self.zlo
        n_atoms = len(atoms[:,0])
        n_bonds = len(bonds[:,0])
       
        e = 0.0 
        Gamma = 0.0
        f =  np.zeros((n_atoms,3), dtype = float)
        for i in range(0, n_bonds):
            ctype=bonds[i,0]
            parameters=self.parameters
            N=parameters[ctype,0]#bonds[i,0] gives the ctype
            b=parameters[ctype, 1]
            K=parameters[ctype, 2]
            fit_param=parameters[ctype, 3]
            E_b=parameters[ctype, 4]
            
            lnk_1 = bonds[i,2]-1
            lnk_2 = bonds[i,3]-1
            delr = atoms[lnk_1,:] - atoms[lnk_2,:]
##            print('lnk_1',lnk_1,'lnk_2',lnk_2)
##            print('delr ',delr, '\n atoms[lnk_1,:]',atoms[lnk_1,:], '\n  atoms[lnk_2,:]', atoms[lnk_2,:])
            delr[0] = delr[0] - int(round(delr[0]/Lx))*Lx
            delr[1] = delr[1] - int(round(delr[1]/Ly))*Ly
            delr[2] = delr[2] - int(round(delr[2]/Lz))*Lz
                 
            r = LA.norm(delr)
            if (r > 0): 
               [fbond, lam_b] = self.get_bondforce(r,i) 
               lam = (r-self.r0)/N
               beta = -fbond*r/K*lam_b # fbond*r/K is fbKT
               e_bond = N*0.5*E_b*math.log(lam_b)**2
##               print('fbond',fbond,"  ",'r',r)
               
               e_stretch = N*( (lam/lam_b)*beta + math.log(beta/math.sinh(beta)))
               e = e + e_bond + e_stretch
               
            else:
               fbond = 0.0
               e = e + 0.0
##            stop 
            Gamma = Gamma + r*r
       
            # apply force to each of 2 atoms        
            if (lnk_1 < n_atoms):
               f[lnk_1,0] = f[lnk_1,0] + delr[0]*fbond
               f[lnk_1,1] = f[lnk_1,1] + delr[1]*fbond
               f[lnk_1,2] = f[lnk_1,2] + delr[2]*fbond
        
            if (lnk_2 < n_atoms):
               f[lnk_2,0] = f[lnk_2,0] - delr[0]*fbond
               f[lnk_2,1] = f[lnk_2,1] - delr[1]*fbond
               f[lnk_2,2] = f[lnk_2,2] - delr[2]*fbond
        
        return f, e, Gamma
  
 
    def fire_iterate(self, ftol, maxiter, write_itr, logfilename):
      
        tstart = time.time()

        ## Optimization parameters:
        eps_energy = 1.0e-8
        delaystep = 5
        dt_grow = 1.1
        dt_shrink = 0.5
        alpha0 = 0.1
        alpha_shrink = 0.99
        tmax = 10.0
        maxmove = 0.1
        last_negative = 0

        dt = 0.005
        dtmax = dt*tmax
        alpha = alpha0
        last_negative = 0       
 
        Lx = self.xhi - self.xlo
        Ly = self.yhi - self.ylo
        Lz = self.zhi - self.zlo
        n_atoms = len(self.atoms[:,0])
        n_bonds = len(self.bonds[:,0])
        v = np.zeros((n_atoms,3), dtype = float)

        n_bonds = len(self.bonds)
        dist = np.zeros((n_bonds,4), dtype=float)


        [f,e,Gamma] = self.get_force()
        dist = self.bondlengths()

 
        fmaxitr = np.max(np.max(np.absolute(f)))
        fnormitr = math.sqrt(np.vdot(f,f))
##        logfile = open(logfilename,'w') 
##        logfile.write('FIRE: iter  Energy  fmax  fnorm  avg(r)/Nb  max(r)/Nb\n')
##        logfile.write('%s: %5d  %9.6f  %9.6f  %9.6f  %9.4f  %9.4f\n' %
##                              ('FIRE', 0, e, fmaxitr, fnormitr, np.mean(dist[:,3])/self.N, np.max(dist[:,3])/self.N))
##        logfile.flush()
        print('FIRE: iter  Energy  fmax  fnorm  ')
        print('%s: %5d  %9.6f  %9.6f  %9.6f' %
                              ('FIRE', 0, e, fmaxitr, fnormitr))#, avg(r)/Nb  max(r)/Nb np.mean(dist[:,3])/self.N, np.max(dist[:,3])/self.N))

        for itr in range (0, maxiter):
         
          vdotf = np.vdot(v,f) 
          if (vdotf > 0.0):
             vdotv = np.vdot(v,v)
             fdotf = np.vdot(f,f) 
             scale1 = 1.0 - alpha
             if (fdotf == 0.0): scale2 = 0.0
             else: scale2 = alpha * math.sqrt(vdotv/fdotf)
             v = scale1*v + scale2*f
              
             if (itr - last_negative > delaystep):
                 dt = min(dt*dt_grow,dtmax)
                 alpha = alpha*alpha_shrink
      
          else:
             last_negative = itr
             dt = dt*dt_shrink
             alpha = alpha0
             v[:] = v[:]*0.0
      
          v = v + dt*f 
          dr = dt*v
          normdr = np.sqrt(np.vdot(dr, dr))
          if (normdr > maxmove):
              dr = maxmove * dr / normdr

          self.atoms = self.atoms + dr
          for i in range(0, n_atoms):
              self.atoms[i,0] = self.atoms[i,0] - math.floor((self.atoms[i,0]-self.xlo)/Lx)*Lx
              self.atoms[i,1] = self.atoms[i,1] - math.floor((self.atoms[i,1]-self.ylo)/Ly)*Ly
              self.atoms[i,2] = self.atoms[i,2] - math.floor((self.atoms[i,2]-self.zlo)/Lz)*Lz
          
  
          [f,e,Gamma] = self.get_force()
          fmaxitr = np.max(np.max(np.absolute(f)))
          fnormitr = math.sqrt(np.vdot(f,f))


          if((itr+1)%write_itr==0):
             dist = self.bondlengths()
##             logfile.write('%s: %5d  %9.6f  %9.6f  %9.6f  %9.4f  %9.4f\n' %
##                                  ('FIRE', itr+1, e, fmaxitr, fnormitr, np.mean(dist[:,3])/self.N, np.max(dist[:,3])/self.N))
##             logfile.flush()

             # Print on screen
             print('%s: %5d  %9.6f  %9.6f  %9.6f' %
                               ('FIRE', itr+1,  e, fmaxitr, fnormitr))#,  np.mean(dist[:,3])/self.N, np.max(dist[:,3])/self.N))
        
   
          # Checking for convergence
          if (fnormitr < ftol):
             dist = self.bondlengths()
             tend = time.time()
##             logfile.write('%s: %5d  %9.6f  %9.6f  %9.6f  %9.4f  %9.4f\n' %
##                                  ('FIRE', itr+1, e, fmaxitr, fnormitr, np.mean(dist[:,3])/self.N, np.max(dist[:,3])/self.N))
##             logfile.flush()
             print('%s: %5d  %9.6f  %9.6f  %9.6f' %
                               ('FIRE', itr+1,  e, fmaxitr, fnormitr))#, np.mean(dist[:,3])/self.N, np.max(dist[:,3])/self.N))
             print('Iterations converged, Time taken: %7.4f' %(tend-tstart))
             break
          elif (itr == maxiter-1):
             print('Maximum iterations reached')
     

##        logfile.close() 
        
        return e, Gamma       
                

    def compute_pressure(self):

##        K = self.K
        bonds=self.bonds
##        ctype=bonds[i,0]
##        parameters=self.parameters
##        N=parameters[ctype,0]#bonds[i,0] gives the ctype
##        b=parameters[ctype, 1]
##        K=parameters[ctype, 2]
##        fit_param=parameters[ctype, 3]
##        E_b=parameters[ctype,4]
        r0 = self.r0
        ftype = self.ftype
        Lx = self.xhi - self.xlo
        Ly = self.yhi - self.ylo
        Lz = self.zhi - self.zlo
        atoms = self.atoms
        bonds = self.bonds
        n_atoms = len(atoms[:,0])
        n_bonds = len(bonds[:,0])
       
        pxx = pyy = pzz = pxy = pyz = pzx = 0.0
        sigma = np.zeros((n_atoms,6), dtype=float)
        inv_volume = 1.0/(Lx*Ly*Lz)
        for i in range(0, n_bonds):

            lnk_1 = bonds[i,2]-1
            lnk_2 = bonds[i,3]-1
            delr = atoms[lnk_1,:] - atoms[lnk_2,:]
            
            delr[0] = delr[0] - int(round(delr[0]/Lx))*Lx
            delr[1] = delr[1] - int(round(delr[1]/Ly))*Ly
            delr[2] = delr[2] - int(round(delr[2]/Lz))*Lz
                 
            r = LA.norm(delr)
            if (r > 0.0):
               if(ftype=='Mao'): [fbond, lam_b] = self.get_bondforce(r,i)
               else: fbond = self.get_bondforce(r) 
            else: fbond = 0.0
            
            # apply pressure to each of the 2 atoms   
            # And for each of the 6 components     
            if (lnk_1 < n_atoms):
               sigma[lnk_1,0] = sigma[lnk_1,0] + 0.5*delr[0]*delr[0]*fbond
               sigma[lnk_1,1] = sigma[lnk_1,1] + 0.5*delr[1]*delr[1]*fbond
               sigma[lnk_1,2] = sigma[lnk_1,2] + 0.5*delr[2]*delr[2]*fbond
               sigma[lnk_1,3] = sigma[lnk_1,3] + 0.5*delr[0]*delr[1]*fbond
               sigma[lnk_1,4] = sigma[lnk_1,4] + 0.5*delr[1]*delr[2]*fbond
               sigma[lnk_1,5] = sigma[lnk_1,5] + 0.5*delr[2]*delr[0]*fbond
        
            if (lnk_2 < n_atoms):
               sigma[lnk_2,0] = sigma[lnk_2,0] + 0.5*delr[0]*delr[0]*fbond
               sigma[lnk_2,1] = sigma[lnk_2,1] + 0.5*delr[1]*delr[1]*fbond
               sigma[lnk_2,2] = sigma[lnk_2,2] + 0.5*delr[2]*delr[2]*fbond
               sigma[lnk_2,3] = sigma[lnk_2,3] + 0.5*delr[0]*delr[1]*fbond
               sigma[lnk_2,4] = sigma[lnk_2,4] + 0.5*delr[1]*delr[2]*fbond
               sigma[lnk_2,5] = sigma[lnk_2,5] + 0.5*delr[2]*delr[0]*fbond


        pxx = np.sum(sigma[:,0])*inv_volume
        pyy = np.sum(sigma[:,1])*inv_volume
        pzz = np.sum(sigma[:,2])*inv_volume
        pxy = np.sum(sigma[:,3])*inv_volume
        pyz = np.sum(sigma[:,4])*inv_volume
        pzx = np.sum(sigma[:,5])*inv_volume

        return pxx, pyy, pzz, pxy, pyz, pzx


    def change_box(self, scale_x, scale_y, scale_z):

        xlo = self.xlo
        xhi = self.xhi
        ylo = self.ylo
        yhi = self.yhi
        zlo = self.zlo
        zhi = self.zhi
        atoms = self.atoms
        bonds = self.bonds
        n_atoms = len(atoms[:,0])
        n_bonds = len(bonds[:,0])

        xmid = (xlo+xhi)/2  
        ymid = (ylo+yhi)/2  
        zmid = (zlo+zhi)/2  

        new_xlo = xmid + scale_x*(xlo-xmid)
        new_ylo = ymid + scale_y*(ylo-ymid)
        new_zlo = zmid + scale_z*(zlo-zmid)

        new_xhi = xmid + scale_x*(xhi-xmid)
        new_yhi = ymid + scale_y*(yhi-ymid)
        new_zhi = zmid + scale_z*(zhi-zmid)
        
        newLx = new_xhi - new_xlo
        newLy = new_yhi - new_ylo
        newLz = new_zhi - new_zlo
        for i in range(0, n_atoms):            
            atoms[i,0] = xmid + scale_x*(atoms[i,0]-xmid)
            atoms[i,1] = ymid + scale_y*(atoms[i,1]-ymid)
            atoms[i,2] = zmid + scale_z*(atoms[i,2]-zmid)

        self.atoms = atoms
        self.xlo = new_xlo
        self.xhi = new_xhi
        self.ylo = new_ylo
        self.yhi = new_yhi
        self.zlo = new_zlo
        self.zhi = new_zhi





    def KMCbondbreak(self, tau, delta_t, pflag, index,frac_weak):
    
        # Material parameters:
        # beta = 1.0 -- All material params, U0 and sigma, are in units of kT. 
        # Main array: Bonds_register = [Activity index, type, index, link1, link2, dist, rate(ri)]
        # All are active at the start (active = 1, break = 0)
   
        def get_link_bonds(link, bonds_register):
        
            conn = {}
            a1 = np.where(bonds_register[:,3]==link)
            a2 = np.where(bonds_register[:,4]==link)
            a = np.concatenate((a1[0],a2[0]))
            a = np.unique(a)
            for i in range(0,len(a)):
                if(bonds_register[a[i],0]==1): 
                  conn.update({a[i] : bonds_register[a[i],5]})
           
            conn = dict(sorted(conn.items(), key=lambda x: x[1]))     

            return conn

        
        ftype = self.ftype
        n_bonds = len(self.bonds[:,0])
        bonds_register = np.zeros((n_bonds,7))
        bonds_register[:,0] = 1   
        bonds_register[:,1:5] = self.bonds
        dist = self.bondlengths()
        bonds_register[:,5] = dist[:,3]   

        # File to write bond broken stats
        if(index%10==0):
           directory = './'+str(int(100*frac_weak))+'/'
           filename = 'bondbroken_%d.txt'%(index)
           file_path = os.path.join(directory, filename)
           if not os.path.isdir(directory):
               os.mkdir(directory)  
           f2 = open(file_path,'w')
           f2.write('#type, atom1, atom2, length, rate(v), t, t_KMC, vmax, active bonds\n') 
       
        # Write probability values in a file (at every KMC call)
        if(pflag==1):
           prob_file = 'prob_%d.txt' %(index)
           directory = './'+str(int(100*frac_weak))+'/'
           filename = prob_file
           file_path = os.path.join(directory, filename)
           if not os.path.isdir(directory):
               os.mkdir(directory)
           fl1 = open(file_path,'w')   
 
        for i in range (0, n_bonds):
            bonds=self.bonds
            ctype=bonds[i,0]
            parameters=self.parameters
            N=parameters[ctype,0]#bonds[i,0] gives the ctype
            b=parameters[ctype, 1]
            K=parameters[ctype, 2]
            fit_param=parameters[ctype, 3]
            E_b=parameters[ctype,4]
            U0_kT=parameters[ctype,5]
            Nb = N*b
            r = bonds_register[i,5]
            if(r > 0):
              [fbond, lam_b] = self.get_bondforce(r,i)
            else: fbond = 0.0

##            fit_param = 1
            fbkT = -fbond*r/K
            bonds_register[i,6] = math.exp(-U0_kT + fbkT*fit_param)
            if(pflag==1): fl1.write('%i %i %i %i %i %6.4f %6.4f\n' %(bonds_register[i,0], 
                               bonds_register[i,1], bonds_register[i,2], bonds_register[i,3], 
                               bonds_register[i,4], bonds_register[i,5], bonds_register[i,6]))
    
        if(pflag==1): fl1.close()
     
        active_bonds = np.where(bonds_register[:,0]==1)
        n_bonds_init = len(active_bonds[0])
        vmax = max(bonds_register[active_bonds[0],6]) #same as rmax in paper
        if(vmax == 0): vmax = 1e-12  
        # if fbkT = 0, vmax = exp(-56). This number below the machine precison.
        # hence, we assign a small detectable number, vmax = 10^{-12}. 
        # Essentially, it implies that bond breaking rate is very low, or 
        # t = 1/(vmax*nbonds) is very high compare to del_t and hence it will not 
        # enter the KMC bond breaking loop 
           
        t = 1/(vmax*len(active_bonds[0])) 
        print('KMC statistics:') 
        print('Max rate, Active bonds, and t_KMC = %6.4E, %5d, %6.4E'%(vmax, len(active_bonds[0]), t))
        if(t < delta_t):
           t = 0
           while(t < delta_t):
                
                t_KMC    = 1/(vmax*len(active_bonds[0])) 
                vmax     = max(bonds_register[active_bonds[0],6])
##                random.seed(10)
                bond_index    = random.randint(0, len(active_bonds[0])-1)
                pot_bond = active_bonds[0][bond_index]
##                random.seed(10)
                rnd_num  = random.uniform(0,1)
                if((bonds_register[pot_bond,6]/vmax) > rnd_num): # bonds_register[pot_bond,6] is chain scission rate r
                   bonds_register[pot_bond,0] = 0   # Bond is broken!
                   t = t + t_KMC
                   if(index%10==0):
                      f2.write('%5d  %5d  %5d  %0.4E  %0.4E  %0.4E  %0.4E  %0.4E  %5d\n'%(bonds_register[pot_bond,2], bonds_register[pot_bond,3], 
                        bonds_register[pot_bond,4], bonds_register[pot_bond,5], bonds_register[pot_bond,6], 
                        t, t_KMC, vmax, len(active_bonds[0])) )
                      f2.flush()
                   # Local Relaxation -- If the bond-broken created a dangling end system
                   # then make the force on the remaining fourth bond
                   link_1 = bonds_register[pot_bond,3]
                   conn = get_link_bonds(link_1, bonds_register)
                   if(len(conn)==3): 
                      if(conn[list(conn)[0]]==0 and conn[list(conn)[1]]==0):
                         bonds_register[list(conn)[2],6]=0
    
                   elif(len(conn)==2):
                      if(conn[list(conn)[0]]==0):
                         bonds_register[list(conn)[1],6]=0

                   else:
                      bonds_register[list(conn)[0],6]=0


                   link_2 = bonds_register[pot_bond,4]
                   conn = get_link_bonds(link_2, bonds_register)
                   if(len(conn)==3): 
                      if(conn[list(conn)[0]]==0 and conn[list(conn)[1]]==0):
                         bonds_register[list(conn)[2],6]=0
    
                   elif(len(conn)==2):
                      if(conn[list(conn)[0]]==0):
                         bonds_register[list(conn)[1],6]=0

                   else:
                      bonds_register[list(conn)[0],6]=0
                       

                else: 
                   t = t + t_KMC

                active_bonds = np.where(bonds_register[:,0]==1)


        if(index%10==0): f2.close()
    
        n_bonds_final = len(active_bonds[0])
        if(n_bonds_final < n_bonds_init):
           bonds_final = np.zeros((n_bonds_final, 4), dtype = int)
           bonds_final[:,0:4] = bonds_register[active_bonds[0],1:5].astype(int)
           self.bonds = bonds_final
 
        print('time, init bonds, final bonds = %6.4E, %5d, %5d'%(t, n_bonds_init, n_bonds_final))
        print('---------------------------------------------------------------')
   
 
        return t, n_bonds_init, n_bonds_final
 
