#!/usr/bin/python2.7

'''
Using a simple Monte Carlo algorithm that spontaneously folds an start structure, the score funtion can be set up.
'''
import numpy as np
import pandas as pd
from pandas import Series, DataFrame
from pyrosetta import *
 
init() # load Rosetta database files

from pyrosetta.teaching import *

pose = pose_from_pdb('MC_start.pdb') # the input starting structure file

#set up score function

scorefxn = ScoreFunction()
scorefxn.set_weight(hbond_lr_bb,1.0) # set long range H-bond backbone weight
scorefxn.set_weight(hbond_sr_bb,1.0) # set short range H-bond backbone weight
scorefxn.set_weight(fa_elec, 0.5) # set distance-dependent dielectric electrostatics weight
scorefxn.set_weight(vdw,0.5) # set van der Waals weight


#set up simulation parameters
ncycles = 100000
kT = 1.0

#define the MC protocol
def MC_mover(E1, E2):
    if E2 <= E1:
        return 1
    else:
        p_control = np.random.random()
        delta_E = E1 - E2
        p = np.exp(delta_E/kT)
        if p <= p_control:
            return 1
        else:
            return 0
        
#collect energy information
s_energy = Series()
s_energy['0'] = scorefxn(pose)
lowest_energy = Series()
lowest_energy['0'] = scorefxn(pose)
a = []  
a.append(scorefxn(pose))
List = ['score']
for i in range(1,pose.total_residue()+1):
    List.append('%02.0f_phi'%i)
    List.append('%02.0f_psi'%i)
    
s_total = pd.DataFrame(index=range(1,ncycles+1), columns=tuple(List)) # collect cycle, energy and (phi,psi) 

#run the MC simulation
for i in range(1,ncycles+1):
    E1 = scorefxn(pose)
    #choose one residue at random 
    R_random = np.random.randint(1,pose.total_residue()+1)
    #choose the phi or psi at random
    new_pose = Pose()
    new_pose.assign(pose)
    Current_phi = pose.phi(R_random)
    Current_psi = pose.psi(R_random)
    if np.random.choice(['phi','psi']) == 'phi':        
        #random.gauss distribution of phi with std of 20
        
        new_pose.set_phi(R_random, np.random.normal(Current_phi, 20))
    else:
        
        new_pose.set_psi(R_random, np.random.normal(Current_psi, 20))

    E2 = scorefxn(new_pose)
    a.append(E2)
    s_energy['%d'%i] = E2
    s_total['score'][i]=E2
    for k in range(1,pose.total_residue()+1):
        s_total['%02.0f_phi'%k][i] = new_pose.phi(k)
        s_total['%02.0f_psi'%k][i] = new_pose.psi(k)


    lowest_energy['%d'%i] = min(a)
    if MC_mover(E1, E2) == 1:
        pose = new_pose
    else:
        del new_pose

#write down final MC conformation
pose.dump_pdb('MC_final.pdb')
s_energy.to_csv('energy.csv')
s_total.to_csv('summary_total.csv')
lowest_energy.to_csv('lowest_energy.csv')
del a

#write lowest energy conformation
b = s_total.sort_values(by='score').head(1)
for i in range(1,pose.total_residue()+1):
    [phi] = b['%02.0f_phi'%i].values.tolist()
    [psi] = b['%02.0f_psi'%i].values.tolist()
    pose.set_phi(i,phi)
    pose.set_psi(i,psi)

pose.dump_pdb('lowest.pdb')
    

    
        
