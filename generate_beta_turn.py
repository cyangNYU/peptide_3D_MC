#! /usr/bin/python2.7

import random
import sys
from pyrosetta import *


def beta_turn_generate(seq):
    '''
    generate beta turn conformation of the short protein sequence
    '''
    init() #load Rosetta database files
    pose = pose_from_sequence(seq,'fa_standard') #fa_standard is residue type of the fasta
    
    #set the (phi, psi) of beta structure as (-150, 150) with 3 deviation. 
    for i in range(1,pose.total_residue()+1):
        pose.set_phi(i,random.gauss(-150,3))
        pose.set_psi(i,random.gauss(150,3))
    #set the (phi, psi) of 2 residues in turn structure as (-60,120), (80,7)
    t1 = len(seq)/2 # first residue number in turn
    t2 = len(seq)/2 +1 # second residue number in turn
    pose.set_phi(t1,random.gauss(-60,3))
    pose.set_psi(t1,random.gauss(120,3))
    pose.set_phi(t2,random.gauss(80,3))
    pose.set_psi(t2,random.gauss(7,3))
        
    return pose
    

def main():
    if len(sys.argv) != 3:
        print 'usage" ./generate_beta_turn.py {--sequence} sequence'
        sys.exit(1)

    option = sys.argv[1]
    if option == '--sequence':
        sequence = sys.argv[2].upper()
        beta_turn_generate(sequence).dump_pdb('beta_turn.pdb')
    else:
        print 'unknown option:' + option
        sys.exit(1)

if __name__ == '__main__':
    main()
