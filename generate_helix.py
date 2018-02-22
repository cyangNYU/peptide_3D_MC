#! /usr/bin/python2.7

import random
import sys
from pyrosetta import *


def helix_generate(seq):
    '''
    generate helix conformation of the short protein sequence
    '''
    init() #load Rosetta database files
    pose = pose_from_sequence(seq,'fa_standard') #fa_standard is residue type of the fasta
    
    #set the (phi, psi) of helix structure as (-60, -43) with 3 deviation. 
    for i in range(1,pose.total_residue()+1):
        pose.set_phi(i,random.gauss(-60,3))
        pose.set_psi(i,random.gauss(-43,3))
        
    return pose
    

def main():
    if len(sys.argv) != 3:
        print 'usage" ./generate_helix.py {--sequence} sequence'
        sys.exit(1)

    option = sys.argv[1]
    if option == '--sequence':
        sequence = sys.argv[2].upper()
        helix_generate(sequence).dump_pdb('helix.pdb')
    else:
        print 'unknown option:' + option
        sys.exit(1)

if __name__ == '__main__':
    main()
        

