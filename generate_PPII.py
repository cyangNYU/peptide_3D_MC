#! /usr/bin/python2.7

import random
import sys
from pyrosetta import *


def PPII_generate(seq):
    '''
    generate PPII conformation of the short protein sequence
    '''
    init() #load Rosetta database files
    pose = pose_from_sequence(seq,'fa_standard') #fa_standard is residue type of the fasta
    
    #set the (phi, psi) of PPII structure as (-70, 150) with 3 deviation. 
    for i in range(1,pose.total_residue()+1):
        pose.set_phi(i,random.gauss(-70,3))
        pose.set_psi(i,random.gauss(150,3))
        
    return pose
    

def main():
    if len(sys.argv) != 3:
        print 'usage" ./generate_PPII.py {--sequence} sequence'
        sys.exit(1)

    option = sys.argv[1]
    if option == '--sequence':
        sequence = sys.argv[2].upper()
        PPII_generate(sequence).dump_pdb('PPII.pdb')
    else:
        print 'unknown option:' + option
        sys.exit(1)

if __name__ == '__main__':
    main()
        
