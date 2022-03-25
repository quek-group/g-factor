#!/usr/bin/env python

import h5py
import numpy as np
from numpy import linalg as LA
#from progress.bar import Bar
import sys
import math
import argparse

# WFN_kkxky WFN_before WFN_tempout WFN_after
parser = argparse.ArgumentParser()
parser.add_argument('fname_in', help='Input WFN.h5 file kkxky')
parser.add_argument('fname_spin2', help='Input WFN.h5 file 12x12')
parser.add_argument('fname_out', help='Output WFN.h5 for processing f1')
parser.add_argument('fname_out2', help='Output WFN.h5 file with dudk stored in 12x12 mesh')
args = parser.parse_args()

f1 = h5py.File(args.fname_in, 'r')
print(list(f1.keys()))
print(list(f1['mf_header'].keys()))
print(list(f1['wfns'].keys()))
f2 = h5py.File(args.fname_spin2, 'r')
fout = h5py.File(args.fname_out, 'w')
fout2 = h5py.File(args.fname_out2, 'w')
assert f1['mf_header/kpoints/nspin'][()] == 1
print('Copying header...')
fout.copy(f1['mf_header'], 'mf_header')
print('Copying wfns...')
fout.copy(f1['wfns'], 'wfns')
#
print('Copying header...')
fout2.copy(f2['mf_header'], 'mf_header')
print('Copying wfns...')
fout2.copy(f2['wfns'], 'wfns')
#
#band index
ib=13
#number of gvector
ngvec= 8676
unitk = 1.154701*2*3.14159/(5.97153455*0.52918)
overlapxR=0
overlapxI=0
overlapyR=0
overlapyI=0
tempoutK=0
tempoutKx=0
tempoutKy=0
#overlapxR= 0.5517574454194285
#overlapxI= 0.8339817670205144
#overlapyR= 0.33262279736121575 
#overlapyI=-0.9430531702004293
qx=0.00001*unitk*math.sqrt(3)
qy=0.00001*unitk
berrycurvR=0
berrycurvI=0
#fout['wfns/coeffs'][band index,spin index,cg index of k1 k2 k3...,Real&imag component]
#calculate overlap (phase)
for ig in range(0,ngvec):
    tempoutK=tempoutK+fout['wfns/coeffs'][ib,0,ig,0]*fout['wfns/coeffs'][ib,0,ig,0]+fout['wfns/coeffs'][ib,0,ig,1]*fout['wfns/coeffs'][ib,0,ig,1]
    igx=ig+ngvec
    tempoutKx=tempoutKx+fout['wfns/coeffs'][ib,0,igx,0]*fout['wfns/coeffs'][ib,0,igx,0]+fout['wfns/coeffs'][ib,0,igx,1]*fout['wfns/coeffs'][ib,0,igx,1]
    overlapxR=overlapxR+fout['wfns/coeffs'][ib,0,ig,0]*fout['wfns/coeffs'][ib,0,igx,0] +fout['wfns/coeffs'][ib,0,ig,1]*fout['wfns/coeffs'][ib,0,igx,1]
    overlapxI=overlapxI+fout['wfns/coeffs'][ib,0,ig,0]*fout['wfns/coeffs'][ib,0,igx,1] -fout['wfns/coeffs'][ib,0,ig,1]*fout['wfns/coeffs'][ib,0,igx,0]
    igy=ig+ngvec+ngvec
    tempoutKy=tempoutKy+fout['wfns/coeffs'][ib,0,igy,0]*fout['wfns/coeffs'][ib,0,igy,0]+fout['wfns/coeffs'][ib,0,igy,1]*fout['wfns/coeffs'][ib,0,igy,1]
    overlapyR=overlapyR+fout['wfns/coeffs'][ib,0,ig,0]*fout['wfns/coeffs'][ib,0,igy,0] +fout['wfns/coeffs'][ib,0,ig,1]*fout['wfns/coeffs'][ib,0,igy,1]
    overlapyI=overlapyI+fout['wfns/coeffs'][ib,0,ig,0]*fout['wfns/coeffs'][ib,0,igy,1] -fout['wfns/coeffs'][ib,0,ig,1]*fout['wfns/coeffs'][ib,0,igy,0]

print(tempoutK,tempoutKx,tempoutKy)
print(overlapxR,overlapxI)
print(overlapyR,overlapyI)
rhox=math.sqrt(overlapxR*overlapxR+overlapxI*overlapxI)
rhoy=math.sqrt(overlapyR*overlapyR+overlapyI*overlapyI)
#calculate partial u partial k
for ig in range(0,ngvec):
    igx=ig+ngvec
    fout['wfns/coeffs'][ib,0,igx,0]=( overlapxR*f1['wfns/coeffs'][ib,0,igx,0]/rhox+overlapxI*f1['wfns/coeffs'][ib,0,igx,1]/rhox-f1['wfns/coeffs'][ib,0,ig,0])/qx
    fout['wfns/coeffs'][ib,0,igx,1]=(-overlapxI*f1['wfns/coeffs'][ib,0,igx,0]/rhox+overlapxR*f1['wfns/coeffs'][ib,0,igx,1]/rhox-f1['wfns/coeffs'][ib,0,ig,1])/qx
    igy=ig+ngvec+ngvec
    fout['wfns/coeffs'][ib,0,igy,0]=( overlapyR*f1['wfns/coeffs'][ib,0,igy,0]/rhoy+overlapyI*f1['wfns/coeffs'][ib,0,igy,1]/rhoy-f1['wfns/coeffs'][ib,0,ig,0])/qy
    fout['wfns/coeffs'][ib,0,igy,1]=(-overlapyI*f1['wfns/coeffs'][ib,0,igy,0]/rhoy+overlapyR*f1['wfns/coeffs'][ib,0,igy,1]/rhoy-f1['wfns/coeffs'][ib,0,ig,1])/qy
#calculate berry curvature
for ig in range(0,ngvec):
    igx=ig+ngvec
    igy=ig+ngvec+ngvec
    berrycurvR = berrycurvR+fout['wfns/coeffs'][ib,0,igx,0]*fout['wfns/coeffs'][ib,0,igy,1]-fout['wfns/coeffs'][ib,0,igx,1]*fout['wfns/coeffs'][ib,0,igy,0]

berrycurvR =-2*berrycurvR
print('berrycurv')
print(berrycurvR)
# store du/dk in fout, which is copied from f2 (12x12)
ntemp = 0
for ig in range(0,30):
    ntemp=ntemp+f2['mf_header/kpoints/ngk'][ig]
    print(ig)
    print(f2['mf_header/kpoints/ngk'][ig])

overlapxR=0.0
overlapyR=0.0
for ig in range(0,ngvec):
    igx=ig+ngvec
    igy=ig+ngvec+ngvec
    overlapxR = overlapxR + fout['wfns/coeffs'][ib,0,igx,0]*fout['wfns/coeffs'][ib,0,igx,0]+fout['wfns/coeffs'][ib,0,igx,1]*fout['wfns/coeffs'][ib,0,igx,1]
    overlapyR = overlapyR + fout['wfns/coeffs'][ib,0,igy,0]*fout['wfns/coeffs'][ib,0,igy,0]+fout['wfns/coeffs'][ib,0,igy,1]*fout['wfns/coeffs'][ib,0,igy,1]

print('normalization')
print(overlapxR,overlapyR,math.sqrt(overlapxR)*math.sqrt(overlapyR))

for ig in range(0,ngvec):
    igt=ig+ntemp
    igx=ig+ngvec
    fout2['wfns/coeffs'][12,0,igt,0]=fout['wfns/coeffs'][ib,0,igx,0]/math.sqrt(overlapxR)
    fout2['wfns/coeffs'][12,0,igt,1]=fout['wfns/coeffs'][ib,0,igx,1]/math.sqrt(overlapxR)
    igy=ig+ngvec+ngvec
    fout2['wfns/coeffs'][13,0,igt,0]=fout['wfns/coeffs'][ib,0,igy,0]/math.sqrt(overlapyR)
    fout2['wfns/coeffs'][13,0,igt,1]=fout['wfns/coeffs'][ib,0,igy,1]/math.sqrt(overlapyR)

print('Finished!')
