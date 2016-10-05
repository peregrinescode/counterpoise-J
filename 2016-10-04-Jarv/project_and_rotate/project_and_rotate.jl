#!/bin/julia
# Forked from https://github.com/jarvist/VASP_LatheOfHeaven/blob/b5f66a49cfd5d43745adc22bb52457b27f0def04/rotate_MA_POSCAR.jl

# What Would Hamilton Do? (WWHD?)
using Quaternions

XYZs = readlines(open("test.xyz"))

# Origin; about which we rotate 
#origin=[0.5 0.5 0.5]' # origin of MA in unit cell, fractional coords

# Read through XYZ file and calculate (atom position average, not mass weighted) centre-of-molecule
sum=[0;0;0]
count=0
for XYZ in XYZs
    r=float(split(XYZ)[2:4]) # "Atom 1.2 3.4 5.6" --> [1.2; 3.4; 5.6]
    sum+=r
    count+=1
end
println(sum)
origin=sum/count
println(STDERR, "Origin: ",origin)

#shift=[5.0*randn() 5.0*randn() 5.0*randn()]'
shift=[0 0 0]' # No shift

# Generate a normalized Quaternion with normally distributed random 4 component
qr=normalize(Quaternion(randn(),randn(),randn(),randn()))
# This should be an evenly distributed rotation matrix
rotate=rotationmatrix(qr)
#rotate=eye(3) #3-dimensional identity matrix; null rotation operation

println(STDERR,"Determinant of 'rotate' matrix (should be 1): ",det(rotate))
println(STDERR,"    Does Q^TQ = I? \n",rotate*rotate')

for XYZ in XYZs
    r=float(split(XYZ)[2:4]) # "Atom 1.2 3.4 5.6" --> [1.2; 3.4; 5.6]
    r=shift+origin+(rotate*(r-origin)) # apply shift + rotate, shift back
    # Left over PBC code for fractional coordinates
    #r=mod(r,1) # Modulo arithmatic so that all {x,y,z} values are on [0,1]
    @printf "%s  %.6f  %.6f  %.6f \n" XYZ[1] r[1] r[2] r[3]
    # 6 digits of precision seems to match VASP output
end

