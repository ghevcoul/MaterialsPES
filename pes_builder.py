#!/usr/bin/python
#
# Written by Gavin Heverly-Coulson
# Email: gavin <at> quantumgeranium.com
# 
# Takes an input cell and uses it to build a potential energy surface following
# the system sheared along the a and b lattice vectors. 
# Creates the Quantum ESPRESSO input files for each point on the PES.
#
# 
# This work is licensed under a Simplified BSD License
# Copyright (c) 2014, Gavin Heverly-Coulson
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met: 
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer. 
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution. 
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import sys
import math
from decimal import Decimal
import os

#####################################################
# Calculates the transpose of a matrix
# Works for arbitrary NxM size
def transpose(mat):
  cols = len(mat) # number of rows in mat
  rows = len(mat[0])    # number of columns in mat
  transMat = [x[:] for x in [[None]*cols]*rows] # cols, rows

  for a in range(rows):
    for b in range(cols):
      transMat[a][b] = mat[b][a]

  return transMat
#####################################################
# Calculates the determinant of a 3x3 matrix, using the 2x2 sub-matrices method
def det3(mat):
  return ( ( mat[0][0]*det2([[mat[1][1], mat[1][2]], [mat[2][1], mat[2][2]]]) ) - ( mat[0][1]*det2([[mat[1][0], mat[1][2]], [mat[2][0], mat[2][2]]]) ) + (mat[0][2]*det2([[mat[1][0], mat[1][1]], [mat[2][0], mat[2][1]]])) )
#####################################################
# Calculates the determinant of a 2x2 matrix
def det2(mat):
  return ((mat[0][0]*mat[1][1]) - (mat[0][1]*mat[1][0]))
#####################################################
# a function that takes the cell parameters, in angstrom, and a list of fractional coordinates
# and returns the structure in Cartesian coordinates
#
# Assumes the cellParam matrix is of the form:
# | ax  ay  az |
# | bx  by  bz |
# | cx  cy  cz |
#
def frac2cart(cellParam, fracCoords):
  cartCoords = []
  for i in fracCoords:
    xPos = i[1]*cellParam[0][0] + i[2]*cellParam[1][0] + i[3]*cellParam[2][0]
    yPos = i[1]*cellParam[0][1] + i[2]*cellParam[1][1] + i[3]*cellParam[2][1]
    zPos = i[1]*cellParam[0][2] + i[2]*cellParam[1][2] + i[3]*cellParam[2][2]
    cartCoords.append([i[0], xPos, yPos, zPos])
  return cartCoords
#####################################################
# a function that takes the cell parameters, in angstrom, and a list of Cartesian coordinates
# and returns the structure in fractional coordinates
#
# Uses Cramer's Rule to solve for the fractional coordinates
#
# Assumes the cellParam matrix is of the form:
# | ax  ay  az |
# | bx  by  bz |
# | cx  cy  cz |
#
# Need to use the transpose of this matrix in calculation, so call transpose function first
#
# Assumes cartCoords are of the form:
# | X  x0  y0  z0 |
# | X  x1  y1  z1 |
# | X  x2  y2  z2 |
# | ............. |
# where X is the element symbol
def cart2frac(cellParam, cartCoords):
  latCnt = transpose(cellParam)

  fracCoords = []
  detLatCnt = det3(latCnt)
  for i in cartCoords:
    aPos = (det3([[i[1], latCnt[0][1], latCnt[0][2]], [i[2], latCnt[1][1], latCnt[1][2]], [i[3], latCnt[2][1], latCnt[2][2]]])) / detLatCnt
    bPos = (det3([[latCnt[0][0], i[1], latCnt[0][2]], [latCnt[1][0], i[2], latCnt[1][2]], [latCnt[2][0], i[3], latCnt[2][2]]])) / detLatCnt
    cPos = (det3([[latCnt[0][0], latCnt[0][1], i[1]], [latCnt[1][0], latCnt[1][1], i[2]], [latCnt[2][0], latCnt[2][1], i[3]]])) / detLatCnt
    fracCoords.append([i[0], aPos, bPos, cPos])
  return fracCoords
#####################################################
# Take a list of strings representing the atomic positions in a system and set it up
# for use in further calculations. Returns a list of lists with each sublist formatted
# as: [AtomicSymbol(str), a/x(float), b/y(float), c/z(float)]
def formatCoords(rawCoords):
  coords = []
  for i in rawCoords:
    temp = i.split()
    temp[1] = float(temp[1])
    temp[2] = float(temp[2])
    temp[3] = float(temp[3])
    coords.append(temp)
  return coords
#####################################################
# Take a list of strings for the atomic positions
# [0] - atom type, keep as string
# [1-3] - fractional coordinates, cast to float
# [4] - freeze 0 = freeze, 1 = move, cast to int
def formatCoordsFrozen(rawCoords):
  coords = []
  for i in rawCoords:
    temp = i.split()
    temp[1] = float(temp[1])
    temp[2] = float(temp[2])
    temp[3] = float(temp[3])
    temp[4] = int(temp[4])
    coords.append(temp)
  return coords
#####################################################

pseudodir = "/home/n/nmosey/nmosey/Espresso/PPs/Database"

bohr2angstrom = 1.8897626

# First read in from input file
reader = open(sys.argv[1], 'r')
sourceFile = reader.readlines()
reader.close()

# Find the positions for the pertinent sections of the file
counter = 0
atomNumPos    = -1
coordPos      = -1
latVectPos    = -1
cutoffPos     = -1
stepSizePos   = -1
minAPos       = -1
maxAPos       = -1
minBPos       = -1
maxBPos       = -1
ecutwfcPos    = -1
ecfixedPos    = -1
ecutrhoPos    = -1
qcutzPos      = -1
q2sigmaPos    = -1
londons6Pos   = -1
londonrcutPos = -1
pseudoPos     = -1
for line in sourceFile:
  if "Number of Atoms" in line:
    atomNumPos = counter
  elif "Atomic Positions" in line:
    coordPos = counter + 1
  elif "Lattice Vectors" in line:
    latVectPos = counter + 1
  elif "Cutoff" in line:
    cutoffPos = counter
  elif "Step Size" in line:
    stepSizePos = counter
  elif "Max A" in line:
    maxAPos = counter
  elif "Max B" in line:
    maxBPos = counter
  elif "Min A" in line:
    minAPos = counter
  elif "Min B" in line:
    minBPos = counter
  elif "ecutwfc" in line:
    ecutwfcPos = counter
  elif "ecfixed" in line:
    ecfixedPos = counter
  elif "ecutrho" in line:
    ecutrhoPos = counter
  elif "qcutz" in line:
    qcutzPos = counter
  elif "q2sigma" in line:
    q2sigmaPos = counter
  elif "Pseudopotentials" in line:
    pseudoPos = counter + 1
  elif "london_s6" in line:
    londons6Pos = counter
  elif "london_rcut" in line:
    londonrcutPos = counter
  counter += 1

numAtoms = int(sourceFile[atomNumPos][sourceFile[atomNumPos].index(':')+1:].strip())
stepSize = Decimal(sourceFile[stepSizePos][sourceFile[stepSizePos].index(':')+1:].strip())
maxA     = Decimal(sourceFile[maxAPos][sourceFile[maxAPos].index(':')+1:].strip())
maxB     = Decimal(sourceFile[maxBPos][sourceFile[maxBPos].index(':')+1:].strip())
minA     = Decimal(sourceFile[minAPos][sourceFile[minAPos].index(':')+1:].strip())
minB     = Decimal(sourceFile[minBPos][sourceFile[minBPos].index(':')+1:].strip())

if cutoffPos == -1:
  coords = formatCoordsFrozen(sourceFile[coordPos:coordPos+numAtoms])
else:
  cutoff = float(sourceFile[cutoffPos][sourceFile[cutoffPos].index(':')+1:].strip())
  coords = formatCoords(sourceFile[coordPos:coordPos+numAtoms])
  for i in range(numAtoms):
    if coords[i][3] >= cutoff:
      coords[i].append(1)
    else:
      coords[i].append(0)

latVectsRaw = sourceFile[latVectPos:latVectPos+3] # Note that this is correct, list slice second number is exclusive
latVects = []
temp1 = latVectsRaw[0].split()
latVects.append([float(temp1[0]), float(temp1[1]), float(temp1[2])])
temp2 = latVectsRaw[1].split()
latVects.append([float(temp2[0]), float(temp2[1]), float(temp2[2])])
temp3 = latVectsRaw[2].split()
latVects.append([float(temp3[0]), float(temp3[1]), float(temp3[2])])

# Use a_x component as celldm(1)
celldmAng = latVects[0][0]
celldmBohr = celldmAng * bohr2angstrom

# determine number of atom types in system
atomTypes = []
for i in coords:
  if i[0] not in atomTypes:
    atomTypes.append(i[0])

# Create a string containing all the non-coordinate & lattice vectors part of the input file
fileStart = "&control\n   calculation=\'scf\',\n   restart_mode=\'from_scratch\',\n   tstress=.false.,\n   pseudo_dir=\'{0}\',\n   wf_collect=.false.\n/\n\n".format(pseudodir)
fileStart += "&system\n   ibrav=0,\n   nosym=.true.,\n"
fileStart += "   celldm(1)={0:.6f},  ! Lattice constant of {1:.6f} angstrom in bohr\n".format(celldmBohr, celldmAng)
fileStart += "   nat={0},\n   ntyp={1},\n   {2},\n   {3},\n   {4},\n   {5},\n".format( numAtoms, len(atomTypes), sourceFile[ecutwfcPos].strip(), sourceFile[ecfixedPos].strip(), sourceFile[qcutzPos].strip(), sourceFile[q2sigmaPos].strip() )
if ecutrhoPos >= 0:
  fileStart += "   " + sourceFile[ecutrhoPos].strip() + ",\n"
if londons6Pos >= 0:
  fileStart += "   london=.true.,\n   " + sourceFile[londons6Pos].strip() + ",\n   " + sourceFile[londonrcutPos].strip() + ",\n"
fileStart += "/\n\n&electrons\n   electron_maxstep=50,\n/\n\n&ions\n/\n\n&cell\n/\n\n"
fileStart += "ATOMIC_SPECIES\n"
for i in range(len(atomTypes)):
  fileStart += sourceFile[pseudoPos + i]
fileStart += "ATOMIC_POSITIONS crystal\n"

# Convert input fractional coordinates to Cartesian coordinates
cartCoords = frac2cart(latVects, coords)

lengthA = math.sqrt((latVects[0][0]**2) + (latVects[0][1]**2) + (latVects[0][2]**2))
lengthB = math.sqrt((latVects[1][0]**2) + (latVects[1][1]**2) + (latVects[1][2]**2))
# Normalize the a and b vectors for use in determining how far along x and y directions atoms are moving
# This will be used for tilting c vector by same amount atoms have moved
unitA = [(latVects[0][0]/lengthA), (latVects[0][1]/lengthA), (latVects[0][2]/lengthA)]
unitB = [(latVects[1][0]/lengthB), (latVects[1][1]/lengthB), (latVects[1][2]/lengthB)]

# Do the meat of the operation
xyzName = sys.argv[1] + ".xyz"
xyzWriter = open(xyzName, 'w')
dirsWriter = open("folders", 'w')

numStepsA = math.ceil((maxA-minA)/stepSize) + 1
numStepsB = math.ceil((maxB-minB)/stepSize) + 1

for a in range(int(numStepsA)):
  stepA = minA + (a * stepSize)
  aDir = "a_" + str(stepA)
  for b in range(int(numStepsB)):
    stepB = minB + (b * stepSize)
    bDir = "b_" + str(stepB)
    
    # Determine movement in x and y from stepA and stepB
    stepX = (float(stepA)*lengthA*unitA[0]) + (float(stepB)*lengthB*unitB[0])
    stepY = (float(stepA)*lengthA*unitA[1]) + (float(stepB)*lengthB*unitB[1])
    
    print "stepA={0:f}  stepB={1:f}\nstepX={2:f}  stepY={3:f}\n".format(stepA, stepB, stepX, stepY)
    
    # move atoms
    cartSlide = []
    for i in range(len(cartCoords)):
      if coords[i][4]:
        newX = cartCoords[i][1] + stepX
        newY = cartCoords[i][2] + stepY
        cartSlide.append([cartCoords[i][0], newX, newY, cartCoords[i][3]])
      else:
        cartSlide.append(cartCoords[i])
    
    # Tilt cell in same direction and magnitude atoms moved
    latVectsTilt = [latVects[0], latVects[1], [(latVects[2][0] + stepX), (latVects[2][1] + stepY), latVects[2][2]]]
    
    # Convert slide coordinates to cartesian, using tilted lattice vectors
    slide = cart2frac(latVectsTilt, cartSlide)
    
    # Build CELL_PARAMETERS section using longest component from above
    cellParam = []
    for i in latVectsTilt:
      temp0 = i[0] / celldmAng
      temp1 = i[1] / celldmAng
      temp2 = i[2] / celldmAng
      cellParam.append([temp0, temp1, temp2])
    
    # print output
    # Write to the xyz movie file
    xyzWriter.write( "{0}\ndelta a={1:f}, delta b={2:f}\n".format(numAtoms, stepA, stepB) )
    for q in cartSlide:
      xyzWriter.write( "{0[0]:<2}  {0[1]:> 10.6f}  {0[2]:> 10.6f}  {0[3]:> 10.6f}\n".format(q) )
    
    # Write the PWSCF input file
    path = "{0}/{1}".format(aDir, bDir)
    os.makedirs(path)
    dirsWriter.write(path)
    dirsWriter.write("\n")
    pwscfFile = "{0}/{1}.in".format(path, sys.argv[1])
    pwscfWriter = open(pwscfFile, 'w')
    pwscfWriter.write(fileStart)
    for i in slide:
      pwscfWriter.write( "{0[0]:<2}  {0[1]:> 10.6f}  {0[2]:> 10.6f}  {0[3]:> 10.6f}\n".format(i) )
    pwscfWriter.write("CELL_PARAMETERS\n")
    pwscfWriter.write(" {0[0][0]:> 10.6f}   {0[0][1]:> 10.6f}   {0[0][2]:>10.6f}\n {0[1][0]:> 10.6f}   {0[1][1]:> 10.6f}   {0[1][2]:> 10.6f}\n {0[2][0]:> 10.6f}   {0[2][1]:> 10.6f}   {0[2][2]:> 10.6f}\n".format(cellParam))
    pwscfWriter.close()
    

xyzWriter.close()
dirsWriter.close()

#eof
