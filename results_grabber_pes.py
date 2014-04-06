#!/usr/bin/python
# 
# Written by Gavin Heverly-Coulson
# Email: gavin <at> quantumgeranium.com
# 
# Reads the energies from the Quantum Espresso output files
# that result from running the jobs built by pes_builder.py
# Prints a file containing the surface energies (in J/m^2) for each
# point on the PES. Also writes a gnuplot input file to make
# an image to visualize the PES.
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

# Find a line containing a particular keyword in a the file list
# created from the reader.readlines() function.
# Returns -1 if the line is not found.
def findLine(lineList, keyword):
  counter = 0
  lineNum = -1
  while (counter < len(lineList)):
    if (keyword in lineList[counter]):
      lineNum = counter
      break
    counter += 1
  return lineNum

filename = sys.argv[1]

# Open the list of folders and read it in
reader = open("folders", 'r')
files = reader.readlines()
reader.close()

# Strip the whitespace/newlines from the folderlist
foo = 0
while foo < len(files):
  files[foo] = files[foo].strip()
  foo += 1

outFilename = filename + ".out"
resultsFilename = filename + "_results.dat"

# Write all energies to this string, will write it to file at end
resultsString = "b\\a"

# Count the number of points on PES in a and b directions
aVals = []
bVals = []
for i in files:
  temp = i.split('/')
  if temp[0][2:] not in aVals:
    aVals.append(temp[0][2:])
  if temp[1][2:] not in bVals:
    bVals.append(temp[1][2:])

#print "A Values:"
#for i in aVals:
#  print i
#print "\nB Values:"
#for i in bVals:
#  print i

for i in aVals:
  resultsString += "   " + i
resultsString += "\n"

# Build the empty matrix for the energies
# Will be indexed as: energies[b][a]
energies = [[None for i in range(len(aVals))] for j in range(len(bVals))]

# Calculate the cross-sectional area (area of AB plane)
reader = open("{0}/{1}.in".format(files[0], filename), 'r')
inFile = reader.readlines()
reader.close()

celldmLine = findLine(inFile, "celldm(1)")
celldmTemp = inFile[celldmLine].split(',')[0]
celldm = float(celldmTemp.split('=')[1]) * 5.2917721e-11 # Store in metre
latVectsLine = findLine(inFile, "CELL_")
aLength = (celldm * float(inFile[latVectsLine+1].split()[0]))
b_x = (celldm * float(inFile[latVectsLine+2].split()[0]))
b_y = (celldm * float(inFile[latVectsLine+2].split()[1]))
bLength = math.sqrt( (b_x**2) + (b_y**2) )
xsArea = aLength * bLength

for a in range(len(aVals)):
  for b in range(len(bVals)):
    print "delta A = {0}, delta B = {1}".format(aVals[a], bVals[b])
    
    filePath = "a_{0}/b_{1}/{2}".format(aVals[a], bVals[b], outFilename)
    reader = open(filePath, 'r')
    outFile = reader.readlines()
    reader.close()
  
    # check if calculation converged
    converged = findLine(outFile, "convergence NOT achieved")
    if (converged > -1):
      energies[b][a] = "XXXX"
      print "   NOT converged"
    else:
      # Working under the assumption that this was an SCF job
      engLine = findLine(outFile, "!")
      energyRy = float(outFile[engLine].split()[4])
      energyJ  = energyRy * 2.1798715e-18
      surfaceEnergy = energyJ / (aLength * bLength) # Energy in J/m^2
      print "   Energy Ry = {0}\n   Surface Energy = {1}".format(energyRy, surfaceEnergy)
      energies[b][a] = surfaceEnergy

# Add the energies we just calculated to the resultsString
for b in range(len(bVals)):
  resultsString += bVals[b]
  for a in range(len(aVals)):
    resultsString += "  " + str(energies[b][a])
  resultsString += "\n"

# Write the resultsString to resultsFilename
writer = open(resultsFilename, 'w')
writer.write(resultsString)
writer.close()

# Calculate relative energies and format them for making the plot in gnuplot
stepSize = float(aVals[1]) - float(aVals[0])

# Ensure that the starting high and low energies aren't "XXXX"
for i in range(len(energies[0])):
  if energies[0][i] != "XXXX":
    lowestEng = energies[0][i]
    highestEng = energies[0][i]
    break

# Find the highest and lowest energies
for i in energies:
  for j in i:
    if (j < lowestEng) and (j != "XXXX"):
      lowestEng = j
    elif (j > highestEng) and (j != "XXXX"):
      highestEng = j

print "Lowest energy  = {0}\nHighest energy = {1}".format(lowestEng, highestEng)

deltaE = []
for i in energies:
  temp = []
  for j in i:
    # If point didn't converge, set its relative energy to the largest difference observed
    if j == "XXXX":
      temp.append((highestEng - lowestEng))
    else:
      temp.append((j - lowestEng))
  deltaE.append(temp)

plotDataName = filename + "_plot.dat"
writer = open(plotDataName, 'w')
for i in deltaE:
  for j in i:
    writer.write("   " + str(j))
  writer.write("\n")
writer.close()

aRange = [float(aVals[0])-(stepSize/2), float(aVals[-1])-(stepSize/2)]
bRange = [float(bVals[0])-(stepSize/2), float(bVals[-1])-(stepSize/2)]

gnuplot = []
gnuplot.append("reset\n\nset terminal png size 700,524 enhanced font 'Verdana,10'\nset output '")
gnuplot.append("{0}.png'\n\nunset key\n\nset style line 11 lc rgb '#808080' lt 1\nset border 3 front ls 11\nset tics nomirror out scale 0.75\n\n".format(filename))
gnuplot.append('set cbtics scale 0\nset palette defined ( 0 "#000090", 1 "#000fff", 2 "#0090ff", 3 "#0fffee", 4 "#90ff70", 5 "#ffee00", 6 "#ff7000", 7 "#ee0000", 8 "#7f0000")\n\n')
# Need to have aMin and bMin minus stepsize/2 and aMax and bMax plus stepsize/2
gnuplot.append("set xrange [{0[0]}:{0[1]}]\nset yrange [{1[0]}:{1[1]}]\n".format(aRange, bRange))
gnuplot.append("set xlabel 'delta A'\nset ylabel 'delta B'\nset cblabel 'Relative energy (J/m^2)'\n\n")
# Need to multiply $1 and $2 by step distance and add minimum a/b values
gnuplot.append("plot '{0}' u (($1*{1})+({2})):(($2*{3})+({4})):($3) matrix with image".format(plotDataName, stepSize, aVals[0], stepSize, bVals[0]))
gnuplotName = sys.argv[1] + "_plot.gp"
gnuWriter = open(gnuplotName, 'w')
gnuWriter.write(''.join(gnuplot))
gnuWriter.close()

# eof
