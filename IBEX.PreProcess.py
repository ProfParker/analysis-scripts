# Python script for formatting IBEX acceptability judgment data for analysis in R.
# Extracts subject, condition, and rating information.
# Dan Parker January 2013


import string
import os
import sys
import fileinput
import re

# Get user information
myFile = raw_input("Name of IBEX data file (enter full path if file is not in current directory): ")
myData = open(myFile)
myData = myData.read().split('\n')

numSubjs = int(input("How many participants did you test? "))
numItems = int(input("How many items were in your experiment  (items + fillers)? "))

out = raw_input("Name of output file (enter full path if you want to save to a different directory): ")
outputFile = open(out, 'w')

# Extract relevant data
prep1 = []
for line in myData:
    #if "NULL,NULL" in line:
	if re.search(r"NULL\,[1-7]\,", line):
		prep1.append(line)
		
# Remove practice data
myData = []
for line in prep1:
    if not "practice" in line:
		myData.append(line)

# Extract condition and rating information
myRatings = []
for line in myData:
	temp = line.split(",")
	item = temp[3]
	cond = temp[5]
	rating = temp[8]
	rt = temp[10]
	#print (cond + "\t" + rating)
	myRatings.append(item + "\t" + cond + "\t" + rating + "\t" + rt)
	
# Add subject information
subj = []
for i in range(1,numSubjs+1):
	for j in range(1,numItems+1):
		subj.append(i)

# Output formatting
output = []
output.append("\t" + "Subj" + "\t" + "Item" + "\t" + "Cond" + "\t" + "Rating" + "\t" + "RT")
for i, (a, b) in enumerate(zip(subj,myRatings)):
	output.append(str(i+1) + "\t" + str(a) + "\t" + b)

# Final Output	
for line in output:
  	outputFile.write(line + "\n")

print "That's all, Folks!"

sys.exit()