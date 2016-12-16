from __future__ import division
import csv
import numpy as np

data1FilePath = "data1.txt"
data2FilePath = "data2.txt"

def readInputFile(filePath):
#    retVal = []
#    with open(filePath, 'rb') as csvfile:
#        filereader = csv.reader(csvfile, delimiter=' ', quotechar='|')
#        for row in filereader:
#            retVal.append([int(row[0]), int(row[1]), int(row[2])])

    retVal=[]
    with open(filePath) as file:
         line=file.readline()
         arr=[float(a) for a in line.split(',')]
 #        retVal.append(file.readline())
         retVal.append(arr)
    return retVal

def transient_removal(x=[]):
    x = readInputFile(data2FilePath)
#    x = np.loadtxt(data1FilePath)
    
    x = np.array(x)
    N = x.shape[1]
    y=[] 
    for kk in np.arange(np.floor(N/2)):
        v = np.var(x[kk+1:])*1.0/(N-kk)
        y.append(v)
    y = np.array(y)
    ind = y.min(axis=0)
    print('index of transient point in the signal:')
    print(ind)
    return ind


x = readInputFile(data2FilePath)
print(x)
index = transient_removal()
#print 
