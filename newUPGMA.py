from Bio import Phylo
from io import StringIO
import copy
import re
def makeLabels(firstLetter, lastLetter):
    labels = []
    for i in range(ord(firstLetter), ord(lastLetter)+1):
        labels.append(chr(i))
    return labels

def findSmallestElement(distanceTable):
    minNumber = float("inf")
    minXCoordinate = -1
    minYCoordinate = -1
    for rowIndex in range(len(distanceTable)):
        for elementIndex in range(len(distanceTable[rowIndex])):
            if distanceTable[rowIndex][elementIndex] < minNumber:
                minXCoordinate = rowIndex
                minYCoordinate = elementIndex
                minNumber = distanceTable[minXCoordinate][minYCoordinate]

    return minXCoordinate, minYCoordinate, minNumber


def theAlgorithm(distanceTable, labels, xCoordinate, yCoordinate, minNumber):
    newDistanceTable = []
    newDistanceTable = copy.deepcopy(distanceTable)
    tableLength = len(distanceTable)
    
    del(newDistanceTable[xCoordinate])
    del(newDistanceTable[yCoordinate])
    
    for j in range(len(newDistanceTable)):
        if len(newDistanceTable[j]) > yCoordinate and len(newDistanceTable[j]) < xCoordinate:
            del(newDistanceTable[j][yCoordinate])
        elif len(newDistanceTable[j]) > yCoordinate and len(newDistanceTable[j]) > xCoordinate:
            del(newDistanceTable[j][xCoordinate])
            del(newDistanceTable[j][yCoordinate])
    
    newRow = []
    for h in range(tableLength):
        if h == xCoordinate or h == yCoordinate:
            continue
        if h < xCoordinate and h < yCoordinate:
            newRow.append((distanceTable[xCoordinate][h] + distanceTable[yCoordinate][h])/2)
        elif h < xCoordinate and h > yCoordinate:
            newRow.append((distanceTable[xCoordinate][h] + distanceTable[h][yCoordinate])/2)
        else:
            newRow.append((distanceTable[h][xCoordinate] + distanceTable[h][yCoordinate])/2)
    newDistanceTable.append(newRow)

    firstElement = labels[xCoordinate]
    secondElement = labels[yCoordinate]

    print(labels)

    newstrX = "0"
    newstrY = "0"
    if firstElement.strip()[-1] == ")":
        newString = firstElement.split(":")
        newstrX = newString[len(newString) - 1].replace(")", "")
    
    if secondElement.strip()[-1] == ")":
        newString = secondElement.split(":")
        newstrY = newString[len(newString) - 1].replace(")", "")
    


    firstElement += ":" + str(minNumber/2-float(newstrX))
    secondElement += ":" + str(minNumber/2-float(newstrY))

    
    del(labels[xCoordinate])
    del(labels[yCoordinate])
    labels.append("(" + firstElement + "," + secondElement + ")")

    return newDistanceTable




def upgma(table, labels):
    newTable = table
    while len(labels) > 1:
        xCoordinate, yCoordinate, minNumber = findSmallestElement(newTable)
        newTable = theAlgorithm(newTable, labels, xCoordinate, yCoordinate, minNumber)
    finalTree = labels[0]
    
    print(finalTree)
    handle = StringIO(finalTree)
    print(handle)
    tree = Phylo.read(handle, "newick")
    print(tree)
    Phylo.draw_ascii(tree)
    Phylo.draw(tree, branch_labels=lambda c: c.branch_length)




labels = makeLabels("A", "E")   #A through G
distanceMatrix = [
    [],                         #A
    [20],                       #B
    [60, 50],                   #C
    [100, 90, 40],              #D
    [90, 80, 50, 30]            #E
    ]




upgma(distanceMatrix, labels)
