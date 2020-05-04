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




def updateLabels(xCoordinate, yCoordinate, minNumber):
    firstElement = labels[xCoordinate]
    secondElement = labels[yCoordinate]

    firstElementTotal = labelsTotal[xCoordinate]
    secondElementTotal = labelsTotal[yCoordinate]

    newstrX = "0"
    newstrY = "0"
    if firstElementTotal.strip()[-1] == ")":
        newString = firstElementTotal.split(":")
        newstrX = newString[len(newString) - 1].replace(")", "")
    
    if secondElementTotal.strip()[-1] == ")":
        newString = secondElementTotal.split(":")
        newstrY = newString[len(newString) - 1].replace(")", "")
    


    firstElement += ":" + str(minNumber/2-float(newstrX))
    secondElement += ":" + str(minNumber/2-float(newstrY))

    firstElementTotal += ":" + str(minNumber/2)
    secondElementTotal += ":" + str(minNumber/2)

    
    del(labels[xCoordinate])
    del(labels[yCoordinate])
    labels.append("(" + firstElement + "," + secondElement + ")")


    del(labelsTotal[xCoordinate])
    del(labelsTotal[yCoordinate])
    labelsTotal.append("(" + firstElementTotal + "," + secondElementTotal + ")")





def theAlgorithm(distanceTable, xCoordinate, yCoordinate):
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

    

    return newDistanceTable




def upgma(table, labels):
    newTable = table
    while len(labels) > 1:
        xCoordinate, yCoordinate, minNumber = findSmallestElement(newTable)
        newTable = theAlgorithm(newTable, xCoordinate, yCoordinate)
        updateLabels(xCoordinate, yCoordinate, minNumber)
    finalTree = labels[0]
    
    print(finalTree)
    handle = StringIO(finalTree)
    print(handle)
    tree = Phylo.read(handle, "newick")
    print(tree)
    Phylo.draw_ascii(tree)
    Phylo.draw(tree, branch_labels=lambda c: c.branch_length)




labels = makeLabels("A", "E")   #A through G
labelsTotal = copy.deepcopy(labels)
distanceMatrix = [
    [],                         #A
    [20],                       #B
    [60, 50],                   #C
    [100, 90, 40],              #D
    [90, 80, 50, 30]            #E
    ]



# labels = makeLabels("A", "H")   #A through G
# labelsTotal = copy.deepcopy(labels)


# distanceMatrix = [
#     [],                         #A
#     [32],                       #B
#     [48, 26],                   #C
#     [51, 34, 42],              #D
#     [50, 29, 44, 44],           #E
#     [48, 33, 44, 38, 24],            #F
#     [98, 84, 92, 86, 89, 90],            #G
#     [148, 136, 152, 142, 142, 142, 148]            #H
#     ]






upgma(distanceMatrix, labels)
