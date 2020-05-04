from Bio import Phylo
from io import StringIO
import copy
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

    # print(xCoordinate)
    # print(yCoordinate)
    # print(minNumber)
    # print(labels)

    firstElement = labels[xCoordinate] + ":" + str(minNumber/2)
    secondElement = labels[yCoordinate] + ":" + str(minNumber/2)
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
    
    # print(finalTree)
    # handle = StringIO(finalTree)
    # print(handle)
    # tree = Phylo.read(handle, "newick")
    # print(tree)
    # Phylo.draw_ascii(tree)


    # stringTree = "(((E:15.0,D:15.0):7.5,C:22.5):13.75,(B:10.0,A:10.0):26.25)"
    # stringHandle = StringIO(stringTree)
    # newStringTree = Phylo.read(stringHandle, "newick")
    # Phylo.draw(newStringTree, branch_labels=lambda c: c.branch_length)



    # Phylo.draw(tree, branch_labels=lambda c: c.branch_length)

    print(finalTree)
    stringTree = "(((E:15.0,D:15.0):7.5,C:22.5):13.75,(B:10.0,A:10.0):26.25)"
    print(stringTree)
    newTestTree = finalTree
    x = newTestTree.split(":")
    for i in x:
        print(i)











# labels = makeLabels("A", "G")   #A through G
# distanceMatrix = [
#     [],                         #A
#     [19],                       #B
#     [27, 31],                   #C
#     [8, 18, 26],                #D
#     [33, 36, 41, 31],           #E
#     [18, 1, 32, 17, 35],        #F
#     [13, 13, 29, 14, 28, 12]    #G
#     ]


labels = makeLabels("A", "E")   #A through G
distanceMatrix = [
    [],                         #A
    [20],                       #B
    [60, 50],                   #C
    [100, 90, 40],              #D
    [90, 80, 50, 30]            #E
    ]


# labels = makeLabels("A", "H")   #A through C
# distanceMatrix = [
#     [],                             #A
#     [10],                           #B
#     [12, 16],                       #C
#     [10, 8, 10],                    #D
#     [12, 14, 6, 10],                #E
#     [12, 14, 16, 18, 20],           #F
#     [10, 12, 14, 16, 18, 20],        #G
#     [20, 22, 24, 26, 30, 32, 40]    #H
#     ]

#upgma(distanceMatrix, labels)


upgma(distanceMatrix, labels)
