#Zucker algorithm
#Goal: minimize free energy

#useful libraries
import numpy as np
import math


RNA_S = str(input("Enter RNA Sequence:"))
#user can test by entering GGACCUUUGGACUC, the given sequence of RNA.
#the input must be in capital letters
L = len(RNA_S)

#initialize the matrix W and E
def initialize_matrix():

    x = np.zeros((L, L))
    y = np.zeros((L, L))
    for i in range(L):
        for j in range(L):
            if j - i < 5:
                x[i][j] = np.inf
                y[i][j] = np.inf
            else:
                x[i][j] = np.nan
                y[i][j] = np.nan
    return x, y


W, E = initialize_matrix()



#pairing energy function
def s(i,j):

    if RNA_S[i]+RNA_S[j] in ["AU","UA","GC","CG",]:
        return -4 #s(i,j) = -4 if A-U, C-G pairings

    if RNA_S[i]+RNA_S[j] in ["GU", "UG"]:
        return 0 #s(i,j) = 0 if G-U pairing
    else:
        return 4

#hairpin loop function
def h(i, j):
    return 2*(j-i+5)


#function to get the minimum between W(k,i) and W(j,k-1)
def get_mini(j, i):

    minimum = float('inf') #temporary minimum
    min_k = -1  #minimum between W(k,i) and W(j,k-1)
    for k in range (j + 2, i):
        result = fillW(k, i) + fillW(j, k - 1)
        if result < minimum :
            minimum = result
            min_k = k
    return minimum

#fill the W matrix
def fill_W(j, i):

    if i >= j + 5:
        W[j][i] = min(
            fillW(j, i -1),
            fillW(j+1, i),
            fill_E(j, i),
            get_mini(j, i)
        )

    return W[j][i]

#energy matrix of optimized folding
def fill_E(j, i):
#for a base pair (i,j)
    E[j][i] = min(s(j,i) + h(j+1, i-1), s(j,i) + fillW(j+1, i-1))
    return E[j][i]


fill_W(0, L-1)
fill_E(0, L-1)

print("Here is the W matrix: \n")
print(W)

#final energy matrix

print("\nHere is the Zuker energy matrix: \n")
print(E)


#backtrack path:  to save directions of path while inspecting values from neighbors
def backtrack():
    backtrack = []
    directions = []
    i, j = len(RNA_S)-1, 0
    backtrack.append((j,i))
    while i-j >= 5:

        if W[j][i] == W[j][i-1]:
            directions.append("left")
            i -= 1
            backtrack.append((j,i))
        elif W[j][i] == W[j+1][i]:
            directions.append("down")
            j += 1
            backtrack.append((j,i))
        elif W[j][i] == E[j][i]:
            if E[j][i] == s(j,i) + h(j+1, i-1):
                directions.append("hairpin")
            elif E[j][i] == s(j,i) + W[j+1, i-1]:
                directions.append("diagonal")
            j += 1
            i -= 1
            backtrack.append((j,i))
    return backtrack

print("\nHere is the traceback path: \n")
print(backtrack())
print("\n Here are the correspondences: \n")
print(directions)

#list of the pairings that occur in the sequence
def pairing_list():
    pair_list = []
    for coord in (backtrack()):
        pos1, pos2 = coord[0], coord[1]

        pairing = RNA_S[pos1], RNA_S[pos2]
       # if pairing_valid(pairing):
        pair_list.append(pairing)
    return pair_list


#print the list of pairings

print("\nHere are the pairings obtained: \n")
print(pairing_list())




