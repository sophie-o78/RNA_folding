#Nussinov algorithm
#Goal: maximize the number of pairings

#useful library
import numpy as np


RNA_S = str(input("Enter RNA Sequence:"))
#user can test by entering GGACCUUUGGACUC, the given sequence of RNA
#the input must be in capital letters
L = len(RNA_S)


def initialize_matrix():

    # create an empty matrix
    N = np.zeros([L, L])
    return N


def pairing_valid(pair):

    valid_bp = {"A": "U", "U": "A", "G": "C", "C": "G"} # or a list of tuples
    if pair in valid_bp.items(): #if pair is in the valid base pairs list
        return True
    return False #otherwise it is not valid

#fill the scoring matrix
def fill_N(N, rna):

    minimal_loop_length = 0
    for k in range(1, len(rna)):
        for i in range(len(rna) - k):
            j = i + k
            if j - i >= minimal_loop_length:
                down = N[i + 1][j] #rule
                left = N[i][j - 1] # 2nd rule
                diag = N[i + 1][j - 1] + pairing_valid((rna[i], rna[j])) # 3rd rule
                rc = max([N[i][t] + N[t + 1][j] for t in range(i, j)]) # 4th rule
                N[i][j] = max(down, left, diag, rc) # max of all
                # returns the score of the optimal pairing between indices i and j

            else:
                N[i][j] = 0
    return N

#trace back the path
def backtrack_pathN(N, rna_seq, fold, i, L):
    j = L
    if i < j:
        if N[i][j] == N[i + 1][j]: #rule 1
            backtrack_pathN(N, rna_seq, fold, i + 1, j) #recursive algorithm
        elif N[i][j] == N[i][j - 1]: #rule 2
            backtrack_pathN(N, rna_seq, fold, i, j - 1)
        elif N[i][j] == N[i + 1][j - 1] + pairing_valid((rna_seq[i], rna_seq[j])): #rule 3
            fold.append((i, j))
            backtrack_pathN(N, rna_seq, fold, i + 1, j - 1)
        else:
            for k in range(i + 1, j - 1):
                if N[i][j] == N[i, k] + N[k + 1][j]: #rule 4
                    backtrack_pathN(N, rna_seq, fold, i, k)
                    backtrack_pathN(N, rna_seq, fold, k + 1, j)
                    break
    return fold


def pairing_list():
    pair_list = []
    for coord in (backtrack_pathN(N, RNA_S, [], 0, L-1)):
        pos1, pos2 = coord[0], coord[1]
        pairing = RNA_S[pos1], RNA_S[pos2]
        if pairing_valid(pairing):
            pair_list.append(pairing)
    return pair_list


#Nussinov's scoring matrix

N = initialize_matrix()

print("\nHere is the Nussinov scoring matrix: \n")
print(fill_N(N, RNA_S))


#traceback path

print("\nHere is the traceback path: \n")
print(backtrack_pathN(N, RNA_S, [], 0, L-1))

#list of pairings

print("\nHere are the pairings obtained: \n")
print(pairing_list())
