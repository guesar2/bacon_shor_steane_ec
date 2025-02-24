import math
import re

def numNonOverlapping11Subs(diffVec):
    """Returns the number of non-overlapping 11 substrings in a given difference vector
    and the index corresponding to the latest syndrome measurments

    Args:
        diff_vec (str): Difference vector
    """
    chunks = re.findall(r'0+|1+', diffVec)
    numSubStr = 0
    latestRound = 0
    for chunk in chunks:
        
        if all(char == '1' for char in chunk):
            numSubStr += math.floor(len(chunk)/2) 
        
        latestRound += len(chunk)
    return numSubStr, latestRound

def find_usable_substring(diff_vec, t):
    """Returns indexes corresponding to usable substrings according to Ink's algorithm and the weight of error that a sta
    bilizer code of distance d can correct.

    Args:
        diff_vec (str): Difference vector
        t (int): Distance
    """
    usables = []
    chunks = re.findall(r'0+|1+', diff_vec)
    index = 0
    for i in range(len(chunks)):
        if all(char == '0' for char in chunks[i]):
            f = 0
            j = 0 
            for j in range(len(chunks)):
                if all(char == '1' for char in chunks[j]):
                    n_ones  = len(chunks[j])
                    if abs(i-j) == 1:
                        n_ones -= 1
                    f += math.ceil(float(n_ones)/2)
            y = len(chunks[i])
            if (y+f) >= t:
                usables.append(index)
        index += len(chunks[i])
    return usables