# Longest Repetition

# Define a procedure, longest_repetition, that takes as input a 
# list, and returns the element in the list that has the most 
# consecutive repetitions. If there are multiple elements that 
# have the same number of longest repetitions, the result should 
# be the one that appears first. If the input list is empty, 
# it should return None.

def longest_repetition(l):
    longest=[None,0] # [element, occurrence]
    current=[None,0]
    for e in l:
        if e != current[0]:
            current=[e, 1]
        else:
            current[1]=current[1]+1
        if current[1]>longest[1]:
            longest=current
    if current[1]>longest[1]:
        longest=current
    return longest[0]
