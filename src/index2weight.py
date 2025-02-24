import sys
j = int(sys.argv[1])
subsets = [[e,f,h] for e in range(11) for f in range(11) for h in range(11) if e+f+h <= 10]
print(sum([i for i in subsets[j]]))