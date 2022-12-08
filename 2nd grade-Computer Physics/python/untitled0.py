import numpy as np
sum = 0
b = np.loadtxt("P13.txt",str)
for i in range(len(b)):
    sum += int(b[i])
number = str(sum)
for j in range(10):
    print(number[j])
