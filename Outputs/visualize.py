import numpy as np 
import matplotlib.pyplot as plt 
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--file", default="khouot.txt")
args = parser.parse_args()

arr = np.loadtxt(args.file)
print (arr.shape)

plt.imshow(arr)
plt.show()