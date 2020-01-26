import numpy as np

def is_sorted(x):
	return np.all(np.diff(x)>0)

def arrayfy(x):
	return [np.array(x) for x in zip(*x)]