import numpy as np

"""This script generates fake data for the example"""

def doit():
	"""Summary
	
	Returns:
	    TYPE: Description
	"""
	a = 10**(10 + 4 * np.random.rand())
	ea = - (10 + 20 * np.random.rand())

	alltemps = np.arange(278, 318+1, 5)
	num = np.random.randint(6, 8)
	temps = np.sort(np.random.choice(alltemps, num, replace=False))

	error1 = (5 + 20 * np.random.rand(temps.shape[0])) / 100.
	error2 = (10 + 70 * np.random.rand(temps.shape[0])) / 100.

	k = a * np.exp(ea / 8.314e-3 / temps) * error1
	dk = k * error2

	t   = ('{:d}\t\t'*num).format(*temps)
	ks  = ('{:.2e}\t'*num).format(*k)
	dks = ('{:.2e}\t'*num).format(*dk)

	for i in range(1, 20):
	    ks = ks.replace('e+{:02d}'.format(i), 'e{:d}'.format(i))
	    dks = dks.replace('e+{:02d}'.format(i), 'e{:d}'.format(i))

	print t
	print ks
	print dks
	print

for i in range(12):
	doit()