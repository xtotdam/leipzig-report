import numpy as np

"""This script generates fake data for the example"""

header = '''\
// ======================================
//            Raw data file
// ======================================

// blank lines are ignored, as well as everything behind '//', just like in C language. No /* */ comment-outs though.
//
// Format:

// full name
// acronym
// pKa  (999 if unknown)
// pH #1
// temperatures #1              (T)
// rate constants #1            (k)
// rate constants errors #1     (dk)
// pH #2
// temperatures #2
// rate constants #2
// rate constants errors #2
// ===  is a delimiter between substances
//
// ============================================================================
// it works only with two sets of data for two different pH
// if only one exists / is needed, use suitable junk data
// minimum two points of data for each pH -- it is needed for regression not to throw exception
// ============================================================================
\n\n'''

acids = [
	['2-chloropropionic acid', '2CPA'],
	['2-fluoropropionic acid', '2FPA'],
	['dichloroacetic acid', 'DCA'],
	['difluoroacetic acid', 'DFA'],
	['chloroacetic acid', 'MCA '],
	['fluoroacetic acid', 'MFA']
]

def generate_tkdk():
	"""
	Returns 3 strings:
	temperatures
	constants
	errors
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

	return '\n'.join((t, ks, dks))

def generate_pka():
	'''
	Returns pKa in bounds [1, 3)
	'''
	return '{:.3f}'.format((np.random.rand(1) * 2 + 1)[0])


if __name__ == '__main__':
	with open('data.txt', 'w') as f:
		f.write(header)
		f.write('\n\n// JUNK DATA BELOW\n\n')

		for (i, (name, acro)) in enumerate(acids):
			f.write(
				name + '\n' +
				acro + '\n' +
				generate_pka() +
				'\n\n1\n' +
				generate_tkdk() +
				'\n\n5\n' +
				generate_tkdk()
				)

			if not i == len(acids) - 1:
				f.write('\n\n===\n\n')
