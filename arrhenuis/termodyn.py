from arrhenius import stringify

"""
This script generates cool table with termodynamic parameters, taken from 'termodyn-acid.data.txt' and 'termodyn-anion.data.txt'
"""

head = '''
\\begin{center}
\\begin{tabular}{ x{1cm}  x{2cm}  x{2cm}  x{2cm}  x{2cm}  x{2cm} } \\hline
Acid & E$_A$ \\newline kJ/mol & A \\newline M$^{-1}$s$^{-1}$ & \
$\\Delta \\text{S} \\ddag$ \\newline J / mol ${\\cdot}$ K & $\\Delta \\text{H} \\ddag$ \\newline kJ / mol & $\\Delta \\text{G} \\ddag$ \\newline kJ / mol \\\\ \\hline
'''

feet = '''
\\hline
\\end{tabular}
\\end{center}
'''

acidstrings, anionstrings = list(), list()

data = open('termodyn-acid.data.txt').readlines()
for line in data:
    acro = line.split()[0]
    ea, a, s, h, g = map(float, line.split()[1:])
    string = '{} & ${:.1f}$ & ${:.2e}$ & ${:.1f}$ & ${:.1f}$ & ${:.1f}$ \\\\'.format(acro, ea, a, s, h, g)
    acidstrings.append(string)

data = open('termodyn-anion.data.txt').readlines()
for line in data:
    acro = line.split()[0]
    ea, a, s, h, g = map(float, line.split()[1:])
    string = '{} & ${:.1f}$ & ${:.2e}$ & ${:.1f}$ & ${:.1f}$ & ${:.1f}$ \\\\'.format(acro, ea, a, s, h, g)
    anionstrings.append(string)

latex = 'Reactions of OH radicals with halogenated acids \n\n'  + \
         head + ' \n '.join(acidstrings)  + '\n' + feet + '\n\n' + \
        'Reactions of OH radicals with haloacetate and halopropiate anions \n\n' + \
         head + ' \n '.join(anionstrings) + '\n' + feet

for i in range(5, 16):
    latex = latex.replace('e+{:02d}'.format(i), ' \\cdot 10^{{{0}}} '.format(i))
latex = latex.replace('$nan$', '---')

with open('termodyn.tex', 'w') as f:
    f.write(latex)
