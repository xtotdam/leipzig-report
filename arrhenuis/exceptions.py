from arrhenius import fractions

"""This script is for compounds, we can only plot fractions for, because their rate constants are too slow"""

# Memo: fractions(acro, name, pKa, pH1, pH2)
acids = [
['TCA', 'trichloroacetic acid',         0.65, 5, 5],
['TFA', 'trifluoroacetic acid',        -0.56, 5, 5],
['DFPA', 'difluoropropionic acid',      0.56, 5, 5],
['TFPA', 'tetrafluoropropionic acid',   1.34, 5, 5]
]

for line in acids:
    acro, name, pKa, pH1, pH2 = line

    print acro,
    latex_fractions = fractions(acro, name, pKa, pH1, pH2, create_table=False)
    with open(acro + '-fractions.tex', 'w') as f:
        f.write(latex_fractions)
    print '...done'
