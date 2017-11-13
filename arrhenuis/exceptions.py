from arrhenius import latex_fractions, draw_fractions, draw_png, draw_legend, draw_bw, create_pH_table, gen_full_pH_table

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

    latex_fracs = latex_fractions(acro, name, pKa, pH1, pH2, create_pH_table=create_pH_table, gen_full_pH_table=gen_full_pH_table)
    with open(acro + '-fractions.tex', 'w') as f:
        f.write(latex_fracs)

    draw_fractions(acro, pKa, pH1, pH2, draw_png=draw_png, draw_legend=draw_legend, draw_bw=draw_bw)

    print '...done'
