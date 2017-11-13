from scipy.stats import linregress
import numpy as np
from pprint import pprint
from matplotlib import pyplot as plt, ticker
import matplotlib.lines as mlines
import sys

import matplotlib as mpl
mpl.style.use('classic')

from hh import calculate_params

plt.rc('font', family='serif')
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15

R = 8.314e-3
temps = np.arange(278, 318+1, 5)

# A bit of settings here
draw_png = False        # should I create also pngs in addition to pdfs?
draw_legend = True      # should I draw legend on plots?
draw_bw = False         # should I draw plots black and white?

def stringify(a, da, digits=2):
    """Converts number and its error to nifty latex string

    Args:
        a (float): number
        da (float): its error
        digits (int, optional): round up to how many digits

    Returns:
        str: latex string of number with its error
    """
    template = '({{:.{0:d}f}} \\pm {{:.{0:d}f}})'.format(digits)    # x +- y, formatted

    alog = np.trunc(np.log10(abs(a)))
    dalog = np.trunc(np.log10(abs(da)))

    b = a / 10**alog
    db = da / 10**alog

    try:
        return (template + '\\cdot 10^{{{:d}}}').format(b, db, int(alog))
    except ValueError:
        return ''

def read_data_file(fn):
    """Obtains data from file with raw data (usually 'data.txt')

    Args:
        fn (str): filename

    Returns:
        list: parsed data from file - datasets
    """
    with open(fn) as f:
        data = f.readlines()

    data = filter(bool, (line.split('//')[0].replace('\t', ' ').strip() for line in data))

    splitby = list()
    for i, v in enumerate(data):
        if '=' in v:
            splitby.append(i)
    splitby.append(len(data))

    datasets = list()
    first = 0
    for i in splitby:
        if not first:
            datasets.append(data[first:i])
        else:
            datasets.append(data[first + 1:i])
        first = i

    # TODO if data.txt ends with === error occurs
    return datasets

def parse(data):
    """Puts data from dataset to right variables

    Args:
        data (list): dataset read from data file

    Returns:
        tuple: many variables containing much information
    """
    name = data[0].strip()
    acro = data[1].strip()
    pKa = float(data[2].strip())
    pKa = pKa if pKa != 999. else None

    pH1 = float(data[3])
    T1 =  np.array(map(float, filter(bool, data[4].split())), dtype=np.float)
    k1 =  np.array(map(float, filter(bool, data[5].split())), dtype=np.float)
    dk1 = np.array(map(float, filter(bool, data[6].split())), dtype=np.float)

    pH2 = float(data[7])
    T2 =  np.array(map(float, filter(bool, data[8].split())), dtype=np.float)
    k2 =  np.array(map(float, filter(bool, data[9].split())), dtype=np.float)
    dk2 = np.array(map(float, filter(bool, data[10].split())), dtype=np.float)

    return name, acro, pKa, pH1, T1, k1, dk1, pH2, T2, k2, dk2

def plain_output(res):
    """Formats arrhenius and thermodynamic data for plain text output

    Args:
        res (dict): arrhenius and thermodynamic data dictionary

    Returns:
        str: formatted data
    """
    return '''
                            :  pH = {ph:.1f}
    ORDINATENABSCHNITT      :  {abschnitt:.5g} +- {abschnittsfehler:.5g}
    REL. FEHLER             :  {rel_fehler_abschnitt:3.1f} %
    STEIGUNG                :  {steigung:.5g} +- {steigungsfehler:.5g}
    REL. FEHLER             :  {rel_fehler_steigung:3.1f} %
    KORRELATIONSKOEFFIZIENT :  {korrelationskoeffizient:.5f}
    AKTIVIERUNGSENERGIE     :  {ea:.5g} +- {eafehler:.5g} KJ/MOL
    ARRHENIUS-VORFAKTOR     :  {a:.5e} +- {afehler:.5e} 1/S
    BIMOLEKULARE REAKTION   :
    AKTIVIERUNGENTHALPIE    :  {h:.5f} +- {hfehler:.5f} KJ/MOL
    AKTIVIERUNGSENTROPIE    :  {s:.5f} +- {sfehler:.5f} J/MOL*K
    GIBBS-AKT.-ENERGIE      :  {g:.5f} +- {gfehler:.5f} KJ/MOL

    FEHLER MIT STUDENT-FAKTOR FUER 95% KONFIDENZINTERVALL.
    '''.format(**res)

def latex_output_params(res):
    """Formats arrhenius and thermodynamic data for LaTeX text output

    Args:
        res (dict): arrhenius and thermodynamic data dictionary

    Returns:
        str: LaTeX formatted data
    """
    string = '''
    \\begin{samepage}
    \\begin{tabular}{r l l}

    ~                                          & pH = {ph:.1f} \\\\
    \\bfseries Ordinatenabschnitt      \\quad~ &  ${abschnitt:.3f} \\pm {abschnittsfehler:.3f}$ \\\\
    \\bfseries Rel. Fehler             \\quad~ &  ${rel_fehler_abschnitt:3.1f}\\ \\% $ \\\\
    \\bfseries Steigung                \\quad~ &  ${steigung:.3f} \\pm {steigungsfehler:.3f}$ \\\\
    \\bfseries Rel. Fehler             \\quad~ &  ${rel_fehler_steigung:3.1f}\\ \\% $ \\\\
    \\bfseries Korrelationskoeffizient \\quad~ &  ${korrelationskoeffizient:.5f}$ \\\\
    \\bfseries Aktivierungsenergie     \\quad~ &  ${ea:.3f} \\pm {eafehler:.3f}$ kJ/mol \\\\
    \\bfseries Arrhenius-Vorfaktor     \\quad~ &  ${a:.3e} \\pm {afehler:.3e}$ 1/s \\\\
    \\bfseries Bimolekulare Reaktion   \\quad~ & ~ \\\\
    \\bfseries Aktivierungenthalpie    \\quad~ &  ${h:.3f} \\pm {hfehler:.3f}$ kJ/mol \\\\
    \\bfseries Aktivierungsentropie    \\quad~ &  ${s:.3f} \\pm {sfehler:.3f}$ J/mol ${{\\cdot}}$ K \\\\
    \\bfseries Gibbs-Akt.-Energie      \\quad~ &  ${g:.3f} \\pm {gfehler:.3f}$ kJ/mol \\\\

    \\end{{tabular}}
    \\end{{samepage}}
    '''.format(**res)

    for i in range(5, 16):
        string = string.replace('e+{:02d}'.format(i), ' \\cdot 10^{{{0}}} '.format(i))

    return string

def latex_double_output_params(res1, res2):
    """Formats two sets of arrhenius and thermodynamic data into one table for LaTeX text output

    Args:
        res1 (dict): first arrhenius and thermodynamic data dictionary
        res2 (dict): second arrhenius and thermodynamic data dictionary

    Returns:
        str: LaTeX formatted data
    """
    res = dict()
    for k in res1.keys():
        res[k + '1'] = res1[k]
    for k in res2.keys():
        res[k + '2'] = res2[k]

    la1 = stringify(res['a1'], res['afehler1'])
    la2 = stringify(res['a2'], res['afehler2'])
    res['la1'], res['la2'] = la1, la2

    string = '''
    \\begin{{samepage}}
    \\begin{{tabular}}{{r | l l}}

    ~                                          & pH = {ph1:.1f}                                   & pH = {ph2:.1f} \\\\ \\hline
    \\bfseries Ordinatenabschnitt      \\quad~ &  ${abschnitt1:.3f} \\pm {abschnittsfehler1:.3f}$ & ${abschnitt2:.3f} \\pm {abschnittsfehler2:.3f}$ \\\\
    \\bfseries Rel. Fehler             \\quad~ &  ${rel_fehler_abschnitt1:3.1f}\\ \\% $           & ${rel_fehler_abschnitt2:3.1f}\\ \\% $ \\\\
    \\bfseries Steigung                \\quad~ &  ${steigung1:.3f} \\pm {steigungsfehler1:.3f}$   & ${steigung2:.3f} \\pm {steigungsfehler2:.3f}$ \\\\
    \\bfseries Rel. Fehler             \\quad~ &  ${rel_fehler_steigung1:3.1f}\\ \\% $            & ${rel_fehler_steigung2:3.1f}\\ \\% $ \\\\
    \\bfseries Korrelationskoeffizient \\quad~ &  ${korrelationskoeffizient1:.5f}$                & ${korrelationskoeffizient2:.5f}$ \\\\
    \\bfseries Aktivierungsenergie     \\quad~ &  ${ea1:.3f} \\pm {eafehler1:.3f}$ kJ/mol         & ${ea2:.3f} \\pm {eafehler2:.3f}$ kJ/mol \\\\
    \\bfseries Arrhenius-Vorfaktor     \\quad~ &  ${la1}$ 1/s                                     & ${la2}$ 1/s \\\\
    \\bfseries Bimolekulare Reaktion   \\quad~ & ~                                                & ~ \\\\
    \\bfseries Aktivierungenthalpie    \\quad~ &  ${h1:.3f} \\pm {hfehler1:.3f}$ kJ/mol           & ${h2:.3f} \\pm {hfehler2:.3f}$ kJ/mol \\\\
    \\bfseries Aktivierungsentropie    \\quad~ &  ${s1:.3f} \\pm {sfehler1:.3f}$ J/mol ${{\\cdot}}$ K & ${s2:.3f} \\pm {sfehler2:.3f}$ J/mol ${{\\cdot}}$ K \\\\
    \\bfseries Gibbs-Akt.-Energie      \\quad~ &  ${g1:.3f} \\pm {gfehler1:.3f}$ kJ/mol           & ${g2:.3f} \\pm {gfehler2:.3f}$ kJ/mol \\\\

    \\end{{tabular}}
    \\end{{samepage}}
    '''.format(**res)

    #     \\bfseries Arrhenius-Vorfaktor     \\quad~ &  ${a1:.3e} \\pm {afehler1:.3e}$ 1/s              & ${a2:.3e} \\pm {afehler2:.3e}$ 1/s \\\\

    for i in range(5, 16):
        string = string.replace('e+{:02d}'.format(i), ' \\cdot 10^{{{0}}} '.format(i))

    return string

def pure_params(resacid, resanion):
    """Forms table of arrhenius and thermodynamic data (without regression data) for two `result` dictionaries

    Args:
        resacid (dict): `result` dictionary #1
        resanion (dict): `result` dictionary #2

    Returns:
        str: formatted LaTeX data string
    """
    res = dict()
    for k in resacid.keys():
        res[k + '1'] = resacid[k]
    for k in resanion.keys():
        res[k + '2'] = resanion[k]

    string = '''
    \\begin{{samepage}}
    \\begin{{tabular}}{{r | l l}}

    ~                                          & acid form                          & anion form \\\\ \\hline
    \\bfseries Aktivierungsenergie     \\quad~ & ${ea1:.3f}$ kJ/mol                 & ${ea2:.3f}$ kJ/mol \\\\
    \\bfseries Arrhenius-Vorfaktor     \\quad~ & ${a1:.3e}$ 1/s                     & ${a2:.3e}$ 1/s \\\\
    \\bfseries Bimolekulare Reaktion   \\quad~ & ~                                  & ~ \\\\
    \\bfseries Aktivierungenthalpie    \\quad~ & ${h1:.3f}$ kJ/mol                  & ${h2:.3f}$ kJ/mol \\\\
    \\bfseries Aktivierungsentropie    \\quad~ & ${s1:.3f}$ J/mol ${{\\cdot}}$ K    & ${s2:.3f}$ J/mol ${{\\cdot}}$ K \\\\
    \\bfseries Gibbs-Akt.-Energie      \\quad~ & ${g1:.3f}$ kJ/mol                  & ${g2:.3f}$ kJ/mol \\\\

    \\end{{tabular}}
    \\end{{samepage}}
    '''.format(**res)

    for i in range(5, 16):
        string = string.replace('e+{:02d}'.format(i), ' \\cdot 10^{{{0}}} '.format(i))

    return string

# TODO divide into 2 functions
def fractions(acro, name, pKa, pH1, pH2, create_table=True, draw_png=False, draw_legend=True, draw_bw=False):
    """Calculates fractions of anion and acid forms at two different pHs,
       creates a LaTeX table from this data and plots it to
       '*.fractions.pdf'

    Args:
        acro (str): acronym for compound
        name (str): its full name
        pKa (float): pKa
        pH1 (float): pH #1
        pH2 (float): pH #2
        create_table (bool, optional): whether we create a table or not.
                                       if False, then return empty string
        draw_png (bool): determines whether to create png image in addition to pdf
        draw_legend (bool): whether to draw legend or not
        draw_bw (bool): use colors or not (True = black'n'white)

    Returns:
        str: LaTeX table with data
    """
    if draw_bw:
        c_acid = c_anion = c_vline = 'k'
    else:
        c_acid, c_anion, c_vline = 'rgb'


    if pKa is None:
        with open(acro + '-fractions.tex', 'w') as f:
            f.write('')
        return

    pH = np.linspace(0, 14, 1000)
    a = 10.**(pH - pKa)
    acid = 1 / (1 + a) * 100.
    anion = a / (1 + a) * 100.

    fig = plt.figure(figsize=(10,4))
    ax = fig.gca()

    ax.plot(pH, acid, c=c_acid, lw=2.5, label='[AH][%]')
    ax.plot(pH, anion, c=c_anion, lw=2.5, label='[A$^{-}$][%]', dashes=(10,10))
    if draw_legend:
        ax.legend(handlelength=2.7)

    ax.axvline(pH1, c=c_vline, dashes=(6,8), lw=1.1)
    ax.axvline(pH2, c=c_vline, dashes=(6,8), lw=1.1)

    ax2 = ax.twiny()
    ax2.set_xlim(0, 14)
    ax2.set_xticks((pH1, pH2))

    ax.set_xlim(0, 14)
    ax.set_xticks(np.arange(0, 14+1))
    ax.set_xlabel('pH')
    # ax.set_ylabel('%')

    ax.yaxis.set_major_locator(ticker.MaxNLocator(10))
    ax.yaxis.set_minor_locator(ticker.MaxNLocator(50))

    fig.savefig(acro + '.fractions.pdf', bbox_inches='tight')
    if draw_png:
        fig.savefig(acro + '.fractions.png', bbox_inches='tight', dpi=250)

    # pH = np.linspace(0, 14, 29)
    pH = np.array([pH1, pH2], dtype=np.float)
    a = 10.**(pH - pKa)
    acid = 1 / (1 + a)
    anion = a / (1 + a)

    if create_table:
        tlines = ['$ {0:.1f} $ & $ {1:.8f} $ & $ {2:.8f} $ \\\\'.format(*item) for item in np.vstack((pH, acid, anion)).T]

        latex = '\\begin{samepage}\n\\begin{tabular}{r | l l}pH & [AH] & [A$^{-}$] \\\\\\hline\n' + \
                '\n'.join(tlines) + '\n\\end{tabular}\n\\end{samepage}'

        for i in range(5, 16):
            latex = latex.replace('e-{:02d}'.format(i), ' \\cdot 10^{{-{0}}} '.format(i))
    else:
        latex = ''

    latex += '''\n\n\n
    \\begin{figure}[H]
        \centering
        \includegraphics[width=0.8\\textwidth]{arrhenuis/{''' + acro + '''.fractions}.pdf}
        \caption{[AH] and [A$^{-}$] in percent for ''' + name + ''' for different pH values}
        \label{fig:''' + acro.lower() + '''-fractions}
    \end{figure}
    '''
    return latex

'''
def draw_plot(T, k, dk, result, pH, draw_png=False, draw_legend=True, draw_bw=False):
    """Draws plot with the data provided

    Args:
        T (1D ndarray): temperatures
        k (1D ndarray): rate constants
        dk (1D ndarray): rate constants errors
        result (dict): arrhenius and thermodynamic data dictionary
        pH (float): pH (for label)
        draw_png (bool): determines whether to create png image in addition to pdf
        draw_legend (bool): whether to draw legend or not
        draw_bw (bool): use colors or not (True = black'n'white)

    Returns:
        None
    """
    x = 1. / T
    y = np.log10(k)
    dylo = abs(y - np.log10(k - dk))
    dyhi = abs(y - np.log10(k + dk))
    yerr = np.stack((dylo, dyhi))

    fig = plt.figure(figsize=(10,8))
    ax = fig.gca()

    xlims = np.min(x) - 5e-5, np.max(x) + 5e-5
    ax.set_xlim(*xlims)
    # ax.set_ylim(np.min(y) - 0.02, np.max(y) + 0.02)

    ax.errorbar(x, y, yerr=yerr, fmt='o', c='k')
    xs = np.linspace(*xlims)
    ys = np.log10(result['a']) - result['ea'] / 8.314 * xs * 1000 / 2.3026
    ax.plot(xs, ys, lw=2, c='k')

    line = mlines.Line2D([], [], color='black', marker='o', ls='-', markersize=8, label='pH = ' + str(pH))
    if draw_legend:
        ax.legend(handles=[line], handlelength=3)

    ax.xaxis.set_major_locator(ticker.MaxNLocator(8))
    ax.xaxis.set_minor_locator(ticker.MaxNLocator(40))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(8))
    ax.yaxis.set_minor_locator(ticker.MaxNLocator(40))

    scale_pow = 3
    def my_formatter_fun(x, p):
        """Formatter for ticks"""
        return "%.1f" % (x * (10 ** scale_pow))
    ax.get_xaxis().set_major_formatter(ticker.FuncFormatter(my_formatter_fun))
    ax.set_xlabel('my label ' + '$10^{{{0:d}}}$'.format(scale_pow))

    ax.set_ylabel('lg (k / [M$^{-1}$s$^{-1}$])')
    ax.set_xlabel('T$^{-1}$ [10$^{-3}$ K$^{-1}]$')
    # ax.ticklabel_format(style = 'sci')

    ax2 = ax.twiny()
    ax2.set_xlim(*xlims)

    ax2.set_xticks(x)

    fig.canvas.draw()
    labels = [item.get_text() for item in ax2.get_xticklabels()]
    newlabels = list()
    for label in labels:
        try:
            newlabels.append(str(int(round(1. / float(label)))))
        except:
            newlabels.append(label)
    ax2.set_xticklabels(newlabels)
    ax2.set_xlabel('T [K]', va='bottom')

    fig.savefig(acro + '.pdf')
    if draw_png:
        fig.savefig(acro + '.png', dpi=300)
'''

def draw_double_plot(T1, k1, dk1, result1, pH1, T2, k2, dk2, result2, pH2, kacid, kanion, draw_png=False, draw_legend=True, draw_bw=False, singleplot=False):
    """Summary

    Args:
        T1 (1D ndarray): temperatures #1
        k1 (1D ndarray): rate constants #1
        dk1 (1D ndarray): rate constants errors #1
        result1 (dict): arrhenius and thermodynamic data dictionary #1
        pH1 (float): pH (for label) #1
        T2 (1D ndarray): temperatures #2
        k2 (1D ndarray): rate constants #2
        dk2 (1D ndarray): rate constants errors #2
        result2 (dict): arrhenius and thermodynamic data dictionary #2
        pH2 (float): pH (for label) #2
        kacid (1D ndarray): ideal rate constants for acid form
        kanion (1D ndarray): ideal rate constants for anion form
        draw_png (bool): determines whether to create png image in addition to pdf
        draw_legend (bool): whether to draw legend or not
        draw_bw (bool): use colors or not (True = black'n'white)
        singleplot (bool): workaround for draw_plot()

    Returns:
        None
    """
    if draw_bw:
        c_ph1 = c_ph2 = 'k'
    else:
        c_ph1, c_ph2 = 'rb'

    x1 = 1. / T1
    y1 = np.log10(k1)
    dylo1 = abs(y1 - np.log10(k1 - dk1))    # low error
    dyhi1 = abs(y1 - np.log10(k1 + dk1))    # high error
    yerr1 = np.stack((dylo1, dyhi1))        # put them together

    x2 = 1. / T2
    y2 = np.log10(k2)
    dylo2 = abs(y2 - np.log10(k2 - dk2))
    dyhi2 = abs(y2 - np.log10(k2 + dk2))
    yerr2 = np.stack((dylo2, dyhi2))

    fig = plt.figure(figsize=(10,8))
    ax = fig.gca()

    xlims = min(np.min(x1), np.min(x2)) - 5e-5, max(np.max(x1), np.max(x2)) + 5e-5
    ax.set_xlim(*xlims)
    ax.set_ylim(np.min(np.concatenate((y1-dylo1, y2-dylo2))) - 0.03, np.max(np.concatenate((y1+dyhi1, y2+dyhi2))) + 0.03)

    ax.errorbar(x1, y1, yerr=yerr1, fmt='s', c=c_ph1, lw=1.2, capthick=1.3)
    ax.errorbar(x2, y2, yerr=yerr2, fmt='o', c=c_ph2, lw=1.2, capthick=1.3, capsize=plt.rcParams['errorbar.capsize']*1.8)
    xs = np.linspace(*xlims)
    ys1 = np.log10(result1['a']) - result1['ea'] / 8.314 * xs * 1000 / 2.3026
    ys2 = np.log10(result2['a']) - result2['ea'] / 8.314 * xs * 1000 / 2.3026
    ax.plot(xs, ys1, lw=2, c=c_ph1, label='pH = ' + str(pH1))
    ax.plot(xs, ys2, lw=2, c=c_ph2, label='pH = ' + str(pH2), ls='--')

    if not singleplot:
        # plotting 'pure' data
        ax.scatter(1. / temps, np.log10(kacid),  c='k', marker='x',     s=120, label='k_acid', zorder=99)
        ax.scatter(1. / temps, np.log10(kanion), c='k', marker=(6,2,0), s=150, label='k_anion', zorder=99)

        # plotting 'pure' regression
        slope, intercept, r_value, p_value, std_err = linregress(1. / temps, np.log10(kacid))
        ax.plot(xs, xs * slope + intercept, c='k', ls=':', lw=1.5)
        slope, intercept, r_value, p_value, std_err = linregress(1. / temps, np.log10(kanion))
        ax.plot(xs, xs * slope + intercept, c='k', ls=':', lw=1.5)

    # legend
    line1 = mlines.Line2D([], [], color=c_ph1, lw=2, marker='s',     ls='-',  markersize=8,  label='pH = ' + str(pH1))
    line2 = mlines.Line2D([], [], color=c_ph2, lw=2, marker='o',     ls='--', markersize=8,  label='pH = ' + str(pH2))
    line3 = mlines.Line2D([], [], color='black',     marker='x',     ls=':',  markersize=11, label='k$_{acid}$')
    line4 = mlines.Line2D([], [], color='black',     marker=(6,2,0), ls=':',  markersize=12, label='k$_{anion}$')
    if draw_legend:
        if not singleplot:
            ax.legend(handles=[line1, line2, line3, line4], handlelength=3)
        else:
            ax.legend(handles=[line2, line4], handlelength=3)

    # messing with ticks
    ax.xaxis.set_major_locator(ticker.MaxNLocator(8))
    ax.xaxis.set_minor_locator(ticker.MaxNLocator(40))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(8))
    ax.yaxis.set_minor_locator(ticker.MaxNLocator(40))

    scale_pow = 3
    def my_formatter_fun(x, p):
        """Formatter for ticks"""
        return "%.1f" % (x * (10 ** scale_pow))
    ax.get_xaxis().set_major_formatter(ticker.FuncFormatter(my_formatter_fun))
    ax.set_xlabel('my label ' + '$10^{{{0:d}}}$'.format(scale_pow))

    ax.set_ylabel('lg (k / [M$^{-1}$s$^{-1}$])')
    ax.set_xlabel('T$^{-1}$ [10$^{-3}$ K$^{-1}]$')
    # ax.ticklabel_format(style = 'sci')

    # add temperatures on top axis
    ax2 = ax.twiny()
    ax2.set_xlim(*xlims)
    ax2.set_xticks(list(set( list(x1) + list(x2) )))

    fig.canvas.draw()
    labels = [item.get_text() for item in ax2.get_xticklabels()]
    newlabels = list()
    for label in labels:
        try:
            newlabels.append(str(int(round(1. / float(label)))))
        except:
            newlabels.append(label)
    ax2.set_xticklabels(newlabels)
    ax2.set_xlabel('T [K]', va='bottom')

    fig.savefig(acro + '.double.pdf')
    if draw_png:
        fig.savefig(acro + '.double.png', dpi=300)

def draw_plot(T, k, dk, result, pH, draw_png=False, draw_legend=True, draw_bw=False):
    """Draws plot with the data provided

    Args:
        T (1D ndarray): temperatures
        k (1D ndarray): rate constants
        dk (1D ndarray): rate constants errors
        result (dict): arrhenius and thermodynamic data dictionary
        pH (float): pH (for label)
        draw_png (bool): determines whether to create png image in addition to pdf
        draw_legend (bool): whether to draw legend or not
        draw_bw (bool): use colors or not (True = black'n'white)

    Returns:
        None
    """
    # !! NB : NOT TESTED AT ALL !!
    return draw_double_plot(T, k, dk, result, pH, T, k, dk, result, pH, None, None, draw_png, draw_legend, draw_bw, singleplot=True)

def split(acro, pKa, pH1, pH2, result1, result2, digits=2):
    """Summary

    Args:
        acro (str): compound acronym
        pKa (float): compound pKa
        pH1 (float): compound pH #1
        pH2 (float): compound pH #2
        result1 (dict): arrhenius and thermodynamic data dictionary #1
        result2 (dict): arrhenius and thermodynamic data dictionary #2
        digits (int, optional): how many digits to leave while rounding up the numbers for the equations

    Returns:
        latex (str): two equations, which determine rate constants for pure acid and anion forms
        kacid (1D ndarray): 'pure' rate constants for acid form
        kanion (1D ndarray): 'pure' rate constants for anion form
    """
    A1 = result1['a']
    A2 = result2['a']

    alpha = 1. - 1. / (1. + 10.**(pKa - pH1))
    beta  = 1. - 1. / (1. + 10.**(pKa - pH2))
    ab = alpha - beta

    ea1, ea2 = result1['ea'], result2['ea']

    # solving this system
    # A1 x1 = a * x1 + b * x1
    # A2 x2 = c * x1 + d * x1

    a =   A1 * (1 - beta)   / ab
    b = - A2 * (1 - alpha)  / ab
    c = - A1 * beta         / ab
    d =   A2 * alpha        / ab

    # print '{:+.4e}   {:+.4e}'.format(a, b)
    # print '{:+.4e}   {:+.4e}'.format(c, d)

    nums = dict()
    template = '{{:+.{}e}}'.format(digits)
    nums['a'] = template.format(a)
    nums['b'] = template.format(b)
    nums['c'] = template.format(c)
    nums['d'] = template.format(d)
    nums['ea1'] = '{:.2f}'.format(ea1 / R)
    nums['ea2'] = '{:.2f}'.format(ea2 / R)

    # forming latex numbers instead of 'Xe+Y'
    for i in range(5, 16):
        nums['a'] = nums['a'].replace('e+{:02d}'.format(i), ' \\cdot 10^{{{0}}} '.format(i))
        nums['b'] = nums['b'].replace('e+{:02d}'.format(i), ' \\cdot 10^{{{0}}} '.format(i))
        nums['c'] = nums['c'].replace('e+{:02d}'.format(i), ' \\cdot 10^{{{0}}} '.format(i))
        nums['d'] = nums['d'].replace('e+{:02d}'.format(i), ' \\cdot 10^{{{0}}} '.format(i))


    latex = '''
    \\begin{{align*}}
        k_\\text{{acid}} =  {a} \\cdot \\exp\\left(-\\frac{{{ea1}}}{{T}}\\right) {b} \\cdot \\exp\\left(-\\frac{{{ea2}}}{{T}}\\right) \\\\
        k_\\text{{anion}} = {c} \\cdot \\exp\\left(-\\frac{{{ea1}}}{{T}}\\right) {d} \\cdot \\exp\\left(-\\frac{{{ea2}}}{{T}}\\right)
    \\end{{align*}}
    '''.format(**nums)

    kacid = a * np.exp(- ea1 / R / temps) + b * np.exp(- ea2 / R / temps)
    kanion = c * np.exp(- ea1 / R / temps) + d * np.exp(- ea2 / R / temps)

    tkacid, tkanion = list(), list()

    for i in kacid:
        x = '{:.2e}'.format(i)
        for i in range(5, 16):
            x = x.replace('e+{:02d}'.format(i), ' \\cdot 10^{{{0}}} '.format(i))
        tkacid.append(x)

    for i in kanion:
        x = '{:.2e}'.format(i)
        for i in range(5, 16):
            x = x.replace('e+{:02d}'.format(i), ' \\cdot 10^{{{0}}} '.format(i))
        tkanion.append(x)

    return latex, kacid, kanion

def create_table(pH1, pH2, T1, k1, dk1, T2, k2, dk2, kacid, kanion):
    """Creates latex table with the input data and 'pure' rate constants

    Args:
        pH1 (float): pH #1
        pH2 (float): pH #2
        T1 (1D ndarray): temperatures #1
        k1 (1D ndarray): rate constants #1
        dk1 (1D ndarray): rate constants errors #1
        T2 (1D ndarray): temperatures #2
        k2 (1D ndarray): rate constants #2
        dk2 (1D ndarray): rate constants errors #2
        kacid (1D ndarray): 'pure' rate constants for acid form
        kanion (1D ndarray): 'pure' rate constants for anion form

    Returns:
        latex (str): latex table with all rate constants
    """
    d1, d2 = dict(), dict()

    for i, t in enumerate(T1):
        try:
            d1[t].append((k1[i], dk1[i]))
        except KeyError:
            d1[t] = [(k1[i], dk1[i])]

    for i, t in enumerate(T2):
        try:
            d2[t].append((k2[i], dk2[i]))
        except KeyError:
            d2[t] = [(k2[i], dk2[i])]

    strings = list()
    for i, t in enumerate(temps):
        if t in d1.keys():
            if len(d1[t]) > 1:
                a = '$ \\newline $'.join(stringify(*x) for x in d1[t])
            else:
                a = stringify(*d1[t][0])
        else:
            a = '~'

        if t in d2.keys():
            if len(d2[t]) > 1:
                b = '$ \\newline $'.join(stringify(*x) for x in d2[t])
            else:
                b = stringify(*d2[t][0])
        else:
            b = '~'

        strings.append('${}$ & ${}$ & ${}$ & ${:.2e}$ & ${:.2e}$ \\\\ \\hline'.format(t, a, b, kacid[i], kanion[i]))

    head = '''
    \\begin{center}
    \\begin{tabular}{| x{1cm} | x{3.5cm} | x{3.5cm} || x{3.5cm} | x{3.5cm} |}    \\hline
    T,\\,K &
    $k^\\text{exp}_\\text{pH='''+ str(pH1) +'''},\\, \\text{M}^{-1}\\text{s}^{-1}$  &
    $k^\\text{exp}_\\text{pH='''+ str(pH2) +'''},\\, \\text{M}^{-1}\\text{s}^{-1}$ &
    $k^\\text{calc}_\\text{acid},\\, \\text{M}^{-1}\\text{s}^{-1}$ &
    $k^\\text{calc}_\\text{anion},\\, \\text{M}^{-1}\\text{s}^{-1}$\\\\ \\hline
    '''

    feet = '''\n
    \\end{tabular}
    \\end{center}
    '''

    latex = head + ' \n '.join(strings) + feet

    for i in range(5, 16):
        latex = latex.replace('e+{:02d}'.format(i), ' \\cdot 10^{{{0}}} '.format(i))

    return latex

if __name__ == '__main__':
    if len(sys.argv) < 2:
        fn = 'data.txt'
    else:
        fn = sys.argv[-1]

    open('termodyn-acid.data.txt', 'w').close()     # clearing files in case they existed
    open('termodyn-anion.data.txt', 'w').close()
    for dataset in read_data_file(fn):

        name, acro, pKa, pH1, T1, k1, dk1, pH2, T2, k2, dk2 = parse(dataset)
        print acro,

        # calculating fractions based on pKa, putting them into a separate file
        latex_fractions = fractions(acro, name, pKa, pH1, pH2,
                                    draw_png=draw_png, draw_legend=draw_legend, draw_bw=draw_bw)
        with open(acro + '-fractions.tex', 'w') as f:
            f.write(latex_fractions)

        # checking for errors
        if not len(T1) == len(k1) == len(dk1):
            print acro, 'pH #1\n length(T) =', len(T1), '\nlength(k) =', len(k1), '\nlength(dk) =', len(dk1), '\nData must be of equal length'
            raise ValueError
        if not len(T2) == len(k2) == len(dk2):
            print acro, 'pH #2\n length(T) =', len(T2), '\nlength(k) =', len(k2), '\nlength(dk) =', len(dk2), '\nData must be of equal length'
            raise ValueError

        # calculating
        result1 = calculate_params(T1, k1)
        result1['ph'] = pH1

        result2 = calculate_params(T2, k2)
        result2['ph'] = pH2

        # generating 'pure' constants and equations, putting them into arrays
        latex_split, kacid, kanion = split(acro, pKa, pH1, pH2, result1, result2)
        with open(acro + '-split.tex', 'w') as f:
            f.write(latex_split)

        # preparing data for final table with all 'pure' arrhenius parameters for `termodyn.py`
        resacid = calculate_params(temps, kacid)
        resanion = calculate_params(temps, kanion)
        with open('termodyn-acid.data.txt', 'a') as f:
            f.write(acro)
            f.write(' {ea} {a} {s} {h} {g}\n'.format(**resacid))

        with open('termodyn-anion.data.txt', 'a') as f:
            f.write(acro)
            f.write(' {ea} {a} {s} {h} {g}\n'.format(**resanion))

        # calculating 'pure' arrhenius parameters, putting them into separate file
        latex_pure = pure_params(resacid, resanion)
        with open(acro + '-pure.tex', 'w') as f:
            f.write(latex_pure)

        # generating table with input data, putting them into separate file
        latex_table = create_table(pH1, pH2, T1, k1, dk1, T2, k2, dk2, kacid, kanion)
        with open(acro + '-table.tex', 'w') as f:
            f.write(latex_table)

        # plain text output
        with open(acro + '-data.txt', 'w') as f:
            f.write(name + '\n' + '='*40 + '\n\n')
            f.write(plain_output(result1))
            f.write('\n\n\n')
            f.write(plain_output(result2))

        # preparing latex parts for main document
        latex_params = latex_double_output_params(result1, result2)
        with open(acro + '-data.tex', 'w') as f:
            f.write(latex_params)

        # creating nifty plots
        draw_double_plot(T1, k1, dk1, result1, pH1, T2, k2, dk2, result2, pH2, kacid, kanion,
                         draw_png=draw_png, draw_legend=draw_legend, draw_bw=draw_bw)

        # creating SigmaPlot tables
        from math import log10

        with open(acro + '-sigma-table.txt', 'w') as f:
            f.write('T\t1/T\tlg k\tlg k -error\tlg k +error\tlg error')

            f.write('\n\npH = ' + str(pH1) + '\n\n')
            for (t, k, dk) in zip(T1, k1, dk1):
                f.write(('{:.6g}\t' * 6 + '\n').format(
                    t, 1./t, log10(k), log10(k-dk), log10(k+dk), log10(dk)
                    ))

            f.write('\n\npH = ' + str(pH2) + '\n\n')
            for (t, k, dk) in zip(T2, k2, dk2):
                f.write(('{:.6g}\t' * 6 + '\n').format(
                    t, 1./t, log10(k), log10(k-dk), log10(k+dk), log10(dk)
                    ))

        print '...done'
