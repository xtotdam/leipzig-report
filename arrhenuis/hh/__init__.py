import numpy as np

def calculate_params(T, k, student=True):
    """
    Calculates arrhenius and thermodynamic data from temperature and rate constants
    Just translated from Prof. Herrmann's pascal code into python

    Args:
        T (1D ndarray): temperatures
        k (1D ndarray): rate constants
        student (bool, optional): whether we use Student coefficient or not

    Returns:
        dict: all the data
    """
    x = 1. / T
    y = np.log10(k)

    anzahl = x.shape[0]

    txi = np.sum(x)
    tyi = np.sum(y)
    xi_quadrat = np.sum(x**2)
    xiyi = np.sum(x * y)
    yi_quadrat = np.sum(y**2)

    x_quer = np.mean(x)
    y_quer = np.mean(y)

    sxx = xi_quadrat - (x_quer * txi)
    sxy = xiyi - (x_quer * tyi)
    syy = yi_quadrat - (y_quer * tyi)

    steigung = sxy / sxx
    b_strich = sxy / syy
    b_gross = steigung * b_strich

    korrelationskoeffizient = np.sqrt(abs(b_gross))

    abschnitt = y_quer - (steigung * x_quer)

    bb = steigung * steigung

    s_quadrat = (syy - bb * sxx) / (anzahl - 2)

    v_y_quer = s_quadrat / anzahl

    varianz_b = s_quadrat / sxx

    varianz_a = v_y_quer + (x_quer * x_quer * varianz_b)

    if student:
        if anzahl < 3: t_strich = 1e9
        elif anzahl == 3:  t_strich = 12.7
        elif anzahl == 4:  t_strich = 4.3
        elif anzahl == 5:  t_strich = 3.2
        elif anzahl == 6:  t_strich = 2.8
        elif anzahl == 7:  t_strich = 2.6
        elif anzahl == 8:  t_strich = 2.4
        elif anzahl == 9:  t_strich = 2.35
        elif anzahl == 10: t_strich = 2.23
        else: t_strich = 2.2
    else:
        t_strich = 1.

    steigungsfehler = t_strich * np.sqrt(abs(varianz_b))
    abschnittsfehler = t_strich * np.sqrt(abs(varianz_a))

    rel_fehler_abschnitt = abs(100 * (abschnittsfehler / abschnitt))
    rel_fehler_steigung = abs(100 * (steigungsfehler / steigung))

    ea = -2.3026 * 8.314 * steigung
    eafehler = 2.3026 * 8.314 * steigungsfehler
    a = np.exp(2.3026 * abschnitt)
    afehler = a * rel_fehler_abschnitt / 100.

    h = ea - (8.314 * 298.15)
    s = 8.314 * (np.log(a) - 30.4576)
    g = h - (298.15 * s)

    hfehler = rel_fehler_steigung * h / 100.
    sfehler = rel_fehler_abschnitt * s / 100.
    gfehler = (rel_fehler_steigung + rel_fehler_abschnitt) * g / 100.

    result = dict()
    result['abschnitt'] = abschnitt
    result['abschnittsfehler'] = abschnittsfehler
    result['rel_fehler_abschnitt'] = rel_fehler_abschnitt
    result['steigung'] = steigung
    result['steigungsfehler'] = steigungsfehler
    result['rel_fehler_steigung'] = rel_fehler_steigung
    result['korrelationskoeffizient'] = korrelationskoeffizient
    result['ea'] = ea / 1000.
    result['eafehler'] = abs(eafehler / 1000.)
    result['a'] = a
    result['afehler'] = abs(afehler)
    result['h'] = h / 1000.
    result['hfehler'] = abs(hfehler / 1000.)
    result['s'] = s
    result['sfehler'] = abs(sfehler)
    result['g'] = g / 1000.
    result['gfehler'] = abs(gfehler / 1000.)

    return result