# Leipzig report

This tool helps to visualize and understand the data, which was obtained during the experiments at Laser Labor at TROPOS in Leipzig.

### Dependencies

* Latex with a number of packages
* Python 2
	- numpy
	- scipy
	- matplotlib

### What it needs

All the original data goes to `data.txt` or how would you like to call it in `arrhenuis/` folder (NB! it is `arrhenUIs`). Name of this file will be given as argument to `arrhenius.py`, who will do almost all the work. Passing the name can be omitted, then it will be default to `data.txt`.

Data format and limitations are explained in `data.txt`.

### Output

For every compound mentioned in `data.txt` following LaTeX parts for the resulting document will be produced. There `*` is acronym for any compound, given in `data.txt`.

* `*-table.tex` A table with original data for both pH and already calculated rate constants for acid and anion forms, this is based on given pKa.
* `*-fractions.tex` A plot of fractions of acid and anion parts at different pH. Also a table with two fractions at used pH is provided, which can be omitted (see `arrhenius.py`, function `fractions()`).
* `*-data.tex` Heart of the tool. A table with Arrhenius  and thermodynamic data for both pH, calculated with another reincarnation of Prof. Dr. H. Herrmann's code.
* `*-split.tex` A set of two equations for acid and anion rate constants, calculated based on pH and Arrhenius data.
* `*-pure.tex` Arrhenius and thermodynamic data for already splitted rate constants.

Also a nice Arrhenius plot with original and calculated data is generated, named `*.double.pdf`. (Generating `png` image is based on bool flag, see `arrhenius.py`, functions `draw_plot()` and `draw_double_plot()`)

You can generate fractions plots for using `exceptions.py`. It was created for compounds with very low rate constant and provides only `*-fractions.tex`.

Final table with various termodynamic data, `termodyn.tex`, is generated by `termodyn.py`, which uses data, generated by `arrhenius.py`, `termodyn-acid.data.txt` and `termodyn-anion.data.txt`.

----

Nothing in `main-leipzig.tex` is generated automatically, it just uses `arrhenius.py`'s outputs, because you never know what data you will have to present in addition to mentioned above and what text you will write for all the numbers and graphs.

### Launching

#### Windows

Just launch `compile.cmd` in the root directory

To delete old files launch `clean.cmd`

#### Linux

`make all` to compile everything

`make clean` to delete temporary files

### License

MIT License, except for `hh.py`, which is published under separate license.

### TODO

 - [x] describe `exceptions.py`
 - [x] describe `termodyn.py`
 - [x] ask HH about license
 - [ ] receive answer
 - [x] code commentary
 - [x] license
 - [ ] rewrite linear regression using scipy / OLS from MR
