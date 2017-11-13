.PHONY: all clean gendata purge

all:
	python arrhenius.py data.txt
	python exceptions.py
	python termodyn.py

clean:
	rm -vf *.tex
	rm -vf *.png
	rm -vf *.pdf
	rm -vf *-data.txt
	rm -vf *-fractions.tex
	rm -vf *-sigma-table.txt
	rm -vf *.fractions.pdf
	rm -vf *.pyc
	rm -vf termodyn*.txt

	cd hh/ && make clean

gendata:
	python junk_data_creator.py

purge: clean
	rm -vf data.txt
