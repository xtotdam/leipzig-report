.PHONY: all clean gendata purge

all:
	cd arrhenuis/ && make all

	cd ../
	pdflatex main-leipzig.tex
	pdflatex main-leipzig.tex
	# pdflatex main-leipzig.tex

latex:
	pdflatex main-leipzig.tex
	pdflatex main-leipzig.tex

clean:
	rm -vf *.aux
	rm -vf *.log
	rm -vf *.out
	rm -vf *.synctex.gz
	rm -vf *.toc

	cd arrhenuis/ && make clean

gendata:
	cd arrhenuis/ && make gendata

purge: clean
	cd arrhenuis/ && make purge
