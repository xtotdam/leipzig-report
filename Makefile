.PHONY: all clean

all:
	cd arrhenuis/ && make all

	cd ../
	pdflatex main-leipzig.tex

clean:
	rm -vf *.aux
	rm -vf *.log
	rm -vf *.out
	rm -vf *.synctex.gz
	rm -vf *.toc

	cd arrhenuis/ && make clean
