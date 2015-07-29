all: all.pdf

all.pdf: all.tex ./tex/*.tex ./conversions.sh

# 	Make the tikz directory if it doesn't exist:
	mkdir -p tikz

# 	Run the LaTeX script:
	latex -shell-escape all.tex
	dvipdf all.dvi
	rm -f *.bcf *.xml *.aux *.auxlock *.log *.out *.dvi

# 	Put all the EPS files in their appropriate place:
	cp tikz/all-figure0.eps ../gnuastro-figures/iandtime.eps
	cp tikz/all-figure1.eps ../gnuastro-figures/samplingfreq.eps

#	Make all the conversions:
	./conversions.sh ../gnuastro-figures/
