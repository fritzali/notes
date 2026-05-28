all: plot report count spell

plot: plot.py matplotlibrc matplotlib.tex
	TEXINPUTS=$$(pwd)/.: MATPLOTLIBRC=./matplotlibrc python plot.py

report: FORCE | build
	TEXINPUTS=.: \
	BIBINPUTS=.: \
	max_print_line=1048576 \
	latexmk \
		--lualatex \
		--output-directory=build \
		--interaction=batchmode \
		--halt-on-error \
	main.tex; cp build/main.pdf report.pdf

count:
	texcount chapter/*.tex *.tex

spell:
	cat chapter/*.tex *.tex | \
	aspell list --mode=tex --master=en_US | \
	aspell list --mode=tex --master=de_DE | sort | uniq

build:
	mkdir -p build

clean:
	rm -rf build

diagnose: FORCE | build
	TEXINPUTS=.: \
	BIBINPUTS=.: \
	max_print_line=1048576 \
	latexmk \
		-f \
		--lualatex \
		--output-directory=build \
		--interaction=batchmode \
		--halt-on-error \
	main.tex; cp build/main.pdf report.pdf

FORCE:

.PHONY: all clean
