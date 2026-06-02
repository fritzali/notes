all: build/plot build/main build/flat

build:
	mkdir -p build

clean:
	rm -rf build

count:
	texcount *.tex

spell:
	cat *.tex | \
	aspell list --mode=tex --master=en_US | \
	aspell list --mode=tex --master=de_DE | sort | uniq

FORCE:

.PHONY: all clean

build/main: FORCE | build
	latexmk \
		--lualatex \
		--output-directory=build \
		--interaction=batchmode \
		--halt-on-error \
	main.tex

build/flat: FORCE | build
	latexmk \
		--lualatex \
		--output-directory=build \
		--interaction=batchmode \
		--halt-on-error \
	flat.tex

build/plot: plot.py matplotlibrc matplotlib.tex | build
	TEXINPUTS=$$(pwd)/.: MATPLOTLIBRC=./matplotlibrc python plot.py
