all: build/plot build/main build/flat

build:
	mkdir -p build

clean:
	rm -rf build

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
