all: plot build/main build/flat

slides: build/main build/flat

plot: plot.py
	python plot.py

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
