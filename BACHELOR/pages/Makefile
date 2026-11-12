all: report

report: FORCE | build
	TEXINPUTS=.: \
	BIBINPUTS=.: \
	max_print_line=1048576 \
	latexmk \
		--lualatex \
		--output-directory=build \
		--interaction=batchmode \
		--halt-on-error \
	main.tex; cp build/main.pdf revision.pdf

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean
