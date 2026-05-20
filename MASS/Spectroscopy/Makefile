all: report count spell

report: FORCE | build
	TEXINPUTS=.: \
	BIBINPUTS=.: \
	max_print_line=1048576 \
	latexmk \
		--lualatex \
		--output-directory=build \
		--interaction=batchmode \
		--halt-on-error \
	main.tex; mv build/main.pdf report.pdf

count:
	texcount main.tex

spell:
	cat main.tex | \
	aspell list --mode=tex --master=en_US | \
	aspell list --mode=tex --master=de_DE | sort | uniq

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean
