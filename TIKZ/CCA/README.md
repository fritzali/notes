## Chaotic Cold Accretion

### Contents

- Feeding & Feedback Cycle: [greyscale](https://cdn.jsdelivr.net/gh/fritzali/notes@main/TIKZ/CCA/greyscaleCCA.pdf) [color](https://cdn.jsdelivr.net/gh/fritzali/notes@main/TIKZ/CCA/colorCCA.pdf)

### Compilation

For a single figure, run:

<pre>pdflatex <i>file</i>.tex</pre>

To compile all serialized or in parallel:

<pre>for f in *.tex; do pdflatex -interaction=nonstopmode "$f"; done</pre>

<pre>ls *.tex | xargs -P4 -n1 pdflatex -interaction=nonstopmode >/dev/null</pre>

And to clean up any auxiliary files afterwards:

<pre>rm -f *.aux *.log</pre>

### Inspiration

- [M. Gaspari, F. Tombesi & M. Cappi 2020](https://doi.org/10.1038/s41550-019-0970-1)
- [M. Gaspari & A. Sądowski 2017](https://doi.org/10.3847/1538-4357/aa61a3)
