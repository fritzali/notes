## Circumgalactic Medium

### Contents

- Inflows & Outflows: [greyscale](https://cdn.jsdelivr.net/gh/fritzali/notes@main/TIKZ/CGM/greyscaleCGM.pdf) [color](https://cdn.jsdelivr.net/gh/fritzali/notes@main/TIKZ/CGM/colorCGM.pdf)

### Compilation

For a single figure, run:

<pre>pdflatex <i>file</i>.tex</pre>

To compile all serialized or in parallel:

<pre>for f in *.tex; do pdflatex -interaction=nonstopmode "$f"; done</pre>

<pre>ls *.tex | xargs -P4 -n1 pdflatex -interaction=nonstopmode >/dev/null</pre>

And to clean up any auxiliary files afterwards:

<pre>rm -f *.aux *.log</pre>

### Inspiration

- [J. Tumlinson, M. S. Peeples & J. K. Werk 2017](https://doi.org/10.1146/annurev-astro-091916-055240)
