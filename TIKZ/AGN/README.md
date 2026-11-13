## Unified Active Galactic Nucleus Model

With inspiration taken from:

### Contents

- Model Sketch: [greyscale](https://cdn.jsdelivr.net/gh/fritzali/notes@main/TIKZ/AGN/greyscaleAGN.pdf) [color](https://cdn.jsdelivr.net/gh/fritzali/notes@main/TIKZ/AGN/colorAGN.pdf)
- Emission Bands: [greyscale](https://cdn.jsdelivr.net/gh/fritzali/notes@main/TIKZ/AGN/greyscaleAGNbands.pdf) [color](https://cdn.jsdelivr.net/gh/fritzali/notes@main/TIKZ/AGN/colorAGNbands.pdf)

### Compilation

For a single figure, run:

<pre>pdflatex <i>file</i>.tex</pre>

To compile all serialized or in parallel:

<pre>for f in *.tex; do pdflatex -interaction=nonstopmode "$f"; done</pre>

<pre>ls *.tex | xargs -P4 -n1 pdflatex -interaction=nonstopmode >/dev/null</pre>

And to clean up any auxiliary files afterwards:

<pre>rm -f *.aux *.log</pre>
