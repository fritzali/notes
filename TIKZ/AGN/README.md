## Unified AGN Model

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

### Inspiration

- [M. L. Menzel, V. Beckmann & C. R. Shrader 2012](https://doi.org/10.22323/1.176.0069)
- [J. Thorne, A. Robotham, L. Davies & S. Bellstedt 2022](https://doi.org/10.5281/zenodo.6381013)
- [C. M. Urry & P. Padovani 1995](https://doi.org/10.1086/133630)
- [E. Alexander 2022](https://emmaalexander.github.io/resources.html)
