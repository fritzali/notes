## Interstellar Turbulence Forcing

### Contents

- Drivers & Stochastic Forcing: [greyscale](https://cdn.jsdelivr.net/gh/fritzali/notes@main/TIKZ/OU/greyscaleOU.pdf) [color](https://cdn.jsdelivr.net/gh/fritzali/notes@main/TIKZ/OU/colorOU.pdf)

### Compilation

For a single figure, run:

<pre>pdflatex <i>file</i>.tex</pre>

To compile all serialized or in parallel:

<pre>for f in *.tex; do pdflatex -interaction=nonstopmode "$f"; done</pre>

<pre>ls *.tex | xargs -P4 -n1 pdflatex -interaction=nonstopmode >/dev/null</pre>

And to clean up any auxiliary files afterwards:

<pre>rm -f *.aux *.log</pre>

### Inspiration

- [C. Federrath, J. Roman-Duval, R. S. Klessen, W. Schmidt & M.-M. Mac Low 2010](https://doi.org/10.1051/0004-6361/200912437)
- [P. Grete, B. W. O'Shea & K. Beckwith 2018](https://doi.org/10.3847/2041-8213/aac0f5)
- [V. Eswaran & S. B. Pope 1988](https://doi.org/10.1016/0045-7930(88)90013-8)
