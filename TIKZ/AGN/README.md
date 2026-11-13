## Unified Active Galactic Nucleus Model

With inspiration taken from:

### Compilation

For a single figure, run:

<pre>pdflatex <i>file</i>.tex</pre>

To compile all serialized or in parallel:

<pre>for f in *.tex; do pdflatex -interaction=nonstopmode "$f"; done</pre>

<pre>ls *.tex | xargs -P4 -n1 pdflatex -interaction=nonstopmode >/dev/null</pre>

And to clean up any auxiliary files afterwards:

<pre>rm -f *.aux *.log</pre>
