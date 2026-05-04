## Cloudy – Spectral Synthesis Code for Astrophysical Plasmas<sup></sup>

### Installation

#### Download

Use version `c25.00` for cloning the current release:

<pre>git clone https://gitlab.nublado.org/cloudy/cloudy.git -b <i>version</i> Cloudy</pre>

#### Compile

Use an appropriate number of cores:

<pre>cd Cloudy/source</pre>

<pre>make -j <i>N</i></pre>

#### <i>Smoke</i> Test

Use the following to verify basic functionality:

<pre>echo test > test.in</pre>

<pre>./cloudy.exe -r test</pre>

<pre>vim test.out</pre>

### Documentation

#### <i>Hazy</i>

Build the `LaTeX` manual and clean auxiliary files:

<pre>cd ../docs/latex</pre>

<pre>perl CompileAll.pl</pre>

<pre>perl cleanAll.pl</pre>
