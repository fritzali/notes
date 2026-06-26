<!-- ## Project Outline — Habitable-Zone Target Selection Using the PLATO Mission Simulator -->

## Tools

### PSLS

The [PLATO Solarlike Lightcurve Simulator](https://sites.lesia.obspm.fr/psls/) is intended for the creation of synthetic stellar lightcurves which model those that will be produced by the PLATO mission once it is operating. It uses mathematical parametrizations for the impacts of different noise sources to enable the fast computation of integrated flux signals at low cost, thereby allowing researchers the creation of large datasets for statistical testing. 

#### Setup

Download the source distribution archive from [PyPI](https://pypi.org/project/psls/#files) and unpack it inside its designated [PSLS](PSLS) directory:

<pre>tar -xvzf psls-<i>version</i>.tar.gz --strip-components=1</pre>

Make the `psls.py` script executable:

<pre>chmod +x psls.py</pre>

After some fixes that are already included in this repository, test the installation via the `psls.yaml` main sequence basic example by running it in an active environment that includes all required packages from its root directory:

<pre>./psls.py -o data --extended-plots -V examples/psls.yaml</pre>

#### Usage

A more comprehensive list of user flags and configuration variables is included with the [PSLS](PSLS) directory. The following is an overview of the most important cases:

### PlatoSim

The [PLATO Simulator](https://ivs-kuleuven.github.io/PlatoSim3/) implements a more granular approach, simulating pixel level data and therefore represents a realistic framework to study the detector itself. Access is public only within the PLATO consortium and requires credentials to gain access to the software.

### Git LFS

The [Git Large File Storage](https://git-lfs.com/) utility replaces specified files with a text pointer to their remote location and manages their versioning like standard Git. Install its package:

<pre>sudo pacman -S git-lfs</pre>

Once installed, initialize it for each user account:

<pre>git lfs install</pre>

Track a specific file or use a wildcard string for a file type:

<pre>git lfs track "<i>file</i>"</pre><pre>git lfs track "*.<i>type</i>"</pre>

This will add a line in `.gitattributes` which needs to be tracked by Git. Other useful commands:

<pre>git lfs status</pre><pre>git lfs track</pre><pre>git lfs pull</pre><pre>git lfs prune</pre><pre>git lfs migrate info</pre>
