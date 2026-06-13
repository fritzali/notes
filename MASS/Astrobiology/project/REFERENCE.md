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

<!-- ### Project Type

Master course project in Computational Astrobiology

 Duration: **3 months**

### Project title

#### Habitable-Zone Target Selection from Synthetic PLATO Observations Using PlatoSim/PLATOnium

### Project description

The aim of this project is to develop a computational pipeline for identifying and ranking potentially habitable exoplanet targets using **synthetic observations generated with the official PLATO mission simulator, PlatoSim, and its user-friendly workflow toolkit, PLATOnium**. PlatoSim is an end-to-end simulator designed to reproduce realistic PLATO observations, including stellar fields, spacecraft systematics, jitter, optics, CCD/electronics, and natural noise sources. PLATOnium provides the standard workflow for setting up and running realistic multi-camera, multi-quarter simulations.

The scientific motivation of the project is directly aligned with the goals of ESA’s PLATO mission, which is designed to detect and characterize terrestrial planets in orbits up to the habitable zone of Sun-like stars, while also characterizing the host stars through high-precision photometry and asteroseismology.

In this project, the student will use PlatoSim/PLATOnium to simulate a small, well-designed sample of planetary systems, recover transiting planets from the simulated data, compute habitability-related quantities such as stellar irradiation and equilibrium temperature, determine which planets lie in or near the habitable zone, and construct a transparent ranking model that identifies the most promising targets for follow-up. The project should explicitly distinguish between a **planet located in the habitable zone** and a **planet that is actually habitable**, since ESA notes that habitable-zone location alone does not guarantee habitability.

### Important technical note

This project must use **PlatoSim / PLATOnium as the main simulator**. Lightweight alternatives such as PSLS may be used only for testing or comparison, but not as the core simulation engine.

PlatoSim is currently available only to members of the **PLATO Mission Consortium** with a signed NDA, according to the official quickstart documentation. Therefore, this project assumes that access to the simulator, credentials, and required auxiliary files are available through the collaboration or through supervised access provided by the instructor.

### Main scientific question

#### Given realistic synthetic PLATO observations, which detected planets are the most promising habitable-zone targets, and how can they be ranked in a transparent, physically motivated, and reproducible way?

### Main objectives

By the end of the project, the student should be able to:

- use the official **PlatoSim / PLATOnium** workflow to generate synthetic PLATO observations
- simulate a small sample of planetary systems with different host-star and orbital properties
- extract or use the resulting synthetic light curves
- recover transiting planets from the simulated data
- compute stellar luminosity, incident flux, and equilibrium temperature
- identify planets located in or near the habitable zone
- rank the detected planets with a transparent habitability scoring model
- discuss the limits of using simple proxies for habitability

### Scope of the project

Because PlatoSim is a realistic end-to-end mission simulator, the project should focus on a **small but carefully chosen sample** rather than a very large population. A good scope for a 3-month course project is:

- **5 to 10 simulated planetary systems**
- at least one case clearly inside or near the habitable zone
- at least one hot inner-orbit control case
- at least one colder outer-orbit control case
- at least one more difficult case affected by stronger noise or less favorable detectability


This is fully sufficient for a strong course project.

### Simulation workflow

PLATOnium’s standard workflow consists of four main steps:

1. `picsim` — generate stellar catalogues
2. `varsim` — generate noiseless stellar/planetary light curves
3. `payload` — generate instrumental systematics
4. `platonium` — run the PlatoSim multi-camera, multi-quarter simulation

The student is expected to use this workflow, at least at a basic level, to set up and run a small number of simulations.

### Project tasks

#### Task 1 — Learn the simulator environment

The student should install or access PlatoSim, understand the required environment variables, test the quickstart example, and inspect the structure of the outputs. The official quickstart shows both command-line and Python usage and recommends PLATOnium for efficient project creation and execution.

#### Task 2 — Define a small synthetic sample

The student should design a small set of planetary systems with varying stellar and orbital parameters. The sample should be intentionally structured to include both promising and non-promising targets. The goal is not to simulate as many systems as possible, but to build a meaningful test sample for habitability-based target selection.

#### Task 3 — Run PlatoSim/PLATOnium simulations

The student should generate synthetic PLATO observations for the selected systems. These simulations should reflect realistic observing conditions rather than idealized toy models.

#### Task 4 — Recover transiting planets

The student should apply a transit search method such as **BLS** or **TLS** to the simulated light curves and compare recovered periods with the injected values. A correct and reproducible detection table is a key part of the project.

#### Task 5 — Compute physical habitability indicators

For each recovered planet candidate, the student should compute:

$\frac{L_*}{L_\odot}=\left(\frac{R_*}{R_\odot}\right)^2\left(\frac{T_*}{T_\odot}\right)^4$

$S=\frac{L_*/L_\odot}{(a/\mathrm{AU})^2}$

$T_{\mathrm{eq}} = T_* \sqrt{\frac{R_*}{2a}}(1-A)^{1/4}$

where $L_∗$ is stellar luminosity, $S$ is incident stellar flux in Earth units, and $T_{\mathrm{eq}}$ is the planetary equilibrium temperature. The student may lso discuss the simplified zero-albedo form if that is the version emphasized in the course.

#### Task 6 — Identify habitable-zone planets

The student should define a working habitable-zone criterion. A practical choice is to use stellar flux as the primary indicator and equilibrium temperature as a secondary diagnostic. The student should clearly explain that these are approximate screening tools, not proof of actual habitability.

#### Task 7 — Construct a habitability ranking model

The student should build a transparent scoring model to rank the recovered planets according to their potential as habitable-zone targets. The score should combine physical relevance and detection reliability.

#### Task 8 — Scientific interpretation

The student should discuss which systems are the most promising, which are less promising, and why. The report should include a clear discussion of limitations, especially the difference between “planet in the habitable zone” and “planet truly habitable.”

### Recommended physics and scoring framework

A reasonable project-level scoring model is:

$H_{\mathrm{HZ}} + 0.2H_{R_p} + 0.15H_{\mathrm{host}} + 0.15H_{\mathrm{det}} + 0.1H_{\mathrm{stab}}$

where:

- $H_{\mathrm{HZ}}$ measures how close the planet is to Earth-like stellar irradiation
- $H_{R_p}$ favors Earth-sized or super-Earth planets
- $tH_{\mathrm{host}}$ rewards suitable host stars
- $tH_{\mathrm{det}}$ measures robustness of detection
- $H_{\mathrm{stab}}​$ is an optional stability term, for example related to orbital eccentricity if available


One possible mathematical form is:
$]H_{\mathrm{HZ}} = \exp\left[-\frac{(S-1)^2}{2\sigma_S^2}\right]$

$H_{R_p} = \exp\left[-\frac{(R_p-1.2)^2}{2\sigma_R^2}\right]$

This is only one suggested model. The student may refine it, but the final score must remain physically interpretable.

### Expected outputs

By the end of the project, the student should produce:

- a small set of PlatoSim/PLATOnium simulation products
- a clean table of injected system parameters
- a table of recovered transit candidates
- computed values of L∗ S, and Teq
- habitable-zone classification flags
- a final ranked list of targets
- several scientific plots
- a short report and oral presentation


### Suggested plots

The final report should include at least the following:
- examples of simulated PLATO light curves
- injected period vs recovered period
- incident flux vs planet radius
- histogram of equilibrium temperatures
- bar chart or table of top-ranked targets



### Required deliverables

The final submission should contain:
- **Code repository or structured project folder**
- **Simulation configuration and workflow notes**
- **Detection and analysis scripts or notebooks**
- **Final results table**
- **Scientific report** of approximately 8–12 pages
- **Presentation** of approximately 20 minutes
-->
