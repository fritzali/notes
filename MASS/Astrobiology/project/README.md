# Exoplanet Habitability with PLATO: An End-to-End Machine Learning Pipeline for Exoplanet Detection and Habitability Classification

This project implements a fully integrated, modular data processing and machine learning pipeline engineered to simulate, clean, validate, and analyze planetary transit signatures from high-precision stellar time-series data. Designed around the instrument capabilities and scientific goals of the European Space Agency's (ESA) upcoming **PLATO** mission, this repository provides researchers with tools to transition from raw stellar parameters to robust exoplanet habitability classifications. By bridging the gap between forward-modeled astrophysics and modern ensemble/deep-learning classification techniques, the pipeline automates the identification of worlds capable of retaining atmospheres and supporting liquid surface water.

---

## 1. Project Overview & Core Concepts

To understand the architecture of this repository, it helps to break down the core astronomical concepts and mission architectures that form its foundation:

* **Exoplanets and Planetary Transits:** Exoplanets are planets orbiting stars outside our solar system. The primary method used here to find them is *transit photometry*. When an exoplanet’s orbit aligns with our line of sight, it periodically passes in front of its host star, blocking a tiny fraction of starlight. This creates a characteristic dip in the observed brightness over time, known as a **light curve**.
* **The PLATO Mission:** Scheduled for launch by ESA, the **PLATO** (PLAnetary Transits and Oscillations of stars) mission is designed to discover and characterize rocky exoplanets orbiting solar-type stars, specifically targeting those within the habitable zone. PLATO utilizes an innovative multi-camera array to achieve unprecedented photometric precision and long, uninterrupted observation baselines.
* **PSLS (Plato Solar-like Light-curve Simulator):** To prepare data channels ahead of the mission launch, the astronomical community relies on high-fidelity synthetic data. PSLS generates ultra-realistic light curves that mimic not only the planetary transit signals but also instrumental noise, photon noise, and stochastic stellar phenomena (such as stellar oscillations and granulation).
* **Exoplanet Habitability:** Determining habitability is a multi-dimensional problem. Beyond merely finding a planet of the right size, this project uses structural parameters and radiation metrics to assess whether an exoplanet can sustain life. We analyze the balance between stable thermal environments and volatile-stripping radiation events to assign a conclusive habitability categorization.

---

## 2. Installation & Repository Setup

### 2.1. Virtual Environment & Dependencies

To ensure isolated dependency resolution and prevent library conflicts, it is highly recommended to deploy a local Python virtual environment.

```bash
# 1. Clone the repository and navigate to its root directory
git clone https://github.com/your-username/astrobio.git
cd astrobio

# 2. Create a localized Python virtual environment (.venv)
python -m venv .venv

# 3. Activate the virtual environment
# On Linux / macOS:
source .venv/bin/activate
# On Windows (Command Prompt):
.venv\Scripts\activate
# On Windows (PowerShell):
.venv\Scripts\Activate.ps1

# 4. Upgrade pip and install the required core packages
pip install --upgrade pip
pip install -r requirements.txt

```

### 2.2. Tracking Large Data Configurations via Git LFS

Because this project interfaces with raw stellar model grids, synthetic light-curve FITS files, and compiled machine learning weights, **Git Large File Storage (Git LFS)** must be initialized. This prevents local repository bloating and guarantees smooth version control histories.

```bash
# Initialize the Git LFS extension on your system
git lfs install

# Configure tracking profiles for large data elements and model checkpoints
git lfs track "*.grid"
git lfs track "*.fits"
git lfs track "*.pkl"
git lfs track "*.csv"
git lfs track "*.h5"

# Ensure tracking attributes are correctly locked into version history
git add .gitattributes
git commit -m "infrastructure: initialize Git LFS data tracking parameters"

```

### 2.3. Data Provenance & Framework Adaptations

* **Stellar Model Grid Acknowledgments:** The underlying physical stellar model grid utilized by this pipeline contains foundational parameters (mass, radius, metallicity, and age) obtained directly via private communication with **Dr. Reza Samadi** (LESIA, Observatoire de Paris). Any external redistribution of these files requires prior explicit authorization.
* **PSLS Customizations:** The default repository version of the **PSLS source code** was debugged, refactored, and directly modified for this project. These modifications enable multi-threaded high-throughput batch execution and allow the simulator to extract intermediate stellar flare energy states and stochastic noise variance directly into the downstream classification pipeline.

---

## 3. Technical Pipeline & Code Workflow

The project's architectural execution flows sequentially through five distinct Jupyter notebooks. Each module represents a core layer of the data-to-prediction pipeline:

```
[ P1explore ] ──> [ P2generate ] ──> [ P3sanitize ] ──> [ P4inspect ] ──> [ P5train ]

```

### 3.1. Explore (`P1explore.ipynb`)

The pipeline begins with Exploratory Data Analysis (EDA) of the underlying physical parameter distributions:

* Parses the custom stellar model grid files provided by Dr. Reza Samadi.
* Maps continuous parameter correlations (e.g., $T_{\text{eff}}$ vs. $\log g$) to build bounds for realistic synthetic systems.
* Calculates baseline theoretical signal-to-noise ratios (SNR) based on different PLATO instrument configurations, factoring in active camera alignments (e.g., 6, 12, 18, or 24 active telescopes).

### 3.2. Generate (`P2generate.ipynb`)

This module functions as our data synthesis engine, utilizing our optimized PSLS variant:

* Executes random grid-sampling algorithms across verified stellar parameter bounds.
* Injects analytical planetary transits using Mandel-Agol limb-darkened geometry models.
* Superimposes multi-component stellar noise (granulation profiles, convective background shifts, and solar-like p-mode oscillations).
* Injects random, high-energy **stellar flares** drawn from power-law distributions. This provides a realistic radiation background crucial for habitability calculations.

### 3.3. Sanitize (`P3sanitize.ipynb`)

Raw stellar time-series are dominated by low-frequency stellar activity and systematic instrumental drift. This notebook isolates high-frequency transit signals using a two-stage filter:

* **Wotan Detrending:** Employs the `wotan` optimization framework to apply time-windowed robust high-pass filters (such as the biweight or Huber estimators). This effectively flattens long-term stellar variability while preserving the sharp, steep boundaries of individual transit events.
* **Box Least Squares (BLS) Periodogram:** Processes the detrended light curves through a periodic grid search using the Astropy BLS algorithm. This extracts key observational features: orbital period ($P$), transit depth ($\delta$), transit duration ($\tau$), and the transit epoch midpoint ($T_0$).

### 3.4. Inspect (`P4inspect.ipynb`)

A structural validation layer designed to eliminate false positives and catch numerical processing issues:

* Generates multi-panel diagnostic plots for each planet candidate, displaying the raw light curve, the flattened detrended curve, the Fourier power spectrum (FFT) for stellar oscillation diagnostics, and the phase-folded transit profile.
* Fits an idealized inverted top-hat model back to the phase-folded data to evaluate goodness-of-fit metrics, filtering out numerical anomalies or severe stellar activity masquerading as planetary transits.

### 3.5. Train (`P5train.ipynb`)

The final phase engineers a robust feature matrix from the physical properties discovered in the previous steps and feeds it into a machine learning classification pipeline.

```
                  ┌──────────────────────┐
                  │ Engineered Features: │
                  │  - Planet Size       │
                  │  - Equilibrium Temp  │
                  │  - Flare Radiation   │
                  └──────────┬───────────┘
                             │
            ┌────────────────┴────────────────┐
            ▼                                 ▼
┌───────────────────────┐         ┌───────────────────────┐
│     Random Forest     │         │   Multi-Layer Neural  │
│ Classifier/Regressor  │         │    Network (Sklearn)  │
└───────────┬───────────┘         └───────────┬───────────┘
            │                                 │
            └────────────────┬────────────────┘
                             ▼
              ┌─────────────────────────────┐
              │ Output Habitability Score   │
              │ (Gaseous vs. Rocky Habitable)│
              └─────────────────────────────┘

```

#### Model Architecture Comparisons

This module implements a competitive benchmarking pipeline that evaluates two separate machine learning frameworks:

1. **Ensemble Random Forest Pipeline:** Evaluates non-linear feature splits via a grid-search optimized `RandomForestClassifier` and `RandomForestRegressor`. This provides feature-importance rankings that help verify which physical traits impact the final classification decision the most.
2. **Multi-Layer Perceptron (MLP) Neural Network:** Implements a dense, fully connected neural network via scikit-learn's `MLPClassifier` (configured with deep structural topologies like `(128, 64, 32)` nodes, ReLU activations, automated early stopping, and an Adam optimization back-end). This maps high-dimensional interactions between radiation profiles and planetary size.

#### Multi-Tiered Habitability Classification Criteria

The models assess habitability by evaluating three core physical constraints:

* **Planet Size & Composition:** Derived from the transit depth ($\delta$) combined with the stellar radius ($R_*$) to yield the precise planetary radius ($R_p$). This stage explicitly separates solid, high-density rocky terrestrial candidates ($R_p \le 1.5 \, R_\oplus$) from volatile-rich, low-density mini-Neptunes or gas giants.
* **Thermodynamic Environment:** Computes the semi-major axis ($a$) via Kepler's Third Law to establish the incident stellar flux and the planet's equilibrium temperature ($T_{\text{eq}}$). This maps the planet's position relative to its host star's conservative and optimistic **Habitable Zone (HZ)** boundaries.
* **Radiation & Flare Exposure Profiles:** Quantifies the cumulative high-energy impact of stellar flares against the background stellar irradiance. The models use this metric to determine if an otherwise thermally stable world is at high risk for atmospheric stripping, ozone depletion, or surface sterilization.

---

## 4. References & Bibliography

1. **Rauer, H., Catala, C., et al.** "The PLATO mission", *Experimental Astronomy*, Vol. 38, Issue 1-2, pp. 249-330 (2014).
ADS Link: [https://ui.adsabs.harvard.edu/abs/2014ExA....38..249R](https://www.google.com/search?q=https://ui.adsabs.harvard.edu/abs/2014ExA....38..249R)
2. **Samadi, R., Belkacem, K., et al.** "The PLATO Solar-like Light-curve Simulator — A tool to generate realistic stellar light-curves with instrumental effects representative of the PLATO mission", *Astronomy & Astrophysics*, Vol. 624, id.A117, 18 pp. (2019).
ADS Link: [https://ui.adsabs.harvard.edu/abs/2019A%26A...624A.117S](https://www.google.com/search?q=https://ui.adsabs.harvard.edu/abs/2019A%2526A...624A.117S)
3. **Hippke, M., David, T. J., et al.** "Wōtan: Comprehensive time-series detrending in Python", *The Astronomical Journal*, Vol. 158, Issue 4, id.143, 21 pp. (2019).
ADS Link: [https://ui.adsabs.harvard.edu/abs/2019AJ....158..143H](https://www.google.com/search?q=https://ui.adsabs.harvard.edu/abs/2019AJ....158..143H)
4. **Kovács, G., Zucker, S., & Mazeh, T.** "A box-fitting algorithm in the search for periodic transits", *Astronomy & Astrophysics*, Vol. 391, pp. 369-377 (2002).
ADS Link: [https://ui.adsabs.harvard.edu/abs/2002A%26A...391..369K](https://www.google.com/search?q=https://ui.adsabs.harvard.edu/abs/2002A%2526A...391..369K)
