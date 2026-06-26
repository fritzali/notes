## Computational Astrobiology – Project: Exoplanet Habitability with PLATO using PSLS

This project implements a modular data processing and machine learning pipeline to simulate, clean, validate, and analyze planetary transit signatures from **PSLS** stellar timeseries data. Designed around the instrument capabilities and scientific goals of the upcoming **PLATO** mission, this repository aims to provide tools to transition from raw stellar parameters to sensible exoplanet habitability classifications. By bridging the gap from forward modeled astrophysics to modern ensemble classification and regression techniques, the pipeline automates the identification of worlds capable of retaining atmospheres and supporting liquid surface water.

---

### 1. Project Overview & Core Concepts

To understand the architecture of this repository, it helps to break down the core astronomical concepts and mission architectures that form its foundation:

* **Exoplanets and Planetary Transits:** Exoplanets are planets orbiting stars outside our solar system. The primary method used here to find them is transit photometry. When an exoplanet orbit aligns with our line of sight, it periodically passes in front of its host star, blocking a tiny fraction of starlight. This creates a characteristic dip in the observed brightness over time, known as a lightcurve.
* **PLATO:** Studying planetary transits and oscillations of stars, this mission is designed to discover and characterize rocky exoplanets orbiting solar type stars, specifically targeting those within the habitable zone. It utilizes an innovative multi camera array to achieve unprecedented photometric precision and long, uninterrupted observation baselines.
* **PSLS:** To prepare data channels ahead of the mission launch, the astronomical community relies on high fidelity synthetic data. This tool generates realistic light curves that mimic not only the planetary transit signals but also instrumental noise, photon noise, as well as stochastic stellar phenomena such as stellar oscillations and granulation.
* **Exoplanet Habitability:** Determining habitability is a multi dimensional problem. Beyond merely finding a planet of the right size, this project uses structural parameters and radiation metrics to assess whether an exoplanet can sustain life, attempting to analyze the balance between stable thermal environments and volatile stripping radiation events to assign a conclusive habitability categorization.

---

### 2. Installation & Repository Setup

#### 2.1. Virtual Environment & Dependencies

To ensure isolated dependency resolution and prevent library conflicts, it is highly recommended to deploy a local Python virtual environment.

```sh
# 1. clone the repository and navigate to its root directory
git clone https://github.com/your-username/astrobio.git
cd astrobio

# 2. create a localized Python virtual environment
python -m venv .venv

# 3. activate the virtual environment
source .venv/bin/activate

# 4. upgrade pip and install the required core packages
pip install --upgrade pip
pip install -r requirements.txt

```

#### 2.2. Tracking Large Data Configurations

Because this project interfaces with raw stellar model grids, synthetic lightcurve files, and compiled machine learning weights, **Git Large File Storage** must be initialized. This prevents local repository bloating and guarantees smooth version control histories.

```sh
# 1. initialize the LFS extension on your system
git lfs install

# 2. ensure tracking attributes are correctly locked into version history
git add .gitattributes
git lfs status

```

#### 2.3. Data Provenance & Framework Adaptations

* **Stellar Model Grid Acknowledgments:** The underlying physical stellar model grid utilized by this pipeline contains foundational parameters obtained directly via private communication with **Reza Samadi**.
* **Customizations:** The default repository version of **PSLS** was debugged for modern package versions, exports for addiational metadata from the model grid was implemented, and notably, the function `AddFlare` in `psls.py` was updated to work with other parameters than the defaults.

---

### 3. Technical Pipeline & Code Workflow

The project execution flows sequentially through five distinct Jupyter notebooks. Each module represents a core layer of the pipeline:

```
[ foundation ] ──> [ generate ] ──> [ harvesting ] ──> [ inspect ] ──> [ jobs ]

```

#### 3.1. Explore (`foundation.ipynb`)

The pipeline begins with data analysis of the underlying physical parameter distributions:

* parses the custom stellar model grid files
* maps continuous parameter correlations to build bounds for realistic synthetic systems

#### 3.2. Generate (`generate.ipynb`)

This module functions as our data synthesis engine:

* executes random grid sampling algorithms across verified stellar parameter bounds
* injects analytical planetary transits using limb darkened geometry models
* superimposes multi component stellar noise
* injects random stellar flares drawn from power law distributions, providing a realistic radiation background crucial for habitability calculations

#### 3.3. Sanitize (`harvesting.ipynb`)

Raw stellar time series are dominated by low frequency stellar activity and systematic instrumental drift. This notebook isolates high frequency transit signals using a staged filter:

* **Wotan Detrending:** Employs the `wotan` optimization framework to apply time windowed robust high pass filters using the biweight estimator. This effectively flattens long term stellar variability while preserving the sharp, steep boundaries of individual transit events.
* **Box Least Squares:** Processes the detrended light curves through a periodic grid search using the BLS algorithm. This extracts key observational features like orbital period, transit depth, transit duration, and the transit epoch midpoint.

#### 3.4. Inspect (`inspect.ipynb`)

A structural validation layer designed to spot false positives and catch numerical processing issues:

* Generates multi panel diagnostic plots for selectable planet candidates, displaying the raw light curve, the flattened detrended curve, the diagnostic periodogram power spectrum, and the phase folded transit profile.
* Fits an idealized inverted top hat model back to the phase folded data to evaluate goodness of fit metrics, allowing the filtering out of numerical anomalies or severe stellar activity masquerading as planetary transits in the next step.

#### 3.5. Train (`jobs.ipynb`)

The final phase engineers a feature matrix from the physical properties discovered in the previous steps and feeds it into a machine learning pipeline.

```
                  ┌──────────────────────┐
                  │ Engineered Features: │
                  │  - Planet Size       │
                  │  - Habitable Zone    │
                  │  - Star Radiation    │
                  └──────────┬───────────┘
                             │
            ┌────────────────┴────────────────┐
            ▼                                 ▼
┌───────────────────────┐         ┌───────────────────────┐
│ Random Forest         │         │   Multi Layer Neural  │
│ Classifier/Regressor  │         │   Network Step        │
└───────────┬───────────┘         └───────────┬───────────┘
            │                                 │
            └────────────────┬────────────────┘
                             ▼
              ┌─────────────────────────────┐
              │ Output Habitability Score   │
              │ (Components/Combined)       │
              └─────────────────────────────┘

```

##### Model Architecture Comparisons

This module implements a benchmarking pipeline that evaluates separate machine learning frameworks:

1. **Ensemble Random Forest Pipeline:** Evaluates nonlinear feature splits via a grid search optimized `RandomForestClassifier` and `RandomForestRegressor`. This provides feature importance rankings that help verify which physical traits impact the final classification decision the most.
2. **Multi Layer Perceptron Neural Network:** Implements a dense, fully connected neural network. This maps previous models and underlying feature connections to give a confidence on each result.

##### Multi Tiered Habitability Classification Criteria

The models assess habitability by evaluating three core physical constraints:

* **Planet Size & Composition:** Derived from the transit depth combined with the stellar radius to yield the precise planetary radius. This stage explicitly separates solid, high density rocky terrestrial candidates from volatile rich, low density gas giants.
* **Thermodynamic Environment:** Computes the semimajor axis to establish the incident stellar flux and the equilibrium temperature. This maps the planet position relative to its host star conservative and optimistic habitable zone boundaries.
* **Radiation & Flare Exposure Profiles:** Quantifies the cumulative high energy impact of stellar flares against the background stellar irradiance. The models use this metric to determine if an otherwise thermally stable world is at high risk for atmospheric stripping, ozone depletion, or surface sterilization.

##### Final Results

* table ranking high confidence candidates
* mediocre performance of planetary metrics due to high noise
* extend in future for more advanced methods

---

### 4. References & Bibliography

1. **Heras, A. M. and Rauer, H.** "The PLATO Mission", *COSPAR Scientific Assembly* (2022) [https://ui.adsabs.harvard.edu/abs/2022cosp...44..585H/abstract](https://ui.adsabs.harvard.edu/abs/2022cosp...44..585H/abstract)
2. **Samadi, R. et al.** "The PLATO Solarlike Lightcurve Simulator", *Astronomy & Astrophysics* (2019) [https://ui.adsabs.harvard.edu/abs/2019A%26A...624A.117S%2F/abstract](https://ui.adsabs.harvard.edu/abs/2019A%26A...624A.117S%2F/abstract)
3. **Hippke, M. et al.** "Wōtan: Comprehensive Timeseries Detrending in Python", *The Astronomical Journal* (2019) [https://ui.adsabs.harvard.edu/abs/2019AJ....158..143H/abstract](https://ui.adsabs.harvard.edu/abs/2019AJ....158..143H/abstract)
4. **Kovács, G. et al.** "A Box Fitting Algorithm in the Search for Periodic Transits", *Astronomy & Astrophysics* (2002) [https://ui.adsabs.harvard.edu/abs/2002A%26A...391..369K/abstract](https://ui.adsabs.harvard.edu/abs/2002A%26A...391..369K/abstract)

---

### 5. Testing & Grading

```sh
# 1. verify working Jupyter notebooks
pytest --nbmake *.ipynb

# 2. score against best practice and code style convention
nbqa pylint --msg-template="" --reports=y *.ipynb
```
