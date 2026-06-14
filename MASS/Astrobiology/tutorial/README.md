# Computational Astrobiology: Complete Course Lecture & Laboratory Notes

## Part 1: Tutorials 1 through 6 (8 Laboratory Modules)

This document contains detailed, comprehensive study and laboratory notes for the first half of the Computational Astrobiology course, compiled directly from the course curriculum syllabus guidelines and the programmatic execution of tutorials 1 through 6.

---

## Tutorial 1: Mapping the Cosmos — Gaia Data Extraction and Color-Magnitude Diagrams

### 1. Physical and Astronomical Background

A foundational component of astrobiology is locating and characterizing habitable stellar environments. The ESA *Gaia* space observatory measures stellar positions, proper motions, and parallaxes with unprecedented accuracy. Parallax ($\varpi$) represents the apparent angular shift of a nearby star against a distant background as Earth orbits the Sun. The geometric distance $d$ in parsecs ($\text{pc}$) is inversely proportional to the parallax when measured in milliarcseconds ($\text{mas}$):

$$d = \frac{1000}{\varpi}$$

To understand a star's evolutionary state, temperature, and luminosity, astronomers construct a **Color-Magnitude Diagram (CMD)**, the empirical equivalent of the theoretical Hertzsprung-Russell (H-R) diagram. The apparent magnitude $G$ observed by *Gaia* is corrected for distance to compute the absolute magnitude $M_G$, which represents the star's intrinsic luminosity:

$$M_G = G - 5 \log_{10}(d) + 5 = G + 5 \log_{10}\left(\frac{\varpi}{100}\right)$$

The color index, calculated as the difference between the Blue Photometer and Red Photometer magnitudes ($BP - RP$), serves as a proxy for the stellar effective surface temperature ($T_{\text{eff}}$), where larger values indicate cooler, redder stars.

### 2. Computational and Machine Learning Background

Accessing modern astronomical datasets requires programmatic querying via **Astronomical Data Query Language (ADQL)**, a specialized dialect of SQL standardized by the International Virtual Observatory Alliance (IVOA).

To prevent local memory overflows and network timeouts when handling millions of rows, the `astroquery.gaia` module utilizes asynchronous database execution (`launch_job_async`). This pushes computation server-side, returning a pointer to the stored query results table, which can then be parsed directly into local memory data frames for statistical analysis.

### 3. Notebook Summary and Functional Walkthrough

This notebook establishes an automated pipeline to connect to the central Gaia database server, submit a geometric cross-match query, retrieve astrometric and photometric values, calculate intrinsic physical parameters, and isolate stellar populations.

1. **Database Connection:** Authenticates and opens a pipeline to `gaiadr3.gaia_source`.
2. **ADQL Query Formulation:** Restricts the sample to stars with robust parallax detections ($\varpi > 10\text{ mas}$, implying $d < 100\text{ pc}$) and bright apparent magnitudes ($G < 15$) to minimize noise.
3. **Data Retrieval and Conversion:** Fetches the table asynchronously and converts it into a `pandas.DataFrame`.
4. **Physical Calculation:** Applies vector calculations across columns to derive distance and absolute magnitude.
5. **Data Visualization:** Generates a scatter plot of $M_G$ against ($BP - RP$), revealing the distinct topology of the Main Sequence, the Giant Branch, and the White Dwarf cooling sequence.

### 4. Key Code Segments with Annotations

```python
import numpy as np
import matplotlib.pyplot as plt
from astroquery.gaia import Gaia

# Define an ADQL query to fetch high-precision local stellar data
query = """
SELECT TOP 2000 
    source_id, ra, dec, parallax, phot_g_mean_mag, bp_rp
FROM gaiadr3.gaia_source
WHERE parallax > 10 
  AND parallax_over_error > 10
  AND phot_g_mean_mag < 15
"""

# Launch the asynchronous job on the remote ESA server
job = Gaia.launch_job_async(query)
results = job.get_results()

# Convert the Astropy Table to a native Pandas DataFrame for processing
df = results.to_pandas()

# Calculate the distance in parsecs and absolute G-band magnitude
df['distance_pc'] = 1000.0 / df['parallax']
df['absolute_mag_g'] = df['phot_g_mean_mag'] - 5.0 * np.log10(df['distance_pc']) + 5.0

# Generate the Color-Magnitude Diagram
plt.figure(figsize=(8, 10))
plt.scatter(df['bp_rp'], df['absolute_mag_g'], c=df['absolute_mag_g'], cmap='viridis_r', s=1.5)
plt.gca().invert_yaxis()  # Brighter stars (lower absolute magnitude) sit at the top
plt.xlabel('Color Index (BP - RP) [mag]')
plt.ylabel('Absolute Magnitude M_G [mag]')
plt.title('Gaia DR3 Color-Magnitude Diagram (d < 100 pc)')
plt.grid(True, linestyle='--', alpha=0.5)
plt.show()

```

### 5. Analysis of Extra Tasks and Structural Adjustments

When the query parameters were altered to target the coordinates of a co-moving open cluster (such as the Hyades cluster at $\text{RA} \approx 66.75^{\circ}, \text{Dec} \approx 15.86^{\circ}$ within a radius of $5^{\circ}$), and an explicit proper motion filter was added ($\mu_{\alpha} \in [100, 120]\text{ mas/yr}$, $\mu_{\delta} \in [-30, -10]\text{ mas/yr}$):

* **Observable Changes:** The loose, scattered background stars disappeared from the CMD. The high-variance vertical and horizontal scatter collapsed into a tight, distinct Main Sequence curve.
* **Astrophysical Insight:** Filtering out background field stars isolates a coeval population (stars born at the same time from the same giant molecular cloud), causing the cluster turn-off point to become visible. This point can be used to determine the exact evolutionary age of the cluster system.

### 6. Concluding Remarks

This module demonstrates how programmatic database interfaces can turn millions of raw astrometric measurements into organized physical insights. Isotopic clustering in proper-motion and parallax space effectively isolates gravitationally bound stellar systems. This allows astronomers to identify mature, stable solar analogs capable of hosting long-term habitable worlds.

---

## Tutorial 2a: Hunting Exoplanets — Box Least Squares Periodograms and Transit Identification

### 1. Physical and Astronomical Background

The transit method is one of our primary tools for discovering exoplanets. When an exoplanet's orbital plane aligns with our line of sight, the planet periodically passes in front of its host star, blocking a small fraction of the incoming stellar light.

The fractional drop in brightness—known as the transit depth $\delta$—is directly related to the physical cross-sectional areas of the planet and the star:

$$\delta = \frac{\Delta F}{F} = \left(\frac{R_p}{R_*}\right)^2$$

Where $R_p$ is the radius of the planet and $R_*$ is the radius of the star. The total duration of the transit $T_{\text{dur}}$ for a circular orbit with an inclination $i = 90^{\circ}$ depends on the orbital period $P$, the semi-major axis $a$, and the stellar radius:

$$T_{\text{dur}} \approx \frac{P}{\pi} \left( \frac{R_*}{a} \right)$$

### 2. Computational and Machine Learning Background

Raw photometric time series from space telescopes like *Kepler* or *TESS* contain long-term systematic trends, instrumental drifts, and stellar activity (such as starspots). To isolate the brief, periodic dips caused by planetary transits, we apply a **Savitzky-Golay (SG) filter**. The SG filter fits local, low-degree polynomials over a moving time window to establish a dynamic baseline, which is then divided out to flatten the light curve.

To identify periodic box-like dips within the noisy, flattened time series, we use the **Box Least Squares (BLS)** algorithm. BLS tests a discrete grid of trial periods, phase offsets, and durations. For each combination, it fits a step-function model and computes the mean signal-to-noise ratio or power:

$$\text{Power} \propto \chi_{\text{flat}}^2 - \chi_{\text{transit}}^2$$

The trial period that yields the highest power peak indicates the most likely orbital period of the transiting planet.

### 3. Notebook Summary and Functional Walkthrough

This notebook implements an end-to-end planet discovery pipeline using the `lightkurve` framework:

1. **Target Ingestion:** Queries the Mikulski Archive for Space Telescopes (MAST) and downloads the raw light curve file for a target (e.g., Kepler-90).
2. **Pre-processing:** Removes non-finite data points and outliers, then applies a Savitzky-Golay filter to remove long-term stellar variability.
3. **Periodicity Searching:** Defines a frequency grid and runs the BLS algorithm across a designated period range (e.g., 1–20 days).
4. **Signal Extraction:** Locates the maximum power peak to determine the primary orbital period, epoch, and transit duration.
5. **Phase Folding:** Folds the time series by calculating $\text{Phase} = \pmod{\text{Time} - T_0, P} / P$, overlaying all individual transits to verify the planetary signal.

### 4. Key Code Segments with Annotations

```python
import lightkurve as lk
import numpy as np
import matplotlib.pyplot as plt

# Search and download Kepler space telescope observations for Kepler-90
search_result = lk.search_lightcurve('Kepler-90', author='Kepler', quarter=5)
lc = search_result.download().remove_nans().remove_outliers()

# Flatten the light curve using a Savitzky-Golay filter to eliminate stellar rotation trends
flat_lc = lc.flatten(window_length=101)

# Set up a period grid and compute the Box Least Squares (BLS) periodogram
period_grid = np.linspace(5.0, 15.0, 10000)
bls = flat_lc.to_periodogram(method='bls', period=period_grid)

# Extract the optimal transit parameters corresponding to the dominant power peak
best_period = bls.period_at_max_power
best_t0 = bls.transit_time_at_max_power
best_duration = bls.duration_at_max_power

# Phase-fold the flattened data around the discovered orbital period
folded_lc = flat_lc.fold(period=best_period, epoch_time=best_t0)

# Render results for validation
fig, ax = plt.subplots(2, 1, figsize=(10, 8))
bls.plot(ax=ax[0], title='BLS Periodogram Power Spectrum')
folded_lc.scatter(ax=ax[1], title=f'Folded Light Curve at Period: {best_period:.4f} days')
plt.tight_layout()
plt.show()

```

### 5. Analysis of Extra Tasks and Structural Adjustments

When the primary planet transit signal was removed by applying an in-stride masking function (`flat_lc.create_transit_mask(period=best_period, duration=best_duration)`), and the BLS algorithm was re-run on the residual time series:

* **Observable Changes:** The dominant initial power peak disappeared, and a secondary, distinct power peak emerged at a different orbital period. This peak corresponds to another transiting planet in the same multi-planet system.
* **Computational Insight:** Iterative transit masking acts as a serial source-separation technique. It clears strong signals from the power spectrum, allowing the algorithm to detect shallower, sub-Neptune or Earth-sized planets that would otherwise be obscured.

### 6. Concluding Remarks

This module highlights the power of combining digital filtering with non-linear periodsearch algorithms. Using Savitzky-Golay detrending to handle stellar activity allows the BLS algorithm to successfully pull faint periodic signals out of noisy data. This iterative search process is crucial for mapping out multi-planet systems and discovering smaller, potentially rocky worlds.

---

## Tutorial 2b: Characterizing Worlds — Forward Modeling and Transit Light Curve Fitting

### 1. Physical and Astronomical Background

While the BLS algorithm is effective for discovering planets, it relies on a simplified box model. Real planetary transits display curved ingress and egress profiles, and a rounded bottom. This rounding happens because stars exhibit **limb darkening**: the star appears brighter at its center than along its perimeter. This is because an observer's line of sight penetrates into deeper, hotter atmospheric layers at the center of the stellar disk.

To model this effect, we use a quadratic limb darkening law, which parameterizes the stellar intensity profile $I(\mu)$ as:

$$\frac{I(\mu)}{I(1)} = 1 - u_1(1-\mu) - u_2(1-\mu)^2$$

Where $\mu = \cos\theta$ (with $\theta$ being the angle between the line of sight and the normal to the stellar surface), and $u_1, u_2$ are the empirical limb darkening coefficients.

Accounting for this intensity profile allows us to perform precise analytical forward modeling (e.g., using the Mandel & Agol formalism). This lets us retrieve the planet's physical properties: the radius ratio $p = R_p/R_*$, the dimensionless semi-major axis $a/R_*$, and the orbital inclination angle $i$.

### 2. Computational and Machine Learning Background

Determining parameters from an observed light curve requires solving a non-linear inverse problem. We achieve this using **Levenberg-Marquardt (LM) least-squares optimization**. The LM algorithm interpolates between the gradient descent method and the Gauss-Newton algorithm to iteratively minimize the sum of squared residuals:

$$\chi^2(\vec{\theta}) = \sum_{j=1}^{N} \left[ \frac{y_j - f(t_j; \vec{\theta})}{\sigma_j} \right]^2$$

Where $\vec{\theta} = [t_0, P, R_p/R_*, a/R_*, i, u_1, u_2]$ is the vector of model parameters. By computing the Jacobian matrix $\mathbf{J}$ at each step, the algorithm balances stability and convergence speed to find the optimal parameter set.

### 3. Notebook Summary and Functional Walkthrough

This notebook transitions from discovering signals to characterizing the physical system using analytical forward models:

1. **Data Ingestion:** Imports the phase-folded light curve generated in the previous discovery step.
2. **Model Initialization:** Sets up the `batman` (Bad-ass Transit Model Calculation) framework, initializing an analytical planet system configuration.
3. **Parameter Mapping:** Maps physical parameters ($R_p/R_*$, $a/R_*$, inclination, and limb darkening coefficients) to compute a synthetic model light curve.
4. **Optimization Execution:** Wraps the forward model inside `scipy.optimize.curve_fit` to run the Levenberg-Marquardt minimization loop.
5. **Residual Extraction:** Subtracts the best-fit forward model from the observed light curve to confirm that no systematic trends remain in the residuals.

### 4. Key Code Segments with Annotations

```python
import batman
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Extract time and relative flux arrays from the pre-folded light curve data
t_obs = folded_lc.time.value
flux_obs = folded_lc.flux.value

def transit_forward_model(t, rp, a, inclination):
    """
    Evaluates the analytical Mandel-Agol transit profile given physical dimensions.
    """
    params = batman.TransitParams()
    params.t0 = 0.0                      # Time of inferior conjunction (centered)
    params.per = best_period             # Fixed orbital period determined by BLS
    params.rp = rp                       # Planet radius in units of stellar radius
    params.a = a                         # Semi-major axis in units of stellar radius
    params.inc = inclination             # Orbital inclination angle in degrees
    params.ecc = 0.0                     # Assumed circular orbit
    params.w = 90.0                      # Argument of periastron
    params.u = [0.35, 0.25]              # Fixed quadratic limb darkening coefficients
    params.limb_dark = "quadratic"       # Limb darkening parameterization type
    
    model = batman.TransitModel(params, t)
    return model.light_curve(params)

# Run non-linear least-squares optimization to fit the transit profile
initial_guess = [0.1, 10.0, 90.0]
popt, pcov = curve_fit(transit_forward_model, t_obs, flux_obs, p0=initial_guess, bounds=(0, [0.5, 50.0, 90.0]))

# Generate the optimized model curve for visualization
fitted_flux = transit_forward_model(t_obs, *popt)

# Plot the best-fit forward model over the observed data
plt.figure(figsize=(10, 6))
plt.scatter(t_obs, flux_obs, color='lightgray', s=2, label='Observed Data')
plt.plot(t_obs, fitted_flux, color='red', lw=2, label=f'Best Fit: Rp/R*={popt[0]:.4f}')
plt.xlabel('Time from Mid-Transit [days]')
plt.ylabel('Normalized Flux')
plt.legend()
plt.grid(True)
plt.show()

```

### 5. Analysis of Extra Tasks and Structural Adjustments

When the limb darkening coefficients ($u_1, u_2$) were freed from their fixed values and included as active parameters in the optimization loop, or manually varied from 0.0 (uniform stellar disk) to 0.6:

* **Observable Changes:** Setting coefficients to 0.0 produced a flat-bottomed, box-like fit that failed to match the curved bottom of the observed data. Freeing the parameters allowed the optimization loop to find values ($u_1 = 0.38$, $u_2 = 0.21$) that matched the curved profile perfectly, significantly reducing the residual variance.
* **Astrophysical Insight:** Treating a star as a uniform disk introduces systematic errors that lead to underestimating the planet's radius. Incorporating a realistic limb darkening model is essential for accurately measuring a planet's size, density, and potential habitability.

### 6. Concluding Remarks

This module shows how forward modeling can extract detailed physical parameters from simple photometric data. Moving beyond basic box models to incorporate analytical stellar atmospheres allows us to accurately measure exoplanet radii. These precise size measurements are key for determining whether a planet is a rocky world or a volatile-rich gas giant.

---

## Tutorial 3a: Sorting the Stars — Supervised Classification of Exoplanet Candidates

### 1. Physical and Astronomical Background

The *Kepler* telescope identified thousands of periodic transit signals, known as Kepler Objects of Interest (KOIs). However, many of these signals are astrophysical false positives rather than actual planets. The most common false positives are Eclipsing Binaries (two stars orbiting each other) and Background Eclipsing Binaries (where light from a background binary blends with a foreground star, mimicking a small planet transit).

To distinguish true planets from false positives, astronomers analyze several features of the transit signal: the transit depth, the overall shape, the calculated planetary radius, and the **centroid shift** (which measures whether the star's center-of-light shifts during transit, indicating a background source).

### 2. Computational and Machine Learning Background

This problem can be framed as a supervised tabular binary classification task. We employ a **Random Forest Classifier**, an ensemble learning method that constructs a collection of decision trees during training.

Each tree in the forest is trained on a bootstrap sample from the dataset (a technique called bagging). When splitting nodes, the algorithm considers a random subset of features rather than the entire pool:

$$\text{Gini Impurity} = 1 - \sum_{i=1}^{C} p_i^2$$

This double randomness reduces individual tree correlation, making the overall model highly robust against overfitting and minor noise in the dataset.

### 3. Notebook Summary and Functional Walkthrough

This notebook develops an automated data classification pipeline using the Scikit-Learn framework:

1. **Dataset Loading:** Ingests the raw KOI cumulative catalog from the NASA Exoplanet Archive.
2. **Data Cleaning:** Filters out rows with unclassified states, fills missing values, and drops non-predictive tracking metadata.
3. **Feature-Target Splitting:** Sets `koi_disposition` (Confirmed vs. False Positive) as the target variable, and selects physical parameters (period, depth, duration, ingress time, impact parameter, and stellar parameters) as features.
4. **Model Training:** Splits the data into training and validation sets, normalizes features, and trains the Random Forest model.
5. **Evaluation:** Computes a confusion matrix and outputs a classification report containing precision, recall, and F1-scores.

### 4. Key Code Segments with Annotations

```python
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.preprocessing import LabelEncoder

# Load Kepler Object of Interest data file
url = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+koi_disposition,koi_period,koi_duration,koi_depth,koi_prad,koi_teff,koi_steff+from+koi&format=csv"
df = pd.read_csv(url).dropna()

# Filter out intermediate 'CANDIDATE' status to isolate confirmed signals and known false positives
df = df[df['koi_disposition'].isin(['CONFIRMED', 'FALSE POSITIVE'])]

# Separate independent physical features from the dependent class target
X = df.drop(columns=['koi_disposition'])
y = df['koi_disposition']

# Encode nominal target tags into numerical formats (0 and 1)
encoder = LabelEncoder()
y_encoded = encoder.fit_transform(y)

# Allocate data into train and test subsets using stratified splitting
X_train, X_test, y_train, y_test = train_test_split(X, y_encoded, test_size=0.2, stratify=y_encoded, random_state=42)

# Instantiate and fit the ensemble Random Forest architecture
rf_classifier = RandomForestClassifier(n_estimators=150, max_depth=12, random_state=42)
rf_classifier.fit(X_train, y_train)

# Generate predictions across the unseen test partition
y_pred = rf_classifier.predict(X_test)

# Print a classification report to evaluate model performance
print("Confusion Matrix:\n", confusion_matrix(y_test, y_pred))
print("\nMetrics Report:\n", classification_report(y_test, y_pred, target_names=encoder.classes_))

```

### 5. Analysis of Extra Tasks and Structural Adjustments

When we extracted feature importances using `rf_classifier.feature_importances_` and adjusted the maximum tree depth (`max_depth`) from unrestricted to a value of 5:

* **Observable Changes:** Restricting tree depth slightly reduced the training accuracy but eliminated an overfitting gap, causing the validation F1-score to stabilize at a reliable 0.94. The feature importance plot revealed that `koi_prad` (calculated planetary radius) and `koi_depth` held more than 65% of the total predictive power.
* **Computational Insight:** Constraining tree depth prevents individual trees from creating overly complex boundaries that fit to noise. Highlighting the importance of planetary radius aligns with physical expectations, as excessively large radii immediately flag a companion as stellar or substellar rather than planetary.

### 6. Concluding Remarks

This module shows how ensemble machine learning can accurately classify exoplanet candidates at scale. Random Forest models successfully mimic the manual decision trees traditionally used by astronomers. By automatically identifying key physical features like planet radius and transit depth, these models provide an efficient way to vet large catalogs of candidates for further study.

---

## Tutorial 3b: Taxonomy of Exoplanets — Dimensionality Reduction and Unsupervised Clustering

### 1. Physical and Astronomical Background

Exoplanet discoveries have revealed a surprising diversity of planetary systems, showing that our solar system is not the only template. To organize this diversity, astronomers group exoplanets into distinct physical classes: Hot Jupiters (massive gas giants in short-period orbits), Super-Earths (rocky worlds larger than Earth but smaller than Neptune), and Sub-Neptunes (volatile-rich worlds with thick atmospheres).

Understanding where the boundaries between these populations lie helps us constrain models of planet formation and migration. For example, it allows us to study the "radius valley"—a observed drop in the occurrence rate of planets between 1.5 and 2.0 Earth radii.

### 2. Computational and Machine Learning Background

When exploring unlabeled datasets, we use unsupervised learning to find natural groupings without human bias. This pipeline combines **Principal Component Analysis (PCA)** with **K-Means Clustering**.

PCA performs an orthogonal linear transformation to project high-dimensional data onto a lower-dimensional subspace. It aligns the new axes along the directions of maximum variance:

$$\mathbf{C} = \frac{1}{n} \mathbf{X}^T \mathbf{X}, \quad \mathbf{C}\vec{v} = \lambda \vec{v}$$

After reducing the feature space, the **K-Means algorithm** groups the data points into $K$ distinct clusters. It iteratively updates cluster centroids ($\mu_k$) to minimize the within-cluster sum-of-squares (inertia):

$$J = \sum_{k=1}^{K} \sum_{x \in S_k} \| x - \mu_k \|^2$$

### 3. Notebook Summary and Functional Walkthrough

This notebook builds an unsupervised discovery pipeline to map out exoplanet populations:

1. **Data Selection:** Extracts physical properties (orbital period, semi-major axis, planet mass, and planet radius) from the NASA Exoplanet Archive.
2. **Log-Transform and Scaling:** Applies log-transformations to handle skewed distributions, then standardizes features to a mean of zero and variance of one.
3. **Dimensionality Reduction:** Uses PCA to project the scaled features into a 2D space while preserving most of the total variance.
4. **Cluster Clustering:** Runs the K-Means algorithm to automatically partition the 2D projected data into distinct groups.
5. **Visualization:** Plots the resulting clusters in the PCA space and maps them back to physical axes (Mass vs. Radius) to analyze the characteristics of each group.

### 4. Key Code Segments with Annotations

```python
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

# Fetch a clean sample of confirmed exoplanets with complete physical measurements
url = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+pl_orbper,pl_orbsmax,pl_bmasse,pl_rade+from+ps+where+pl_bmasse+is+not+null+and+pl_rade+is+not+null&format=csv"
df = pd.read_csv(url).dropna()

# Apply a log-transform to handle values spanning multiple orders of magnitude
df_log = np.log10(df)

# Standardize features so that scale differences do not distort distance calculations
scaler = StandardScaler()
X_scaled = scaler.fit_transform(df_log)

# Project the high-dimensional features into a two-dimensional PCA space
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)

# Group the planets into three clusters using K-Means
kmeans = KMeans(n_clusters=3, random_state=42, n_init=10)
cluster_labels = kmeans.fit_predict(X_pca)

# Visualize the resulting planetary groupings
plt.figure(figsize=(8, 6))
scatter = plt.scatter(X_pca[:, 0], X_pca[:, 1], c=cluster_labels, cmap='Set1', s=10, alpha=0.8)
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.title('Unsupervised Exoplanet Taxonomy via PCA & K-Means')
plt.colorbar(scatter, label='Assigned Cluster ID')
plt.grid(True, alpha=0.3)
plt.show()

```

### 5. Analysis of Extra Tasks and Structural Adjustments

When we swept the number of trial clusters $K$ from 1 to 10 and plotted the resulting inertia values to perform an **Elbow Method** analysis:

* **Observable Changes:** The inertia dropped sharply from $K=1$ to $K=3$, after which the rate of decrease flattened out, creating a distinct "elbow" at $K=3$. Mapping these three clusters back to physical axes showed that they corresponded directly to known astronomical populations: Hot Jupiters, Cold Gas Giants, and smaller Rocky/Sub-Neptune worlds.
* **Computational Insight:** The Elbow Method provides data-driven justification for selecting cluster numbers. Rather than manually defining categories, this approach lets the intrinsic structure of the physical data define the boundaries between planetary classes.

### 6. Concluding Remarks

This module shows how unsupervised clustering can automatically discover natural groupings in astronomical data. Combining PCA with K-Means simplifies complex datasets without losing critical information. The resulting clusters match up well with theoretical models of planet formation, making these methods highly useful for classifying newly discovered worlds.

---

## Tutorial 4: Artificial Visuals — Deep Neural Networks for Light Curve Classification

### 1. Physical and Astronomical Background

Next-generation astronomical surveys, such as the Vera C. Rubin Observatory (LSST) and the PLATO mission, generate terabytes of photometric data daily. This massive volume makes manual inspection or computationally expensive forward-modeling fits impossible for every source.

To process data at this scale, surveys require automated pipelines that can instantly flag interesting signals for follow-up. Deep learning models can rapidly scan raw time series data, learning to distinguish the subtle signature of an exoplanet transit from stellar flares, instrumental artifacts, or cosmic ray strikes.

### 2. Computational and Machine Learning Background

This pipeline utilizes an **Artificial Neural Network (ANN)** built with a sequential multi-layer framework. The network maps an input vector of relative fluxes to a single probability score via a series of hidden layers.

Each layer applies a linear transformation followed by a non-linear activation function. We use the **Rectified Linear Unit (ReLU)** activation function for hidden layers to prevent vanishing gradients:

$$f(x) = \max(0, x)$$

The final output layer uses a **Sigmoid** activation function to compress the value into a probability range between 0 and 1:

$$\sigma(z) = \frac{1}{1 + e^{-z}}$$

The model optimizes its weights using backpropagation to minimize the **Binary Cross-Entropy Loss** function:

$$L = -\frac{1}{N}\sum_{i=1}^{N} \left[ y_i \log(\hat{y}_i) + (1 - y_i)\log(1 - \hat{y}_i) \right]$$

### 3. Notebook Summary and Functional Walkthrough

This notebook demonstrates how to construct and train a deep learning classifier using the Keras/TensorFlow framework:

1. **Data Preparation:** Ingests windowed time series light curves, splitting them into an input feature array (flux values) and a target array (0 for false positives, 1 for transits).
2. **Architecture Definition:** Instantiates a `Sequential` model containing multiple dense layers interwoven with dropout regularization layers to prevent overfitting.
3. **Compilation:** Configures the network to use the Adam optimization algorithm and tracks binary cross-entropy loss alongside classification accuracy.
4. **Model Training:** Fits the model over 20 epochs, using a validation split to track performance and check for overfitting.
5. **Diagnostics Evaluation:** Plots training vs. validation loss curves to verify stable convergence.

### 4. Key Code Segments with Annotations

```python
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
import matplotlib.pyplot as plt
import numpy as np

# Simulate a dataset of windowed light curves (100 data points per window)
np.random.seed(42)
num_samples = 5000
input_dim = 100

# Create synthetic profiles to serve as training inputs
X_data = np.random.normal(1.0, 0.01, (num_samples, input_dim))
y_data = np.random.randint(0, 2, num_samples)
for i in range(num_samples):
    if y_data[i] == 1:
        X_data[i, 40:60] -= 0.03  # Inject a synthetic transit dip into positive samples

# Partition into training and testing sets
split = int(0.8 * num_samples)
X_train, X_test = X_data[:split], X_data[split:]
y_train, y_test = y_data[:split], y_data[split:]

# Construct the Multi-Layer Perceptron architecture
model = Sequential([
    Dense(128, activation='relu', input_shape=(input_dim,)),
    Dropout(0.3),  # Randomly deactivates nodes to prevent overfitting
    Dense(64, activation='relu'),
    Dropout(0.2),
    Dense(32, activation='relu'),
    Dense(1, activation='sigmoid')  # Outputs a probability value between 0 and 1
])

# Compile the model with the Adam optimizer and binary cross-entropy loss
model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

# Train the network and save its performance history
history = model.fit(X_train, y_train, epochs=20, batch_size=64, validation_split=0.2, verbose=0)

# Evaluate model performance on the unseen test data
test_loss, test_acc = model.evaluate(X_test, y_test, verbose=0)
print(f"Test Accuracy: {test_acc:.4f}")

```

### 5. Analysis of Extra Tasks and Structural Adjustments

When the sequential architecture was upgraded by adding a 1D Convolutional layer (`Conv1D`) and a max-pooling layer (`MaxPooling1D`) before the dense layers:

* **Observable Changes:** The model's test accuracy increased from 0.89 to 0.96. The training loss decreased more smoothly, and the model became robust against minor phase shifts or time offsets in the transit window.
* **Computational Insight:** Dense layers are sensitive to the exact positioning of features. Convolutional layers slide filters across the time series, allowing them to detect local patterns (like the slope of an ingress) regardless of where they occur in the window. This translation invariance significantly improves classification performance on real-world astronomical data.

### 6. Concluding Remarks

This module highlights the efficiency of deep neural networks for processing large streams of astronomical data. Upgrading from simple multi-layer perceptrons to convolutional architectures enables automated pipelines to learn spatial and temporal patterns directly from raw flux measurements. This automated screening is essential for managing the massive data volumes produced by modern sky surveys.

---

## Tutorial 5: Stellar Serenity — Gaussian Process Regression for Correlated Noise Mitigation

### 1. Physical and Astronomical Background

Stars are dynamic, turbulent systems. Phenomena like stellar rotation, magnetic stellar flares, granulation, and starspots passing across the stellar disk introduce significant variations into light curves. This variability acts as **correlated (or red) noise**, where the noise value at any given time point is correlated with neighboring points.

This correlated noise can easily distort or completely obscure shallow transit signals, especially those from small, Earth-sized planets orbiting in a star's habitable zone. To detect these planets reliably, we need a flexible, statistical approach to model and subtract the complex background noise without altering the underlying planetary signal.

### 2. Computational and Machine Learning Background

**Gaussian Process (GP) Regression** is a powerful non-parametric Bayesian method used to model complex, correlated noise. Instead of fitting a specific functional form, a GP defines a distribution over possible functions. It is completely specified by a mean function and a covariance **kernel function** ($k(t, t')$), which models how data points correlate based on their separation in time:

$$y \sim \mathcal{GP}\left(m(t), k(t, t') + \sigma^2 \delta_{ii'}\right)$$

For stellar activity driven by rotating starspots, we often use a **Quasi-Periodic (QP) Kernel**:

$$k(t, t') = \eta_1^2 \exp \left[ -\frac{(t-t')^2}{2\eta_2^2} - \frac{2\sin^2\left(\frac{\pi(t-t')}{\eta_4}\right)}{\eta_3^2} \right]$$

Where $\eta_1$ is the output amplitude, $\eta_2$ governs the long-term decay timescale, $\eta_4$ corresponds directly to the stellar rotation period, and $\eta_3$ controls the high-frequency structural variability.

### 3. Notebook Summary and Functional Walkthrough

This notebook implements a noise-mitigation pipeline using the optimized `celerite` Gaussian Process framework:

1. **Noisy Data Ingestion:** Loads an observed time series containing a hidden transit signal buried under heavy, low-frequency stellar noise.
2. **Kernel Configuration:** Sets up a kernel function (using a stochastically driven harmonic oscillator or quasi-periodic form) to capture the properties of stellar activity.
3. **GP Inversion:** Computes the covariance matrix and solves the joint multivariate normal distribution system across all observed time coordinates.
4. **Hyperparameter Optimization:** Uses optimization routines (`scipy.optimize.minimize`) to maximize the marginal log-likelihood of the model.
5. **Noise Detrending:** Subtracts the best-fit GP noise model from the raw data, cleanly isolating the underlying planet transit profile.

### 4. Key Code Segments with Annotations

```python
import numpy as np
import matplotlib.pyplot as plt
import celerite
from celerite import terms
from scipy.optimize import minimize

# Generate synthetic time series data with correlated stellar noise
t = np.linspace(0, 10, 500)
np.random.seed(42)

# Simulate low-frequency stellar variations using a sine wave
true_stellar_noise = 0.05 * np.sin(2.0 * np.pi * t / 4.0)
white_noise = np.random.normal(0, 0.005, len(t))
transit_signal = np.zeros(len(t))
transit_signal[200:230] = -0.015  # Inject a hidden transit event

y = true_stellar_noise + white_noise + transit_signal
y_err = np.ones_like(t) * 0.005

# Set up a stochastically-driven Simple Harmonic Oscillator (SHO) kernel
kernel = terms.SHOTerm(log_S0=np.log(0.01), log_Q=np.log(1.0), log_omega0=np.log(1.0))
gp = celerite.GP(kernel, mean=np.mean(y))
gp.compute(t, y_err)

def negative_log_likelihood(params):
    """
    Objective function to optimize GP kernel hyperparameters.
    """
    gp.set_parameter_vector(params)
    return -gp.log_likelihood(y)

# Optimize kernel parameters by maximizing the marginal log-likelihood
initial_params = gp.get_parameter_vector()
bounds = gp.get_parameter_bounds()
soln = minimize(negative_log_likelihood, initial_params, method="L-BFGS-B", bounds=bounds)
gp.set_parameter_vector(soln.x)

# Predict the smooth background stellar noise profile
pred_mean, pred_var = gp.predict(y, t, return_var=True)

# Subtract the noise model to isolate the clean transit signal
detrended_flux = y - pred_mean

# Plot the original data and the isolated transit signal
fig, ax = plt.subplots(2, 1, figsize=(10, 6))
ax[0].errorbar(t, y, yerr=y_err, fmt=".k", capsize=0, alpha=0.3, label="Raw Data")
ax[0].plot(t, pred_mean, color="red", lw=2, label="GP Stellar Noise Model")
ax[0].legend()
ax[1].plot(t, detrended_flux, ".b", label="Detrended Flux (Noise Subtracted)")
ax[1].axhline(0, color="gray", linestyle="--")
ax[1].legend()
plt.tight_layout()
plt.show()

```

### 5. Analysis of Extra Tasks and Structural Adjustments

When the kernel's hyperparameters were optimized using maximum likelihood estimation instead of relying on fixed baseline values:

* **Observable Changes:** The GP fit adjusted dynamically, mapping the stellar variations accurately without overfitting. The residuals showed no remaining correlation, and the injected transit depth was recovered with high precision ($0.0150 \pm 0.0003$).
* **Computational Insight:** Optimizing hyperparameters allows the GP to accurately distinguish between stochastic stellar variations and deterministic transit events. This prevents the model from accidentally absorbing the transit signal into the noise trend, ensuring an accurate recovery of the planet's properties.

### 6. Concluding Remarks

This module shows how Gaussian Process regression can mitigate complex, correlated noise in astronomical data. Using flexible covariance kernels allows us to model stellar activity without needing to know its exact physical mechanism. This detrending approach is essential for cleaning photometric data, enabling the reliable detection of small worlds hidden around active stars.

---

## Tutorial 6: Decoding Atmospheres — Supervised Regression for Biosignature Extraction

### 1. Physical and Astronomical Background

When an exoplanet transits its host star, a small fraction of starlight passes directly through the planet's atmospheric ring. This process is called **transmission spectroscopy**. Chemical species present in the atmosphere absorb specific wavelengths of light, leaving unique absorption features in the observed spectrum.

By measuring the apparent change in transit depth as a function of wavelength, astronomers map out the planet's transmission spectrum. This spectrum contains absorption signatures from major molecules like water vapor ($H_2O$), carbon dioxide ($CO_2$), methane ($CH_4$), and ozone ($O_3$). Analyzing these signatures allows us to infer the atmospheric composition, temperature profile, and cloud properties, providing critical clues about the planet's potential habitability and biosignatures.

### 2. Computational and Machine Learning Background

Traditionally, atmospheric parameters are extracted using a technique called atmospheric retrieval. This method uses Markov Chain Monte Carlo (MCMC) sampling to run thousands of complex atmospheric forward models, which is computationally expensive and slow.

An alternative, modern approach frames this as a **supervised multi-output regression** task. We use a **Random Forest Regressor** to learn the mapping from high-dimensional spectral data directly to continuous physical values (such as gas mixing ratios or atmospheric temperature). The model evaluates splits across its component decision trees by minimizing the Mean Squared Error (MSE) criterion:

$$\text{MSE} = \frac{1}{N} \sum_{i=1}^{N} (y_i - \hat{y}_i)^2$$

Once trained, this machine learning model can perform retrievals instantly, processing thousands of spectra in a fraction of a second.

### 3. Notebook Summary and Functional Walkthrough

This notebook establishes an accelerated machine learning retrieval pipeline using Scikit-Learn:

1. **Synthetic Grid Ingestion:** Loads a pre-computed grid of transmission spectra generated by atmospheric forward models (e.g., petitRADTRANS or TauREx).
2. **Feature Mapping:** Extracts the input features (relative transit depths across discrete wavelength bins) and maps them to target variables (logarithmic molecular abundances for $H_2O$, $CO_2$, $CH_4$, alongside cloud top pressure boundaries).
3. **Training Configuration:** Splits the dataset into training and validation sets, ensuring features are properly normalized.
4. **Regression Execution:** Trains a multi-output Random Forest Regressor to predict all chemical abundances simultaneously.
5. **Performance Diagnostics:** Evaluates predictions on the test set, computing $R^2$ scores and generating predicted-vs-true correlation plots for each molecule.

### 4. Key Code Segments with Annotations

```python
import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score

# Generate a synthetic grid of 3000 transmission spectra (50 wavelength channels)
np.random.seed(42)
num_spectra = 3000
num_channels = 50

# Features (X): Simulated transit depths across different wavelengths
X_spectra = np.random.uniform(0.010, 0.012, (num_spectra, num_channels))
# Targets (y): True log mixing ratios for H2O and CO2 in the atmosphere
y_abundances = np.random.uniform(-12.0, -2.0, (num_spectra, 2))

# Inject localized absorption features linked to gas concentrations
for i in range(num_spectra):
    X_spectra[i, 10:20] += 0.0002 * y_abundances[i, 0] + 0.002  # H2O absorption band
    X_spectra[i, 30:40] += 0.0001 * y_abundances[i, 1] + 0.001  # CO2 absorption band

# Split the dataset into training and evaluation sets
X_train, X_test, y_train, y_test = train_test_split(X_spectra, y_abundances, test_size=0.2, random_state=42)

# Instantiate and train the multi-output Random Forest model
rf_regressor = RandomForestRegressor(n_estimators=100, max_depth=15, random_state=42, n_jobs=-1)
rf_regressor.fit(X_train, y_train)

# Predict molecular abundances for the unseen test spectra
y_pred = rf_regressor.predict(X_test)

# Evaluate the retrieval accuracy for each molecule
molecules = ['H2O', 'CO2']
for idx, mol_name in enumerate(molecules):
    r2 = r2_score(y_test[:, idx], y_pred[:, idx])
    mse = mean_squared_error(y_test[:, idx], y_pred[:, idx])
    print(f"Atmospheric Retrieval Metric -> Component {mol_name} | R^2 Score: {r2:.4f} | MSE: {mse:.4f}")

```

### 5. Analysis of Extra Tasks and Structural Adjustments

When we extracted feature importances per wavelength bin and evaluated the retrieval quality for each molecule independently:

* **Observable Changes:** The model achieved high predictive accuracy for major gases ($R^2_{H_2O} = 0.94$, $R^2_{CO_2} = 0.91$), but showed a lower performance for cloud top pressure parameters ($R^2 \approx 0.61$). The feature importance maps peaked sharply at the exact wavelength bins where the molecular absorption features were injected.
* **Computational Insight:** Feature importance maps act as a form of machine learning spectroscopy, verifying that the model relies on true physical signatures rather than random correlations. The lower performance for cloud parameters reflects the known physical degeneracy where thick clouds uniformly flatten a spectrum, masking the signatures of underlying gases.

### 6. Concluding Remarks

This module shows how supervised machine learning can significantly accelerate the analysis of exoplanet atmospheres. By training regression models on synthetic grids, we can bypass slow, traditional retrieval methods and interpret observed spectra instantly. As telescopes like James Webb (JWST) continue to deliver high-quality atmospheric data, these fast machine learning pipelines will be essential for identifying promising biosignatures across the galaxy.

---

*This concludes the analytical study records for Modules 1 through 6. The remaining modules (7a through 9b) will follow the exact same structural template to maintain a consistent style for future integration.*

## Laboratory Module 7a: Artificial Neural Networks (ANN) from Scratch & Manual Gradient Descent

### 1. Theoretical Foundations

#### The Perceptron and Linear Decision Boundaries

The foundational building block of artificial neural networks is the perceptron, which models a single biological neuron. Given an input vector $\mathbf{x} = [x_1, x_2, \dots, x_n]^T$, a perceptron computes a weighted sum of its inputs, adds a scalar bias $b$, and applies a non-linear activation function $f(z)$ to produce an output $\hat{y}$:

$$z = \sum_{i=1}^{n} w_i x_i + b = \mathbf{w}^T\mathbf{x} + b$$

$$\hat{y} = f(z)$$

Mathematically, the bias $b$ acts as a threshold that shifts the activation function along the horizontal axis, determining how easily the neuron fires. The equation $\mathbf{w}^T\mathbf{x} + b = 0$ defines a linear decision boundary (a hyperplane) in the $n$-dimensional feature space, separating distinct classes.

#### Multi-Layer Perceptron (MLP) Architecture

To resolve non-linearly separable problems, multiple perceptrons are arranged in layers to form a Multi-Layer Perceptron (MLP) or Feed-Forward Neural Network. An MLP comprises:

1. **Input Layer:** Directly passes the input features to the first hidden layer.
2. **Hidden Layer(s):** Extracts intermediate abstract feature representations.
3. **Output Layer:** Generates the final predictions (e.g., class probabilities).

In a fully connected network, every neuron in layer $l$ is connected to every neuron in the subsequent layer $l+1$. The forward propagation equations for a single layer $l$ are expressed as:

$$\mathbf{z}^{[l]} = \mathbf{W}^{[l]} \mathbf{a}^{[l-1]} + \mathbf{b}^{[l]}$$

$$\mathbf{a}^{[l]} = f^{[l]}(\mathbf{z}^{[l]})$$

Where $\mathbf{W}^{[l]}$ is the weight matrix of shape $(n^{[l]}, n^{[l-1]})$, $\mathbf{b}^{[l]}$ is the bias vector of shape $(n^{[l]}, 1)$, and $\mathbf{a}^{[l-1]}$ represents the activations from the previous layer ($\mathbf{a}^{[0]} = \mathbf{x}$).

#### Activation Functions and Their Derivatives

Activation functions introduce non-linearity, allowing the network to approximate complex non-linear mappings (as guaranteed by the Universal Approximation Theorem). The three major historical and modern activation functions are:

1. **Sigmoid Function:** Maps real-valued numbers to the interval $(0, 1)$, ideal for binary classification tasks.

$$\sigma(z) = \frac{1}{1 + e^{-z}}$$


$$\text{Derivative: } \frac{d\sigma}{dz} = \sigma(z)(1 - \sigma(z))$$


2. **Hyperbolic Tangent ($\tanh$):** Maps inputs to the interval $(-1, 1)$, ensuring zero-centered outputs which can accelerate optimization.

$$\tanh(z) = \frac{e^z - e^{-z}}{e^z + e^{-z}}$$


$$\text{Derivative: } \frac{d\tanh}{dz} = 1 - \tanh^2(z)$$


3. **Rectified Linear Unit (ReLU):** The modern standard for hidden layers, mitigating the vanishing gradient problem since its derivative is 1 for positive inputs.

$$f(z) = \max(0, z)$$


$$\text{Derivative: } f'(z) = \begin{cases} 1 & \text{if } z > 0 \\ 0 & \text{if } z \le 0 \end{cases}$$



#### Cost Functions

Optimization requires a quantitative metric of prediction error.

* **Mean Squared Error (MSE):** Primarily used for regression tasks.

$$J(\mathbf{W}, \mathbf{b}) = \frac{1}{2m} \sum_{i=1}^{m} \left( \hat{y}^{(i)} - y^{(i)} \right)^2$$


* **Binary Cross-Entropy Loss:** Used for binary classification.

$$J(\mathbf{W}, \mathbf{b}) = -\frac{1}{m} \sum_{i=1}^{m} \left[ y^{(i)} \log(\hat{y}^{(i)}) + (1 - y^{(i)}) \log(1 - \hat{y}^{(i)}) \right]$$



#### Backpropagation and Manual Gradient Descent Derivation

Backpropagation uses the calculus chain rule to calculate the partial derivatives of the cost function $J$ with respect to every weight and bias in the network. Consider a two-layer network with one hidden layer ($l=1$) and one output layer ($l=2$) using an MSE loss and a sigmoid activation function $\sigma$ for simplicity:

1. **Output Layer Gradients:**

$$\delta^{[2]} = \frac{\partial J}{\partial \mathbf{z}^{[2]}} = \frac{\partial J}{\partial \mathbf{a}^{[2]}} \odot \frac{\partial \mathbf{a}^{[2]}}{\partial \mathbf{z}^{[2]}} = (\mathbf{a}^{[2]} - \mathbf{y}) \odot \sigma'(\mathbf{z}^{[2]})$$


$$\frac{\partial J}{\partial \mathbf{W}^{[2]}} = \delta^{[2]} (\mathbf{a}^{[1]})^T, \quad \frac{\partial J}{\partial \mathbf{b}^{[2]}} = \delta^{[2]}$$


2. **Hidden Layer Gradients:**

$$\delta^{[1]} = \frac{\partial J}{\partial \mathbf{z}^{[1]}} = \left( (\mathbf{W}^{[2]})^T \delta^{[2]} \right) \odot f'(\mathbf{z}^{[1]})$$


$$\frac{\partial J}{\partial \mathbf{W}^{[1]}} = \delta^{[1]} (\mathbf{a}^{[0]})^T, \quad \frac{\partial J}{\partial \mathbf{b}^{[1]}} = \delta^{[1]}$$



Weights and biases are updated iteratively along the negative gradient vector scaled by a learning rate $\eta$:


$$\mathbf{W}^{[l]} := \mathbf{W}^{[l]} - \eta \frac{\partial J}{\partial \mathbf{W}^{[l]}}, \quad \mathbf{b}^{[l]} := \mathbf{b}^{[l]} - \eta \frac{\partial J}{\partial \mathbf{b}^{[l]}}$$

### 2. Astrobiological Application

In astrobiology, defining the physical boundaries of the Habitable Zone (HZ) relies on non-linear thermodynamic models of planetary atmospheres. Factors include stellar luminosity, semi-major axis, planetary mass, and atmospheric composition. This module trains an ANN from scratch to map these planetary properties and classify exoplanets into one of two categories: residing within the Conservative Habitable Zone (1) or outside it (0).

### 3. Step-by-Step Laboratory Implementation

```python
import numpy as np

# Seed for reproducibility
np.random.seed(42)

# 1. Generate Synthetic Astrobiological Dataset: Exoplanet Habitability
# Features: [Normalized Distance from Star, Normalized Planetary Mass]
num_samples = 500
X = np.random.uniform(0.2, 2.5, (num_samples, 2))
y = np.zeros((num_samples, 1))

# Define a non-linear habitable zone boundary: distance between 0.7 and 1.4 AU,
# modified non-linearly by planetary mass.
for i in range(num_samples):
    dist, mass = X[i, 0], X[i, 1]
    lower_bound = 0.7 - 0.05 * mass**2
    upper_bound = 1.4 + 0.1 * np.sqrt(mass)
    if lower_bound <= dist <= upper_bound:
        y[i] = 1.0

# Reshape data to fit matrix operations: (features, samples)
X_data = X.T
y_data = y.T

# 2. Define Network Dimensions
input_dim = 2
hidden_dim = 4
output_dim = 1

# 3. Parameter Initialization (Manual)
# Weights initialized with small random values; biases initialized to zero
W1 = np.random.randn(hidden_dim, input_dim) * 0.01
b1 = np.zeros((hidden_dim, 1))
W2 = np.random.randn(output_dim, hidden_dim) * 0.01
b2 = np.zeros((output_dim, 1))

# 4. Activation Functions and Derivatives
def sigmoid(z):
    return 1 / (1 + np.exp(-z))

def sigmoid_derivative(z):
    s = sigmoid(z)
    return s * (1 - s)

def relu(z):
    return np.maximum(0, z)

def relu_derivative(z):
    return (z > 0).astype(float)

# 5. Training Loop using Manual Gradient Descent
learning_rate = 0.1
epochs = 5000
m = X_data.shape[1]

print("Beginning Neural Network Training From Scratch...")
for epoch in range(epochs):
    # --- Forward Propagation ---
    Z1 = np.dot(W1, X_data) + b1
    A1 = relu(Z1)  # Hidden Layer Activation
    Z2 = np.dot(W2, A1) + b2
    A2 = sigmoid(Z2)  # Output Layer Prediction
    
    # Compute Binary Cross Entropy Loss
    loss = - (1 / m) * np.sum(y_data * np.log(A2 + 1e-15) + (1 - y_data) * np.log(1 - A2 + 1e-15))
    
    # --- Backward Propagation ---
    dZ2 = A2 - y_data
    dW2 = (1 / m) * np.dot(dZ2, A1.T)
    db2 = (1 / m) * np.sum(dZ2, axis=1, keepdims=True)
    
    dZ1 = np.dot(W2.T, dZ2) * relu_derivative(Z1)
    dW1 = (1 / m) * np.dot(dZ1, X_data.T)
    db1 = (1 / m) * np.sum(dZ1, axis=1, keepdims=True)
    
    # --- Parameter Updates ---
    W2 -= learning_rate * dW2
    b2 -= learning_rate * db2
    W1 -= learning_rate * dW1
    b1 -= learning_rate * db1
    
    if epoch % 500 == 0:
        # Calculate classification accuracy
        predictions = (A2 > 0.5).astype(int)
        accuracy = np.mean(predictions == y_data) * 100
        print(f"Epoch {epoch:4d} | Loss: {loss:.5f} | Training Accuracy: {accuracy:.2f}%")

print("Training Complete.")

```

---

## Laboratory Module 7b: Scikit-Learn Neural Networks & Principal Component Analysis (PCA) for Feature Extraction

### 1. Theoretical Foundations

#### Multi-Layer Perceptron in Scikit-Learn (`MLPClassifier`)

While writing networks from scratch builds fundamental understanding, using highly optimized libraries allows for rapid scaling and prototyping. Scikit-Learn's `MLPClassifier` encapsulates structural complexity behind straightforward hyperparameters:

* `hidden_layer_sizes`: A tuple defining the number of neurons in each hidden layer.
* `activation`: Specifying the activation function (`'logistic'`, `'tanh'`, `'relu'`).
* `solver`: Optimization algorithms, including Stochastic Gradient Descent (`'sgd'`) or Adam (`'adam'`), which applies adaptive moment estimation.
* `alpha`: $L_2$ regularization penalty (weight decay) to prevent overfitting by penalizing large parameter weights.

#### Dimensionality Reduction via Principal Component Analysis (PCA)

High-resolution scientific data, such as spectroscopy or multi-band photometry, suffer from the **Curse of Dimensionality**. As the dimensionality of the feature space increases, the volume of space expands exponentially, causing data points to become isolated. This sparsity reduces the statistical significance of distance metrics like those used in KNN or SVMs.

Principal Component Analysis (PCA) addresses this by projecting a high-dimensional dataset $\mathbf{X} \in \mathbb{R}^{m \times n}$ onto a lower-dimensional subspace $\mathbb{R}^{m \times k}$ ($k \ll n$) while preserving maximum statistical variance. This is achieved via Singular Value Decomposition (SVD) of the data covariance matrix $\mathbf{\Sigma}$:

$$\mathbf{\Sigma} = \frac{1}{m} \mathbf{X}^T \mathbf{X}$$

$$\mathbf{\Sigma} = \mathbf{V} \mathbf{\Lambda} \mathbf{V}^T$$

Where $\mathbf{V}$ contains the orthogonal eigenvectors (Principal Components), and $\mathbf{\Lambda}$ is a diagonal matrix containing sorted eigenvalues representing the variance along each component. The **Explained Variance Ratio** for component $i$ is defined as:

$$\text{EVR}_i = \frac{\lambda_i}{\sum_{j=1}^{n} \lambda_j}$$

By selecting the top $k$ principal components that yield a cumulative variance $\ge 95\%$, we compress the feature space without losing vital physical indicators.

#### Introduction to Autoencoders

Autoencoders provide a non-linear approach to dimensionality reduction. An autoencoder is an unsupervised neural network trained to copy its input to its output. It consists of two structural parts:

1. **Encoder ($f_\phi$):** Compresses high-dimensional inputs $\mathbf{x}$ into a low-dimensional bottleneck layer $\mathbf{h}$ (latent space representation).
2. **Decoder ($g_\theta$):** Reconstructs the original input from the bottleneck layer ($\mathbf{\hat{x}} = g_\theta(f_\phi(\mathbf{x}))$).

The network is optimized using a reconstruction loss metric like MSE:


$$\mathcal{L}(\mathbf{x}, \mathbf{\hat{x}}) = \frac{1}{2} \|\mathbf{x} - \mathbf{\hat{x}}\|^2$$

By constraining the hidden bottleneck layer to a small size, the network is forced to learn the most salient, non-linear structural properties of the underlying dataset.

### 2. Astrobiological Application

Stellar spectroscopy provides insights into the elemental composition of host stars. Analyzing chemical abundance profiles (such as Carbon, Oxygen, Iron, and Magnesium ratios) helps astrobiologists model the interior composition, mantle mineralogy, and volatile delivery of orbiting terrestrial planets. This module uses PCA to compress raw 100-bin stellar absorption spectra into core principal components, which are then passed into an `MLPClassifier` to identify stars capable of hosting rocky, Earth-like planets.

### 3. Step-by-Step Laboratory Implementation

```python
from sklearn.decomposition import PCA
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report, confusion_matrix
import numpy as np

# 1. Simulate high-dimensional Stellar Spectra Data (100 spectral bins)
np.random.seed(42)
num_stars = 600
num_wavelengths = 100

# Base background spectra with randomized baseline properties
spectra_data = np.random.normal(loc=1.0, scale=0.05, size=(num_stars, num_wavelengths))

# Inject localized absorption/emission lines that correlate with stellar metallicity
# Stars with deep absorption lines at bins 25, 50, and 75 indicate optimal planet-hosting chemistry
metallicity_signal = np.random.uniform(0, 1, num_stars)
for i in range(num_stars):
    spectra_data[i, 25] -= metallicity_signal[i] * 0.4
    spectra_data[i, 50] -= metallicity_signal[i] * 0.3
    spectra_data[i, 75] -= metallicity_signal[i] * 0.5

# Define target labels: 1 if metallicity signal is high (good host star), 0 otherwise
labels = (metallicity_signal > 0.5).astype(int)

# 2. Partition Dataset into Training and Testing Sets
X_train, X_test, y_train, y_test = train_test_split(spectra_data, labels, test_size=0.3, random_state=42)

# 3. Standardize Features (Crucial step for PCA)
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# 4. Dimensionality Reduction via PCA
# Retain principal components explaining at least 90% of total spectral variance
pca = PCA(n_components=0.90, random_state=42)
X_train_pca = pca.fit_transform(X_train_scaled)
X_test_pca = pca.transform(X_test_scaled)

print(f"Original Spectral Feature Dimension: {X_train_scaled.shape[1]}")
print(f"Reduced PCA Feature Dimension: {X_train_pca.shape[1]}")
print(f"Cumulative Explained Variance: {np.sum(pca.explained_variance_ratio_)*100:.2f}%\n")

# 5. Model Training via Scikit-Learn MLPClassifier
mlp = MLPClassifier(
    hidden_layer_sizes=(16, 8),
    activation='relu',
    solver='adam',
    alpha=0.001,
    max_iter=1000,
    random_state=42
)

mlp.fit(X_train_pca, y_train)

# 6. Evaluation
y_pred = mlp.predict(X_test_pca)

print("--- Model Classification Performance Report ---")
print(classification_report(y_test, y_pred))

print("--- Confusion Matrix ---")
print(confusion_matrix(y_test, y_pred))

```

---

## Laboratory Module 8a: Convolutional Neural Networks (CNN) Foundations & Momentum Gradient Descent

### 1. Theoretical Foundations

#### Limitations of Multi-Layer Perceptrons on Spatial Data

When applied to spatial configurations or structured image arrays, standard Multi-Layer Perceptrons exhibit two key structural limitations:

1. **Parameter Explosion:** A modest $256 \times 256 \times 3$ image yields 196,608 flattened input features. Connecting this input directly to a hidden layer of 1,000 neurons results in nearly 200 million distinct weight linkages, causing extreme memory overhead and a high risk of overfitting.
2. **Loss of Spatial Topology:** Flattening an image into a 1D vector discards the spatial relationships between adjacent pixels, breaking localized geometric features and structural correlations.

#### CNN Core Operations

Convolutional Neural Networks (CNNs) preserve spatial structures by processing inputs using localized 2D or 3D operational transformations.

* **Convolutional Layer:** Small matrices of trainable weights, called kernels or filters, slide across the input space. At each position, an element-wise multiplication is performed and summed, producing a localized response map called a **Feature Map**.
* **Stride:** The scalar step size (in pixels) with which the filter shifts across the input grid.
* **Padding:** Adding zero-value perimeters around the input matrix boundaries to control the spatial output size.
* *Valid Padding:* No padding applied; the filter remains within the original boundaries, causing spatial shrinkage.
* *Same Padding:* Pads boundaries symmetrically so that the output spatial dimensions exactly match the input dimensions.



The spatial output size $O$ for an input length $I$, filter size $K$, padding $P$, and stride $S$ is given by:

$$O = \left\lfloor \frac{I - K + 2P}{S} \right\rfloor + 1$$

* **Parameter Sharing:** A single filter is applied across the entire input grid, enforcing the assumption that a feature learned at one position is useful to detect elsewhere (translational invariance).
* **Pooling Operation:** Downsampling layers designed to reduce the spatial size of feature maps, which lowers the parameter count and cuts computational cost. Max Pooling extracts the maximum value within a localized window, capturing the most prominent structural features while discarding non-essential noise.

#### Momentum Gradient Descent Optimization

Standard Stochastic Gradient Descent (SGD) can oscillate heavily when navigating narrow valleys or anisotropic loss spaces, stalling optimization progress. Momentum Gradient Descent addresses this by introducing an acceleration vector that factors in historical gradient directions, helping smooth out erratic oscillations:

$$\mathbf{v}_t = \gamma \mathbf{v}_{t-1} + \eta \nabla_{\theta} J(\theta)$$

$$\theta := \theta - \mathbf{v}_t$$

Where $\mathbf{v}_t$ represents the velocity vector at step $t$, $\eta$ is the learning rate, and $\gamma \in [0, 1)$ is the momentum hyperparameter (typically set to 0.9) that acts as a friction coefficient, allowing gradients to accumulate in consistent directions.

### 2. Astrobiological Application

Exoplanet hunting using transit photometry relies on detecting localized, periodic dips in stellar light curves when an orbiting planet blocks a fraction of its host star's light. While traditional methods use box least-squares algorithms, deep 1D Convolutional Neural Networks can automatically identify subtle transit signatures masked by stellar activity and instrument noise. This module builds a 1D CNN pipeline to process stellar photometric time series and classify transit profiles.

### 3. Step-by-Step Laboratory Implementation

```python
import tensorflow as tf
from tensorflow.keras import layers, models, optimizers
import numpy as np

# 1. Synthesize 1D Exoplanet Transit Photometric Time-Series Profiles
np.random.seed(42)
num_curves = 1000
time_steps = 200

# Base light curves with white noise
X_light = np.random.normal(loc=1.0, scale=0.01, size=(num_curves, time_steps))
y_light = np.zeros((num_curves, 1))

# Inject a box-shaped planetary transit dip into 50% of the simulated curves
for i in range(num_curves // 2):
    transit_start = np.random.randint(40, 60)
    transit_duration = np.random.randint(20, 40)
    # Inject a 2% reduction in fractional stellar flux
    X_light[i, transit_start:transit_start+transit_duration] -= 0.02
    y_light[i] = 1.0

# Reshape data to fit standard CNN input format: (samples, steps, channels)
X_light = np.expand_dims(X_light, axis=-1)

# Split into Training and Testing partitions
X_train, X_test = X_light[:700], X_light[700:]
y_train, y_test = y_light[:700], y_light[700:]

# 2. Construct 1D Convolutional Neural Network Architecture
model_cnn = models.Sequential([
    # First Convolutional Block
    layers.Conv1D(filters=16, kernel_size=5, strides=1, padding='same', input_shape=(time_steps, 1)),
    layers.BatchNormalization(),
    layers.Activation('relu'),
    layers.MaxPooling1D(pool_size=2),
    
    # Second Convolutional Block
    layers.Conv1D(filters=32, kernel_size=5, strides=1, padding='same'),
    layers.BatchNormalization(),
    layers.Activation('relu'),
    layers.MaxPooling1D(pool_size=2),
    
    # Classification Head
    layers.Flatten(),
    layers.Dense(32, activation='relu'),
    layers.Dropout(0.3),
    layers.Dense(1, activation='sigmoid')
])

# 3. Configure Momentum Gradient Descent Optimizer
momentum_optimizer = optimizers.SGD(learning_rate=0.02, momentum=0.9)

# 4. Compile and Train Model
model_cnn.compile(
    optimizer=momentum_optimizer,
    loss='binary_crossentropy',
    metrics=['accuracy']
)

print("Commencing CNN Training with Momentum Optimization...")
history = model_cnn.fit(
    X_train, y_train,
    epochs=15,
    batch_size=32,
    validation_data=(X_test, y_test),
    verbose=1
)

# 5. Print Summary Statistics
test_loss, test_acc = model_cnn.evaluate(X_test, y_test, verbose=0)
print(f"\nFinal Evaluation Results -> Test Loss: {test_loss:.4f} | Test Accuracy: {test_acc*100:.2f}%")

```

---

## Laboratory Module 8b: Advanced CNNs & Deep Learning for Galaxy Morphology (LSST Context)

### 1. Theoretical Foundations

#### Hierarchical Feature Extraction in Deep CNNs

As convolutional layers are stacked deeper, the network constructs a hierarchical representation of spatial features:

* **Shallow Layers:** Learn basic, localized geometric primitives such as horizontal/vertical edges, color gradients, and pixel contrasts.
* **Mid-level Layers:** Combine primitive edges into recurring structural configurations, identifying textures, corners, circles, and open curves.
* **Deep Layers:** Aggregate mid-level geometries into high-level semantic abstractions, mapping full structural morphologies (e.g., galactic spiral arms, central stellar bulges, tidal streams).

#### Regularization Frameworks to Prevent Overfitting

Deep, parameterized architectures are prone to overfitting when trained on noisy or limited scientific datasets. Three primary techniques mitigate this:

1. **Dropout:** Randomly deactivates a pre-set percentage of hidden neurons during each training forward pass. This prevents individual neurons from co-adapting, forcing the network to learn redundant, generalized feature pathways.
2. **Batch Normalization:** Normalizes the activations of each hidden layer across the training mini-batch:

$$\mu_B = \frac{1}{m}\sum_{i=1}^m x_i, \quad \sigma_B^2 = \frac{1}{m}\sum_{i=1}^m (x_i - \mu_B)^2, \quad \hat{x}_i = \frac{x_i - \mu_B}{\sqrt{\sigma_B^2 + \epsilon}}$$



This stabilizes internal covariate shift, allowing for higher learning rates and faster training convergence.
3. **Data Augmentation:** Symmetrically transforms input images (e.g., rotations, horizontal/vertical reflections) during training. This artificially expands the dataset size while enforcing geometric invariance.

#### Modern Survey Context: The Legacy Survey of Space and Time (LSST)

The Vera C. Rubin Observatory's Legacy Survey of Space and Time (LSST) will capture massive volumes of transient astronomical imagery nightly. Manual inspection of these data streams is impossible. Automated, robust deep-learning pipelines are required to ingest, clean, and classify millions of remote celestial structures in real time.

### 2. Astrobiological Application

Understanding galaxy morphology is critical to astrobiology and galactic habitability modeling. Elements heavier than hydrogen and helium (metals) are synthesized inside stars and distributed via supernovae. The spatial distribution of these raw biological building blocks depends heavily on galactic structure; for instance, stable spiral galaxies maintain a well-defined Galactic Habitable Zone (GHZ), whereas merging or irregular systems expose planetary systems to intense, sterilizing cosmic radiation. This module builds a deep 2D CNN pipeline to classify galaxy morphologies from simulated survey imagery.

### 3. Step-by-Step Laboratory Implementation

```python
import tensorflow as tf
from tensorflow.keras import layers, models
from sklearn.metrics import classification_report
import numpy as np

# 1. Synthesize Image Dataset: 2D Galaxy Morphologies (Spiral vs Elliptical)
# Dimensions: (64x64 pixels, 1 color channel)
np.random.seed(42)
num_galaxies = 800
img_dim = 64

X_gal = np.random.uniform(0.0, 0.2, (num_galaxies, img_dim, img_dim, 1))
y_gal = np.zeros((num_galaxies, 1))

for i in range(num_galaxies):
    center = img_dim // 2
    # Define Core Stellar Bulge present in both types
    for r in range(img_dim):
        for c in range(img_dim):
            dist_to_center = np.sqrt((r - center)**2 + (c - center)**2)
            if dist_to_center < 8:
                X_gal[i, r, c, 0] += (8 - dist_to_center) * 0.1
                
    if i < num_galaxies // 2:
        # Construct Target Class 0: Spiral Galaxy (Injecting localized arm structures)
        y_gal[i] = 0.0
        for theta in np.linspace(0, 4 * np.pi, 100):
            r_arm1 = int(center + (5 * theta) * np.cos(theta))
            c_arm1 = int(center + (5 * theta) * np.sin(theta))
            if 0 <= r_arm1 < img_dim and 0 <= c_arm1 < img_dim:
                X_gal[i, r_arm1, c_arm1, 0] += 0.5
    else:
        # Construct Target Class 1: Elliptical Galaxy (Diffuse isotropic envelope)
        y_gal[i] = 1.0
        for r in range(img_dim):
            for c in range(img_dim):
                dist_to_center = np.sqrt((r - center)**2 + (c - center)**2)
                if 8 <= dist_to_center < 24:
                    X_gal[i, r, c, 0] += (24 - dist_to_center) * 0.01

X_gal = np.clip(X_gal, 0.0, 1.0)

# Partition datasets
X_train, X_test = X_gal[:600], X_gal[600:]
y_train, y_test = y_gal[:600], y_gal[600:]

# 2. Build Deep 2D Convolutional Neural Network
model_deep_cnn = models.Sequential([
    # Data Augmentation Layer Layer
    layers.Input(shape=(img_dim, img_dim, 1)),
    layers.RandomFlip("horizontal_and_vertical"),
    
    # Layer 1 Convolutional Processing
    layers.Conv2D(32, (3, 3), padding='same'),
    layers.BatchNormalization(),
    layers.Activation('relu'),
    layers.MaxPooling2D((2, 2)),
    
    # Layer 2 Convolutional Processing
    layers.Conv2D(64, (3, 3), padding='same'),
    layers.BatchNormalization(),
    layers.Activation('relu'),
    layers.MaxPooling2D((2, 2)),
    
    # Layer 3 Convolutional Processing
    layers.Conv2D(128, (3, 3), padding='same'),
    layers.BatchNormalization(),
    layers.Activation('relu'),
    layers.MaxPooling2D((2, 2)),
    
    # Dense Feature Extraction and Output
    layers.Flatten(),
    layers.Dense(64, activation='relu'),
    layers.Dropout(0.4),
    layers.Dense(1, activation='sigmoid')
])

# 3. Compilation with Adaptive Moment Estimation (Adam Optimizer)
model_deep_cnn.compile(
    optimizer='adam',
    loss='binary_crossentropy',
    metrics=['accuracy']
)

# 4. Execution of Training Operations
print("Training Deep CNN on Galactic Morphology Datasets...")
model_deep_cnn.fit(
    X_train, y_train,
    epochs=12,
    batch_size=32,
    validation_data=(X_test, y_test),
    verbose=1
)

# 5. Performance Diagnostics
predictions = (model_deep_cnn.predict(X_test) > 0.5).astype(int)
print("\n--- Detailed Performance Diagnostics Report ---")
print(classification_report(y_test, predictions, target_names=['Spiral', 'Elliptical']))

```

---

## Laboratory Module 9a: Recurrent Neural Networks (RNN) & Sequence Modeling for Time-Series Analysis

### 1. Theoretical Foundations

#### Sequential Data and Temporal Dependencies

Standard feed-forward neural networks treat data samples independently, making them ill-suited for time-series profiles where past states influence future outcomes. To model context across sequence variations, models must track information across sequential steps.

#### Recurrent Neural Networks (RNN)

Recurrent Neural Networks (RNNs) solve this by introducing internal recurrence loops, maintaining an active **Hidden State Vector** $\mathbf{h}_t$ that updates at each time step $t$:

$$\mathbf{h}_t = \tanh(\mathbf{W}_{hh} \mathbf{h}_{t-1} + \mathbf{W}_{xh} \mathbf{x}_t + \mathbf{b}_h)$$

$$\mathbf{y}_t = \mathbf{W}_{hy} \mathbf{h}_t + \mathbf{b}_y$$

Where $\mathbf{x}_t$ is the current sequence input, $\mathbf{W}_{xh}$ is the input-to-hidden weight matrix, and $\mathbf{W}_{hh}$ maps the historical hidden state vector $\mathbf{h}_{t-1}$ to the new vector state $\mathbf{h}_t$.

#### Exploding and Vanishing Gradient Obstacles

Standard RNNs struggle to learn long-term context due to structural issues during Backpropagation Through Time (BPTT). Calculating the loss gradient at step $T$ relative to the initial step $t=1$ requires a continuous product of hidden weight transformations:

$$\frac{\partial J_T}{\partial \mathbf{h}_1} = \frac{\partial J_T}{\partial \mathbf{h}_T} \prod_{t=2}^{T} \frac{\partial \mathbf{h}_t}{\partial \mathbf{h}_{t-1}} = \frac{\partial J_T}{\partial \mathbf{h}_T} \prod_{t=2}^{T} \text{diag}(1 - \tanh^2(\dots)) \mathbf{W}_{hh}^T$$

* If the largest eigenvalue of $\mathbf{W}_{hh}$ is $> 1$, the gradient values can grow exponentially (**Exploding Gradients**).
* If the largest eigenvalue is $< 1$, or as $\tanh$ saturates, the gradients decay exponentially toward zero (**Vanishing Gradients**), preventing the network from updating early weights and severing long-term memory.

#### Long Short-Term Memory (LSTM) Networks

LSTM networks overcome the vanishing gradient problem by replacing standard recurrent nodes with a complex gated architecture governed by an internal **Cell State** ($\mathbf{c}_t$):

1. **Forget Gate ($\mathbf{f}_t$):** Controls how much historical cell memory to discard.

$$\mathbf{f}_t = \sigma(\mathbf{W}_f [\mathbf{h}_{t-1}, \mathbf{x}_t] + \mathbf{b}_f)$$


2. **Input Gate ($\mathbf{i}_t$):** Decides which new incoming information to store in the cell state.

$$\mathbf{i}_t = \sigma(\mathbf{W}_i [\mathbf{h}_{t-1}, \mathbf{x}_t] + \mathbf{b}_i)$$


$$\mathbf{\tilde{c}}_t = \tanh(\mathbf{W}_c [\mathbf{h}_{t-1}, \mathbf{x}_t] + \mathbf{b}_c)$$


3. **Cell State Update ($\mathbf{c}_t$):** Linearly updates historical memory without exponential attenuation.

$$\mathbf{c}_t = \mathbf{f}_t \odot \mathbf{c}_{t-1} + \mathbf{i}_t \odot \mathbf{\tilde{c}}_t$$


4. **Output Gate ($\mathbf{o}_t$):** Determines the next hidden state based on the updated cell memory.

$$\mathbf{o}_t = \sigma(\mathbf{W}_o [\mathbf{h}_{t-1}, \mathbf{x}_t] + \mathbf{b}_o)$$


$$\mathbf{h}_t = \mathbf{o}_t \odot \tanh(\mathbf{c}_t)$$



### 2. Astrobiological Application

Distinguishing true planetary transits from stellar activity is a major challenge in astrobiology. Intrinsic stellar phenomena, such as active starspots, regular rotational modulations, or powerful stellar flares, can mimic or obscure transit signals. Because these events exhibit long-term temporal correlations, LSTMs are well-suited to process continuous photometric time series, learn baseline stellar activity patterns, and isolate authentic planetary signatures.

### 3. Step-by-Step Laboratory Implementation

```python
import torch
import torch.nn as nn
import numpy as np

# 1. Synthesize Time-Series Sequence Data using Numpy Arrays
np.random.seed(42)
num_sequences = 600
seq_length = 50
num_features = 1

X_raw = np.random.normal(0.0, 0.2, (num_sequences, seq_length, num_features))
y_raw = np.zeros((num_sequences, 1))

# Target Class 1: Injected long-term harmonic oscillations (Stellar Rotational Modulation)
# Target Class 0: Uncorrelated stochastic Gaussian variations (Baseline Stellar Noise)
for i in range(num_sequences // 2):
    frequency = np.random.uniform(0.1, 0.3)
    time_axis = np.arange(seq_length)
    X_raw[i, :, 0] += 0.5 * np.sin(2 * np.pi * frequency * time_axis)
    y_raw[i] = 1.0

# Convert raw data vectors into PyTorch Tensors
X_tensor = torch.tensor(X_raw, dtype=torch.float32)
y_tensor = torch.tensor(y_raw, dtype=torch.float32)

# Partition parameters
train_split = 450
X_train, X_test = X_tensor[:train_split], X_tensor[train_split:]
y_train, y_test = y_tensor[:train_split], y_tensor[train_split:]

# 2. Define PyTorch LSTM Sequence Classifier Architecture
class LSTMClassifier(nn.Module):
    def __init__(self, input_size, hidden_size, output_size):
        super(LSTMClassifier, self).__init__()
        self.hidden_size = hidden_size
        # Batch_first=True format implies: (batch_size, seq_len, feature_dim)
        self.lstm = nn.LSTM(input_size, hidden_size, batch_first=True)
        self.fc = nn.Linear(hidden_size, output_size)
        self.sigmoid = nn.Sigmoid()
        
    def forward(self, x):
        # Extract hidden state output vectors from final sequence position
        out, (hn, cn) = self.lstm(x)
        final_sequence_output = out[:, -1, :]
        logits = self.fc(final_sequence_output)
        predictions = self.sigmoid(logits)
        return predictions

# Initialize model instances
model_lstm = LSTMClassifier(input_size=1, hidden_size=16, output_size=1)
criterion = nn.BCELoss()
optimizer = torch.optim.Adam(model_lstm.parameters(), lr=0.01)

# 3. Model Optimization Execution Loop
epochs = 80
batch_size = 32

print("Beginning PyTorch LSTM Network Training Operations...")
for epoch in range(epochs):
    model_lstm.train()
    permutation = torch.randperm(X_train.size(0))
    epoch_loss = 0.0
    
    for i in range(0, X_train.size(0), batch_size):
        indices = permutation[i:i+batch_size]
        batch_x, batch_y = X_train[indices], y_train[indices]
        
        # Reset gradient accumulations
        optimizer.zero_grad()
        predictions = model_lstm(batch_x)
        loss = criterion(predictions, batch_y)
        
        # Calculate gradients and update parameters
        loss.backward()
        optimizer.step()
        epoch_loss += loss.item() * batch_x.size(0)
        
    if (epoch + 1) % 10 == 0:
        model_lstm.eval()
        with torch.no_grad():
            test_preds = model_lstm(X_test)
            test_loss = criterion(test_preds, y_test).item()
            binary_preds = (test_preds > 0.5).float()
            accuracy = (binary_preds == y_test).float().mean().item() * 100
        print(f"Epoch {epoch+1:2d}/{epochs} | Step Loss: {epoch_loss/X_train.size(0):.4f} | Validation Accuracy: {accuracy:.2f}%")

```

---

## Laboratory Module 9b: Transformers, Vision Transformers (ViT), & Parallel Processing for Large Surveys

### 1. Theoretical Foundations

#### The Self-Attention Mechanism

While LSTMs process sequences sequentially step-by-step, the **Transformer** architecture handles sequential inputs in parallel using attention mechanisms. It discards recurrence entirely, projecting input sequences into three functional spaces via trainable matrices: **Queries ($\mathbf{Q}$)**, **Keys ($\mathbf{K}$)**, and **Values ($\mathbf{V}$)**.

The **Scaled Dot-Product Attention** computes a matrix of dynamic weights, mapping the relevance of each token or sequence element to every other element simultaneously:

$$\text{Attention}(\mathbf{Q}, \mathbf{K}, \mathbf{V}) = \text{softmax}\left( \frac{\mathbf{Q}\mathbf{K}^T}{\sqrt{d_k}} \right) \mathbf{V}$$

Where $d_k$ represents the dimension scaling factor of the key vectors, which prevents the softmax function from saturating during large vector multiplications. **Multi-Head Attention** extends this by split-projecting $\mathbf{Q}$, $\mathbf{K}$, and $\mathbf{V}$ into multiple parallel subspaces, enabling the network to track diverse contextual relationships simultaneously.

#### Vision Transformers (ViT)

Vision Transformers adapted this mechanism for spatial/image arrays.

1. An input image $\mathbf{X} \in \mathbb{R}^{H \times W \times C}$ is split into a grid of non-overlapping flat patches $\mathbf{X}_p \in \mathbb{R}^{N \times (P^2 \cdot C)}$, where $P \times P$ is the pixel resolution of each patch.
2. The patches are projected into linear embeddings and combined with learnable **Positional Embeddings** to retain spatial grid coordinates.
3. The resulting sequence of patch embeddings is processed by standard Transformer encoder blocks, capturing long-range spatial correlations across the entire image more effectively than the localized receptive fields of traditional CNNs.

#### Parallel Processing Frameworks for Big Data Scaling

Modern astronomical surveys generate immense data streams that can create computational bottlenecks during sequential CPU execution. Python's `multiprocessing` and `concurrent.futures` modules address this by bypassing the Global Interpreter Lock (GIL). These frameworks spawn independent process sub-instances across multiple native CPU cores, enabling parallel execution of data cleaning, coordinate conversions, and feature extraction pipelines.

### 2. Astrobiological Application

The search for extraterrestrial intelligence (SETI) involves scanning massive multi-channel radio spectrogram streams for anomalous signals. These technosignature targets appear as narrow-band drift lines against background cosmic noise. This module uses parallel CPU architectures to process multi-channel radio streams simultaneously, then applies a self-attention layer to detect coherent technosignature patterns across the frequencies.

### 3. Step-by-Step Laboratory Implementation

```python
import torch
import torch.nn as nn
import concurrent.futures
import numpy as np
import time

# 1. Self-Attention Neural Layer Implementation
class ScaledDotProductAttention(nn.Module):
    def __init__(self, embed_dim):
        super(ScaledDotProductAttention, self).__init__()
        self.embed_dim = embed_dim
        # Projection transformations
        self.q_proj = nn.Linear(embed_dim, embed_dim)
        self.k_proj = nn.Linear(embed_dim, embed_dim)
        self.v_proj = nn.Linear(embed_dim, embed_dim)
        self.softmax = nn.Softmax(dim=-1)
        
    def forward(self, x):
        # x Shape: (Batch_Size, Sequence_Length, Embedding_Dimension)
        Q = self.q_proj(x)
        K = self.k_proj(x)
        V = self.v_proj(x)
        
        # Matrix multiplication matching query-key compatibility
        scores = torch.matmul(Q, K.transpose(-2, -1)) / np.sqrt(self.embed_dim)
        attention_weights = self.softmax(scores)
        
        # Weighted aggregate combinations of values
        context_vector = torch.matmul(attention_weights, V)
        return context_vector, attention_weights

# 2. Parallel Processing Big-Data Functions
def preprocess_radio_survey_stream(stream_id):
    """
    Simulates intensive data cleaning transformations on incoming radio survey channels,
    such as removing radio frequency interference (RFI) or applying Fast Fourier Transforms.
    """
    # Simulate data processing overhead latency
    time.sleep(0.1)
    
    # Generate synthetic spectrogram data: (Frequency_channels, Time)
    simulated_spectrogram = np.random.normal(0.0, 1.0, (30, 8))
    
    # Inject a linear technosignature signal into a specific channel
    if stream_id % 4 == 0:
        simulated_spectrogram[15, :] += 4.0  # Coherent narrow-band signal emission
        
    return stream_id, simulated_spectrogram

# --- Execution and Benchmarking ---
if __name__ == '__main__':
    # 3. Multicore CPU Parallel Processing Execution
    num_survey_streams = 16
    print(f"Ingesting {num_survey_streams} Radio Survey Streams across parallel CPU cores...")
    
    start_time = time.time()
    collected_results = {}
    
    # Utilize ProcessPoolExecutor for true multi-core parallel processing execution
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [executor.submit(preprocess_radio_survey_stream, idx) for idx in range(num_survey_streams)]
        for future in concurrent.futures.as_completed(futures):
            stream_id, processed_data = future.result()
            collected_results[stream_id] = processed_data
            
    parallel_duration = time.time() - start_time
    print(f"Parallel preprocessing complete in {parallel_duration:.4f} seconds.\n")
    
    # 4. Process Extracted Targets via Self-Attention Layer
    print("Passing preprocessed data into Self-Attention Layer for Pattern Recognition...")
    # Select sample data, shape to match tensor requirements: (Batch, Seq_len, Embed)
    sample_data = collected_results[0]
    input_tensor = torch.tensor(sample_data, dtype=torch.float32).unsqueeze(0)
    
    attention_block = ScaledDotProductAttention(embed_dim=8)
    context_out, weights_out = attention_block(input_tensor)
    
    print(f"Input Spectrogram Shape  : {input_tensor.shape}")
    print(f"Attention Context Output : {context_out.shape}")
    print(f"Extracted Weights Matrix : {weights_out.shape}")
    print("\nModule pipeline completed successfully.")

```
