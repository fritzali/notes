## PLATO Solarlike Lightcurve Simulator

### Capabilities & Modeling

PSLS simulates solarlike oscillators, defined as stars with outer convective zones that through turbulence excite acoustic oscillations, representative of PLATO mission targets. Stars fulfilling these characteristics make up most of the lower main sequence including dwarfs, as well as subgiants and red giant stars, and also some pre main sequence candidates, while classic high amplitude pulsators such as cepheids, hot massive stars with radiative envelopes, or remnants like white dwarfs are excluded. Due to uncertainty in the pulsation mechanism, the simulator is not suited for the modeling of main sequence dwarf stars.

This tool models such stochastic oscillations, includes planetary transits, spot modulation, flares and granulation, as well as instrumental errors and random noises expected from the detector. These are implemented as follows:
- **Oscillation Spectra:** [Samadi (2019)](https://www.aanda.org/articles/aa/full_html/2019/04/aa34822-18/aa34822-18.html) in `sls.py`
  - **MS & SGB:** [Kjeldsen & Bedding (1995)](https://arxiv.org/abs/astro-ph/9403015) in `sls.py` 
  - **RGB:** [Mosser (2011)](https://www.aanda.org/articles/aa/abs/2011/01/aa15440-10/aa15440-10.html) in `universal_pattern.py`
- **Planetary Transits:** [Mandel & Agol (2002)](https://ui.adsabs.harvard.edu/abs/2002ApJ...580L.171M/abstract) in `transit.py` 
- **Spot Modulation:** [Dorren (1987)](https://ui.adsabs.harvard.edu/abs/1987ApJ...320..756D/abstract) in `spotintime.py`
- **Flares:** [Baudin (2025)](https://github.com/fritzali/notes/blob/main/MASS/Astrobiology/project/PSLS/PKG-INFO) in `flares.py`
- **Granulation:** [Kallinger (2014)](https://www.aanda.org/articles/aa/full_html/2014/10/aa24313-14/aa24313-14.html) in `sls.py`
- **Systematic Influences:** [Marchiori (2019)](https://www.aanda.org/articles/aa/abs/2019/07/aa35269-19/aa35269-19.html) in `systematics/`

### Command Flags

- `-v` prints the program version.
- `-V` makes the output verbose.
- `-P` outputs the power spectral density and lightcurve as plots.
- `-f` saves individual lightcurves for each camera instead of default averaging over all.
- `-m` averages camera groups and then merges interlaced lightcurve while taking into account temporal offset to increase time resolution in exoplanet transits,
as opposed to default averaging of all sensors for noise suppression in astroseismology.
- <code>-M <i>number</i></code> sets amount of performed simulations.
- <code>-o <i>path</i></code> specfies output directory instead of default working directory.
- `--extended-plots` displays an extended set of plots.
- `--psd` saves the power spectral density associated with the lightcurve averaged over all cameras.
- `--pdf` saves plots as `.pdf` instead of `.png` default format.
- `--hdf5` saves averaged lightcurve and simulation components in `.hdf5` file.
- `--proto-sas` formats data saved in `.hdf5` file to be compatible with prototype PLATO SAS pipeline.

### Configuration Variables

The parameters used in the configuration file are explained below:

- **Observation:**
  - **Duration:** simulation duration [days]
  - **MasterSeed:** master seed of the pseudorandom number generator (integer number)
  - **Gaps:**
    - **Enable:** include  [1] or not [0] gaps
    - **Seed:** Seed of the pseudo-random number generator used to generate the gaps ; negative value if controlled by MasterSeed
    - **InterQuarterGapDuration:**  Duration of the inter-quarter interruptions [days] ; inter-quarter gaps are ignored if zero or negative value ; nominal value: 3
    - **RandomGapDuration:** Duration of the random interruptions [minutes] ; random gaps are ignored if zero or negative value ; nominal value: 0
    - **RandomGapTimeFraction:** fraction [in %] of the total time lost by the random gaps ; random gaps are ignored if zero or negative value ; nominal value: 0.5
    - **RandomGapStep:** Random gap step in % ; nominal value: 0
    - **PeriodicGapCadence:** Cadence of the periodic interruptions [days] ; periodic gaps are ignored if zero or negative value ; nominal value: 5.
    - **PeriodicGapDuration:** Duration of each periodic interruption [minutes] ; periodic gaps are ignored if zero or negative value ; nominal value: 20
    - **PeriodicGapJitter:**  Jitter in the time instants the periodic gaps occur [hours] ; normal distribution assumed ; nominal value: 2.
    - **PeriodicGapStep:** Periodic gap step in % ; nominal value: 0
- **Instrument:**
  - **Sampling:** Sampling cadence [s] (nominal value: 25s)
  - **IntegrationTime:** Integration time [s] (nominal value: 21s)
  - **GroupID:** IDs of the camera groups included in the simulation, example [2,3] to simulate group 2 and group 3.  Nominal value: [1,2,3,4].
  - **NCamera:** Number of cameras per group (1 -> 6) (nominal value: 6)
  - **TimeShift:** Time shift between camera groups [s] (nominal value: 6.25s)
  - **RandomNoise:**
    - **Enable:** active [1] or dis-activate [0]  the random noise
    - **Type:** type of random noise, either ‘User’ or ‘PLATO_SCALING’ or ‘PLATO_SIMU’.
      - **‘User’:** the NSR value is specified by the user (see below)
      - **‘PLATO_SCALING’:** the NSR value is obtained by interpolating, at the given magnitude, the NSR scaling relation expected for PLATO
      - **‘PLATO_SIMU’:** the NSR is taken from realistic simulated ligth-curves (stored in the systematics error table, see below) and vary with the mask shapes and thus then the latter are updated
    - **NSR:** User-specified Noise to Signal Ratio [ppm in one hour] for a single camera. This value takes into account all random noises but does not include systematic errors.
    - **Systematics:**
      - **Enable:** active [1] or dis-activate [0]  the  systematic errors
      - **Table:** name of the binary file containing the parameters for the systematic errors
      - **Version:** table version  (the latest version is recommended2)
      - **DriftLevel:** Amplitude of the drift. Can be either  ‘min’, ‘low’, ‘medium’, ‘high’, ‘max’ or ‘any’. Applicable only for Version>0
      - **Seed:** Seed of the pseudo-random number generator used for the systematic errors ; negative value if controlled by MasterSeed
      - **Note:** When systematic errors are enabled, PSLS picks from the systematic error table the target with magnitude close to the magnitude specified by the user (within a +/- 0.25 around the magnitude specified below by the user) and with a drift amplitude in a given range of amplitude (low: 0-0.4 px/90days, medium: 0.4-0.8 px/90days, and high: >0.8 pix/90days). When several targets fulfil the criteria (magnitude and drift level), PSLS randomly selects one of those.
- **Star:**
  - **Mag:** V magnitude (John V magnitude). The V magnitude is converted into the PLATO P magnitude using the star effective temperature following Marchiori et al (2019, A&A)’s Teff-magnitude relation.
  - **ID:** star ID (an arbitrary integer number)
  - **ModelDir:** Directory containing the pulsation models (a single ADIPLS file, a grid of ADIPLS files, or a simple TEXT file, see below)
  - **ModelType:** Type of pulsation model, this can be either ‘UP’, ‘grid’, ‘grid-old’ ‘single’, or ‘text’.  The option ‘grid-old’ shall be used with the old type of grid, the new one being stored in HDF5 files
  - **ModelName:**  Name of the input pulsation model, to be specified when ModelType = ‘single’ or ‘text’.  The input pulsation model can be either a .gsm file (generated by ADIPLS),  a simple TEXT file,   or an HDF5 file storing a grid of CESAM2K models. In the cas of a TEXT file, the latter shall provide the mode properties (in 3 columns: frequency, width and height). Mode frequencies and mode widths are in muHz and mode heights are in ppm^2/muHz. For the mode heights a « single-sided » spectrum is assumed.
  - **ES:** Evolutionary status,  ‘ms’ for the main-sequence phase, ‘sg’ for the sub-giant phase, ‘rg’ for redgiants (Red Giant Branch or clump stars)
  - **Teff:** Star effective temperature [K]
  - **Logg:** Surface gravity, ignored when ModelType = ‘UP’
  - **SurfaceRotationPeriod:** Surface rotation period [days], not used when ModelType = ‘UP’
  - **CoreRotationFreq:** Core rotation frequency [muHz], this is by definition Omega/2pi*1e6 where Omega is the angular rate [rad/s], used only when ModelType = ‘UP’
  - **Inclination:** Inclination angle [deg.]
  - **Seed:** Seed of the pseudo-random number generator used to stellar signal (activity, granulation, and oscillations ; spot excluded) ; negative value if controlled by MasterSeed
- **Oscillations:**
  - **Enable:** include  [1] or not [0]  the solar-like oscillations
  - **numax:** frequency at maximum power [muHZ], used only when ModelType = ‘UP’
  - **delta_nu:** Mean large separation [muHz], used only when ModelType = ‘UP’ , -1 if you want this parameter to be derived from a scaling relation
  - **DPI:** Asymptotic values of the gravity mode period spacing [s], used only when ModelType = ‘UP’, -1 if you don’t want mixed modes to be included
  - **q:** Mixed mode coupling factor, used only when ModelType = ‘UP’
  - **SurfaceEffects:** 1 Include near-surface effects in mode frequencies, not implemented when ModelType = ‘UP’
  - **Seed:** Seed of the pseudo-random number generator used to simulate the oscillations component ; negative value if controlled by MasterSeed
- **Activity:**
  - **Enable:** include  [1] or not [0]  stochastic activity component (Lorentzian component). Should not be used when spot modulation is enabled (see below). Not operating when    - ModelType = « UP ».
  - **Sigma:** Amplitude of the activity component [ppm]
  - **Tau:** Time-scale of the activity component [days]
  - **Seed:** Seed of the pseudo-random number generator used to simulate the (stochastic) activity component ; negative value if controlled by MasterSeed
  - **Spot:**
    - **Enable:** include  [1] or not [0] stellar spots
    - **dOmega:** Differential rotation, dimensionless
    - **MuStar:** Limb darkening coefficient of the star (a linear limb darkening law is assumed)
    - **MuSpot:** Limb darkening coefficient of the spot
    - **Radius:** Spots radii, in degrees, as many values as spots modelled
    - **Latitude:** Spots latitudes, in degrees
    - **Longitude:** Spots latitudes, in degrees
    - **Lifetime:** lifetime of the spot in days
    - **TimeMax:** The time of maximum spot contrast, in days. Negative value if you want to be drawn randomly
    - **Contrast:** Maximum contrast of the spot (flux of the spot in units of unspotted stellar flux), dimensionless
    - **Modulation:** Modulation period of the spot radii [days], ignored if negative or zero value  ; default value: 0
    - **Seed:** Seed of the pseudo-random number generator used to generate TimeMax ; negative value if controlled by MasterSeed
  - **Flare:**
    - **Enable:** 0 # include  [1] or not [0] flares
    - **MeanPeriod:** 2 # mean period btw 2 flares (days)
    - **Amplitude:** 2500. # mean flares amplitudes (ppm) ; the flares amplitudes are drawn from normal distribution centred at FlareAmplitude and with a dispersion of FlareAmplitude / 10
    - **UpDown:** 0.1 # Ratio of the time taken for the flow to rise to the time taken for the flow to fall
    - **MeanDuration:** -1 # mean flare duration (days). If negative, MeanDuration = MeanPeriod/5
    - **DurationDispersion:** -1 # dispersion in the flare duration (days). If negative,  DurationDispersion = MeanPeriod/20
    - **Seed:** -1 # Seed of the pseudo-random number generator used to generate the flares. Negative value if controlled by MasterSeed
- **Granulation:**
  - **Enable:** include  [1] or not [0] stellar granulation
  - **Type:**  Model type. 0-> single Lorentzian component ; 1-> Kallinger et al(2014)’s empirical model
  - **Seed:** Seed of the pseudo-random number generator used to simulate the granulation component ; negative value if controlled by MasterSeed
- **Transit:**
  - **Enable:** include  [1] or not [0] planetary transits
  - **PlanetRadius:** planet raidus [jupiter radii]
  - **OrbitalPeriod:** orbilat period [days]
  - **PlanetSemiMajorAxis:** semi major axis [A.U.]
  - **OrbitalAngle:** orbital angle [deg]
  - **LimbDarkeningCoefficients:** limb darkening coefficients (2 or 4 coefficients). 2 for quadratic law and 4 for a non-linear low
