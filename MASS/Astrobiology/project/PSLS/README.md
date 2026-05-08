## PLATO Solarlike Lightcurve Simulator

### Capabilities & Modeling

PSLS simulates solarlike oscillators, defined as stars with outer convective zones that through turbulence excite acoustic oscillations, representative of PLATO mission targets. Stars fulfilling these characteristics make up most of the lower main sequence including dwarfs, as well as subgiants and red giant stars, and also some pre main sequence candidates, while classic high amplitude pulsators such as cepheids, hot massive stars with radiative envelopes, or remnants like white dwarfs are excluded.

> Note: Due to uncertainty in the pulsation mechanism, the simulator is not suited for the modeling of main sequence dwarf stars which are therefore excluded from any useful analysis.

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
  - **Duration:** simulation duration in days
  - **MasterSeed:** integer number master seed of the pseudorandom number generator
  - **Gaps:**
    - **Enable:** include *1* or exclude *0* gaps
    - **Seed:** seed of the pseudorandom number generator used to generate gaps, use negative value if controlled by **MasterSeed**
    - **InterQuarterGapDuration:** duration of inter quarter interruptions in days, inter quarter gaps are ignored if zero or negative value, nominal value *3*
    - **RandomGapDuration:** duration of random interruptions in minutes, random gaps are ignored if zero or negative value, nominal value *0*
    - **RandomGapTimeFraction:** fraction in percent of total time lost by random gaps, random gaps are ignored if zero or negative value, nominal value *0.5*
    - **RandomGapStep:** random gap step in percent, nominal value *0*
    - **PeriodicGapCadence:** cadence of periodic interruptions in days, periodic gaps are ignored if zero or negative value, nominal value *5*
    - **PeriodicGapDuration:** duration of each periodic interruption in minutes, periodic gaps are ignored if zero or negative value, nominal value *20*
    - **PeriodicGapJitter:**  jitter of time instants with which the periodic gaps occur in hours, normal distribution assumed, nominal value *2*
    - **PeriodicGapStep:** periodic gap step in percent, nominal value *0*
- **Instrument:**
  - **Sampling:** sampling cadence in seconds, nominal value *25*
  - **IntegrationTime:** integration time in seconds, nominal value *21*
  - **GroupID:** identifiers of the camera groups included in simulation, example *[2,3]* to simulate groups two and three, nominal value *[1,2,3,4]*
  - **NCamera:** number of cameras per group, between one and six, nominal value *6*
  - **TimeShift:** time shift between camera groups in seconds, nominal value *6.25*
  - **RandomNoise:**
    - **Enable:** enable *1* or disable *0* the random noise
    - **Type:** type of random noise, either ***User*** or ***PLATO_SCALING*** or ***PLATO_SIMU***
      - **User:** the **NSR** value is specified by the user
      - **PLATO_SCALING:** the **NSR** value is obtained by interpolating at a given magnitude the **NSR** scaling relation expected for PLATO
      - **PLATO_SIMU:** the **NSR** is taken from realistic simulated ligth curves stored in the systematics error table, varied with the mask shapes and then updated
    - **NSR:** user specified noise to signal ratio in parts per million over one hour for a single camera, takes into account all random noises but does not include systematic errors
    - **Systematics:**
      - **Enable:** enable *1* or disable *0* systematic errors
      - **Table:** name of the binary file containing the parameters for the systematic errors
      - **Version:** table version, latest being recommended
      - **DriftLevel:** amplitude of the drift, can be either  *min* or *low* or *medium* or *high* or *max* or *any*
      - **Seed:** seed of the pseudorandom number generator used for the systematic errors, negative value if controlled by **MasterSeed**
      <br>
      
      > when systematic errors are enabled, PSLS picks from the systematic error table a target with magnitude closest to the one specified by the user, within *0.25* around the specification, and with a drift amplitude in a given range of amplitude
      > <br>
      > - low: *0.0* to *0.4* pixels per *90* days
      > - medium: *0.4* to *0.8* pixels per *90* days
      > - high: above *0.8* pixels per *90* days
      > <br>
      > when several targets fulfill the criteria for magnitude and drift level, PSLS randomly selects one

- **Star:**
  - **Mag:** specify John V magnitude, then converted into the PLATO P magnitude using the star effective temperature
  - **ID:** star identifier, arbitrary integer number
  - **ModelDir:** directory containing the pulsation models consisting of a single ADIPLS file or a grid of ADIPLS files or a simple text file
  - **ModelType:** type of pulsation model, this can be either *UP* or *grid* or *grid-old* or *single* or *text*
  - **ModelName:** name of the input pulsation model to be specified when *ModelType* is *single* or *text* loaded from either a *.gsm* file generated by ADIPLS or a simple *.txt* file or an *.hdf5* file storing a grid of CESAM2K models
    <br>

    > in case of a text file, the latter shall provide the mode properties in three columns of frequency and width and height, with mode frequencies and mode widths in micro hertz as well as mode heights are in squared parts per million over micro hertz, where for the mode heights a single sided spectrum is assumed

  - **ES:** evolutionary status, *ms* for the main sequence phase or *sg* for the subgiant phase or *rg* for red giant branch or clump stars
  - **Teff:** star effective temperature in kelvin
  - **Logg:** surface gravity, ignored when **ModelType** is *UP*
  - **SurfaceRotationPeriod:** surface rotation period in days, not used when **ModelType** is *UP*
  - **CoreRotationFreq:** core rotation frequency in micro hertz, by definition connected to the angular rate in radiants per second, used only when **ModelType** is *UP*
  - **Inclination:** inclination angle in degrees
  - **Seed:** seed of the pseudorandom number generator used for stellar signals of activity and granulation with oscillations but excluding spots, negative value if controlled by **MasterSeed**
- **Oscillations:**
  - **Enable:** include *1* or exclude *0* the solarlike oscillations
  - **numax:** frequency at maximum power in micro hertz, used only when **ModelType** is *UP*
  - **delta_nu:** mean large separation in micro hertz, used only when **ModelType** is *UP*, set *-1* if you want this parameter to be derived from a scaling relation
  - **DPI:** asymptotic values of the gravity mode period spacing in seconds, used only when **ModelType** is *UP*, set *-1* if you do not want mixed modes to be included
  - **q:** mixed mode coupling factor, used only when **ModelType** is *UP*
  - **SurfaceEffects:** set *1* to include near surface effects in mode frequencies, not implemented when **ModelType** is *UP*
  - **Seed:** seed of the pseudorandom number generator used to simulate the oscillations component, negative value if controlled by **MasterSeed**
- **Activity:**
  - **Enable:** include *1* or exclude *0* lorentzian stochastic activity component, should not be used when spot modulation is enabled, not operating when **ModelType** is *UP*
  - **Sigma:** amplitude of the activity component in parts per million
  - **Tau:** timescale of the activity component in days
  - **Seed:** seed of the pseudorandom number generator used to simulate the stochastic activity component, negative value if controlled by **MasterSeed**
  - **Spot:**
    - **Enable:** include *1* or exclude *0* stellar spots
    - **dOmega:** differential rotation, dimensionless
    - **MuStar:** limb darkening coefficient of the star, a linear limb darkening law is assumed
    - **MuSpot:** limb darkening coefficient of the spot
    - **Radius:** spots radii in degrees, as many values as spots modelled
    - **Latitude:** spots latitudes in degrees
    - **Longitude:** spots longitudes in degrees
    - **Lifetime:** lifetime of the spot in days
    - **TimeMax:** time of maximum spot contrast in days, negative value if you want to be drawn randomly
    - **Contrast:** maximum contrast of the spot, flux of the spot in units of unspotted stellar flux, dimensionless
    - **Modulation:** modulation period of the spot radii in days, ignored if negative or zero value, default value *0*
    - **Seed:** seed of the pseudorandom number generator used to generate **TimeMax**, negative value if controlled by **MasterSeed**
  - **Flare:**
    - **Enable:** include *1* or exclude *0* flares, default value *1*
    - **MeanPeriod:** mean period between two flares in days, default value *2*
    - **Amplitude:** mean flare amplitude in parts per million, amplitudes are drawn from normal distribution centred at **FlareAmplitude** with dispersion one tenth of **FlareAmplitude**, default value *2500*
    - **UpDown:** ratio of the time it takes for the flow to rise versus the time for it to fall, default value *0.1*
    - **MeanDuration:** mean flare duration in days, if negative **MeanDuration** is one fifth of **MeanPeriod**, default value *-1*
    - **DurationDispersion:** dispersion in the flare duration in days, if negative **DurationDispersion** is one twentieth of **MeanPeriod**, default value *-1*
    - **Seed:** seed of the pseudorandom number generator used to generate the flares, negative value if controlled by **MasterSeed**. default value *-1*
- **Granulation:**
  - **Enable:** include *1* or exclude *0* stellar granulation
  - **Type:**  model type of *0* for single lorentzian component or *1* for more complex empirical model
  - **Seed:** seed of the pseudorandom number generator used to simulate the granulation component, negative value if controlled by **MasterSeed**
- **Transit:**
  - **Enable:** include *1* or exclude *0* planetary transits
  - **PlanetRadius:** planet radius in jupiter radii
  - **OrbitalPeriod:** orbital period in days
  - **PlanetSemiMajorAxis:** semi major axis in astronomical units
  - **OrbitalAngle:** orbital angle in degrees
  - **LimbDarkeningCoefficients:** limb darkening coefficients, two for quadratic and four for a nonlinear law
