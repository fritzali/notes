## PSLS - PLATO Solarlike Lightcurve Simulator

### Capabilities & Modeling

### Command Flags

- `-v` prints the program version.
- `-V` makes the output verbose.
- `-P` outputs the power spectral density and lightcurve as plots.
- `-f` saves individual lightcurves for each camera instead of default averaging over all.
- `-m` averages camera groups and then merges interlaced lightcurve while taking into account temporal offset to increase time resolution in exoplanet transits,
as opposed to default averaging all sensors for noise suppression in astroseismology.
- <code>-M <i>number</i></code> sets amount of performed simulations.
- <code>-o <i>path</i></code> specfies output directory instead of default working directory.
- `--extended-plots` displays an extended set of plots.
- `--psd` saves the power spectral density associated with the lightcurve averaged over all cameras.
- `--pdf` saves plots as `.pdf` instead of `.png` default format.
- `--hdf5` saves averaged lightcurve and simulation components in `.hdf5` file.
- `--proto-sas` formats data saved in `.hdf5` file to be compatible with prototype PLATO SAS pipeline.

### Configuration Variables
