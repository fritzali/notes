## PLUTO

This directory documents the installation and basic usage of the PLUTO code.

### Setup

1. Download the current release here: [The PLUTO Code for Astrophysical GasDynamics](https://plutocode.ph.unito.it/)
2. Unpack the archive in the desired location: <code>tar -xvf <i>patch</i>.tar.gz</code>
3. Change the environment variable in your shell configuration file like `~/.bashrc` or `~/.zshrc` to point at your PLUTO version directory: <code>export PLUTO_DIR=/<i>install</i>/<i>location</i>/PLUTO</i></code>
4. Modify the same file to include some useful shortcuts for running the setup script and multicore runs:
   
   ```alias plutosetup='python $PLUTO_DIR/setup.py'```
   
   ```
      plutorun() {
          touch pluto.0.log
          mpirun -np 6 ./pluto "$@" > /dev/null 2>&1 &
          tail -f --pid=$! pluto.0.log
      }
   ```

> *The former of these compiles your run, while latter uses six parallel processes for the run specified in the `pluto` executable located in the current directory. This is a good default for a processor with eight physical cores on consumer grade hardware.*

### Run

1. For a basic simulation, you must copy the following files to your working directory: `init.c`, `pluto.ini`, `definitions.h`
2. Only `definitions.h` is ever modified in the working directory, other changes are made directly in the PLUTO source.
3. Once complete, we can now modify settings via the interactive `plutosetup` script and generate our makefile by choosing: `Linux.gcc.defs`
4. To compile the code, in the working directory run: `make`

> *I ran into a problem here, where `drand` could not be called and the compiler suggested `srand` instead. To fix this, I inserted `#define _XOPEN_SOURCE 700` into `$PLUTO_DIR/Src/Math_Tools/math_random.c*` at the very top before the `#include "pluto.h"` statement.*

5. Finally, run the executable: `plutorun`

### Environment

1. Create: `mamba env create -f environment.yaml`
2. Update: `mamba env update -f environment.yaml --prune`
3. Activate: `mamba activate pluto`
