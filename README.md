## Notebooks

#### Cloning

<pre>git clone https://fritzali:<i>PAT</i>@github.com/fritzali/notes.git notebooks</pre>

#### Installing

This repository includes `Jupyter Notebook` documents relating to class notes and general experimentation.
The following describes how to install an all purpose `conda` environment.

1. Change into the `~/.local` directory:

   <pre>cd */.local</pre>

2. Download the system appropriate installer script:

   <pre>curl -LO "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"</pre>

3. Run the installation into the specified location:

   <pre><i>shell</i> Miniforge3-$(uname)-$(uname -m).sh -p ~/.local/conda</pre>

   Here, use `zsh` or `bash` as a shell, purely based on preference.
   
4. After this, expand the license agreement by entering `↵` and type `yes` to confirm. Proceed with `↵` as
   the installation location has been chosen in the previous command.

5. Next, automatically update your shell profile for `conda` initialization:

   <pre>mamba shell init</pre>

   If `mamba` cannot be found, export `~/.local/miniforge3/bin/mamba` and `~/.local/miniforge3` as the `MAMBA_EXE` and `MAMBA_ROOT_PREFIX`
   variables in your shell, respectively, or try the legacy option

   <pre>. $HOME/.local/miniforge3/etc/profile.d/mamba.sh</pre>

   as a line in your shell configuration file.

6. To initialize the `mamba` interface, also include

   <pre>eval "$(mamba shell hook --shell <i>shell</i>)"</pre>

   in your shell configuration.

7. Test the installation by creating and activating an environment, such as from the included file,

   <pre>mamba create -f environment.yml<br>mamba activate <i>environment</i></pre>

   and running any command of interest inside an `ipython` prompt.

8. After confirming success, remove the installer script:

   <pre>rm ~/.local/Miniforge3-*.sh</pre>

9. Keep the environment up to date:

   <pre>mamba update -n notes --all</pre>
