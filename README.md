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

10. Add new packages by editing the environment file and running:

	<pre>mamba env update -f environment.yml --prune</pre>

*Adapted from the [Toolbox Workshop](https://toolbox.pep-dortmund.org/install/linux/).*

> Mamba is a drop in replacement for Conda that is significantly faster at resolving dependencies
> and installing packages due to utilizing `C++` and parallel processing.

To compile `LaTeX` documents, a current `TeX Live` installation is required. For this case, a full setup is used.

1. Change directories into the default installation location:

   <pre>cd ~/.local</pre>

2. Download the installer into the directory and unpack it:

	<pre>curl -L http://mirror.ctan.org/systems/texlive/tlnet/install-tl-unx.tar.gz | tar xz</pre>

3. Run the latest install script and enter `I` to proceed:

   <pre>TEXLIVE_INSTALL_PREFIX=~/.local/texlive ./install-tl-<i>version</i>/install-tl</pre>

4. Append to the path in `.bashrc` or `.zshrc` shell configuration files:

   <pre>echo 'export PATH="$HOME/.local/texlive/<i>year</i>/bin/x86_64-linux:$PATH"' >> <i>~/.shellrc</i></pre>

5. Adjust the new environment manager to keep previous packages, use dynamic mirror for updates, and reinitialize fonts:

   <pre>tlmgr option autobackup -- -1</pre>
   <pre>tlmgr option repository https://mirror.ctan.org/systems/texlive/tlnet</pre>
   <pre>luaotfload-tool --update --force</pre>

*Adapted from the [Toolbox Workshop](https://toolbox.pep-dortmund.org/install/linux/).*

### Surrender

To keep up with the demanding volume of coding tasks today, I bow to the rule of *Large Language Model* assistants.

#### Backend

1. Create an account on [*OpenRouter*](openrouter.ai) as an LLM API aggregator.
2. Buy usage credits and ensure automatic payments are off to only rely on manual deposits.
3. Generate a key with an optional hard limit for added safeguarding against overcharges.

#### Frontend

1. Install *Docker* dependencies:

   <pre>sudo pacman -Syu docker docker-compose</pre>

   <pre>sudo systemctl enable --now docker.service</pre>

   <pre>sudo usermod -aG docker $USER</pre>

   <pre>newgrp docker</pre>

2. Clone *LibreChat* Repository:

   <pre>git clone https://github.com/danny-avila/LibreChat.git</pre>

   <pre>cd LibreChat</pre>

   <pre>cp .env.example .env</pre>

   <pre>cp docker-compose.override.yml.example docker-compose.override.yml</pre>

3. Open `.env` in a text editor and add the *OpenRouter* key:

   <pre>OPENROUTER_KEY=<i>key</i></pre>

4. Next, edit `docker-compose.override.yml` and check that the `librechat.yaml` volume mount is active:

   <pre>
   services:
      api:
         volumes:
           - type: bind
             source: ./librechat.yaml
             target: /app/librechat.yaml
   </pre>

5. Create `librechat.yaml` and define your models:



6. Launch the containerized application stack:

   <pre>docker compose up -d</pre>
