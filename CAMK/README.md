## Centrum Astronomiczne imienia Mikołaja Kopernika Polskiej Akademii Nauk<br>Nicolaus Copernicus Astronomical Center of the Polish Academy of Sciences

### Account

1. To access the email account, visit [this](https://webmail.camk.edu.pl) webmail.
2. Enter <code><i>username</i></code> and <code><i>password</i></code> when prompted.

### Cluster

1. Run the `openssh` command <code>ssh <i>username</i>@ssh.camk.edu.pl</code> and provide your <code><i>password</i></code> to enter the <code><i>username</i>@gatekeeper</code> shell. Alternatively and more securely, one can use `ssh-keygen` and <code>ssh-copy-id <i>username</i>@ssh.camk.edu.pl</code> for access to the server as [this](https://www.camk.edu.pl/en/camknet/access/) site explains in more detail.
2. From here, enter the actual compute command prompt with <code>ssh <i>host</i></code> and your <code><i>password</i></code> for setting up your environment.
3. In <code><i>username</i>@<i>host</i></code> you are in your own separate *Ubuntu* partition and can set up installations and runs however you like.
4. Check `lscpu` for available hardware. Use `logout` to exit the connection.

#### Example Installation

1. Setting up the terminal emulator, run <code>infocmp -a | ssh <i>username</i>@monster "tic -x -o ~/.terminfo -"</code> on your local machine.
2. Since packages such as `gcc` and `python` should come preinstalled, you only need to get your tools, like `conda` by running <br><code>curl -L -O "htt<span>ps://</span>github.com/conda-forge/miniforge/releases/latest/download/Miniforge-<i>version</i>-Linux-x86_64.sh"</code> and executing <code>bash Miniforge-<i>version</i>-Linux-x86_64.sh</code> in your home directory. Follow the installation script, reload `bash` and paste `conda config --set auto_activate_base false` to prevent automatically starting in the base environment every time.
3. To copy files and directories from the local to the remote machine, run these on your computer for simple transfers:
   
   <pre>scp /<i>local</i>/<i>path</i>/<i>file</i> <i>username</i>@<i>host</i>:/<i>remote</i>/<i>destination</i>/</pre>

   <pre>scp -r /<i>local</i>/<i>path</i>/<i>folder</i> <i>username</i>@<i>host</i>:/<i>remote</i>/<i>destination</i>/</pre>

   For larger sizes, follow this syntax instead:

   <pre>rsync -avzP /<i>local</i>/<i>path</i>/<i>file</i> <i>username</i>@<i>host</i>:/<i>remote</i>/<i>destination</i>/</pre>

   <pre>rsync -avzP /<i>local</i>/<i>path</i>/<i>folder</i>/ <i>username</i>@<i>host</i>:/<i>remote</i>/<i>destination</i>/</pre>

   When transferring from remote to local storage, simply swap the origin and destination order.

#### Screen Tool

1. If trying to leave processes like simulations running in the background without being logged into remote, use the `screen` package.
2. Type <code>screen -S <i>identifier</i></code> to open a session, use `CTRL + A` followed by `CTRL + D` to exit.
3. Any process started in the session will continue to run even after logging out. Show active screens with the `screen -ls` flag.
4. To rejoin the session, run <code>screen -r <i>identifier</i></code> with the same identity as specified before.
5. And to terminate a detached session, <code>screen -X -S <i>identifier</i> quit</code> can be used.
