# gLinear

## Install libraries that gLinear depends on

First the necessary libraries that `gLinear` depends on. This guide assumes you are running the latest version of one the supporting operating systems, i.e. Ubuntu 14.04 (Trusty Tahr) or OS X 10.8 (Mountain Lion).

For Ubuntu:

    sudo apt-get install build-essential liblapack-dev libblas-dev libexpat1-dev zlib1g-dev cmake

For OS X download and install the latest release of CMake for OS X [here](http://cmake.org/cmake/resources/software.html). Let the installer create symbolic links to the CMake command line tools in /usr/bin or /usr/local/bin.

Alternatively, you can install CMake via one of the package managers available for OS X, [Homebrew](http://mxcl.github.com/homebrew/) or [MacPorts](http://www.macports.org/).

## Clone gLinear from bitbucket

Use the bitbucket clone button to get a local copy of the source code -- follow the instructions given there. This will create a `glinear` directory in your local directory. Amongst others it contains the source code as well as a git repository. In the following we are going to assume that your initial local directory is `~/devel` in which case bitbucket would have created a `~/devel/glinear` directory. If you have done differently you will have to replace this with your actual gLinear directory in the instructions below. Lets now first fix bitbucket's case confusion (which also has to cater for windoze's case-blindness):

    mv glinear gLinear
    cd gLinear

## Building the library

First you need to configure `gLinear` using CMake via the `configure` script in the cmake subfolder:

    ./cmake/configure -g Release

If the `configure` script is not found it may be necessary to initialise and update the git submodule first:

    git submodule init
    git submodule update

The `configure` script will print the necessary commands required to continue the build process.

For Ubuntu the output might look something like this:

    cd build/Linux-Release-x86_64
    make -j3
    sudo make install

And on OS X:

    cd build/Darwin-Release-x86_64
    make -j3
    make install

To compile it, copy and paste the first two of these printed lines. It is not necessary to perform the installation as well, i.e. `sudo make install` or `make install`, but it may be convenient for you to do so. If at this stage your system still complains about packages that it can not find, use your favourite package manager to install them.

## Add the library to the search path

If you chose not the install `gLinear` you will need to add a path to the directory containing it to your `PATH` environment variable in order to allow applications to find and use the library. A flexible way to do this, is to first create a `bin` directory in your home directory and include that in your `PATH`:

    cd ~/
    mkdir bin
    cd ~/bin
    ln -sf $HOME/devel/gLinear/build/$(uname -s)-Release-$(uname -m)/libgLinear.so
    
If this `~/bin` directory is not already in your search `PATH`, add it by including the following line in your ~/.bashrc file:

    export PATH=$HOME/bin:$PATH

New terminal windows will automatically be aware of this updated path. Existing terminal windows can be made aware of this change by doing a:

    source ~/.bashrc

The above should work on both Ubuntu and OS X provided you did not deviate from the instructions above, and if so, change the path accordingly.

## Test it:

Open a new terminal window and test whether `gLinear` has been successfully installed or set up by running the following command:

    cd $HOME/devel/gLinear/build/$(uname -s)-Release-$(uname -m)/
    ./testnew
    
### If your system for some reason still can not find the library

Do the following ONLY if at this stage your system complains that it still can't find `libgLinear.so`  (horrible hack)

    cd /usr/local/lib
    sudo ln -sf <full_path_to_the_above_bin_directory>/libgLinear.so
    sudo ldconfig -v 


## Debug builds

The instructions above will produce binaries built using the *Release* configuration, in other words the binaries will be optimised for speed. However, if you need to debug serious problems in `patrecII` or your own code that uses `patrecII`, you might want to build the libraries and executables in the *Debug* configuration. In the *Debug* configuration optimisations will be turned off and symbols will be preserved. This will be handy if you plan to step through your code and inspect variables using a debugger like `gdb`.

To build `gLinear` and `patrecII` in *Debug* mode, simply omit the `-g Release` option from the configure commands above:

    ./cmake/configure

In addition this, you can also turn on index checking, which will detect and notify you of most array indexing errors, which is a very common cause of otherwise hard to debug crashes. To turn on index checking simply specify the `-i` option:

    ./cmake/configure -i

Turning on index checking will dramatically slow down your code, but it is advisable you do so if you are struggling to resolve an issue. It is important to note that you can also turn on index checking for *Release* builds:

    ./cmake/configure -ig Release
    
## Editor settings

Please set your editor so that your code will:
  contain no trailing whitespace,
  uses spaces instead of tabs,
  uses unix lf and not windows crlf. 
These things complicate file comparisons.   
     