# FAST++

This is a C++ version of the popular SED fitting code [FAST](http://w.astro.berkeley.edu/~mariska/FAST.html) (Kriek et al. 2009). Below is a list of the major differences and selling points:

 - FAST++ is free software and does not require an IDL license
 - FAST++ is faster
 - FAST++ can handle larger parameter grids
 - FAST++ can work with FITS tables
 - FAST++ can be used as a drop-in replacement for FAST (same input and output formats)

If you use this code for your own work, please cite this repository and the 2009 paper of M. Kriek where FAST was first introduced.

# Install instructions

To install FAST++, you first need to build it. For this, you will need a recent C++ compiler and [CMake](https://cmake.org/). From there you ave two choices: the eazy "do it for me" road, and the more expert "I know how to install and build stuff" road.

## Eazy road

1. Download the script [install.sh](wip) and save it anywhere (it does not matter where).
2. If you do not have root access rights on your computer, open the script file and find the line:
    ```
    INSTALL_ROOT_DIR=""
    ```

   Then replace the ```""``` with a path to a directory of your choice where the program will be installed. For example:
    ```
    INSTALL_ROOT_DIR="/home/cschreib/programming/fastpp"
    ```
   
   In this case, the FAST++ executable will be placed in ```/home/cschreib/programming/fastpp/bin```. Then save the file.
   
3. Make the script file executable and run it:
    ```
    chmod +x install.sh
    ./install.sh
    ```
4. If necessary (i.e., if you have done step 2 above), do not forget to add the path to the FAST++ executable to your ```$PATH``` environment variable.

## Expert road
### With root access rights

1. Install the C library [cfitsio](http://heasarc.gsfc.nasa.gov/fitsio/fitsio.html) from your distribution's repositories.
2. Manually install the [phy++](http://cschreib.github.io/phypp/) library. You can disable all optional components, as they will not be used by FAST++.
3. Download the latest FAST++ release from [here](https://github.com/cschreib/fastpp/releases).
4. Extract the content of this archive to a directory of your choice.
5. Navigate inside the extracted directory with a terminal and run the following commands:

    ```
    mkdir build && cd build
    cmake ../
    make
    sudo make install
    ```

### Without root access rights

1. Download the libraries [cfitsio](http://heasarc.gsfc.nasa.gov/fitsio/fitsio.html) and [phy++](http://cschreib.github.io/phypp/) and configure them to be installed in a directory of your choice, which we will refer to as ```$DIR``` in the following. You can disable all optional components of the phy++ library, as they will not be used by FAST++.
2. Download the latest FAST++ release from [here](https://github.com/cschreib/fastpp/releases).
3. Extract the content of this archive to a directory of your choice.
4. Navigate inside the extracted directory with a terminal and run the following commands:

    ```
    mkdir build && cd build
    cmake ../ -DCMAKE_INSTALL_PREFIX=$DIR -DPHYPP_ROOT_DIR=$DIR -DCFITSIO_ROOT_DIR=$DIR
    make install
    ```
