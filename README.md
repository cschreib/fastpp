# FAST++

This is a C++ version of the popular SED fitting code [FAST](http://w.astro.berkeley.edu/~mariska/FAST.html) (Kriek et al. 2009). Below is a list of the main selling points:

 - FAST++ is free software and does not require an IDL license.
 - FAST++ can be used as a drop-in replacement for FAST (same input and output formats).
 - FAST++ is on average 5 times faster, and up to 30 times with multi-threading.
 - FAST++ uses 5 to 600 times less memory.
 - FAST++ can handle *much* larger parameter grids.

There are two main differences with the original FAST: in dealing with redshifts from EAzY, and in estimating uncertainties. First, the EAzY redshifts and confidence intervals are used differently: by default the photo-z will not be enforced for the best-fitting solution, and only the confidence intervals are used to restrict the parameter space. The Monte Carlo simulations always use the same constraints on the redshift as the true observations. Second, uncertainties are derived directly from the scatter of the best-fitting values in the Monte Carlo simulations, rather than from their chi2 distribution in the observed grid.

If you use this code for your own work, please cite this repository, as well as Kriek et al. (2009) where FAST was first introduced.

# Install instructions

To install FAST++, you first need to build it. For this, you will need a recent C++ compiler and [CMake](https://cmake.org/). From there you ave two choices: the easy "do it for me" road, and the more expert "I know how to install and build stuff" road.

## Easy road

1. Download the script [install.sh](wip) and save it anywhere (it does not matter where).
2. By default the script will install FAST++ in the standard directories of your operating system. This will probably require you to have root permission. If you do not have this permission, or if you want to install the program elsewhere, open the script file and find the line:
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


# Benchmarks
## Hardware/software

* 16 GB of RAM
* 3.4 TB hard drive disk
* 8 core 64bit Intel CPU
* Fedora
* C++ compiler: gcc 6.3.1 (-O3 -std=c++11)
* IDL: 8.4

## Run 1: a catalog of galaxies with broadband fluxes
### Parameters

The point of this run is to test the performance of FAST++ when applied to catalogs of galaxies with a number of broadband fluxes. For the purpose of the test, the redshifts derived by EAzY were not used, and the entire grid was explored. This is done to have relatively long computing time (to avoid being dominated by run-to-run fluctuations), and also to avoid the complex behavior of FAST-IDL when coupled with EAzY output (a behavior which is only partly emulated in FAST++). The memory consumption will be small in all cases, so the main test will be to measure the execution speed.

* SFH: delayed
* Stellar population model: BC03
* Metallicity: 0.02 (solar)
* tau: 8.5 to 10, step 0.5 (4 values)
* age: 8 to 10, step 0.2 (11 values)
* Av: 0 to 3, step 0.1 (31 values)
* z: 0 to 6, step 0.05 (121 values)
* 165044 models to fit (total grid size: 6.6MB)
* 1067 galaxies to fit, 7 broadband fluxes per galaxy
* 100 Monte Carlo simulations per galaxy

### Runs

Both FAST-IDL and FAST++ were ran from the same parameter file to fit a catalog of about a thousand galaxies. To test different scenarios, the codes were ran once without any cache, then a second time using the cached model fluxes, and a third time with Monte Carlo simulations enabled.

### Recorded times

For each run, execution times are given in seconds, as measured by ```/usr/bin/time```. This includes all steps of the computation, including the startup time of IDL (3 seconds).

| Run name                    | FAST-IDL      | FAST++ (single thread)       | FAST++ (8 threads)          |
| --------------------------- | ------------- | ---------------------------- | --------------------------- |
| no_z, no_sim, no_cache      | 117           | 24.9 (**4.7** times faster)  | 17.0 (**6.8** times faster) |
| no_z, no_sim, with_cache    | 99.4          | 11.8 (**8.4** times faster)  | 4.2 (**24** times faster)   |
| no_z, with_sim, with_cache  | 1992          | 209 (**9.5** times faster)   | 58.6 (**34** times faster)  |

FAST++ is from **5** to **10** times faster in a single thread, and **7** to **34** times faster with multi-threading enabled.

### Memory consumption

For each run, peak memory consumption of the process is given in MB, as measured by ```/usr/bin/time```. For reference, a fresh IDL session uses 16 MB of memory.

| Run name                    | FAST-IDL      | FAST++ (single thread)    | FAST++ (8 threads)            |
| --------------------------- | ------------- | ------------------------- | ----------------------------- |
| no_z, no_sim, no_cache      | 195           | 22.9 (**8.5** times less) | 23.0 (**8.5** times less)     |
| no_z, no_sim, with_cache    | 43.6          | 9.9 (**4.4** times less)  | 10.1 (**4.3** times less)     |
| no_z, with_sim, with_cache  | 48.2          | 11.6 (**4.2** times less) | 12.1 (**4.0** times less)     |

FAST++ consumes from **4** to **9** times less memory.


## Run 2: one galaxy with a high resolution spectrum
### Parameters

This run is meant to test the memory consumption using a large model flux grid. To do so we use a template grid a bit bigger than in the previous run, but also a much larger number of observed data points (about a thousand) using a high-resolution spectrum.

* SFH: delayed
* Stellar population model: BC03
* Metallicity: 0.008, 0.02 (solar), 0.05
* tau: 6.5 to 9, step 0.1 (14 values)
* age: 6.1 to 9.5, step 0.1 (35 values)
* Av: 0 to 2, step 0.1 (14 values)
* z: 3.710 to 3.720, step 0.001 (11 values)
* 226380 models to fit (total grid size: 1.1GB)
* 1 galaxy to fit, 35 broadband fluxes and 1136 spectral chanels
* 1000 Monte Carlo simulations

### Recorded times

| Run name                    | FAST-IDL      | FAST++ (single thread)       | FAST++ (8 threads)          |
| --------------------------- | ------------- | ---------------------------- | --------------------------- |
| no_z, no_sim, no_cache      | 208           | 95.2 (**2.2** times faster)  | 140 (**1.5** times faster)  |
| no_z, no_sim, with_cache    | 15.0          | 3.1 (**4.8** times faster)   | 0.9 (**17** times faster)   |
| no_z, with_sim, with_cache  | 11733         | 432 (**27.2** times faster)  | 74.4 (**158** times faster) |

FAST++ is from **2** to **30** times faster in a single thread, and **1.5** to **160** times faster with multi-threading enabled. When the cache is not created yet, the multi-threaded version is actually *slower* than the single-threaded version because the main performance bottleneck is generating models, rather than fitting them. When a thread is given a model to fit, it finishes to do so before the gridder is able to provide the next model to fit, and the thread thus has to wait. If the cache is already created (and/or if Monte Carlo simulations are enabled), the fitting stage becomes the most time-consuming process, and the advantage of the multi-threaded version becomes clear.

### Memory consumption

| Run name                    | FAST-IDL      | FAST++ (single thread)    | FAST++ (8 threads)            |
| --------------------------- | ------------- | ------------------------- | ----------------------------- |
| no_z, no_sim, no_cache      | 6216          | 16.9 (**368** times less) | 17.7 (**351** times less)     |
| no_z, no_sim, with_cache    | 6081          | 10.2 (**596** times less) | 15.7 (**387** times less)     |
| no_z, with_sim, with_cache  | 6089          | 19.1 (**319** times less) | 25.3 (**241** times less)     |

FAST++ consumes from **240** to **600** times less memory. The amount required is actually almost the same as for the other run, about 10 to 30 MB, while FAST-IDL requires 6 GB! If we had asked for a finer grid, FAST-IDL would not have been able to run.

# What is the trick?
## Why is it faster?

There are a number of reasons why FAST++ outperforms FAST-IDL in terms of speed. The fact that FAST-IDL is only single-threaded is the cause for a significant part of the difference: enabling multiple concurrent threads can reduce execution times by a factor of up to 6 (with 8 threads). Apart from this, the algorithms employed by FAST++ are mostly identical to that used in FAST-IDL (the main difference is in the way uncertainties are computed from the Monte Carlo simulations). The main improvement in execution time is thus due to a more efficient implementation of these algorithms, by reusing computations and memory and having overall faster computations.

FAST++ has the advantage of being compiled from C++ code: the program is executed directly from assembly instructions rather than being interpreted by the IDL virtual machine. FAST-IDL was written to make the best out of IDL, using vectorized instructions as much as possible (where the bulk of the work is executed by the IDL virtual machine using pre-compiled routines written, e.g., in Fortran). But there is only so much one can do, and some loops still have to be written explicitly. On the other hand, to avoid explicit loops the IDL code needs to perform some redundant calculations, create large arrays with repeated information, or unnecessary copies and temporaries. Since C++ code can execute loops very efficiently, these issues are avoided. This also makes the code more straightforward to read, since operations are executed explicitly, without having to resort to tricks.

It should be noted that FAST++ was compiled with the default options of GCC, i.e., not using architecture-dependent optimizations (SSE instructions mainly). This is a potential for even further improvements.

## Why does it use so little memory?

Reducing memory usage was the main reason behind the creation of FAST++. The architecture of the program was therefore designed from the start to avoid having to store all models in memory at once: models are generated and adjusted to the observed data one by one (or in parallel), and are discarded immediately after, when they are no longer needed. For this reason the program also does not keep the entire chi2 grid in memory.

The end result is that, in both runs described above, FAST++ never uses more than 30MB. This is less than any of the IDL runs. FAST-IDL is limited mostly by the model flux grid size, namely, the product of the number of observations (observed fluxes) by the number of templates in the library. In the case of FAST++, the main limitation will be the number of galaxies in the input catalog. For a single-threaded run on a 64bit machine, without Monte Carlo simulations, and a 30 band catalog with photo-z from EAzY, the required memory will be about 360 bytes per galaxy. If you have 8 GB of RAM to spare, this corresponds to 24 million galaxies. Actually, most of this space (66%) is used by the input fluxes and uncertainties, so as a rule of thumb if your entire photometry catalog in binary format would occupy less than half your available RAM memory (only fluxes and uncertainties), you can fit it with FAST++. If you also ask for "N" Monte Carlo simulations, the memory per galaxy is increased by 20 times "N" bytes, so 100 simulations will roughly double the required memory per galaxy.

## There must be a price to pay... did the code become more complex in C++?

FAST++ is about 104 kB (including comments and empty lines), while FAST-IDL weights 82 kB. The size of the C++ code is thus mildly larger (27% more characters). In fact, 33% of the C++ code is dedicated to read the inputs (flux catalogs, redshifts and filter database), while this makes up only 16% of the IDL version. The reason why is that the C++ code includes many consistency checks and detailed error messages so the user can solve issues with their input; most of these checks are absent in the IDL version. Taking out the input reading code, the code of FAST++ is actually 2% *smaller* than that of FAST-IDL.

At similar code size, deciding if a code is more complex than another is a subjective question, especially if the languages are not the same. While C++ has suffered from a reputation of being an obscure language for many years, the situation has vastly improved with the 2011 version of the C++ standard, on which FAST++ relies greatly. These changes made C++ much more expressive and clear, and as always without compromising on performances. It is therefore my opinion that the code of FAST++ is *not* complex nor obscure, at least not more than that of FAST-IDL, but this is the opinion of one who has always programmed in C++.
