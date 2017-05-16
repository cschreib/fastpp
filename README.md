# FAST++

<!-- MarkdownTOC depth=0 -->

- [Description](#description)
- [Install instructions](#install-instructions)
- [Benchmarks](#benchmarks)
    - [Hardware/software](#hardwaresoftware)
    - [Run 1: a catalog of galaxies with broadband fluxes](#run-1-a-catalog-of-galaxies-with-broadband-fluxes)
        - [Parameters](#parameters)
        - [Runs](#runs)
        - [Recorded times](#recorded-times)
        - [Memory consumption](#memory-consumption)
    - [Run 2: one galaxy with a high resolution spectrum](#run-2-one-galaxy-with-a-high-resolution-spectrum)
        - [Parameters](#parameters-1)
        - [Recorded times](#recorded-times-1)
        - [Memory consumption](#memory-consumption-1)
- [What is the trick?](#what-is-the-trick)
    - [Why is it faster?](#why-is-it-faster)
    - [Why does it use so little memory?](#why-does-it-use-so-little-memory)
    - [There must be a price to pay... did the code become more complex in C++?](#there-must-be-a-price-to-pay-did-the-code-become-more-complex-in-c)
- [Additional features](#additional-features)
    - [Controlling output to the terminal](#controlling-output-to-the-terminal)
    - [Multithreading](#multithreading)
    - [Photometric redshifts from EAzY](#photometric-redshifts-from-eazy)
    - [Monte Carlo simulations](#monte-carlo-simulations)
    - [Controlling the cache](#controlling-the-cache)
    - [Better treatment of spectra](#better-treatment-of-spectra)

<!-- /MarkdownTOC -->

# Description

This is a C++ version of the popular SED fitting code [FAST](http://w.astro.berkeley.edu/~mariska/FAST.html) (Kriek et al. 2009). Below is a list of the main selling points:

 - FAST++ is free software and does not require an IDL license.
 - FAST++ can be used as a drop-in replacement for FAST (same input and output formats).
 - FAST++ is on average 5 times faster, and up to 30 times with multi-threading.
 - FAST++ uses 5 to 600 times less memory.
 - FAST++ can handle *much* larger parameter grids.

There are two main differences with the original FAST: in dealing with redshifts from EAzY, and in estimating uncertainties. First, the EAzY redshifts and confidence intervals are used differently: by default the photo-z will not be enforced for the best-fitting solution, and only the confidence intervals are used to restrict the parameter space. The Monte Carlo simulations always use the same constraints on the redshift as the true observations. Second, uncertainties are derived directly from the scatter of the best-fitting values in the Monte Carlo simulations, rather than from their chi2 distribution in the observed grid.

If you use this code for your own work, please cite this repository, as well as Kriek et al. (2009) where FAST was first introduced.


# Install instructions

To install FAST++, you first need to build it. For this, you will need a recent C++ compiler, [git](https://git-scm.com/) and [CMake](https://cmake.org/). All of these tools are available by default on most modern linux distributions, and for MacOS you may need to install them with MacPorts. From there installing FAST++ is very easy. Navigate to a directory of your choosing where the code will be downloaded (it will be placed in a subdirectory called ```fastpp```), then execute the following commands:
```
# Download the code and dependencies
git clone --recursive https://github.com/cschreib/fastpp.git

# Compile
mkdir build && cd build
cmake ../
make
```

This will create an executable called ```fast++``` in the ```fastpp/build``` directory, which you can use immediately to replace FAST.


# Benchmarks
## Hardware/software

* 16 GB of RAM
* 3.4 TB hard drive disk
* 8 core 64bit Intel CPU
* Fedora
* C++ compiler: gcc 6.3.1 (-O3 -std=c++11)
* IDL: 8.4
* When run multithreaded, FAST++ is set to ```PARALLEL='models'```, ```MAX_QUEUED_FITS=1000``` and ```N_THREAD=8```.

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
| no_sim, no_cache            | 117           | 24.9 (**4.7** times faster)  | 17.0 (**6.8** times faster) |
| no_sim, with_cache          | 99.4          | 11.8 (**8.4** times faster)  | 4.2 (**24** times faster)   |
| with_sim, with_cache        | 1992          | 209 (**9.5** times faster)   | 58.6 (**34** times faster)  |

FAST++ is from **5** to **10** times faster in a single thread, and **7** to **34** times faster with multi-threading enabled.

### Memory consumption

For each run, peak memory consumption of the process is given in MB, as measured by ```/usr/bin/time```. For reference, a fresh IDL session uses 16 MB of memory.

| Run name                    | FAST-IDL      | FAST++ (single thread)    | FAST++ (8 threads)            |
| --------------------------- | ------------- | ------------------------- | ----------------------------- |
| no_sim, no_cache            | 195           | 22.9 (**8.5** times less) | 23.0 (**8.5** times less)     |
| no_sim, with_cache          | 43.6          | 9.9 (**4.4** times less)  | 10.1 (**4.3** times less)     |
| with_sim, with_cache        | 48.2          | 11.6 (**4.2** times less) | 12.1 (**4.0** times less)     |

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
* 1 galaxy to fit, 35 broadband fluxes and 1136 spectral channels
* 1000 Monte Carlo simulations

### Recorded times

| Run name                    | FAST-IDL      | FAST++ (single thread)       | FAST++ (8 threads)          |
| --------------------------- | ------------- | ---------------------------- | --------------------------- |
| no_sim, no_cache            | 208           | 95.2 (**2.2** times faster)  | 140 (**1.5** times faster)  |
| no_sim, with_cache          | 15.0          | 3.1 (**4.8** times faster)   | 0.9 (**17** times faster)   |
| with_sim, with_cache        | 11733         | 432 (**27.2** times faster)  | 74.4 (**158** times faster) |

FAST++ is from **2** to **30** times faster in a single thread, and **1.5** to **160** times faster with multi-threading enabled. When the cache is not created yet, the multi-threaded version is actually *slower* than the single-threaded version because the main performance bottleneck is generating models, rather than fitting them. When a thread is given a model to fit, it finishes to do so before the gridder is able to provide the next model to fit, and the thread thus has to wait. If the cache is already created (and/or if Monte Carlo simulations are enabled), the fitting stage becomes the most time-consuming process, and the advantage of the multi-threaded version becomes clear.

### Memory consumption

| Run name                    | FAST-IDL      | FAST++ (single thread)    | FAST++ (8 threads)            |
| --------------------------- | ------------- | ------------------------- | ----------------------------- |
| no_sim, no_cache            | 6216          | 16.9 (**368** times less) | 17.7 (**351** times less)     |
| no_sim, with_cache          | 6081          | 10.2 (**596** times less) | 15.7 (**387** times less)     |
| with_sim, with_cache        | 6089          | 19.1 (**319** times less) | 25.3 (**241** times less)     |

FAST++ consumes from **240** to **600** times less memory. The amount required is actually almost the same as for the other run, about 10 to 30 MB, while FAST-IDL requires 6 GB! If we had asked for a finer grid, FAST-IDL would not have been able to run.

# What is the trick?
## Why is it faster?

There are a number of reasons why FAST++ outperforms FAST-IDL in terms of speed. The fact that FAST-IDL is only single-threaded is the cause for a significant part of the difference: enabling multiple concurrent threads can reduce execution times by a factor of up to 6 (with 8 threads). Apart from this, the algorithms employed by FAST++ are mostly identical to that used in FAST-IDL (the main difference is in the way uncertainties are computed from the Monte Carlo simulations). The main improvement in execution time is thus due to a more efficient implementation of these algorithms, by reusing computations and memory and having overall faster computations.

FAST++ has the advantage of being compiled from C++ code: the program is executed directly from assembly instructions rather than being interpreted by the IDL virtual machine. FAST-IDL was written to make the best out of IDL, using vectorized instructions as much as possible (where the bulk of the work is executed by the IDL virtual machine using pre-compiled routines written, e.g., in Fortran). But there is only so much one can do, and some loops still have to be written explicitly. On the other hand, to avoid explicit loops the IDL code needs to perform some redundant calculations, create large arrays with repeated information, or unnecessary copies and temporaries. Since C++ code can execute loops very efficiently, these issues are avoided. This also makes the code more straightforward to read, since operations are executed explicitly, without having to resort to tricks.

It should be noted that FAST++ was compiled with the default options of GCC, i.e., not using architecture-dependent optimizations (SSE instructions mainly). This is a potential for even further improvements.

## Why does it use so little memory?

Reducing memory usage was the main reason behind the creation of FAST++. The architecture of the program was therefore designed from the start to avoid having to store all models in memory at once: models are generated and adjusted to the observed data one by one (or in parallel), and are discarded immediately after, when they are no longer needed. For this reason the program also does not keep the entire chi2 grid in memory.

The end result is that, in both runs described above, FAST++ never uses more than 30MB. This is less than any of the IDL runs. FAST-IDL is limited mostly by the model flux grid size, namely, the product of the number of observations (observed fluxes) by the number of templates in the library. In the case of FAST++, the main limitation will be the number of galaxies in the input catalog. For a single-threaded run on a 64bit machine, without Monte Carlo simulations, and a 30 band catalog with photo-z from EAzY, the required memory will be about 360 bytes per galaxy. If you have 8 GB of RAM to spare, this corresponds to 24 million galaxies. Actually, most of this space (66%) is used by the input fluxes and uncertainties, so as a rule of thumb if your entire photometry catalog in binary format would occupy less than half your available RAM memory (only fluxes and uncertainties), you can fit it with FAST++. If you also ask for ```N``` Monte Carlo simulations, the memory per galaxy is increased by 20 times ```N``` bytes, so 100 simulations will roughly double the required memory per galaxy.

## There must be a price to pay... did the code become more complex in C++?

FAST++ is about 104 kB (including comments and empty lines), while FAST-IDL weights 82 kB. The size of the C++ code is thus mildly larger (27% more characters). In fact, 33% of the C++ code is dedicated to read the inputs (flux catalogs, redshifts and filter database), while this makes up only 16% of the IDL version. The reason why is that the C++ code includes many consistency checks and detailed error messages so the user can solve issues with their input; most of these checks are absent in the IDL version. Taking out the input reading code, the code of FAST++ is actually 2% *smaller* than that of FAST-IDL.

At similar code size, deciding if a code is more complex than another is a subjective question, especially if the languages are not the same. While C++ has suffered from a reputation of being an obscure language for many years, the situation has vastly improved with the 2011 version of the C++ standard, on which FAST++ relies greatly. These changes made C++ much more expressive and clear, and as always without compromising on performances. It is therefore my opinion that the code of FAST++ is *not* complex nor obscure, at least not more than that of FAST-IDL, but this is the opinion of one who has always programmed in C++.


# Additional features

In addition to the base feature set provided by FAST-IDL, FAST++ has a number of additional feature, some of which you can enable in the parameter file. These features are all listed below.

## Controlling output to the terminal
 * ```VERBOSE```: possible values are ```1``` or ```0```. The default is ```0```, and the program will not print anything in the terminal. To enable terminal output this variable must be set to ```1```.

## Multithreading
 * ```N_THREAD```: possible values are ```0``` or any positive number. The default is ```0```. This determines the number of concurrent threads that the program can use to speed up calculations. The best value to choose depends on a number of parameters, but as a rule of thumb you should not set it to a number larger than the number of independent CPU cores available on your machine (e.g., ```4``` for a quad-core CPU), and it should be at least ```2``` to start seeing significant improvements. Using a value of ```1``` will still enable parallel execution for some of the code, but the overheads generated by the use of threads will probably make it slower than using no thread at all.
 * ```PARALLEL```: possible values are ```'none'```, ```'sources'``` or ```'models'```. The default is ```'none'```. This determines which part of the code to parallelize (i.e., execute in multiple threads to go faster). Using ```'none'``` will disable parallel execution. Setting the value to ```'models'``` will use the available threads (see ```N_THREAD``` above) to fit multiple models from the grid simultaneously. This is the optimal setup if you have many models in your grid but little computation to do per model (e.g., if you have very few sources, or no Monte Carlo simulations). If you have a large input catalog (more than a few hundred sources) and especially if you have enabled Monte Carlo simulations, you can set this value to ```'sources'```, in which case the code will divide the input catalog in equal parts that will be fit simultaneously.
 * ```MAX_QUEUED_FITS```: possible values are ```0``` or any positive number. The default value is ```1000```. This defines the maximum number of models that are produced and waiting to be fit at any given instant. It is only used if multithreading is enabled, since single-threaded execution will always have a single model in memory at a time. Setting this to ```0``` will remove the restriction. The goal of this parameter is to limit the amount of consumed memory: the higher the value, the more models can be present in memory at once, waiting to be processed. The default value of ```1000``` has a *very* slight impact on performances (less than 10%), so you can most often ignore this parameter. Else, you can disable it if you know your model grid has modest size and memory usage will not be an issue, on on the contrary decrease the value if your models are very large.

## Photometric redshifts from EAzY
 * ```FORCE_ZPHOT```: possible values are ```0``` or ```1```. The default is ```0```, and FAST++ will ignore the photometric redshifts obtained by EAzY. This is different from the behavior of FAST-IDL, which always forces the redshift to that derived by EAzY (except for Monte Carlo simulations). You can recover the FAST-IDL behavior by setting this value to ```1```, but note that contrary to FAST-IDL this will also affect the Monte Carlo simulations. In practice, this option amounts to treating photometric redshifts as spectroscopic redshifts. It makes more sense than the FAST-IDL behavior if you really want to enforce the EAzY redshifts.
 * ```BEST_AT_ZPHOT```: possible values are ```0``` or ```1```. The default is ```0```. Setting this to ```1``` will get you closer to the FAST-IDL behavior. This will force the best-fitting solution to be taken at the photometric redshift determined by EAzY, but the Monte Carlo simulations will still use the entire redshift range allowed (or the one constrained by the EAzY confidence intervals). This will give the closest output compared to FAST-IDL, however note that the derived confidence intervals will be determined on the entire redshift grid, and will be centered on a solution which may not be at the redshift derived by EAzY.
 * ```ZPHOT_CONF```: possible values are ```0```, ```68```, ```95``` or ```99```. The default is ```0```. If set to anything else than ```0```, this will read the corresponding confidence interval on the photometric redshift determined by EAzY, and restrict the redshift grid to this interval for each source. So if the entire redshift grid ranges from ```0.1``` to ```8.5```, but a source has ```l68 = 2.2``` and ```u68 = 3.5```, then FAST++ will only consider redshifts between ```2.2``` and ```3.5``` for this source.

To get the closest behavior to that of FAST-IDL, you should set ```C_INTERVAL=68```, ```BEST_AT_ZPHOT=1``` and ```ZPHOT_CONF=68``` (you can replace ```68``` by ```95``` or ```99```).

## Monte Carlo simulations
 * ```SAVE_SIM```: possible values are ```0``` or ```1```. The default is ```0```. If set to ```1```, the program will save the best fitting parameters in the Monte Carlo simulations of each galaxy in a FITS table located at ```best_fits/[catalog]_[source].sim.fits```. This table contains the values of all fitting parameters and the chi2. This will consume some more disk space, but will not slow down the program significantly. It can be useful to identify covariances and degeneracies that are not apparent from the confidence intervals printed in the output catalog.
 * ```BEST_FROM_SIM```: possible values are ```0``` or ```1```. The default is ```0```, and the best fitting solution will be chosen as the one providing the smallest chi2 value in the grid. If set to ```1```, the program will instead determine the best solution from the median of all the Monte Carlo simulations. This will ensure that the "best fit" values are more consistent with the confidence intervals (i.e., usually more centered), and erases large fluctuations when multiple solutions with very different fit parameters lead to very close chi2 values. This typically happens for galaxies with poor photometry: there are a large number of models which give similarly good chi2, but one of them has a chi2 better by a very small amount (say 0.001) and it thus picked as the "best fit".

## Controlling the cache
 * ```NO_CACHE```: possible values are ```0``` or ```1```. The default is ```0```, and the program will read and/or create a cache file, storing the pre-computed model fluxes for reuse. If you are changing your grid often or if the grid is very large and you do not want to store it on the disk, you can set this value to ```1``` and the program will neither read from nor write to the cache. Because it avoids some IO operations, it may make the program faster when the grid has to be rebuilt.

## Better treatment of spectra

The file format for spectra in FAST-IDL is not well defined. The spectra must be given on a wavelength grid where only the central wavelength of each spectral element is specified. FAST-IDL then assumes that this wavelength grid is contiguous (no gap) and uniform (all spectral elements are defined on the same delta_lambda). No check is made to ensure this is true.

For the sake of backward compatibility, this format is also available in FAST++, with the same assumptions (this time with explicit checks). However you can also use another, more explicit and flexible syntax where the lower and upper range covered by a spectral element are provided in place of the central wavelength.

A picture is worth a thousand words:
```
Wavelength axis:   [-----|-----|-----|-----]
Value (um):        2    2.1   2.2   2.3   2.4
Spectral element:     0     1     2     3

Spectrum file (FAST-IDL format):
# bin    wl     tr     F1     E1
  0     2.05    1      1.0    0.5
  1     2.15    1      2.5    0.6
  2     2.25    1      0.8    0.5
  3     3.35    1      3.2    0.4

Spectrum file (FAST++ format):
# bin  wl_low  wl_up   F1     E1
  0     2.0     2.1    1.0    0.5
  1     2.1     2.2    2.5    0.6
  2     2.2     2.3    0.8    0.5
  3     2.3     2.4    3.2    0.4
```

In this example the information in both catalogs in the same. But the new syntax allows more possibilities, for example adaptive binning, or combining spectra from different instruments (or passbands) with large gaps in between or different spectral resolutions. The ```tr``` (transmission) column, which is just a binary "use/don't use" flag becomes useless since the grid does not need to be uniform anymore.

Even when using the old format, the treatment of these spectrum files is also more correct in FAST++. The ```bin``` column correctly combines multiple data points into a single measurement (using inverse variance weighting) rather than simply using the first value of a bin (why this is implemented in this way in FAST-IDL, I do not know).
