# FAST++

[![Build Status](https://travis-ci.com/cschreib/fastpp.svg?branch=master)](https://travis-ci.com/cschreib/fastpp)

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
- [Differences with FAST](#differences-with-fast)
    - [Photometry](#photometry)
    - [Spectra](#spectra)
    - [Photo-z](#photo-z)
    - [Monte Carlo uncertainties](#monte-carlo-uncertainties)
    - [Chi2 grid](#chi2-grid)
- [Additional features](#additional-features)
    - [Controlling output to the terminal](#controlling-output-to-the-terminal)
    - [Multithreading](#multithreading)
    - [Photometric redshifts from EAzY](#photometric-redshifts-from-eazy)
    - [Monte Carlo simulations](#monte-carlo-simulations)
    - [Confidence intervals from chi2 grid](#confidence-intervals-from-chi2-grid)
    - [Controlling the cache](#controlling-the-cache)
    - [More output options](#more-output-options)
    - [Non-parametric SFH quantities](#non-parametric-sfh-quantities)
    - [Custom star formation histories](#custom-star-formation-histories)
    - [Using priors on the infrared luminosity](#using-priors-on-the-infrared-luminosity)
    - [Better treatment of spectra](#better-treatment-of-spectra)
    - [Velocity broadening](#velocity-broadening)
    - [Continuum indices](#continuum-indices)
- [Additional documentation](#additional-documentation)
    - [Adding new filters](#adding-new-filters)

<!-- /MarkdownTOC -->

# Description

This is a C++ version of the popular SED fitting code [FAST](http://w.astro.berkeley.edu/~mariska/FAST.html) (Kriek et al. 2009). Below is a list of the main selling points:

 - FAST++ is free software and does not require an IDL license.
 - FAST++ can be used as a drop-in replacement for FAST (same input and output formats).
 - FAST++ is at least 4 times faster, and up to 60 times with multi-threading.
 - FAST++ uses no more than 30 MB of RAM memory.
 - FAST++ can handle *much* larger parameter grids.
 - FAST++ can generate models with arbitrary star formation histories.
 - FAST++ can use observational constraints on dust emission.
 - ... and more! See [Additional features](#additional-features).

There are number of small differences between FAST++ and the original FAST, these are listed below in [Differences with FAST](#differences-with-fast).

If you use this code for your own work, please cite this repository, as well as Kriek et al. (2009) where FAST was first introduced.


# Install instructions

To install FAST++, you first need to build it. For this, you will need a recent C++ compiler, [git](https://git-scm.com/) and [CMake](https://cmake.org/). All of these tools are available by default on most modern linux distributions, and for MacOS you may need to install them with MacPorts. From there installing FAST++ is very easy. Navigate to a directory of your choosing where the code will be downloaded (it will be placed in a subdirectory called ```fastpp```), then execute the following commands:
```
# Download the code and dependencies
git clone https://github.com/cschreib/fastpp.git
cd fastpp

# Compile
mkdir build
cd build
cmake ../
make install
```

This will create an executable called ```fast++``` in the ```fastpp/bin``` directory, which you can use immediately to replace FAST. If you have not installed FAST, you will have to download some template libraries from [the FAST website](http://w.astro.berkeley.edu/~mariska/FAST_Download.html) before you can start fitting galaxies (note that some libraries can only be obtained by downloading FAST itself, such as ```ised_exp.pr```). The latest FAST template error function and EAzY filter response database are provided with FAST++ in the ```fastpp/share``` directory.

If you want to use arbitrary star formation histories beyond what the original FAST supports, you will have to download the single stellar populations libraries from Bruzual & Charlot (2003). To simplify this task for you, a script is provided in the ```fastpp/share/libraries``` folder. This script will download the libraries and rename the files to FAST++ convention for you, all you have to do is run this script and it will take care of the rest. You will need about 400 MB of free disk space.


# Benchmarks
## Hardware/software

* 16 GB of RAM
* 3.4 TB free hard drive disk space
* 12 core 64bit Intel Xeon CPU (3.5 GHz)
* Fedora 26
* C++ compiler: gcc 6.4.1 (-O3 -std=c++11)
* IDL: 8.4
* When run multithreaded, FAST++ is set to ```PARALLEL='generators'```, ```MAX_QUEUED_FITS=1000``` and ```N_THREAD=8```.
* These benchmarks were executed with FAST++ v1.2.
* Timings and memory usage were recorded with ```/usr/bin/time -v```

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

| Run name                    | FAST-IDL      | FAST++ (single thread)       | FAST++ (8 threads)           |
| --------------------------- | ------------- | ---------------------------- | ---------------------------- |
|            no_sim, no_cache |           131 |    23.1 (**5.7** times less) |      5.5 (**24** times less) |
|          no_sim, with_cache |           106 |    14.4 (**7.3** times less) |      3.5 (**30** times less) |
|        with_sim, with_cache |          3031 |      303 (**10** times less) |      131 (**23** times less) |

FAST++ is from **6** to **10** times faster, and **23** to **30** times faster with multi-threading enabled.

### Memory consumption

For each run, peak memory consumption of the process is given in MB, as measured by ```/usr/bin/time```. For reference, a fresh IDL session uses 16 MB of memory.

| Run name                    | FAST-IDL      | FAST++ (single thread)       | FAST++ (8 threads)           |
| --------------------------- | ------------- | ---------------------------- | ---------------------------- |
|            no_sim, no_cache |           191 |    21.7 (**8.8** times less) |    24.6 (**7.7** times less) |
|          no_sim, with_cache |          41.9 |     7.5 (**5.6** times less) |     9.9 (**4.2** times less) |
|        with_sim, with_cache |          46.1 |    12.5 (**3.7** times less) |      15 (**3.1** times less) |

FAST++ consumes from **3** to **9** times less memory.


## Run 2: one galaxy with a high resolution spectrum
### Parameters

This run is meant to test the memory consumption using a large model flux grid. To do so we use a template grid a bit bigger than in the previous run, but also a much larger number of observed data points (about 700) using a high-resolution spectrum.

* SFH: delayed
* Stellar population model: BC03
* Metallicity: 0.008, 0.02 (solar), 0.05
* tau: 6.5 to 9, step 0.1 (14 values)
* age: 6.1 to 9.5, step 0.1 (35 values)
* Av: 0 to 2, step 0.1 (14 values)
* z: 3.710 to 3.720, step 0.001 (11 values)
* 226380 models to fit (total grid size: 1.1GB)
* 1 galaxy to fit, 35 broadband fluxes and 644 spectral channels
* 1000 Monte Carlo simulations

### Recorded times

| Run name                    | FAST-IDL      | FAST++ (single thread)       | FAST++ (8 threads)           |
| --------------------------- | ------------- | ---------------------------- | ---------------------------- |
|            no_sim, no_cache |           160 |    39.2 (**4.1** times less) |      9.7 (**17** times less) |
|          no_sim, with_cache |           8.3 |     1.7 (**4.9** times less) |      0.5 (**18** times less) |
|        with_sim, with_cache |          5538 |      343 (**16** times less) |     87.9 (**63** times less) |

FAST++ is from **4** to **16** times faster, and **17** to **63** times faster with multi-threading enabled.

### Memory consumption

| Run name                    | FAST-IDL      | FAST++ (single thread)       | FAST++ (8 threads)           |
| --------------------------- | ------------- | ---------------------------- | ---------------------------- |
|            no_sim, no_cache |          3575 |    15.6 (**230** times less) |    19.3 (**186** times less) |
|          no_sim, with_cache |          3530 |       8 (**443** times less) |    12.3 (**286** times less) |
|        with_sim, with_cache |          3539 |    13.3 (**266** times less) |    19.7 (**180** times less) |

FAST++ consumes from **200** to **450** times less memory. The amount required is actually almost the same as for the other run, about 10 to 30 MB, while FAST-IDL requires 3.5 GB! If we had asked for a finer grid, FAST-IDL would not have been able to run on this machine.

# What is the trick?
## Why is it faster?

There are a number of reasons why FAST++ outperforms FAST-IDL in terms of speed. The fact that FAST-IDL is only single-threaded is the cause for a significant part of the difference: enabling multiple concurrent threads can reduce execution times by a factor of up to 6 (with 8 threads). Apart from this, the algorithms employed by FAST++ are mostly identical to that used in FAST-IDL (the main difference is in the way uncertainties are computed from the Monte Carlo simulations). The main improvement in execution time is thus due to a more efficient implementation of these algorithms, by reusing computations and memory and having overall faster computations.

FAST++ has the advantage of being compiled from C++ code: the program is executed directly from assembly instructions rather than being interpreted by the IDL virtual machine. FAST-IDL was written to make the best out of IDL, using vectorized instructions as much as possible (where the bulk of the work is executed by the IDL virtual machine using pre-compiled routines written, e.g., in Fortran). But there is only so much one can do, and some loops still have to be written explicitly. On the other hand, to avoid explicit loops the IDL code needs to perform some redundant calculations, create large arrays with repeated information, or unnecessary copies and temporaries. Since C++ code can execute loops very efficiently, these issues are avoided. This also makes the code more straightforward to read, since operations are executed explicitly, without having to resort to tricks.

It should be noted that FAST++ was compiled with the default options of GCC, i.e., not using architecture-dependent optimizations (SSE instructions mainly). This is a potential for even further improvements.

## Why does it use so little memory?

Reducing memory usage was the main reason behind the creation of FAST++. The architecture of the program was therefore designed from the start to avoid having to store all models in memory at once: models are generated and adjusted to the observed data one by one (or in parallel), and are discarded immediately after, when they are no longer needed. For this reason the program also does not keep the entire chi2 grid in memory.

The end result is that, in both runs described above, FAST++ never uses more than 30MB. This is less than any of the IDL runs. FAST-IDL is limited mostly by the model flux grid size, namely, the product of the number of observations (observed fluxes) by the number of templates in the library. In the case of FAST++, the main limitation will be the number of galaxies in the input catalog. For a single-threaded run on a 64bit machine, without Monte Carlo simulations, and a 30 band catalog with photo-z from EAzY, the required memory will be about 360 bytes per galaxy. If you have 8 GB of RAM to spare, this corresponds to 24 million galaxies. Actually, most of this space (66%) is used by the input fluxes and uncertainties, so as a rule of thumb if your entire photometry catalog in binary format would occupy less than half your available RAM memory (only fluxes and uncertainties), you can fit it with FAST++. If you also ask for ```N``` Monte Carlo simulations, the memory per galaxy is increased by 20 times ```N``` bytes, so 100 simulations will roughly double the required memory per galaxy.

## There must be a price to pay... did the code become more complex in C++?

The first version of FAST++ which reproduced all features of FAST-IDL was weighting about 104 kB (including comments and empty lines), while FAST-IDL weights 82 kB. Note that FAST++ has since grown larger than this, to support new features which were not available in FAST-IDL. For a similar feature level, the size of the C++ code is thus mildly larger (27% more characters). In fact, 33% of the C++ code is dedicated to read the inputs (flux catalogs, redshifts and filter database), while this makes up only 16% of the IDL version. The reason why is that the C++ code includes many consistency checks and detailed error messages so the user can solve issues with their input; most of these checks are absent in the IDL version. Taking out the input reading code, the code of FAST++ is actually 2% *smaller* than that of FAST-IDL.

At similar code size, deciding if a code is more complex than another is a subjective question, especially if the languages are not the same. While C++ has suffered from a reputation of being an obscure language for many years, the situation has vastly improved with the 2011 version of the C++ standard, on which FAST++ relies greatly. These changes made C++ much more expressive and clear, and as always without compromising on performances. It is therefore my opinion that the code of FAST++ is *not* complex nor obscure, at least not more than that of FAST-IDL, but this is the opinion of one who has always programmed in C++.


# Differences with FAST

While the implementation of FAST++ was designed to resemble FAST-IDL as much as possible, some aspects of FAST-IDL were not preserved. These are listed here.

## Photometry
 * The default "zero point" of the fluxes is 23.9 instead of 25. This means that your input fluxes are expected to be given in micro Janskys by default. You can of course specify your own "zero point" if this is not the case, as in FAST-IDL.
 * If a galaxy has no measured flux in a band (i.e., we have no information about it), the _uncertainty_ must be set to a negative value or "NaN". Note that there is a difference between "no measured flux" and "non-detection": the latter provides information that the fit can use. In FAST-IDL the _flux_ had to be set to -99, which was problematic because -99 could be a valid flux measurement (i.e., a non-detection).
 * The input photometry catalog can contain columns starting with "F" or "E" which are not flux columns, provided the "F" or "E" is not followed by numbers. In FAST-IDL a catalog containing such columns was rejected.

## Spectra
 * The template error function is now also applied to the spectra (except if ```AUTO_SCALE``` is set to ```1```, see below). In FAST-IDL this was not the case.
 * Galaxies with a spectrum can still use the input photo-z. In FAST-IDL their photo-z was ignored.
 * An alternative format for the spectrum file is allowed (see [Better treatment of spectra](#better-treatment-of-spectra) below), but the old format is still supported.
 * The ```bin``` column is used in FAST++ to combine multiple spectral elements into one before the fit, using inverse variance weighting, which can speed up the fit considerably if the spectral resolution is very high. In FAST-IDL only the first element in a bin was used, for some reason.
 * The ```AUTO_SCALE``` option behaves very differently. In FAST-IDL, enabling this option will automatically rescale the spectrum to match the broad band photometry. It will also perform the rescaling during the Monte Carlo simulations, so the confidence intervals will account for the uncertainty in the rescaling. In FAST++, enabling the option will automatically rescale the observed spectrum to the _model_ spectrum, independently of the rest of the photometry. Therefore, the spectrum does not influence the best-fitting (and Monte Carlo) value for the amplitude of the model, but it will still be counted in the chi2. This amounts to assuming that the absolute flux calibration of the spectrum is less well known than the _shape_ of this spectrum, i.e., there is an absorption feature here and a spectral break there. When this option is enabled, the template error function is not applied to the spectrum because the uncertainty on the template normalization at this wavelength has become irrelevant, only its detailed shape.

## Photo-z
 * If an input photo-z is negative, FAST++ will simply ignore the photo-z and let the redshift vary freely. In FAST-IDL a galaxy with a negative photo-z was ignored completely.
 * The EAzY redshifts and confidence intervals are used differently: by default the photo-z will not be enforced for the best-fitting solution, and only the confidence intervals are used to restrict the parameter space. This behavior can be modified to reproduce the original behavior of FAST-IDL. See [Photometric redshifts from EAzY](#photometric-redshifts-from-eazy) below.

## Monte Carlo uncertainties
 * The way FAST++ computes uncertainties using the Monte Carlo simulations is different from FAST-IDL, and this can lead to small differences in the resulting confidence intervals. In particular, the Monte Carlo simulations always use the same constraints on the redshift as the best-fit solution (unless ```BEST_AT_ZPHOT``` is set). Second, uncertainties are derived directly from the scatter of the best-fitting values in the Monte Carlo simulations, rather than from their chi2 distribution in the observed grid. This should not be significant, and the accuracy of the uncertainties computed by FAST++ was verified using mock catalogs with known physical parameters.

## Chi2 grid
 * FAST-IDL uses the IDL "save" format to write the chi2 grid on the disk. This format is proprietary and can only be opened by IDL itself. Instead, FAST++ uses a custom, open binary format, described below. In addition, a convenience tool called ``fast++-grid2fits`` is provided to convert the grid file into a FITS table.

The binary grid format is the following. The file starts with a header, and then lists the chi2 data for each model of the grid. In addition to the chi2, the model's best-fitting properties are also saved.
```
# Begin header
# ------------
[uint32]: number of bytes in header
[uchar]:  'C' (identifies the file type)
[uint32]: number of galaxies
[uint32]: number of properties (mass, sfr, etc)
  # For each property:
  [char*]: name of the property ("lmass", "lsfr", etc)
[uint32]: number of grid parameter
  # For each grid parameter:
  [char*]: name of the grid parameter ("z", "lage", "metal", etc)
  [uint32]: number of elements in the grid
  [float*]: values of the grid for this parameter
# ----------
# End header

# Begin data
# ----------
# For each model in the grid
[float*]: chi2 of each galaxy
  # For each property
  [float*]: values of the property for each galaxy
# --------
# End data
```

Strings (```[char*]```) are null-terminated (the string ends when the null character is found). Galaxies are ordered as in the input catalog. Models are ordered along the grid: the first model in the file corresponds to the first value of all the grid parameters, as defined in the header. The second model corresponds to the next value of the *last* grid parameter, and so on and so forth. For example, if the grid had only three parameters ```z=[1.0,1.2,1.4]```, ```lage=[8,9,10]``` and ```metal=[0.005,0.02,0.05]```, then models in the file would be ordered as:

```
# modelID         z      lage     metal
  0               1.0    8        0.005
  1               1.0    8        0.02
  2               1.0    8        0.05
  3               1.0    9        0.005
  4               1.0    9        0.02
  5               1.0    9        0.05
  6               1.0    10       0.005
  7               1.0    10       0.02
  8               1.0    10       0.05
  9               1.2    8        0.005
  10              1.2    8        0.02
  11              1.2    8        0.05
  ...             ...    ...      ...
```


# Additional features

In addition to the base feature set provided by FAST-IDL, FAST++ has a number of additional feature, some of which you can enable in the parameter file. These features are all listed below.

## Controlling output to the terminal
 * ```VERBOSE```: possible values are ```1``` or ```0```. The default is ```1```, and the program will print the progress of the fit in the terminal. To disable terminal output this variable should be set to ```0```.

## Multithreading
 * ```N_THREAD```: possible values are ```0``` or any positive number. The default is ```0```. This determines the number of concurrent threads that the program can use to speed up calculations. The best value to choose depends on a number of parameters, but as a rule of thumb you should not set it to a number larger than the number of independent CPU cores available on your machine (e.g., ```4``` for a quad-core CPU), and it should be at least ```2``` to start seeing significant improvements. Using a value of ```1``` will still enable parallel execution for some of the code, but the overheads generated by the use of threads will probably make it slower than using no thread at all.
 * ```PARALLEL```: possible values are ```'none'```, ```'sources'```, ```'models'```, or ```'generators'```. The default is ```'none'```. This determines which part of the code to parallelize (i.e., execute in multiple threads to go faster). Using ```'none'``` will disable parallel execution. Setting the value to ```'generators'``` will use the available threads (see ```N_THREAD``` above) to generate and fit multiple models from the grid simultaneously. This is the optimal setup if you have many models in your grid but little computation to do per model (e.g., if you have very few sources to fit, or no Monte Carlo simulations). If you have a large input catalog (more than a few hundred sources) and especially if you have enabled Monte Carlo simulations, you can set this value to ```'sources'```, in which case the code will divide the input catalog in equal parts that will be fit simultaneously. The ```'models'``` option is a compromise between the two other options: models are generated (or read from the cache) by the main thread, but are adjusted to the photometry in parallel. If the model cache exists, ```'generators'``` will fallback to ```'models'``` automatically, so you should not have to choose this option explicitly. Ultimately, the best choice depends on what is the main performance bottleneck. Are there few models, but many fits to do for each model? Then pick ```'sources'```. Are there many models to fit for each source, but few sources? Then pick ```'generators'```.
 * ```MAX_QUEUED_FITS```: possible values are ```0``` or any positive number. The default value is ```1000```. This defines the maximum number of models that are produced and waiting to be fit at any given instant. It is only used if multithreading is enabled, since single-threaded execution will always have a single model in memory at a time. Setting this to ```0``` will remove the restriction. The goal of this parameter is to limit the amount of consumed memory: the higher the value, the more models can be present in memory at once, waiting to be processed. The default value of ```1000``` has a *very* slight impact on performances (less than 10%), so you can most often ignore this parameter. Else, you can disable it if you know your model grid has modest size and memory usage will not be an issue, on on the contrary decrease the value if your models are very large.

## Photometric redshifts from EAzY
 * ```FORCE_ZPHOT```: possible values are ```0``` or ```1```. The default is ```0```, and FAST++ will ignore the photometric redshifts obtained by EAzY. This is different from the behavior of FAST-IDL, which always forces the redshift to that derived by EAzY (except for Monte Carlo simulations). You can recover the FAST-IDL behavior by setting this value to ```1```, but note that contrary to FAST-IDL this will also affect the Monte Carlo simulations. In practice, this option amounts to treating photometric redshifts as spectroscopic redshifts. It makes more sense than the FAST-IDL behavior if you really want to enforce the EAzY redshifts.
 * ```BEST_AT_ZPHOT```: possible values are ```0``` or ```1```. The default is ```0```. Setting this to ```1``` will get you closer to the FAST-IDL behavior. This will force the best-fitting solution to be taken at the photometric redshift determined by EAzY, but the Monte Carlo simulations will still use the entire redshift range allowed (or the one constrained by the EAzY confidence intervals). This will give the closest output compared to FAST-IDL, however note that the derived confidence intervals will be determined on the entire redshift grid, and will be centered on a solution which may not be at the redshift derived by EAzY.
 * ```ZPHOT_CONF```: possible values are ```0```, ```68```, ```95``` or ```99```. The default is ```0```. If set to anything else than ```0```, this will read the corresponding confidence interval on the photometric redshift determined by EAzY, and restrict the redshift grid to this interval for each source. So if the entire redshift grid ranges from ```0.1``` to ```8.5```, but a source has ```l68 = 2.2``` and ```u68 = 3.5```, then FAST++ will only consider redshifts between ```2.2``` and ```3.5``` for this source.

To get the closest behavior to that of FAST-IDL, you should set ```C_INTERVAL=68```, ```BEST_AT_ZPHOT=1``` and ```ZPHOT_CONF=68``` (you can replace ```68``` by ```95``` or ```99```).

## Monte Carlo simulations
 * ```SAVE_SIM```: possible values are ```0``` or ```1```. The default is ```0```. If set to ```1```, the program will save the best fitting parameters in the Monte Carlo simulations of each galaxy in a FITS table located at ```best_fits/[catalog]_[source].sim.fits```. This table contains the values of all fitting parameters and the chi2. This will consume some more disk space, but will not slow down the program significantly. It can be useful to identify covariances and degeneracies that are not apparent from the confidence intervals printed in the output catalog.
 * ```BEST_FROM_SIM```: possible values are ```0``` or ```1```. The default is ```0```, and the best fitting solution will be chosen as the one providing the smallest chi2 value in the grid. If set to ```1```, the program will instead determine the best solution from the median of all the Monte Carlo simulations. This will ensure that the "best fit" values are more consistent with the confidence intervals (i.e., usually more centered), and erases large fluctuations when multiple solutions with very different fit parameters lead to very close chi2 values. This typically happens for galaxies with poor photometry: there are a large number of models which give similarly good chi2, but one of them has a chi2 better by a very small amount (say 0.001) and it thus picked as the "best fit".

## Confidence intervals from chi2 grid
 * ```INTERVAL_FROM_CHI2```: possible values are ```0``` or ```1```. The default is ```0```, in which case confidence intervals (error bars) are computed using a series of Monte Carlo simulations of the data. An alternative, much faster way to obtain similar results is to read off the range of parameter values covered by models in the grid, selecting only the models that have a chi2 within a certain threshold from the best value (see Avni 1976). This can be done by setting this value to ```1```. While in most cases the resulting confidence intervals will be almost the same as with the MC simulations, there are several advantages to using this approach. The first is that, by construction, the resulting confidence intervals will always include the best fit value; this is not guaranteed with MC simulation. The second is that this method behaves better when the best fit solution is close to the edge of the grid; intervals from MC simulations are computed from percentiles and therefore may under-represent edge values. The latter is particularly important for parameters that have a natural upper/lower bound, such as parameters that cannot be negative but for which zero is a perfectly fine value (e.g., Av, or some SFH parameters). There are some drawbacks however. To obtain precise confidence intervals, the grid must be relatively fine: if the true confidence interval of a parameter is 1 +/- 0.01, but the fitting grid for that parameter only has a step of 0.1, then the code will report 1 +/- 0, which is too optimistic. Finally, in terms of performance, for this option to work the program has to write part of the chi2 grid on the disk, and the amount of space required can be large if you have a very fine grid.

## Controlling the cache
 * ```NO_CACHE```: possible values are ```0``` or ```1```. The default is ```0```, and the program will read and/or create a cache file, storing the pre-computed model fluxes for reuse. If you are changing your grid often or if the grid is very large and you do not want to store it on the disk, you can set this value to ```1``` and the program will neither read from nor write to the cache. Because it avoids some IO operations, it may make the program faster when the grid has to be rebuilt.

## More output options
 * ```SFR_AVG```: possible values are any positive number, which define the averaging time for the output SFR (in Myr). The default is ```0```, and the output SFR is the "instantaneous" SFR at the chosen age of the corresponding template, as in FAST. This is not necessarily a good choice, because photometry alone is mostly unable to distinguish variations of SFR on timescales lower than a hundred million years. For this reason in FAST++ you have the option to average the SFRs over an arbitrary interval of time prior to observation. This has no impact on the chosen best-fit SED, and only affects the value of the best-fit SFR (and its error bar).
 * ```REST_MAG```: array of integers corresponding to IDs in the filter database. Much like EAzY, FAST++ can output rest-frame magnitudes (absolute magnitudes, at 10pc), which can be used to classify galaxies into broad classes (blue vs. red, or quiescent vs. star-forming). If you list some filter IDs in ```REST_MAG```, FAST++ will add new columns to the output catalog, listing the absolute magnitudes in these filters for each galaxy of the input catalog. The magnitudes are given in the AB system, using the zero point you specified in ```AB_ZEROPOINT```.
 * ```OUTPUT_COLUMNS```: array of strings determining which columns to show in the output file. The default is to show the same columns as FAST-IDL, namely, ```['id','z','ltau','metal','lage','Av','lmass','lsfr','lssfr','la2t','chi2']```. Possible values are:
     * ```'id'```: the ID of each galaxy, as in the photometric catalog
     * ```'chi2'```: the reduced chi squared of the best fit template
     * ```'z'```: the redshift
     * ```'metal'```: the metallicity (solar value is about 0.02)
     * ```'lage'```: the log10 of the template's age (in years)
     * ```'Av'```: the dust attenuation (in mags)
     * ```'lmass'```: the log10 of the stellar mass (in units of solar mass)
     * ```'lsfr'```: the log10 of the star formation rate (in untis of solar mass per year)
     * ```'lssfr'```: the log10 of the specific star formation rate (in units of 1/year)
     * ```'ldust'```: the log10 of the dust luminosity (in units of total solar luminosity)
     * ```'ltau'```: (only for gridded SFH) the log10 of the star formation timescale (in years)
     * ```'la2t'```: (only for gridded SFH) the log10 of the ratio of age to SF timescale
     * ... plus any of the rest-frame magnitudes listed in ```REST_MAG``` (see above)
     * ... plus any of the custom SFH parameters (see ```CUSTOM_SFH``` below)
 * ```INTRINSIC_BEST_FIT```: possible values are ```0``` or ```1```. The default is ```0```. If set to ```1``` and ```BEST_FIT``` is also set to ```1```, the program will output the intrinsic best-fit SED of a galaxy (i.e., the SED of the galaxy prior to attenuation by dust) alongside the best-fitting template.
 * ```BEST_SFHS```: possible values are ```0``` or ```1```. The default is ```0```. If set to ```1```, the program will output the best fit star formation history (SFH) to a file, in the ```best_fits``` directory (as for the best fit SEDs). If Monte Carlo simulations are enabled, the program will also output confidence intervals on the SFH for each time step, as well as the median SFH among all Monte Carlo simulations. This median may not correspond to any analytical form allowed by your chosen SFH model.
 * ```SFH_OUTPUT_STEP```: possible values are any strictly positive number, which defines the size of a time step in the output SFH (in Myr). The default is ```10``` Myr.
 * ```SFH_OUTPUT```: possible values are ```'sfr'``` or ```'mass'```. The default is ```'sfr'```, and the program outputs as "SFH" the evolution of the instantaneous SFR of each galaxy with time. If set to ```'mass'```, the program will output instead the evolution of the stellar mass with time (which is usually better behaved, see Glazebrook et al. 2017). Note that the evolution of the mass accounts for mass loss, so the mass slowly _decreases_ with time after a galaxy has quenched.
 * ```SAVE_BESTCHI```: FAST++ can save the entire chi2 grid on the disk with the ```SAVE_CHI_GRID``` option. However, if you have *huge* grids, this can require too much disk space (I have been in situations where the chi2 grid would be as large as several TB!). Usually, one is not interested in the chi2 of *all* models, but only those that match the data within some tolerance threshold. This option allows you to only save on the disk the models that are worst than the best chi2 by some amount ```chi2 - best_chi2 < SAVE_BESTCHI``` (where ```SAVE_BESTCHI=1``` if you are interested in standard 68% confidence intervals, or ```2.71``` for 90% confidence, etc., see Avni 1976). If you set the ```INTERVAL_FROM_CHI2``` option, the program will use these saved grids to automatically compute parameter confidence intervals. The parameters of all these "good" models are saved in a separate ".grid" file for each galaxy of the input catalog, inside the ```best_chi2``` folder. The format is similar to the ".grid" file for the full chi2 grid (which is described above), but not identical. The ``fast++-grid2fits`` tool can also convert these files into FITS tables. The binary format is the following:
```
# Begin header
# ------------
[uint32]: number of bytes in header
[uchar]:  'B' (identifies the file type)
[uint32]: number of properties (mass, sfr, etc)
  # For each property:
  [char*]: name of the property ("lmass", "lsfr", etc)
[uint32]: number of grid parameter
  # For each grid parameter:
  [char*]: name of the grid parameter ("z", "lage", "metal", etc)
  [uint32]: number of elements in the grid
  [float*]: values of the grid for this parameter
# ----------
# End header

# Begin data
# ----------
# For each good model
[uint32]: ID of the model in the grid (= modelID, see description of SAVE_CHI_GRID)
[float]: chi2 of the model
[float*]: values of the properties of this model
# --------
# End data
```

## Custom star formation histories
In the original FAST, one has access to three star formation histories: the tau model (exponentially declining), the delayed tau model (delayed exponentially declining) and the constant truncated model (constant, then zero). All three are parametrized with the star formation timescale ```tau``` (either the exponential timescale for the first two SFH, or the duration of the star formation episode for the last one), which can be adjusted in the fit. These star formation histories are distributed as pre-gridded template libraries at various ages, which are then interpolated during the fit to obtain arbitrary ages.

In addition, one could use the ```MY_SFH``` option to supply an arbitrary star formation history to replace the three SFHs listed above, again, as a pre-gridded template library. The downside is that this must be a single SFH, not a family of SFHs (such as the declining model, parametrized with ```tau```), so it is not possible to fit for your custom SFH's parameters this way.

For this reason in FAST++ you have the option to build the template library on the fly, using an analytical formula for the SFH with any number of free parameters. These parameters will be varied on a grid, and participate in determining the best fit and confidence intervals of the other parameters, such as the mass or the SFR.

Using custom SFHs is easy. First you have write down the analytical formula in the ```CUSTOM_SFH``` parameter. This formula should return the star formation rate (in arbitrary units) as a function of time since onset of star formation ```t``` (given in years), where ```t=0``` is the birth of the galaxy. The formula can involve any math function, including: ```abs```, ```acos```, ```asin```, ```atan```, ```atan2```, ```ceil```, ```cos```, ```cosh```, ```exp```, ```floor```, ```log``` (natural logarithm), ```log10``` (base-10 logarithm), ```pow```, ```sin```, ```sinh```, ```sqrt```, ```tan```, ```tanh```, ```min```, ```max```, ```step``` (returns ```1``` if argument is zero or positive, and ```0``` otherwise). See the documentation of [tinyexpr](https://github.com/codeplea/tinyexpr) for more detail.

In addition, if your SFH is parametrized by one or more parameters, the formula can reference these parameters by name. For this to work, your parameters must be listed in the ```CUSTOM_PARAMS``` parameter as an array of strings, and for each parameter ```x``` you must provide the grid parameters as ```X_MIN```, ```X_MAX``` and ```X_STEP``` (with ```X``` in upper case by convention, but lower case works too).

For example you can replicate FAST's delayed SFH using:
```
CUSTOM_SFH = '(t/10^log_tau)*exp(-(t/10^log_tau))'
CUSTOM_PARAMS = ['log_tau']
LOG_TAU_MIN = 6.0
LOG_TAU_MAX = 10.0
LOG_TAU_STEP = 0.2
```

Then if you want to add a second burst with the same ```tau``` and intensity, but with a time delay ```tburst```, you could write:
```
CUSTOM_SFH = '(t/10^log_tau)*exp(-t/10^log_tau) + max(0, ((t - 10^log_tburst)/10^log_tau)*exp(-(t - 10^log_tburst)/10^log_tau))'
CUSTOM_PARAMS = ['log_tau','log_tburst']
LOG_TAU_MIN = 6.0
LOG_TAU_MAX = 10.0
LOG_TAU_STEP = 0.2
LOG_TBURST_MIN = 6.0
LOG_TBURST_MAX = 9.0
LOG_TBURST_STEP = 0.3
```

As you can see the analytical formula can be quite complex, and this has only a minimal impact on performances. However you should pay attention to the size of your grid: this feature allows you to have as many parameter as you wish, and sampling all of them properly may require a very large grid. FAST++ will deal with any grid size, but it may take some time.

Depending on which SFH expression you choose, there may be some parts of the parameter space that you do not want to probe, but that you cannot exclude simply from reducing the extents of the grid. For example with the custom SFH above, you could want to exclude cases where the burst happens before one e-folding time, namely, ```log_tburst < log_tau```. To specify this, you can use the option ```GRID_EXCLUDE```. This must be a custom function of the grid parameters that should return ```1``` ("true") if you want to exclude a particular combination, and ```0``` ("false") otherwise. The above scenario would then be written as:
```
GRID_EXCLUDE = 'log_tburst < log_tau'
```

Lastly, by default the variable ```t``` in the SFH expression is the "cosmic" time, with ```t=0``` being the instant when the galaxy was born, and where larger values of ```t``` correspond to later times, in the future. An alternative parametrization is to work with the "lookback" time, where ```t=0``` is the instant when the galaxy is observed, and larger values of ```t``` correspond to earlier times, back towards the Big Bang. The relation between the two is simply ```t_lookback = 10^lage - t_cosmic```, so this is a simple transformation you can perform in the SFH expression. But to make it more convenient, you can enable the option ```CUSTOM_SFH_LOOKBACK=1```, in which case ```t``` will be defined as the lookback time.


## Non-parametric SFH quantities
Generally, the grid parameters which are used to define the SFH are not particularly meaningful. For example, the 'age' parameter is the elapsed time since the birth of the very first star. This quantity is, in theory, impossible to measure; the only reason it looks like a well constrained quantity with standard fit setups is because most SFH parametrizations are quite rigid and will not explore a large enough variety of SFHs. Likewise, the 'tau' parameter in exponentially-declining SFHs is only meaningful in the case where the galaxy's true SFH is indeed an exponential. If not, interpreting this value becomes difficult.

For this reason, it is often advocated to not give any credit to these grid parameters, but instead compute non-parametric quantities. The simplest example would be to compute the elapsed time since half of the stars were born (which is essentially the mass-weighted age). FAST++ provides a number of such non-parametric quantities for each SFH, which you can ask the program to compute by listing them in the ```SFH_QUANTITIES``` option.

List of available quantities:

 - ```tsf```: shortest time interval over which 68% of SFR took place (= duration of star formation).
 - ```past_sfr```: average SFR during the above time interval (= mean past SFR).
 - ```sfr```X : average SFR over the last ```X``` Myr.
 - ```brate```X : ```sfr```X / ```past_sfr``` (= birth-rate parameter).
 - ```tquench```X : elapsed time since the SFR dropped below a factor ```X``` of the ```past_sfr``` (= time since quenching).
 - ```tform```X : elapsed time since the galaxy had formed ```X```% of its current mass (= time since formation).

Some of these parameters are defined and used in [Schreiber et al. (2018)](http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2018arXiv180702523S&db_key=PRE&link_type=ABSTRACT&high=5b4379bf1108433). Look at this paper for more information.


## Using priors on the infrared luminosity
One of the main degeneracy that arises when fitting UV-to-NIR data is that of dust versus age. When a galaxy has a red SED, unless the signal to noise and the wavelength sampling are high, it is very difficult to say if this is caused by a large amount of dust, or by an older stellar population, or a combination of both. It is for this reason that SFRs obtained from such fits are very uncertain, and that SFRs determined from the far-IR are preferred.

The typical approach is thus to use the UV-to-NIR data to determine the photometric redshift and the stellar mass, and the FIR data for the SFR. But this is inconsistent: by not using all the available information in a single fit, we may be using a stellar mass derived assuming the galaxy is very old, while it was actually dusty. It turns out that the stellar mass does not change dramatically under the age-dust degeneracy, but other parameters do (such as the star formation history, but also the redshift!).

Therefore in FAST++ you have the option to use the observed FIR luminosity to constrain the fit, and break the age-dust degeneracy. This is done by adding the FIR luminosity as an additional "data point" in the SED, and comparing it to the predicted dust luminosity. The dust luminosity is defined as the total energy that was absorbed by dust, and it should be directly comparable to the observed infrared luminosity if the chosen dust law and the assumed dust geometry are accurate. The observed FIR luminosity will participate in determining both the normalization of a model and its chi squared.

To do this you must enable the ```USE_LIR``` option, and provide a luminosity catalog file with extension ```.lir```. This file must contain three columns: ```ID``` as in the photometric catalog (sources must be in the same order), ```LIR``` is the 8-to-1000um luminosity expressed in units of [total solar luminosity](https://en.wikipedia.org/wiki/Solar_luminosity) (Lsun), and ```ELIR``` is the associated uncertainty, which by default is assumed to be Gaussian (see below). If you could not measure the luminosity for a galaxy, simply set both ```LIR``` and ```ELIR``` to zero or any invalid value (a negative number such as -99, or NaN). If you use this option, it is advisable to display ```'ldust'``` in the output catalog (see ```OUTPUT_COLUMNS``` above).

Since you have to provide the luminosity, this implicitly assumes that you know the redshift of the galaxy in advance... In the future FAST++ may be able to fit the FIR photometry directly, so this restriction will eventually be removed.

**Warning.** The method described above assumes that the uncertainty on the FIR luminosity is Gaussian. This is a good assumption only if the S/N of the FIR luminosity measurement is very low (i.e., because the source is not detected on the FIR images), and/or if the dust temperature is well known (i.e., within +/- 1 or 2 Kelvins). A common situation that does *not* satisfy these criteria is when only a single sub-millimeter flux measurement is available to derive the FIR luminosity. This is typical with ALMA observations, for example. In such cases, the dominant source of uncertainty is the unknown dust temperature, and the uncertainty on the FIR luminosity is highly non-Gaussian. Instead, a better assumption is that the noise is log-normally distributed, or in other words, that the uncertainty on log(LIR) is Gaussian.

FAST++ can deal with log-normal uncertainties too. To enable this feature, you need to specify an extra column in the ```.lir``` catalog. This column must be called ```LOG```, and must contain values that are either ```0``` or ```1``` for each source in the catalog. A value of ```0``` means that the FIR luminosity is given in natural units, with Gaussian uncertainties. A value of ```1``` means that the FIR luminosity is given in log10 units, with Gaussian uncertainties on the log. Note that, in this case, the luminosity does not participate in determining the model normalization (because the minimization problem would become non-linear), and instead only contributes to the chi squared. Here is an example ```.lir``` file which displays the two methods:

```
# ID    LIR    ELIR LOG
  1  1.2e12  0.5e12 0
  2   12.15    0.23 1
  3  3.5e12  1.2e12 0
```

In this example, the fit for source ```2``` will assume Gaussian errors on the log(luminosity). Note how the *values* in the ```LIR``` and ```ELIR``` columns are different, since now the luminosity must be given in log10 units, and the uncertainty is the uncertainty on the log (in dex).

Which value of ```LOG``` you should use depends on what probability distribution you estimate for the FIR luminosity. I give here a few guidelines to help you decide. If you are limited by the S/N on your flux measurements (i.e., for non-detections or weak detections) then ```0``` is preferable. If you are limited by the uncertainty in the conversion from flux to luminosity (i.e., when Tdust is not known) then ```1``` is preferable. Lastly, if you have high S/N fluxes and a well constrained Tdust, then either ```0``` or ```1``` will be fine.

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

Even when using the old format, the treatment of these spectrum files is also more correct in FAST++. The ```bin``` column correctly combines multiple data points into a single measurement (using inverse variance weighting) rather than simply using the first value of a bin (why this is implemented in this way in FAST-IDL, I do not know). The order of the columns in the file do not matter, while FAST-IDL assumes a fixed format (but does not tell you).

## Velocity broadening
By default, the template spectra produced by FAST++ assume no velocity dispersion of the stars inside the galaxy. This is a good approximation when only fitting broadband photometry, but it becomes incorrect when fitting high spectral resolution data, such as narrow band photometry or spectra, which can sample the spectral profile of the various absorption lines that compose the spectrum (i.e., mostly the Balmer series).

The solution implemented in FAST++ is to broaden all the spectral templates with a fixed velocity dispersion, specified in ```APPLY_VDISP``` (in km/s). Note that this is the value of the velocity *dispersion* (i.e., the "sigma" of the Gaussian velocity profile), not the FWHM. Unfortunately, because of the architecture of the code, it is not possible to specify different velocity dispersions for each galaxy of the input catalog. If you need to do this, you will have to fit each galaxy separately, in different FAST++ runs.

NB: applying the velocity dispersion needs to be done for each SED library that is used in the fit. This incurs a small performance penalty at the beginning of the fit (and whenever the library is changed, e.g., when switching to a new metallicity). If your grid is small, this step can actually dominate the total computation time.


## Continuum indices
Traditionally, absorption line strengths are measured as "equivalent widths" (in wavelength unit). These are part of a large set of continuum features call "spectral indices", see for example [Balogh et al. (1999)](http://adsabs.harvard.edu/abs/1999ApJ...527...54B), which includes other commonly used quantities such as the Dn4000 index. These features have been used extensively in the past to characterize the SEDs of galaxies and infer their star formation histories and/or metal abundances.

When FAST++ is given a spectrum to fit, it implicitly uses these various continuum features to constrain the models. This does not require any additional input on your part. Yet, one may wish to inspect how well each of these features is reproduced by the models, as this may signal a shortcoming in the set of allowed models (e.g., requiring models of different metallicity), limitations of the adopted stellar library, or inconsistencies in the data. Otherwise, one may wish to make predictions for features that are not observed, to see how well they are constrained by the presently available data (e.g., if you only have photometry, or only UV spectroscopy), or even to subtract the absorption component from emission line measurements.

For this reason, FAST++ allows you to define a set of continuum indices that will be computed for each model SED, the values of which you can output either in the final catalog, or in the various binary grid files. These indices must be listed in a file, and you must specify the path to this file in the variable ```CONTINUUM_INDICES```.

The format of the file is the following. Each feature is placed on a separate line, formatted as ```name = type, params...```. The line must start with the name of the feature, which cannot contain blank spaces (e.g., ```hdelta```), and must be followed by an equal sign. Then, you must specify the "type" of the feature, which must be either ```abs``` for absorption lines, or ```ratio``` for a flux ratio (such as Dn4000). Following the type, you must specify the parameters required to define a feature of this type.

When the type is ```ratio```, the parameters must be, in order: ```low1, up1, low2, up2```, where ```low```/```up``` is the lower/upper wavelength (in Angstrom) of each continuum window. The flux ratio is then computed as the average flux between ```low2``` and ```up2```, divided by the average flux between ```low1``` and ```up1```. For example, the Dn4000 defined in [Balogh et al. (1999)](http://adsabs.harvard.edu/abs/1999ApJ...527...54B) would be defined in FAST++ as:

    dn4000 = ratio, 3850, 3950, 4000, 4100
    #               <-------->  <-------->
    #                  blue        red

When the type is ```abs```, the parameters must be, in order: ```line_low, line_up, cont1_low, cont1_up, ...```. The first two values specify the lower and upper wavelengths (in Angstrom) over which the line will be integrated. The next two values specify the lower and upper wavelengths of the first continuum window. The line equivalent width (in Angstrom) will be computed as the integral of the spectrum between ```line_low``` and ```line_up```, divided by the average flux between ```cont1_low``` and ```cont1_up```, minus ```line_up - line_low```. Following the decades-old convention, these values are always *positive* for absorption, and *negative* for emission.

Optionally, you may specify *two* continuum windows by adding one more pair of values at the end of the line. In this case, the continuum flux at the location of the line will be modeled as a straight line passing through the flux of the two continuum windows. For example, the Hdelta equivalent width defined in [Balogh et al. (1999)](http://adsabs.harvard.edu/abs/1999ApJ...527...54B) would be defined in FAST++ as:

    hdelta = abs, 4082, 4122, 4030, 4082, 4122, 4170
    #             <-------->  <-------->  <-------->
    #                line        blue        red

Pre-defined sets of continuum indices are provided in the ```share``` directory, including that of [Balogh et al. (1999)](http://adsabs.harvard.edu/abs/1999ApJ...527...54B) and the [Lick indices](http://astro.wsu.edu/worthey/html/system.html); feel free to create a new file and copy from these lists the indices you need.


# Additional documentation

## Adding new filters

For compatibility reasons, the filter database format is the same as that of EAzY and FAST. This is an ASCII/text file which contains all the filters, listed one after the other. The format must be:
```
num_pts [optional extra information...]
id lam trans
id lam trans
id lam trans
...
num_pts [optional extra information...]
id lam trans
id lam trans
id lam trans
...
```

In this example ```num_pts``` must be the number of data points in the filter response curve. Then for each data point, ```id``` is the identifier of that point. This ID is unused by FAST++, but it is recommended to make it start at one and increment by one for each data point (consistency checks may be implemented in the future). Then, ```lam``` is the wavelength in Angstrom, and ```trans``` is the filter transmission at that wavelength.

The overall normalization factor of the filter transmission does not matter, as the filters are automatically re-normalized to unit integral before the fit. If ```FILTERS_FORMAT=0```, the integral of ```lam*trans``` is normalized to unity, and if ```FILTERS_FORMAT=1``` then the integral of ```lam^2*trans``` is set to unity. The flux in a band is then computed as the integral of ```lam*trans*flx``` (for ```FILTERS_FORMAT=0```) or ```lam^2*trans*flx``` (for ```FILTERS_FORMAT=1```), where ```flx``` is the flux density (```Slambda```) in ```erg/s/cm^2/A```. This is exactly the same behavior as in FAST and EAzY.

To add new filters, simply append them to the end of the ```FILTER.RES``` file following the format above. For example, if your filter has 11 data points you would write:
```
...
11 my custom favorite filter (lambda0=5000 A)
1    4000.0   0.00
2    4200.0   0.10
3    4400.0   0.20
4    4600.0   0.40
5    4800.0   0.50
6    5000.0   0.55
7    5200.0   0.60
8    5400.0   0.40
9    5600.0   0.20
10   5800.0   0.10
11   6000.0   0.00
```

The wavelength step of the tabulated filter does not need to be constant, and the filter is assumed to have zero transmission below the minimum and above the maximum tabulated wavelength.


# Acknowledgments

Thanks to Francesco Valentino, Themiya Nanayakkara, Leindert Boobart for reporting issues with the code or documentation. Particular thanks to Hugh Dickinson for contributing changes to the code.

Particular thanks go to the authors of the original FAST implementation: Mariska Kriek and Ivo Labbé.
