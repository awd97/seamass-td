seaMass-TD
=======

seaMass-TD is the first method to deconvolute top down proteomics spectra that infers high resolution output on isotopically unresolved input.

See [here](http://www.biospi.org/research/ms/seamass-td/) for more information.

dependencies
-------
Depends on _CMake_, _HDF5_ and _Boost_.
Windows developers can get these from the [dependency repository](https://github.com/biospi/seamass-windeps).
Also currently requires _Intel MKL_ with the _Intel C++ Compiler 14.0 or greater_.

building
-------
Run the cmake.sh or cmake.bat file for your platform, which will create make files for compilation in
build/debug and build/release subfolders.
