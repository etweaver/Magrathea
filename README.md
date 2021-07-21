Magrathea is a C++ library for fitting and modeling millimeter 
and submillimeter emission from protoplanetary disks. This 
preliminary version contains only the modeling sections.

PREREQUISITS
------------

C++ compiler with support for C++17 features

CFITSIO:heasarc.gsfc.nasa.gov/fitsio/

FFTW(>3.3):www.fftw.org

Threadpool:github.com/progschj/ThreadPool/

VCL:github.com/vectorclass


INSTALLATION
------------

Installation is done with the usual sequence:
./configure
make 
make install

The install location can be set using the --prefix command, and
libraries not in the search path can be included using 
--with-[library]. For a full list of option, use --help or -h, 
or consult Section 1.2 of the documentation.
