PREFIX=/usr/local/
OS_NAME=`uname -s`

check_pkgconfig(){
	if [ "$CHECKED_PKGCONFIG" ]; then return; fi
	echo "Looking for pkg-config..."
	which pkg-config 2>&1 > /dev/null
	if [ "$?" -ne 0 ]; then
		echo "Error: pkg-config not found; you will need to specify library locations manually" 1>&2
		exit 1
	fi
	CHECKED_PKGCONFIG=1
}

find_package(){
	PKG=$1
	VAR_PREFIX=`echo $PKG | tr [:lower:] [:upper:]`
	TMP_FOUND=`eval echo "$"${VAR_PREFIX}_FOUND`
	if [ "$TMP_FOUND" ]; then return; fi
	check_pkgconfig
	echo "Looking for $PKG..."

	pkg-config --exists $PKG
	if [ "$?" -ne 0 ]; then
		echo " $PKG not found with pkg-config"
		return
	fi
	if [ $# -ge 2 ]; then
		MIN_VERSION=$2
		pkg-config --atleast-version $MIN_VERSION $PKG
		if [ "$?" -ne 0 ]; then
			echo "Error: installed $PKG version ("`pkg-config --modversion $PKG`") is too old; version >=$MIN_VERSION is required" 1>&2
			exit 1
		fi
	fi
	echo " Found $PKG version `pkg-config --modversion $PKG`"
	eval ${VAR_PREFIX}_FOUND=1
	eval ${VAR_PREFIX}_VERSION=\"`pkg-config --modversion $PKG`\"
	eval ${VAR_PREFIX}_CFLAGS=\"`pkg-config --cflags $PKG`\"
	eval ${VAR_PREFIX}_LDFLAGS=\"`pkg-config --libs $PKG`\"
	eval ${VAR_PREFIX}_INCDIR=\"`pkg-config --variable=includedir $PKG`\"
	eval ${VAR_PREFIX}_LIBDIR=\"`pkg-config --variable=libdir $PKG`\"
}

ensure_found(){
	PKG=$1
	VAR_PREFIX=`echo $PKG | tr [:lower:] [:upper:]`
	TMP_FOUND=`eval echo "$"${VAR_PREFIX}_FOUND`
	if [ "$TMP_FOUND" ]; then return; fi
	#not found
	lowername=`echo $PKG | tr [A-Z] [a-z]`

	TMP_INCDIR=`eval echo "$"${VAR_PREFIX}_INCDIR`
	TMP_LIBDIR=`eval echo "$"${VAR_PREFIX}_LIBDIR`
	if [ "$TMP_INCDIR" -a "$TMP_LIBDIR" ]; then
		echo "Error: $PKG not found in $TMP_INCDIR and $TMP_LIBDIR or with pkg-config" 1>&2
		echo "Please verify that the path given to --with-${lowername} is correct" 1>&2
	else
		echo "Error: $PKG not installed or not registered with pkg-config" 1>&2
		echo "Please specify location using the --with-${lowername} flag" 1>&2
	fi
	unset TMP_INCDIR
	unset TMP_LIBDIR
	exit 1
}


find_cfitsio(){
	PKG=cfitsio
	VAR_PREFIX=`echo $PKG | tr [:lower:] [:upper:]`
	TMP_FOUND=`eval echo "$"${VAR_PREFIX}_FOUND`
	if [ "$TMP_FOUND" ]; then return; fi
	echo "Looking for $PKG..."
	POSSIBLE_PREFIXES="/usr /usr/local"
	POSSIBLE_LIBDIRS="/lib /lib64 /lib/x86_64-linux-gnu"
	for PREFIX in $POSSIBLE_PREFIXES; do
		for LIBDIR in $POSSIBLE_LIBDIRS; do
			CFITSIO_LIBDIR="${PREFIX}${LIBDIR}"
			if [ -d "$CFITSIO_LIBDIR" \
				-a \( -e "$CFITSIO_LIBDIR/libcfitsio.a" -o -e "$CFITSIO_LIBDIR/libcfitsio.so" \) ]; then
				CFITSIO_LIB_FOUND=1
				break
			fi;
		done
		if [ "$CFITSIO_LIB_FOUND" ]; then break; fi
	done
	for PREFIX in $POSSIBLE_PREFIXES; do
		for SUFFIX in "" "/cfitsio"; do 
			CFITSIO_INCDIR="${PREFIX}/include${SUFFIX}"
			if [ -d "$CFITSIO_INCDIR" -a -e "$CFITSIO_INCDIR/fitsio.h" ]; then
				CFITSIO_INC_FOUND=1
				break
			fi
		done
		if [ "$CFITSIO_INC_FOUND" ]; then break; fi
	done
	if [ "$CFITSIO_LIB_FOUND" -a "$CFITSIO_INC_FOUND" ]; then
		CFITSIO_FOUND=1
		echo " Found cfitsio headers at $CFITSIO_INCDIR"
		echo " Found cfitsio libraries at $CFITSIO_LIBDIR"
		CFITSIO_CFLAGS="-I${CFITSIO_INCDIR}"
		CFITSIO_LDFLAGS="-L${CFITSIO_LIBDIR} -lcfitsio"
	fi
}

#find vectorclass.h, vectormath_trig.h, and vectormath_exp.h"
find_vcl(){
	PKG=vcl
	VAR_PREFIX=`echo $PKG | tr [:lower:] [:upper:]`
	TMP_FOUND=`eval echo "$"${VAR_PREFIX}_FOUND`
	if [ "$TMP_FOUND" ]; then return; fi
	echo "Looking for $PKG..."
	POSSIBLE_VCL_INCDIRS="/usr/include /usr/local/include /usr/include/vcl /usr/local/include/vcl"
	for VCL_INCDIR in $POSSIBLE_VCL_INCDIRS; do
		if [ -d "$VCL_INCDIR" \
			-a -e "$VCL_INCDIR/vectorclass.h" ]; then
			VCL_FOUND=1
			#we don't have a libdir for VCL, but ensure_found() expects one
			VCL_LIBDIR="NotActuallyUsed"
			VCL_CFLAGS="-I${VCL_INCDIR}"
			echo " Found VCL headers at $VCL_INCDIR"
			return
		fi;
	done
}

GUESS_CC=gcc
GUESS_CXX=g++
GUESS_AR=ar
GUESS_LD=ld
if [ "$OS_NAME" = Linux ]; then
	DYN_SUFFIX=.so
	DYN_OPT='-shared'
fi
if [ "$OS_NAME" = Darwin ]; then
	GUESS_CC=clang
	GUESS_CXX=clang++
	GUESS_LD=clang++
	DYN_SUFFIX=.dylib
	DYN_OPT='-dynamiclib'
fi

CC=${CC-$GUESS_CC}
CXX=${CXX-$GUESS_CXX}
AR=${AR-$GUESS_AR}
LD=${LD-$GUESS_LD}

find_package fftw3

HELP="Usage: ./config.sh [OPTION]... 
Installation directories:
  --prefix=PREFIX         install files in PREFIX
                          [$PREFIX]
By default, \`make install' will install all the files in
\`$PREFIXbin', \`$PREFIXlib' etc.  You can specify
an installation prefix other than \`$PREFIX' using \`--prefix'.
	--with-cfitsio=DIR          use the copy of cfitsio in DIR
	                            assuming headers are in DIR/include
	                            and libraries in DIR/lib
	--with-fftw=DIR             use the copy of fftw in DIR
	--with-vcl=DIR              use the copy of vcl in DIR
"

for arg in "$@"
do
	if [ "$arg" = "--help" -o "$arg" = "-h" ]; then
		echo "$HELP"
		exit 0
	fi
	
	TMP=`echo "$arg" | sed -n 's/^--prefix=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then 
		PREFIX="$TMP"; 
	continue; fi
		
	TMP=`echo "$arg" | sed -n 's/^--with-cfitsio=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then
		CFITSIO_INCDIR="${TMP}/include";
		CFITSIO_LIBDIR="${TMP}/lib";
	continue; fi
	
	TMP=`echo "$arg" | sed -n 's/^--with-fftw=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then
		FFTW_INCDIR="${TMP}/include";
		FFTW_LIBDIR="${TMP}/lib";
	continue; fi
	
	TMP=`echo "$arg" | sed -n 's/^--with-vcl=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then
		VCL_INCDIR="${TMP}";
	continue; fi
		
	echo "config.sh: Unknown or malformed option '$arg'" 1>&2
	exit 1
done

if [ "$CFITSIO_INCDIR" -a "$CFITSIO_LIBDIR" ]; then
	echo "Checking manually specified cfitsio directories"
	CFITSIO_CFLAGS="-I${CFITSIO_INCDIR}"
	CFITSIO_LDFLAGS="-L${CFITSIO_LIBDIR}"
	if [ -d "$CFITSIO_LIBDIR" \
		-a \( -e "$CFITSIO_LIBDIR/libcfitsio.a" -o -e "$CFITSIO_LIBDIR/libcfitsio.so" \) \
		-a -d "$CFITSIO_INCDIR" -a -e "$CFITSIO_INCDIR/fitsio.h" ]; then
		CFITSIO_FOUND=1
		CFITSIO_CFLAGS="-I${CFITSIO_INCDIR}"
		CFITSIO_LDFLAGS="-L${CFITSIO_LIBDIR} -lcfitsio"
		echo " Found cfitsio"
	fi
fi
find_cfitsio

if [ "$FFTW_INCDIR" -a "$FFTW_LIBDIR" ]; then
	echo "Checking manually specified fftw3 directories"
	FFTW3_CFLAGS="-I${FFTW_INCDIR}"
	FFTW3_LDFLAGS="-L${FFTW_LIBDIR} -lfftw3"
	if [ -d "$FFTW_LIBDIR" \
		-a \( -e "$FFTW_LIBDIR/libfftw3.a" -o -e "$FFTW_LIBDIR/libfftw3.so" \) \
		-a -d "$FFTW_INCDIR" -a -e "$FFTW_INCDIR/fftw3.h" ]; then
		FFTW3_FOUND=1
		FFTW3_CFLAGS="-I${FFTW_INCDIR}"
		FFTW3_LDFLAGS="-L${FFTW_LIBDIR} -lfftw3"
		echo " Found FFTW"
	fi
fi;

if [ "$VCL_INCDIR" ]; then
	echo "Checking manually specified VCL directory"
	if [ -d "$VCL_INCDIR" \
		-a -e "$VCL_INCDIR/vectorclass.h" ]; then
		VCL_FOUND=1
		VCL_CFLAGS="-I${VCL_INCDIR}"
		echo " Found VCL"
	fi;
fi;
find_vcl

if [ ! -d ./build/ ]; then
    mkdir build;
fi

if [ ! -d ./lib/ ]; then
    mkdir lib;
fi

ensure_found fftw3
ensure_found vcl
ensure_found cfitsio

echo "Generating makefile..."
echo "
CC=$CC
CXX=$CXX
AR=$AR
LD=$LD
DYN_SUFFIX=$DYN_SUFFIX
DYN_OPT=$DYN_OPT
PREFIX=$PREFIX
CXXFLAGS+= -fPIC -O3 -mavx -mfma -std=c++17 -pthread
LDFLAGS+= ${CFITSIO_LDFLAGS} ${FFTW3_LDFLAGS}
INCFLAGS+= ${CFITSIO_CFLAGS} ${FFTW3_CFLAGS} ${VCL_CFLAGS}
EXAMPLES := examples/PowerLawDisk \
	examples/GapDisk \
	examples/PowerLawDiskCO\
	examples/Fourier \
	examples/Fit
" > ./Makefile

echo '
.PHONY: all clean
all : lib/libmagrathea$(DYN_SUFFIX)

clean : 
	rm -rf build/*.o
	rm -rf lib/libmagrathea$(DYN_SUFFIX)
	rm -rf examples/PowerLawDisk
	rm -rf examples/GapDisk
	rm -rf examples/PowerLawDiskCO
	rm -rf examples/Fourier
	
lib/libmagrathea$(DYN_SUFFIX) : build/magrathea.o build/geometry.o build/diskPhysics.o build/ParameterSet.o
	$(CXX) $(INCFLAGS) -fPIC $(DYN_OPT) -o lib/libmagrathea$(DYN_SUFFIX) $^ $(LDFLAGS)
	
build/magrathea.o : src/magrathea.cpp include/magrathea/magrathea.h include/magrathea/grid.h include/magrathea/diskPhysics.h include/magrathea/diskStructures.h include/magrathea/image.h include/magrathea/mcmc.h include/magrathea/ParameterSet.h
	$(CXX) $(CXXFLAGS) $(INCFLAGS) src/magrathea.cpp -c -o build/magrathea.o
build/geometry.o : src/geometry.cpp include/magrathea/diskPhysics.h include/magrathea/geometry.h
	$(CXX) $(CXXFLAGS) $(INCFLAGS) src/geometry.cpp -c -o build/geometry.o
build/image.o : src/image.cpp include/magrathea/image.h include/magrathea/diskPhysics.h
	$(CXX) $(CXXFLAGS) $(INCFLAGS) src/image.cpp -c -o build/image.o
build/diskPhysics.o : src/diskPhysics.cpp include/magrathea/diskPhysics.h
	$(CXX) $(CXXFLAGS) $(INCFLAGS) src/diskPhysics.cpp -c -o build/diskPhysics.o
build/ParameterSet.o : src/ParameterSet.cpp include/magrathea/ParameterSet.h 
	$(CXX) $(CXXFLAGS) $(INCFLAGS) src/ParameterSet.cpp -c -o build/ParameterSet.o
	
install :
	cp lib/libmagrathea$(DYN_SUFFIX) $(PREFIX)lib/
	cp -r include/magrathea $(PREFIX)include/
	
examples : $(EXAMPLES)

examples/PowerLawDisk : examples/PowerLawDisk.cpp
	$(CXX) $(CXXFLAGS) -Iinclude -Llib $(INCFLAGS) examples/PowerLawDisk.cpp  -o examples/PowerLawDisk -lmagrathea $(LDFLAGS)
examples/GapDisk : examples/GapDisk.cpp
	$(CXX) $(CXXFLAGS) -Iinclude -Llib $(INCFLAGS) examples/GapDisk.cpp -o examples/GapDisk -lmagrathea $(LDFLAGS)
examples/PowerLawDiskCO : examples/PowerLawDiskCO.cpp
	$(CXX) $(CXXFLAGS) -Iinclude -Llib $(INCFLAGS) examples/PowerLawDiskCO.cpp -o examples/PowerLawDiskCO -lmagrathea $(LDFLAGS)
examples/Fourier : examples/Fourier.cpp
	$(CXX) $(CXXFLAGS) -Iinclude -Llib $(INCFLAGS) examples/Fourier.cpp -o examples/Fourier -lmagrathea $(LDFLAGS)
examples/Fit : examples/Fit.cpp
	$(CXX) $(CXXFLAGS) -Iinclude -Llib $(INCFLAGS) examples/Fit.cpp -o examples/Fit -lmagrathea $(LDFLAGS)
	
uninstall : 
	rm -rf $(PREFIX)lib/libmagrathea$(DYN_SUFFIX)
	rm -rf $(PREFIX)include/magrathea

' >> ./Makefile
