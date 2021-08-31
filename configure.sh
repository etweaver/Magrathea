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
#find_package libcfitsio

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
		VCL_INCDIR="${TMP}/vcl";
	continue; fi
		
	echo "config.sh: Unknown or malformed option '$arg'" 1>&2
	exit 1
done

if [ "$CFITSIO_INCDIR" -a "$CFITSIO_LIBDIR" ]; then
	CFITSIO_CFLAGS="-I${CFITSIO_INCDIR}"
	CFITSIO_LDFLAGS="-L${CFITSIO_LIBDIR}"
fi;

if [ "$FFTW_INCDIR" -a "$FFTW_LIBDIR" ]; then
	CFITSIO_CFLAGS="-I${FFTW_INCDIR}"
	CFITSIO_LDFLAGS="-L${FFTW_LIBDIR}"
fi;

if [ "$VCL_INCDIR" ]; then
	VCL_CFLAGS="-I${VCL_INCDIR}"
fi;

if [ ! -d ./build/ ]; then
    mkdir build;
fi

if [ ! -d ./lib/ ]; then
    mkdir lib;
fi

echo "Generating makefile..."
echo "
CC=$CC
CXX=$CXX
AR=$AR
LD=$LD
DYN_SUFFIX=$DYN_SUFFIX
DYN_OPT=$DYN_OPT
PREFIX=$PREFIX
CXXFLAGS+= -fPIC -O3 -mavx -std=c++17 -U__STRICT_ANSI__
LDFLAGS+= ${CFITSIO_CFLAGS} -lcfitsio ${FFTW_CFLAGS} -lfftw3
INCFLAGS+= ${CFITSIO_CFLAGS} ${FFTW_CFLAGS} ${VCL_CFLAGS}
EXAMPLES := examples/PowerLawDisk \
	examples/GapDisk \
	examples/PowerLawDiskCO
" > ./Makefile

echo '
.PHONY: all clean
all : lib/libmagrathea$(DYN_SUFFIX)

clean : 
	rm -rf build/disk.o
	rm -rf build/geometry.o
	rm -rf build/image.o
	rm -rf build/diskPhysics.o
	rm -rf build/magrathea.o
	rm -rf lib/libmagrathea$(DYN_SUFFIX)
	rm -rf examples/PowerLawDisk
	rm -rf examples/GapDisk
	rm -rf examples/PowerLawDiskCO
	
lib/libmagrathea$(DYN_SUFFIX) : build/magrathea.o build/geometry.o build/diskPhysics.o
	$(CXX) $(LDFLAGS) $(INCFLAGS) -fPIC $(DYN_OPT) -o lib/libmagrathea$(DYN_SUFFIX) $^
	
build/magrathea.o : src/magrathea.cpp include/magrathea/magrathea.h include/magrathea/grid.h include/magrathea/diskPhysics.h
	$(CXX) $(CXXFLAGS) $(INCFLAGS) src/magrathea.cpp -c -o build/magrathea.o
build/geometry.o : src/geometry.cpp include/magrathea/diskPhysics.h include/magrathea/geometry.h
	$(CXX) $(CXXFLAGS) $(INCFLAGS) src/geometry.cpp -c -o build/geometry.o
build/image.o : src/image.cpp include/magrathea/image.h include/magrathea/diskPhysics.h
	$(CXX) $(CXXFLAGS) $(INCFLAGS) src/image.cpp -c -o build/image.o
build/diskPhysics.o : src/diskPhysics.cpp include/magrathea/diskPhysics.h
	$(CXX) $(CXXFLAGS) $(INCFLAGS) src/diskPhysics.cpp -c -o build/diskPhysics.o
	
install :
	cp lib/libmagrathea$(DYN_SUFFIX) $(PREFIX)lib/
	cp -r include/magrathea $(PREFIX)include/
	
examples : $(EXAMPLES)

examples/PowerLawDisk : examples/PowerLawDisk.cpp
	$(CXX) $(CXXFLAGS) -Iinclude -Llib $(INCFLAGS) -lmagrathea $(LDFLAGS) examples/PowerLawDisk.cpp  -o examples/PowerLawDisk
examples/GapDisk : examples/GapDisk.cpp
	$(CXX) $(CXXFLAGS) -Iinclude -Llib $(INCFLAGS) -lmagrathea $(LDFLAGS) examples/GapDisk.cpp -o examples/GapDisk
examples/PowerLawDiskCO : examples/PowerLawDiskCO.cpp
	$(CXX) $(CXXFLAGS) -Iinclude -Llib $(INCFLAGS) -lmagrathea $(LDFLAGS) examples/PowerLawDiskCO.cpp -o examples/PowerLawDiskCO
	
uninstall : 
	rm -rf $(PREFIX)lib/libmagrathea$(DYN_SUFFIX)
	rm -rf $(PREFIX)include/magrathea

' >> ./Makefile
