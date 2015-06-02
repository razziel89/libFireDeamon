#!/bin/bash

if [ ! "$PREFIX" ]
then
    echo "PREFIX undefined. You need to specify the environmental variable PREFIX to the prefix of your CGAL and SWIG installations. aborting." >&2
    exit 100
else
    if [ ! -d $PREFIX/include -o ! -d $PREFIX/lib ]
    then
        echo "There is no include or lib subdirectory in the given prefix $PREFIX. aborting." >&2
        exit 101
    fi
fi

thisdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

INCLUDEFLAGS="-I${PREFIX}/include -I${PREFIX}/include/python2.7"
LIBFLAGS="-L${PREFIX}/lib"
OBJDIR=$thisdir/objdir
LIBDIR=$thisdir/lib
INCDIR=$thisdir/include
SRCDIR=$thisdir/src
BINDDIR=$thisdir/bindings

#generate python bindings
echo "Generating python bindings."
swig -c++ -I$INCDIR -I$(ls -d -1 $PREFIX/share/swig/*/python | tail -n1) -python $BINDDIR/FireDeamon.i
cp $BINDDIR/FireDeamon_wrap.cxx $SRCDIR/FireDeamon_wrap.cpp

#compile source files for library
echo "Compiling source files."
for f in $SRCDIR/*.cpp
do
    g++ -fPIC $INCLUDEFLAGS -I${INCDIR} $LIBFLAGS -lCGAL -c $f -o $OBJDIR/$(basename $f .cpp).o
done

#link shared library
echo "Linking shared library."
g++ -shared -fPIC $INCLUDEFLAGS -I${INCDIR} $LIBFLAGS -lCGAL -o $LIBDIR/libFireDeamon.so $OBJDIR/*.o

#copy it with the appropriate name to the bindings directory
echo "Creating properly named copy of library for use with python"
cp $LIBDIR/libFireDeamon.so $BINDDIR/_FireDeamon.so

#link static library
echo "Linking static library."
ar rcs $LIBDIR/libFireDeamon.a $OBJDIR/*.o
