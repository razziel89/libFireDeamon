#!/bin/bash

if [ ! "$PREFIX" ]
then
    echo "PREFIX undefined. You need to specify the environmental variable PREFIX to the prefix of your CGAL installation. aborting." >&2
    exit 100
else
    if [ ! -d $PREFIX/include -o ! -d $PREFIX/lib ]
    then
        echo "There is no include or lib subdirectory in the given prefix $PREFIX. aborting." >&2
        exit 101
    fi
fi

thisdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

INCLUDEFLAGS="-I${PREFIX}/include"
LIBFLAGS="-L${PREFIX}/lib"
OBJDIR=$thisdir/objdir
LIBDIR=$thisdir/lib
INCDIR=$thisdir/include
SRCDIR=$thisdir/src

#compile source files for library
echo "Compiling source files."
for f in $SRCDIR/*.cpp
do
    g++ -fPIC $INCLUDEFLAGS -I${INCDIR} $LIBFLAGS -lCGAL -c $f -o $OBJDIR/$(basename $f .cpp).o
done

#link shared library
echo "Linking shared library."
g++ -shared -fPIC $INCLUDEFLAGS -I${INCDIR} $LIBFLAGS -lCGAL -o $LIBDIR/libFireDeamon.so $OBJDIR/*.o

#link static library
echo "Linking static library."
ar rcs $LIBDIR/libFireDeamon.a $OBJDIR/*.o
