#!/bin/bash

thisdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

OBJDIR=$thisdir/objdir
LIBDIR=$thisdir/lib
INCDIR=$thisdir/include
TESTDIR=$thisdir/test
BINDIR=$thisdir/bin

#compile test example
echo "Compiling test example."
g++ -I$INCDIR -L${LIBDIR} -lFireDeamon ${TESTDIR}/main.cpp -o ${BINDIR}/testprog.exe

export LD_LIBRARY_PATH=$LIBDIR:$LD_LIBRARY_PATH
echo "Running test programme."
$BINDIR/testprog.exe

code=$?
echo "Return value was $code"
