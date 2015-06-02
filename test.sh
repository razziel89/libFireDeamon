#!/bin/bash

thisdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#
#INCLUDEFLAGS="-I${PREFIX}/include -I${PREFIX}/include/python2.7"
#
#OBJDIR=$thisdir/objdir
#LIBDIR=$thisdir/lib
#INCDIR=$thisdir/include
#TESTDIR=$thisdir/test
#BINDIR=$thisdir/bin
#
##compile test example
#echo "Compiling test example."
#g++ $INCLUDEFLAGS -I$INCDIR -L${LIBDIR} -lFireDeamon ${TESTDIR}/main.cpp -o ${BINDIR}/testprog.exe
#
#export LD_LIBRARY_PATH=$LIBDIR:$LD_LIBRARY_PATH
#echo "Running test programme."
#$BINDIR/testprog.exe

export PYTHONPATH=$thisdir/bindings:$PYTHONPATH
export LD_LIBRARY_PATH=$thisdir/bindings:$LD_LIBRARY_PATH
echo "Running test script"
python $thisdir/test/test.py

code=$?
echo "Return value was $code"
