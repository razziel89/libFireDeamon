#!/bin/bash

thisdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

OBJDIR=$thisdir/objdir
LIBDIR=$thisdir/lib
BINDIR=$thisdir/bin
BINDDIR=$thisdir/bindings

cat << EOF
cleaning:
$OBJDIR/*.o
$LIBDIR/*.a
$LIBDIR/*.so
$BINDIR/*.exe
$BINDDIR/*.so
$BINDDIR/*.cxx
$BINDDIR/*.py
$BINDDIR/*.pyc
EOF
rm $OBJDIR/*.o $LIBDIR/*.a $LIBDIR/*.so $BINDIR/*.exe $BINDDIR/*.so $BINDDIR/*.cxx $BINDDIR/*.py $BINDDIR/*.pyc 2>/dev/null
