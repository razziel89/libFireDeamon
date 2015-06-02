#!/bin/bash

thisdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

OBJDIR=$thisdir/objdir
LIBDIR=$thisdir/lib
BINDIR=$thisdir/bin

echo "cleaning $OBJDIR, $LIBDIR and $BINDIR"
rm $OBJDIR/*.o $LIBDIR/*.a $LIBDIR/*.so $BINDIR/*.exe
