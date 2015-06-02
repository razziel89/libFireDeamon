#!/bin/bash

if [ ! "$INSTALLPREFIX" ]
then
    echo "INSTALLPREFIX undefined. You need to specify the environmental variable INSTALLPREFIX to the prefix where you want this installed. aborting." >&2
    exit 100
fi

thisdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

LIBDIR=$thisdir/lib
INCDIR=$thisdir/include

echo "Installing "
[ ! -d $INSTALLPREFIX/lib ] && mkdir -p $INSTALLPREFIX/lib
[ ! -d $INSTALLPREFIX/include ] && mkdir -p $INSTALLPREFIX/include
cp $LIBDIR/* $INSTALLPREFIX/lib
cp $INCDIR/* $INSTALLPREFIX/include
