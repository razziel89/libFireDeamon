#!/bin/bash
set -x
libdir="$1"
python="$2"
export PYTHONPATH=$(readlink -f $libdir):$PYTHONPATH
thisdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
if [ ! "$libdir" -o ! "$python" ]
then
    echo "First parameter needs to be the libdir and the second a python executable, aborting." >&2
    exit 27
fi
if [ ! -x $python ]
then
    echo "Given python executable $python is not executable, aborting." >&2
    exit 28
fi
version=$($python -c "import sys; version=sys.version_info; print str(version[0])+'.'+str(version[1])")
if [ ! $? -eq 0 ]
then
    echo "Given python executable $python could not determine it's version, aborting." >&2
    exit 29
fi
cd $libdir
$python -c "import FireDeamon as fd" &>/dev/null
if [ ! $? -eq 0 ]
then
    echo "Python in $python could not import the FireDeamon module" >&2
    if [ ! -s FireDeamon.py -o ! -x _FireDeamon.so ]
    then
        echo "and the files FireDeamon.py and _FireDeamon.so are NOT present" >&2
        echo "in $libdir or _FireDeamon.so is not executable. Please ivestigate." >&2
    else
        echo "although the files FireDeamon.py and _FireDeamon.so are present" >&2
        echo "in $libdir Please ivestigate." >&2
    fi
    exit 31
fi
$python $thisdir/test.py >&/dev/null
exit $?
set +x
