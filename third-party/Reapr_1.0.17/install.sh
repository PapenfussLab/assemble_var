#!/usr/bin/env bash
set -e
rootdir=$PWD

unamestr=$(uname)

# Check if the OS seems to be Linux. This should be informative to
# Windows/Mac users that are trying to install this, despite it being
# developed for Linux. There is a VM avilable for Mac/Windows.
if [ "$unamestr" != "Linux" ]
then
    echo "Operating system does not appear to be Linux: uname says it's '$unamestr'"

    if [ "$1" = "force" ]
    then
        echo "'force' option used, carrying on anyway. Best of luck to you..."
    else
        echo "
If you *really* want to try installing anyway, run this:
./install.sh force

If you're using a Mac or Windows, then the recommended way to run REAPR is
to use the virtual machine."
        exit 1
    fi
fi


echo "------------------------------------------------------------------------------
                           Checking prerequisites
------------------------------------------------------------------------------
"
echo "Checking Perl modules..."

modules_ok=1

for module in File::Basename File::Copy File::Spec File::Spec::Link Getopt::Long List::Util
do
    set +e
    perl -M$module -e 1 2>/dev/null

    if [ $? -eq 0 ]
    then
        echo "  OK $module"
    else
        echo "  NOT FOUND: $module"
        modules_ok=0
    fi
    set -e
done

if [ $modules_ok -ne 1 ]
then
    echo "Some Perl modules were not found - please install them. Cannot continue"
    exit 1
else
    echo "... Perl modules all OK"
fi


echo
echo "Looking for R..."

if type -P R
then
    echo "... found R OK"
else
    echo "Didn't find R. It needs to be installed and in your path. Cannot continue"
fi


cd third_party

echo "
------------------------------------------------------------------------------
                           Compiling cmake
------------------------------------------------------------------------------
"
cd cmake
./bootstrap --prefix $PWD
make
cd ..
echo "
------------------------------------------------------------------------------
                           cmake compiled
------------------------------------------------------------------------------
"

echo "
------------------------------------------------------------------------------
                           Compiling Bamtools
------------------------------------------------------------------------------
"
cd bamtools
mkdir build
cd build
$rootdir/third_party/cmake/bin/cmake ..
make
cd $rootdir
echo "
------------------------------------------------------------------------------
                            Bamtools compiled
------------------------------------------------------------------------------
"

echo "
------------------------------------------------------------------------------
                              Compiling Tabix
------------------------------------------------------------------------------
"
cd third_party/tabix
make
cd ..
echo "
------------------------------------------------------------------------------
                              Tabix compiled
------------------------------------------------------------------------------
"

echo "
------------------------------------------------------------------------------
                             Compiling snpomatic
------------------------------------------------------------------------------
"
cd snpomatic
make
cd ..
echo "
------------------------------------------------------------------------------
                              snpomatic compiled
------------------------------------------------------------------------------
"

echo "
------------------------------------------------------------------------------
                             Compiling samtools
------------------------------------------------------------------------------
"
cd samtools
make
echo "
------------------------------------------------------------------------------
                              samtools compiled
------------------------------------------------------------------------------
"

echo "
------------------------------------------------------------------------------
                               Compiling Reapr
------------------------------------------------------------------------------
"
cd $rootdir/src
make
cd ..
ln -s src/reapr.pl reapr
echo "
Reapr compiled

All done!

Run
./reapr
for usage.

Read the manual
manual.pdf
for full instructions
"

