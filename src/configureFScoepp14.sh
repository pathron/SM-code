#!/bin/sh

# Wrapper script for calling FlexibleSUSY's configure script on coepp14.
# So that you don't have to remember all of the library and
# include paths etc. The script takes all of the same 
# command line arguments as configure does, and automatically
# provides the paths (though these can be overwritten if you
# choose to do so). Place this script in your FlexibleSUSY directory
# (i.e. the same directory as configure is located in).

# Directory of this script
BASEDIR=$(dirname $0)

# Absolute path to this script
ABSBASEDIR=$(cd $BASEDIR; pwd)

# Optionally provide the path to configure, otherwise
# assume that the script is located in the same
# directory as configure.
configure_path=${ABSBASEDIR}

# Relevant library and include paths
boost_lib_path="/home/software/boostbuild/lib"
boost_inc_path="/home/software/boostbuild/include"
eigen_inc_path="/home/software/sharedincludes/eigen3"
lapack_lib_path="/home/software/sharedlibs"

trap "exit 1" INT QUIT TERM

# Parse command line arguments.
arg_list=""

if test $# -gt 0 ; then
    while test ! "x$1" = "x" ; do
        case "$1" in
            -*=*) optarg=`echo "$1" | sed 's/[-_a-zA-Z0-9]*=//'` ;;
            *) optarg= ;;
        esac

	# For all FlexibleSUSY arguments, just append to argument list.
	# Note this isn't a particularly efficient way of doing this,
	# but it works.
        case $1 in
            --with-boost-libdir=*)     boost_lib_path=${optarg}    ;;
            --with-boost-incdir=*)     boost_inc_path=${optarg}    ;;
            --with-eigen-incdir=*)     eigen_inc_path=${optarg}    ;;
            --with-lapack-libdir=*)    lapack_lib_path=${optarg}   ;;
            --with-configure-dir=*)    configure_path=${optarg}    ;;  
            *)  arg_list="${arg_list} $1" ;;
        esac
        shift
    done
fi

if test ! "x${boost_lib_path}" = "x"  ; then
    arg_list="${arg_list} --with-boost-libdir=${boost_lib_path}"
fi

if test ! "x${boost_inc_path}" = "x"  ; then
    arg_list="${arg_list} --with-boost-incdir=${boost_inc_path}"
fi

if test ! "x${eigen_inc_path}" = "x"  ; then
    arg_list="${arg_list} --with-eigen-incdir=${eigen_inc_path}"
fi

if test ! "x${lapack_lib_path}" = "x"  ; then
    arg_list="${arg_list} --with-lapack-libdir=${lapack_lib_path}"
fi


# Check that we can find the configure script
if [ ! -f "${configure_path}/configure" ] ; then
    echo "ERROR: configure not found in '${configure_path}'."
    echo "       Check that you have the correct location and try again."
    exit 1;
fi

# Call FlexibleSUSY's configure script with the
# arguments provided
${configure_path}/configure ${arg_list}
