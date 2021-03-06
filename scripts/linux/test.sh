#!/bin/bash
set -e

scriptdir=$(cd $(dirname $0); pwd -P)
sourcedir=$(cd $scriptdir/../..; pwd -P)
. $scriptdir/shlibs/check-env.sh

builddir=$CRYPTOTEMPLATE_BUILD_DIR

cd "$builddir"
ctest --verbose
