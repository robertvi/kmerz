#!/usr/bin/env bash

# run all tests, exit with non-zero code on the first failure
# run from within the test folder

set -eu

trap 'echo "Errorcode $? on line $LINENO" ; exit 1' ERR

export KMERZDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && cd .. && pwd )
export PATH=${KMERZDIR}/scripts:${PATH}
export PATH=${KMERZDIR}/test:${PATH}

#test python prototype pipeline
test_assemble_10k_random.sh
