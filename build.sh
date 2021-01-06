#!/bin/bash

pushd build

cmake ..

rA=$?
if [ $rA != 0 ]; then
	popd
	exit $rA
fi

make

rB=$?
if [ $rB != 0 ]; then
	popd
	exit $rB
fi

popd
exit 0
