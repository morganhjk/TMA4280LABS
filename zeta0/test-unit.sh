#!/bin/sh

prog="zeta0"

[ ! -f "./$prog" ] && exit 1

expected="2.85773803324704145"
computed="$(./$prog 3)"

message="$prog utest: expected=$expected, computed=$computed, test="

if [ "$computed" = "$expected" ]
then
	message="${message}OK"
else
	message="${message}FAIL"
fi

echo "$message"
