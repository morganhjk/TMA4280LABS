#!/bin/sh

prog="mach0"

[ ! -f "./$prog" ] && exit 1

expected="3.14162102932503462"
computed="$(./$prog 3)"

message="$prog utest: expected=$expected, computed=$computed, test="

if [ "$computed" = "$expected" ]
then
	message="${message}OK"
else
	message="${message}FAIL"
fi

echo "$message"
