#!/bin/bash

# Test script for lssolve
# Usage: ./test_lssolve.sh

# Check if lssolve executable exists

if [ ! -f lssolve ]; then
    echo "Error: lssolve executable not found. Run 'make' first."
    exit 1
fi

# Run tests

echo "Testing with full-rank matrix..."
./lssolve data/A1.txt data/b1.txt data/x1.txt
if [ $? -ne 0 ]; then
    echo "  ***Test failed: lssolve returned unexpected exit status"
    exit 1
else
    echo "  Passed."
fi

echo "Testing with rank-deficient matrix..."
./lssolve data/A2.txt data/b2.txt data/x2.txt
if [ $? -ne 1 ]; then
    echo "  ***Test failed: lssolve returned unexpected exit status"
    exit 1
else
    echo "  Passed."
fi

echo "Testing with m < n ..."
./lssolve data/A3.txt data/b3.txt data/x3.txt
if [ $? -ne 1 ]; then
    echo "  ***Test failed: lssolve returned unexpected exit status"
    exit 1
else
    echo "  Passed."
fi

echo "Testing with non-existent input files..."
./lssolve data/missing.txt data/missing.txt data/x4.txt
if [ $? -ne 1 ]; then
    echo "  ***Test failed: lssolve returned unexpected exit status"
    exit 1
else
    echo "  Passed."
fi

exit 0