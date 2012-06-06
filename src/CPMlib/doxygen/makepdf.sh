#!/bin/sh

rm -rf latex html

doxygen
cd latex
make
dvipdfmx refman.dvi
cp refman.pdf ..


