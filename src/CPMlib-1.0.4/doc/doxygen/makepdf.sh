#!/bin/sh

rm -rf latex html

doxygen
cd latex
make
dvipdfmx refman.dvi
#cp refman.pdf ..
cp refman.pdf ../../reference.pdf
cd ..

rm -rf latex html

