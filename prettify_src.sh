#!/bin/bash

for i in src/*.F90
do
  fprettify --indent 4 -l 1000 --enable-decl $i
done

