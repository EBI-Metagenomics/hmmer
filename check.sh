#!/bin/bash

make
src/hmmc2 < ross.txt
md5sum -c files.md5
