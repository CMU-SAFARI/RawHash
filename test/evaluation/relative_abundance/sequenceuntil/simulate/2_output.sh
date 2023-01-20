#!/bin/bash

for i in `echo *.abundance`; do echo $i; cat $i; done;
