#!/bin/bash

while read p; do ./check_prokka_output_quality.sh $p; done<$1
