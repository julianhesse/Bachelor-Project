#!/bin/bash

awk '!/NN+/ && !a[$0]++' RS='>' ORS='>' $1 | head -n -1 > $2
