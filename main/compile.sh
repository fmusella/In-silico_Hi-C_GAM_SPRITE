#!/bin/bash

in_name=$1

out_name=${in_name%".c"}
out_name="${out_name}.out"

gcc $in_name ../settings.c ../src/utilities.c ../src/io.c ../src/clustering.c ../src/insilico_gam_utilities.c ../src/insilico_hic.c ../src/insilico_sprite.c ../src/insilico_gam.c -o $out_name -lm -std=gnu99

chmod +x $out_name
