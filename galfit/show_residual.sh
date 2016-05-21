#!/bin/bash 

img1=$1
png_name="${1%.*}.png"

ds9 -nan black -view info no -view panner no -view magnifier no -view buttons no  -view colorbar no -geometry 800x800 -scale asinh -scale mode 95.0 -zoom to fit $img1[3] -cmap hsv -zoom to fit -export png $png_name -exit 

