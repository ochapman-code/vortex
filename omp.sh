#!/usr/bin/env bash

echo Creating directory omp_img
mkdir omp_img &> /dev/null

gcc lbm_omp.c -fopenmp -lm -O3 -march=native -std=c99 -o omp.out
./omp.out 4  #enter thread num

python bin_to_png.py omp    # This might only work on windows

# Optional: install ffmpeg with 'sudo apt-get install ffmpeg' to create movie

echo Creating omp_img/omp_rho.h264
ffmpeg -i omp_img/omp_rho_%05d.png omp_img/omp_rho.h264  &> /dev/null
echo Creating omp_img/omp_vel.h264
ffmpeg -i omp_img/omp_vel_%05d.png omp_img/omp_vel.h264  &> /dev/null

# Output movie can be played via VLC media player
# Optional: install VLC with 'sudo apt-get install vlc' to create movie



