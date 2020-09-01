#!/usr/bin/env bash

echo Creating directory mpi_img
mkdir mpi_img &> /dev/null

mpicc lbm_mpi.c -lm -O3 -march=native -std=c99 -o mpi.out
mpiexec ./mpi.out -n 2  #enter rank num

python gantt.py             # This might only work on windows
python bin_to_png.py mpi    # This might only work on windows

# Optional: install ffmpeg with 'sudo apt get install ffmpeg' to create movie

echo Creating mpi_img/mpi_rho.h264
ffmpeg -i mpi_img/mpi_rho_%05d.png mpi_img/mpi_rho.h264  &> /dev/null
echo Creating mpi_img/mpi_vel.h264
ffmpeg -i mpi_img/omp_vel_%05d.png mpi_img/omp_vel.h264  &> /dev/null

# Output movie can be played via VLC media player



