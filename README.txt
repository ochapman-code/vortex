I have submitted an OpenMP code and an MPI code.

Each can be run by its bash file which compiles and runs it, and also makes a folder for the image outputs.
More detail is given in comments in the bash files.

The bin_to_png.py converts the output from the codes to nice graphs is configured to be run from
 the command line as
python bin_to_png.py <omp or mpi>   ->    e.g. python bin_to_png.py mpi

I've had a little issue getting bin_to_png to run from Linux, but when run from windows (using Anaconda prompt)
or from Spyder, it works well.

gantt.py produces a graph of the timings for the MPI code. It will popup a graph and it can then be magnified.
The files are quite big, but I thought I'd include this since I use it in my report.

I have also included the option to create a movie from the files. If the software mentioned in the bash files
is installed then movies of the pressure and flow speed will be created. It's worth a go - they're quite hypnotic!

I have set the codes up to do 5,000 iterations. Vortices should be observed by 50,000.