import struct
import numpy as np
import matplotlib.pyplot as plt
import glob
import sys

imgs = glob.glob("mpi_img/*.bin")   #Folder containing .bin files
pngs = glob.glob("mpi_img/*.png")   #Folder containing .png files

print("Converting .bin files to .png files")
print("Total:", len(imgs))

imgs.sort()

for img in imgs:

	print(img+"\t\t\t", end="\r")
	
	with open(img, mode='rb') as file:      # Open file
		fileContent = file.read()
	
	# Read binary file
	#   (binary files used for memory efficiency and data preservation)
	NX, NY, N = struct.unpack("d" * 3, fileContent[:24])
	NX = int(NX)
	NY = int(NY)
	N = int(N)
	frames = struct.unpack("d" * (NX * NY)  , fileContent[24:])

	f = np.empty([NX, NY])
	
	for x in range(NX):
		for y in range(NY):
			f[x,y] = frames[NY *  x + y]    # Extract data into array

	img = img.replace(".bin", ".png")

	plt.figure(figsize=(8,3))
	
	if "rho" in img:
		plt.imshow(f.T, cmap="RdBu_r", origin='lower', vmin=0.975, vmax=1.025)
		plt.colorbar(shrink=0.7, label='\nRelative Density')
	elif 'vel' in img:
		plt.imshow(f.T, cmap="RdBu_r", origin='lower', vmin=0, vmax=0.2)
		plt.colorbar(shrink=0.7, label='\nFluid Speed (m/s)')
	
	plt.tight_layout()
	plt.savefig(img, dpi=150)   #Save file
	plt.close()

print("\nFile conversion complete\n")
