import struct
import numpy as np
import matplotlib.pyplot as plt
import glob
import sys

# Note this might only work on Windows. My Linux disro. and BC3 seem to have issue with matplotlib

folder = str(sys.argv[1])

print(folder)

imgs = glob.glob(folder+"_img/*.bin")   #Folder containing .bin files
pngs = glob.glob(folder+"_img/*.png")   #Folder containing .png files

print("\nConverting .bin files to .png files")
print("Total:", len(imgs))

imgs.sort()

for img in imgs:
    
    if img.replace(".bin", ".png")  in pngs:    # Skip if already converted
        continue
    else:
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
    
        plt.figure(figsize=(12,4))
        
        if "rho" in img:
            plt.imshow(f.T, cmap="coolwarm", origin='lower')
        elif "vel" in img:
            plt.imshow(f.T, cmap="coolwarm", origin='lower')
        else:
            plt.imshow(f.T, origin='lower')
  
        plt.colorbar(shrink=0.8)
        plt.tight_layout()
        plt.savefig(img, dpi=250)   #Save file
        plt.close()

print("\nFile conversion complete\n")
