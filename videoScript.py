import os
import sys

numImages = int(sys.argv[1]) # First CLA is total number of images
dt = float(sys.argv[2]) # Second CLA is time step
time = 0

if (not os.path.isdir("images")):
    os.system("mkdir images")
else:
    os.system("rm -f ./images/fr*.jpeg") # Remove old images if they exist

os.system("rm -f out.mp4") # Remove old video if it exists

print("--Generating Images--")
for i in range(1,numImages+1):

    # Need to pad image name with zeros for correct ordering
    outStr = str(i)
    if i < 10:
        outStr = "00" + outStr
    elif i < 100:
        outStr = "0" + outStr

    # Making plots
    cmd = "gnuplot -e \"set terminal jpeg; set view map; set title 't="+str(round(time,1)) + "'; splot 'out' matrix index " + str(int(i)-1) + " with pm3d \"> ./images/fr"+ outStr + "\.jpeg"
    os.system(cmd)
    time = time + dt
    
print("--Generating Movie out.mp4--")
# Making the movie 
cmd = "ffmpeg -hide_banner -loglevel error -framerate 20 -pattern_type glob -i \'images/*.jpeg\' -c:v libx264 -pix_fmt yuv420p out.mp4"
os.system(cmd)
