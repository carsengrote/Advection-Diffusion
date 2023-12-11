import os
import sys 
import numpy as np

def worker(startIndex,endIndex,dt, tId):
    
    time = startIndex * dt
    for i in range(startIndex, endIndex + 1):
        outStr = str(i)
        if i < 10: 
            outStr = "000" + outStr
        elif i < 100:
            outStr = "00" + outStr
        elif i < 1000:
            outStr = "0" + outStr

        # Making plots
        cmd = "gnuplot -e \"set terminal jpeg; set view map; set cbrange [0:0.005]; set title 't="+str(round(time,1)) + "'; splot 'out' matrix index " + str(int(i)) + " with pm3d \"> ./images/fr"+ outStr + "\.jpeg"
        
        os.system(cmd)
        time = time + dt

    print("Thread " + str(tId) + " finished")

numImages = int(sys.argv[1]) # First CLA is total number of images
dt = float(sys.argv[2]) # Second CLA is time step

if (not os.path.isdir("images")):
    os.system("mkdir images")
else:
    os.system("rm -f ./images/fr*.jpeg") # Remove old images if they exist

os.system("rm -f out.mp4") # Remove old video if it exists

numThreads = 10
imagesPerThread = int(np.floor(numImages/numThreads))
leftOver = np.mod(numImages,numThreads)

curIndex = 0
for t in range(0,numThreads): 
    
    j = os.fork()
    if j == 0: 
        if t == numThreads-1:
            worker(curIndex, numImages, dt, t)
        else:
            worker(curIndex, curIndex+imagesPerThread - 1, dt,t)
        os._exit(0)
    
    curIndex = curIndex + imagesPerThread

print("All " +  str(numThreads) + " threads launched")

for t in range(0, numThreads):
    os.wait()

print("--Generating Movie out.mp4--")
# Making the movie 
cmd = "ffmpeg -hide_banner -loglevel error -framerate 30 -pattern_type glob -i \'images/*.jpeg\' -c:v libx264 -pix_fmt yuv420p out.mp4"
os.system(cmd)
