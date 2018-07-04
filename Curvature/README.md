# Curvature
This is a small script for visualizing the curvature in an image.
This was not part of the paper, but I found it useful.

## How to Build
Please install the following software
- ITK v4.12
- CMake

Then run the following commands to biuld:
```bash
mkdir build
cd build
cmake ..
make -j
```

Run the script to get the required outputs
```bash
./Curv input output [dt]
```

Finally, this code has only been tested on OSX.
