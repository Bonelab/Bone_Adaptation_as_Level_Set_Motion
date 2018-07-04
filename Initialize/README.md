# Initialize
The source code used to initialize the distance transform.
This performs the Maurer signed distance transform and runs the reinitialization code.
However, it was added in the initialization step so I could visualize the actual data going into the algorithm.
That is, I wanted to visualize the phi between reinitialization in Algorithm 1 and the update.
Reinitialization shouldn't change the signed distance transform (by definition) so this isn't a problem.

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

Run the script as such:
```bash
./Initialize input_seg output_phi
```

## Known problems
Central differences were used for calculting the sign of the embedding function in `itkCurvatureBasedBoneResorptionFunction`.
On closer inspection of the original work of Ping et al, this should be a forward difference.
This error hasn't been corrected so the original work can be recreated.

Finally, this code has only been tested on OSX.

