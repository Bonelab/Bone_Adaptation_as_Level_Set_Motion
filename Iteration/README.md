# Iteration
This algorithm performs one iteration of the algorithm as defined in Algorithm 1.
Importantly, the while loop may be ran multiple times in the algorithm.
Say you select a time of 1 year and the CFL condition is 0.5 years.
Two loops will occur.
This is unnatural if you are accustomed to simulated bone atrophy.

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
./Iteration input_phi output_phi output_seg curvature_weight propogation_weight total_time [input_mask]
```

To recreate the paper, run with the following parameters:

| Parameter          | Value  | Units     |
|:------------------ |:------ |:--------- |
| curvature_weight   | 0.0001 | mm^2/year |
| propogation_weight | -0.001 | mm/year   |
| total_time         | 10     | year      |

You would then run the algorithm three times.
Note that since image spacing is in mm, the parameters must be passed in mm.

## Known problems
Central differences were used for calculting the sign of the embedding function in `itkCurvatureBasedBoneResorptionFunction`.
On closer inspection of the original work of Ping et al, this should be a forward difference.
This error hasn't been corrected so the original work can be recreated.

Finally, this code has only been tested on OSX.
