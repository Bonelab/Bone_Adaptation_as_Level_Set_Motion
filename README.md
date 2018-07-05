# Bone Adaptation as Level Set Motion
Try curvature based bone adaptation for yourself!
![Loss over time][C1-6]

<sub>*These images are created by incrementing 50 years of loss one year at a time. They are rendered directly from the embedding function using Marching Cubes. From left to right the surfaces are: Variable thickness rod, constant thickness rod, resorbed rod, resorbing rod, resorbed plate, constant plate.</sub>

If you find any of this code useful, please cite the following work:
```
@inproceedings{besler2018adaptation,
	title={Bone Adaptation as Level Set Motion},
	author={Besler, Bryce A and Gabel, Leigh and Burt, Lauren A and Forkert, Nils D and Boyd, Steven K},
 	booktitle={International Workshop and Challenge on Computational Methods and Clinical Applications in Musculoskeletal Imaging},
	year={2018},
	organization={Springer}
}
```

## Requirements
- ITK v4.12
- CMake
- Conda (Python only)

## File Structure
There exists four directories
- **Curvature** Handy tool for visualizing curvature
- **Initialize** Given a segmented dataset, use this to initialize phi
- **Iteration** Given a phi, use this to run one iteration
- **Python** Various python scripts that go along with the paper

`Initialize` and `Iteration` form the majority of the work.
The Python scripts are useful for preprocessing, visualization, and implicit functions.

[C1-6]: Figures/C1-6.gif
