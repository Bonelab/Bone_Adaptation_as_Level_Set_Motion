# Bone Adaptation as Level Set Motion
Try curvature based bone adaptation for yourself!

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
- Python

## File Structure
There exists four directories
- **Curvature** Handy tool for visualizing curvature
- **Initialize** Given a segmented dataset, use this to initialize phi
- **Iteration** Given a phi, use this to run one iteration
- **Python** Various python scripts that go along with the paper

