SPACING="0.061"

# Create cone in
python rod.py --spacing ${SPACING} ${SPACING} ${SPACING} --alpha 1.75 --length 1.0 --radius_1 0.20 --radius_2 0.15 C1.nii > CREATE_C1.LOG

# Create cylinder
python cylinder.py --spacing ${SPACING} ${SPACING} ${SPACING} --alpha 1.75  --length 1.0 --radius 0.20 C2.nii > CREATE_C2.LOG

# Create resorbed rod
python hyperbola.py --spacing ${SPACING} ${SPACING} ${SPACING} --alpha 1.75 --length 1.0 --radius 0.20 --distance 0.20 C3.nii > CREATE_C3.LOG

# Create partially resorbed rod
python one_sheet_hyperbola.py --spacing ${SPACING}  ${SPACING} ${SPACING} --alpha 1.75 --length 1.0 --small_radius 0.10 --large_radius 0.20 C4.nii > CREATE_C4.LOG

# Create resorbed plate
python torus.py --spacing ${SPACING} ${SPACING} ${SPACING} --alpha 1.75 --diameter 1.0 --thickness 0.20 --distance 0.1 C5.nii > CREATE_C5.LOG

# Create solid plate
python torus.py --spacing ${SPACING} ${SPACING} ${SPACING} --alpha 1.75 --diameter 1.0 --thickness 0.10 --distance -0.2 C6.nii > CREATE_C6.LOG

