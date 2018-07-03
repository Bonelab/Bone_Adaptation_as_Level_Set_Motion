# B.A. Besler 2018
# Bone Imaging Laboratory

import SimpleITK as sitk
import argparse
import numpy as np
import os

parser = argparse.ArgumentParser(
    description='Generate a one sheet hyperbola (partially resorbed rod).',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument(
    'output_file_name',
    help='Output file name'
)
parser.add_argument(
    '--spacing', default=[0.0607 for _ in range(3)],
    type=float, nargs='+',
    help='Resolution in x,y,z [mm]'
)
parser.add_argument(
    '--length', default=2.0,
    type=float,
    help='Length of the rod [mm]'
)
parser.add_argument(
    '--small_radius', default=0.05,
    type=float,
    help='Radius at the center of the rod [mm]'
)
parser.add_argument(
    '--large_radius', default=0.20,
    type=float,
    help='Radius at the end of the rod [mm]'
)
parser.add_argument(
    '--alpha', default=1.75,
    type=float,
    help='Padding around as a multiple of radius.'
)
args = parser.parse_args()

# Print Arguments
print('Arguments')
for arg in vars(args):
    print('  {}: {}'.format(arg, getattr(args, arg)))
print('')

# Test inputs
for space in args.spacing:
    if space <= 0:
        os.sys.exit('[ERROR] Spacing cannot be less than or equal to zero ({} <= 0)'.format(space))

if args.small_radius <= 0:
    os.sys.exit('[ERROR] Small radius cannot be less than or equal to zero ({} <= 0)'.format(args.small_radius))

if args.large_radius <= 0:
    os.sys.exit('[ERROR] Large radius cannot be less than or equal to zero ({} <= 0)'.format(args.large_radius))

if args.large_radius < args.small_radius:
    os.sys.exit('[ERROR] Large radius cannot be smaller than small radius ({} < {})'.format(args.large_radius, args.small_radius))

if args.length <= 0:
    os.sys.exit('[ERROR] Length cannot be less than or equal to zero ({} <= 0)'.format(args.length))

print('Computing a and b')
a = 2.0 / float(args.length) * np.sqrt(float(args.large_radius)**2 / float(args.small_radius)**2 -1)
b = 1.0 / float(args.small_radius)
print('  a :    {:12.3e}'.format(a))
print('  b:     {:12.3e}'.format(b))
print('')

print('Computing bounds')
spacing = np.array(args.spacing)
z_lo = 0
z_up = float(args.length)
x_lo = y_lo = float(-1 * args.alpha * args.large_radius)
x_up = y_up = float(args.alpha * args.large_radius)
low = np.array([x_lo, y_lo, z_lo])
up = np.array([x_up, y_up, z_up])
print('  lower: {}'.format(low))
print('  upper: {}'.format(up))

grid = np.array([int(np.round((up - low)/spacing)) for up, low, spacing in zip(up, low, args.spacing)])
print('  grid:  {}'.format(grid))
print('')

# Define our implicit function
def hyperbola(point, origin, a, b):
    '''Implicit function of a hyperbola.

    Tests the following equation:
        (x-x_0)^2 + (y-y_0)^2 - a^2 * (z-z_0)^2 <= 1
    Note the weak inequality. Only implemented for three
    dimensions.

    Args:
        point: The point to test for membership.
        origin: The cone origin.
        a: The slope scaling parameter.

    Returns:
        An integer where one represents inside, zero outside.
    '''
    # Shift
    res = np.array(point) - np.array(origin)
    test = (-1 * b**2 * res[0]**2 - b**2 * res[1]**2 + a**2 * res[2]**2 >= -1)
    return int(test)

print('Generating image...')
# Generate hyperbola origin.
origin = low
center = (grid[0]-1)/2*spacing[0] + origin[0]
center_z = (grid[2]-1)/2*spacing[2] + origin[2]
cone_origin = np.array([center, center, center_z])

print('  Center: {}'.format(center))
print('  Origin: {}'.format(cone_origin))

# Numpy indexes (z,y,x)
grid = grid[::-1]
data = np.zeros(grid, dtype=np.int8)
for z in range(grid[0]):
    for y in range(grid[1]):
        for x in range(grid[2]):
            # Create point in world coordinates
            point = np.array([x,y,z])*spacing + origin

            # Test implicit function
            data[z,y,x] = 127 * hyperbola(point, cone_origin, a, b)

# Create ITK image
out = sitk.GetImageFromArray(data)
out.SetSpacing(args.spacing)
out.SetOrigin(origin)
print('  Size:         {}'.format(out.GetSize()))
print('  Spacing:      {}'.format(out.GetSpacing()))
print('  Origin:       {}'.format(out.GetOrigin()))
print('')

print('Writing to {}'.format(args.output_file_name))
sitk.WriteImage(out, args.output_file_name)

