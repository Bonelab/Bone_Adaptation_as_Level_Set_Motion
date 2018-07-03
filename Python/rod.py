# B.A. Besler 2018
# Bone Imaging Laboratory

import SimpleITK as sitk
import argparse
import numpy as np

parser = argparse.ArgumentParser(
    description='Generate a rod.',
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
    '--radius_1', default=0.1,
    type=float,
    help='Radius of the first side [mm]'
)
parser.add_argument(
    '--radius_2', default=0.15,
    type=float,
    help='Radius of the second side [mm]'
)
parser.add_argument(
    '--alpha', default=1.75,
    type=float,
    help='Padding around as a multiple of largest radius.'
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

if args.radius_1 <= 0:
    os.sys.exit('[ERROR] First radius cannot be less than or equal to zero ({} <= 0)'.format(args.radius_1))

if args.radius_2 <= 0:
    os.sys.exit('[ERROR] Second radius cannot be less than or equal to zero ({} <= 0)'.format(args.radius_2))

if args.length <= 0:
    os.sys.exit('[ERROR] Length cannot be less than or equal to zero ({} <= 0)'.format(args.length))

print('Computing a and theta')
a = np.abs(float(args.radius_2 - args.radius_1)) / float(args.length)
theta = np.degrees(np.arctan(a))
print('  a :    {:12.3e}'.format(a))
print('  theta: {:12.3e}'.format(theta))
print('')

print('Computing bounds')
spacing = np.array(args.spacing)
r_max = np.max([args.radius_1, args.radius_2])
r_min = np.min([args.radius_1, args.radius_2])
z_lo = 0
z_up = float(args.length)
x_lo = y_lo = float(-1 * args.alpha * r_max)
x_up = y_up = float(args.alpha * r_max)
low = np.array([x_lo, y_lo, z_lo])
up = np.array([x_up, y_up, z_up])
print('  lower: {}'.format(low))
print('  upper: {}'.format(up))

grid = np.array([int(np.round((up - low)/spacing)) for up, low, spacing in zip(up, low, args.spacing)])
print('  grid:  {}'.format(grid))
print('')

# Define our implicit function
def cone(point, origin, a):
    '''Implicit function of a cone.

    Tests the following equation:
        (x-x_0)^2 + (y-y_0)^2 - a^2 * (z-z_0)^2 <= 0
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
    test = (res[0]**2 + res[1]**2 - a**2 * res[2]**2 <= 0)
    return int(test)

def cone_test(point, origin, a):
    res = np.array(point) - np.array(origin)
    test = (res[0]**2 + res[1]**2 - a**2 * res[2]**2)
    print('{} - {} {} = {}'.format(point, origin, a, test))

print('Generating image...')
# Generate cone origin. Need to flip if the origin changed.
origin = low
center = (grid[0]-1)/2*spacing[0] + origin[0]
cone_origin = np.array([center, center, -float(r_min)/a])
if args.radius_1 > args.radius_2:
    cone_origin[2] = -1*cone_origin[2] + float(args.length)

print('  Center:      {}'.format(center))
print('  Cone Origin: {}'.format(cone_origin))

# Numpy indexes (z,y,x)
grid = grid[::-1]
data = np.zeros(grid, dtype=np.int8)
for z in range(grid[0]):
    for y in range(grid[1]):
        for x in range(grid[2]):
            # Create point in world coordinates
            point = np.array([x,y,z])*spacing + origin

            # Test implicit function
            data[z,y,x] = 127 * cone(point, cone_origin, a)

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

