# B.A. Besler 2018
# Bone Imaging Laboratory

import SimpleITK as sitk
import argparse
import numpy as np

parser = argparse.ArgumentParser(
    description='Generate a torus (resorbed plate).',
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
    '--diameter', default=1.0,
    type=float,
    help='Diameter of the plate [mm]'
)
parser.add_argument(
    '--distance', default=0.10,
    type=float,
    help='Resorption distance [mm]'
)
parser.add_argument(
    '--thickness', default=0.10,
    type=float,
    help='Plate thickness [mm]'
)
parser.add_argument(
    '--alpha', default=1.75,
    type=float,
    help='Padding around as a multiple of diameter.'
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

if args.diameter <= 0:
    os.sys.exit('[ERROR] Diameter cannot be less than or equal to zero ({} <= 0)'.format(args.diameter))

if args.thickness <= 0:
    os.sys.exit('[ERROR] Thickness cannot be less than or equal to zero ({} <= 0)'.format(args.thickness))

print('Computing inner and outer radius')
inner_radius = float(args.diameter - args.distance)/4.0
outer_radius = float(args.diameter + args.distance)/4.0
b = float(args.diameter - args.distance)/(2.0 * float(args.thickness))
print('  Inner: {}'.format(inner_radius))
print('  Outer: {}'.format(outer_radius))
print('  b: {}'.format(b))
print('')

print('Computing bounds')
spacing = np.array(args.spacing)
z_lo = float(-1 * args.alpha * args.thickness)
z_up = float(args.alpha * args.thickness)
x_lo = y_lo = float(-1 * args.alpha * args.diameter)
x_up = y_up = float(args.alpha * args.diameter)
low = np.array([x_lo, y_lo, z_lo])
up = np.array([x_up, y_up, z_up])
print('  lower: {}'.format(low))
print('  upper: {}'.format(up))

grid = np.array([int(np.round((up - low)/spacing)) for up, low, spacing in zip(up, low, args.spacing)])
print('  grid:  {}'.format(grid))
print('')

# Define our implicit function
def torus(point, origin, inner_radius, outer_radius, b):
    '''Implicit function of a torus.

    Tests the following equation:
        (sqrt( (x-x_0)^2 + (y-y_0)^2) ) - R )^2 + b^2(z-z_0)^2 <= r^2
    Note the weak inequality. Only implemented for three
    dimensions.

    Args:
        point: The point to test for membership.
        origin: The cone origin.
        inner_radius: The torus inner radius, r (xy radius)
        outer_radius: The torus outer radius, R (offset from z)
        b: Thickness scaling constant (ratio of xy radius to z radius)

    Returns:
        An integer where one represents inside, zero outside.
    '''
    # Shift
    res = np.array(point) - np.array(origin)
    first = np.sqrt(res[0]**2 + res[1]**2) - outer_radius
    test = (first**2 + (b**2)*(res[2]**2) <= inner_radius**2)
    return int(test)

print('Generating image...')
# Generate torus origin.
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
            data[z,y,x] = 127 * torus(point, cone_origin, inner_radius, outer_radius, b)

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

