
import SimpleITK as sitk
import argparse
import numpy as np

parser = argparse.ArgumentParser(
    description='Generate a cylinder.',
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
    '--radius', default=0.15,
    type=float,
    help='Radius of the cylinder [mm]'
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

if args.radius <= 0:
    os.sys.exit('[ERROR] Radius cannot be less than or equal to zero ({} <= 0)'.format(args.radius))

if args.length <= 0:
    os.sys.exit('[ERROR] Length cannot be less than or equal to zero ({} <= 0)'.format(args.length))

print('Computing bounds')
spacing = np.array(args.spacing)
z_lo = 0
z_up = float(args.length)
x_lo = y_lo = float(-1 * args.alpha * args.radius)
x_up = y_up = float(args.alpha * args.radius)
low = np.array([x_lo, y_lo, z_lo])
up = np.array([x_up, y_up, z_up])
print('  lower: {}'.format(low))
print('  upper: {}'.format(up))

grid = np.array([int(np.round((up - low)/spacing)) for up, low, spacing in zip(up, low, args.spacing)])
print('  grid:  {}'.format(grid))
print('')

# Define our implicit function
def cylinder(point, origin, radius):
    '''Implicit function of a cylinder.

    Tests the following equation:
        (x-x_0)^2 + (y-y_0)^2 - r^2 <= 0
    Note the weak inequality. Only implemented for three
    dimensions.

    Args:
        point: The point to test for membership.
        origin: The cone origin.
        radius: The cylinder radius

    Returns:
        An integer where one represents inside, zero outside.
    '''
    # Shift
    res = np.array(point) - np.array(origin)
    test = (res[0]**2 + res[1]**2 - radius**2 <= 0)
    return int(test)

print('Generating image...')
# Generate cylinder origin. Need to flip if the origin changed.
origin = low
center = (grid[0]-1)/2*spacing[0] + origin[0]
cone_origin = np.array([center, center, 0])

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
            data[z,y,x] = 127 * cylinder(point, cone_origin, float(args.radius))

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

