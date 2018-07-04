#! /bin/env/python
#
#   History:
#   2018.05.06  Besler  Created
#
#   Description:
#       Convert Aim to NII
#
#   Notes:
#       - python convertAim.py ../MODELS/* --output=../DOC/CSV/convertAim.csv

# Imports
import argparse
import os
import vtkbone
import vtk
import re
import logging

# Parse arguments
parser = argparse.ArgumentParser(description='Convert AIM to NII')
parser.add_argument(
    'input_files', nargs='+',
    help='Input files')
parser.add_argument(
    '--aim_ending', default='.AIM',
    help='Ending of AIM files')
parser.add_argument(
    '--nii_ending', default='.nii',
    help='Ending of NII files')
parser.add_argument(
    '--output', '-o',
    help='The output file to store results.')
parser.add_argument(
    '--delimiter', '-d', default=',',
    help='The output file to store results.')
args = parser.parse_args()

# Setup logger
logging.basicConfig(
    format='[%(asctime)s] - %(name)s - %(levelname)s - %(message)s'
    ,datefmt='%Y.%M.%d %I:%M:%S %p'
    ,level=logging.INFO # DEBUG for debug output
    ,filename='../LOG/convertAim.log'
)
log = logging.getLogger(os.path.splitext(os.path.basename(__file__))[0])

log.info('Arguments:')
for arg in vars(args):
    log.info('  {0: <16} {1}'.format(arg, getattr(args, arg)))
log.info('')

def write(entry):
    '''Yeah, I know it's strange'''
    # Determine what we're writting to
    if args.output == None:
        out = os.sys.stdout
    else:
        out = open(args.output, "a")

    # Write
    out.write(args.delimiter.join([str(x) for x in entry]))
    out.write(os.linesep)

    # Close
    if not out == os.sys.stdout:
        out.close()

# Do conversion
log.info('Converting files')
if args.output is None or not os.path.isfile(args.output):
    write(['Filename', 'aim_ending', 'nii_ending', 'mu_scaling', 'hu_mu_water', 'hu_mu_air', 'density_slope', 'density_intercept', 'm', 'b'])
n_files = 0
for input_file in args.input_files:
    output_file = input_file.replace(args.aim_ending, args.nii_ending)
    if output_file == input_file:
        log.info('   Skipping {} since it does not match ending {}'.format(input_file, args.aim_ending))
        continue
    else:
        log.info('  Converting {} to {}'.format(input_file, output_file))

    # Read AIM
    log.info('    Reading {}'.format(input_file))
    reader = vtkbone.vtkboneAIMReader()
    reader.SetFileName(input_file)
    reader.DataOnCellsOff()
    reader.Update()

    # Parse header
    proclog = reader.GetProcessingLog()
    mu_scaling_match = re.search(r'Mu_Scaling\s+(\d+)', proclog)
    hu_mu_water_match = re.search(r'HU: mu water\s+(\d+.\d+)', proclog)
    mu_scaling = int(mu_scaling_match.group(1))
    hu_mu_water = float(hu_mu_water_match.group(1))
    hu_mu_air = 0
    log.info('    Found the following calibration data in the header:')
    log.info('      mu_scaling: {}'.format(mu_scaling))
    log.info('      hu_mu_water: {}'.format(hu_mu_water))

    m = 1000.0 / mu_scaling / (hu_mu_water - hu_mu_air)
    b = -1000.0 * hu_mu_water / (hu_mu_water - hu_mu_air)

    den_slope_match = re.search(r'Density: slope[\s=]+([+-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+))', proclog)
    den_intercept_match = re.search(r'Density: intercept[\s=]+([+-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+))', proclog)
    den_slope = float(den_slope_match.group(1))
    den_intercept = float(den_intercept_match.group(1))
    log.info('      density_slope: {}'.format(den_slope))
    log.info('      density_intercept: {}'.format(den_intercept))

    m = 1.0 / mu_scaling * den_slope
    b = den_intercept

    # Apply shift
    log.info('    Converting native to HU')
    caster = vtk.vtkImageShiftScale()
    caster.SetOutputScalarTypeToFloat()
    caster.ClampOverflowOn()
    caster.SetInputConnection(reader.GetOutputPort())
    caster.SetShift(b/m)
    caster.SetScale(m)
    caster.Update()

    log.info('    Writing to {}'.format(output_file))
    writer = vtk.vtkNIFTIImageWriter()
    writer.SetInputConnection(caster.GetOutputPort())
    writer.SetFileName(output_file)
    writer.Update()

    entry = [input_file, args.aim_ending, args.nii_ending, mu_scaling, hu_mu_water, hu_mu_air, den_slope, den_intercept, m, b]
    write(entry)

    n_files += 1

log.info('')
log.info('Converted {} files'.format(n_files))
