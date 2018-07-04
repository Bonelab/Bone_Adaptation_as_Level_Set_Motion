#! /bin/env/python
#
#   History:
#   2018.05.06  Besler  Created
#
#   Description:
#       Convert Aim to NII
#
#   Notes:
#       - python convertBinaryAim.py ../MODELS/*

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
args = parser.parse_args()

# Setup logger
logging.basicConfig(
    format='[%(asctime)s] - %(name)s - %(levelname)s - %(message)s'
    ,datefmt='%Y.%M.%d %I:%M:%S %p'
    ,level=logging.INFO # DEBUG for debug output
    ,filename='../LOG/convertBinaryAim.log'
)
log = logging.getLogger(os.path.splitext(os.path.basename(__file__))[0])

log.info('Arguments:')
for arg in vars(args):
    log.info('  {0: <16} {1}'.format(arg, getattr(args, arg)))
log.info('')

# Do conversion
log.info('Converting files')
n_files = 0
for input_file in args.input_files:
    output_file = input_file.replace(args.aim_ending, args.nii_ending)
    if output_file == input_file:
        log.info('   Skipping {} since it does not match ending {}'.format(input_file, args.aim_ending))
        continue
    else:
        log.info('  Converting {} to {}'.format(input_file, output_file))

    # Read AIM
    reader = vtkbone.vtkboneAIMReader()
    reader.SetFileName(input_file)
    reader.DataOnCellsOff()

    writer = vtk.vtkNIFTIImageWriter()
    writer.SetInputConnection(reader.GetOutputPort())
    writer.SetFileName(output_file)
    writer.Update()

    n_files += 1

log.info('')
log.info('Converted {} files'.format(n_files))
