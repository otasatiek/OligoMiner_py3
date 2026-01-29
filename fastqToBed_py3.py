#!/usr/bin/env python
# --------------------------------------------------------------------------
# OligoMiner
# fastqToBed.py
#
# (c) 2017 Molecular Systems Lab
#
# Wyss Institute for Biologically-Inspired Engineering
# Harvard University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# --------------------------------------------------------------------------
#
# This file is a modified version of "fastqToBed.py" originally part of
# OligoMiner. Modified for Python 3 and modern Biopython compatibility by
# Keita Sato.

# Specific script name.
scriptName = 'fastqToBed'

# Specify script version.
Version = '1.7-py3-fix'

# Import module for handling input arguments.
import argparse
import sys

# Import Biopython mt module.
try:
    from Bio.SeqUtils import MeltingTemp as mt
except ImportError:
    print("Error: Biopython is not installed. Please install it with 'pip install biopython'.")
    sys.exit(1)

def probeTm(seq1, saltConc, formConc):
    """Calculates the melting temperature of a given sequence under the
    specified salt and formamide conditions."""
    
    # Modern Biopython Tm_NN returns a float directly.
    # Na parameter expects mM concentration (e.g., 390).
    try:
        tmval = mt.Tm_NN(seq1, Na=saltConc)
        
        # Apply chemical correction for formamide.
        # fmd expects formamide % v/v (e.g., 50).
        fcorrected = mt.chem_correction(tmval, fmd=formConc)
        
        # Return as string formatted to 2 decimal places
        return '%0.2f' % fcorrected
    except ValueError as e:
        # Handle cases where sequence might contain invalid characters (e.g. 'N')
        # In strict mode, Tm_NN might fail.
        return 'NaN'


def convertFastqToBed(inputFile, saltConc, formConc, outNameVal):
    """Converts a given .fastq file to a .bed file."""

    # Determine the stem of the input filename.
    fileName = str(inputFile).split('.')[0]

    try:
        # Open input file for reading with explicit encoding.
        with open(inputFile, 'r', encoding='utf-8') as f:
            # Read all lines (Note: For very large files, line-by-line processing is better,
            # but we keep list logic to match original structure for compatibility)
            file_read = [line.strip() for line in f]
    except UnicodeDecodeError:
        # Fallback for systems/files with legacy encoding
        with open(inputFile, 'r', encoding='latin-1') as f:
            file_read = [line.strip() for line in f]

    # Create list to hold output.
    outList = []

    # Parse .fastq and extract probe information.
    # Checks added to ensure file format is respected.
    for i in range(0, len(file_read), 4):
        # Ensure we don't go out of bounds
        if i + 1 >= len(file_read):
            break
            
        header = file_read[i]
        seq = file_read[i+1]
        
        # Check if header matches expected format @chr:start-stop
        if not header.startswith('@') or ':' not in header or '-' not in header:
            continue

        try:
            # Extract chrom, start, stop
            # Expected format: @chr:start-stop
            # Split by ':' first -> ['@chr', 'start-stop']
            chrom_part = header.split(':')[0]
            coords_part = header.split(':')[1]
            
            chrom = chrom_part.split('@')[1]
            start = coords_part.split('-')[0]
            stop = coords_part.split('-')[1]
            
            Tm = probeTm(seq, saltConc, formConc)
            outList.append('%s\t%s\t%s\t%s\t%s' % (chrom, start, stop, seq, Tm))
            
        except IndexError:
            # Skip malformed lines
            continue

    # Determine the name of the output file.
    if outNameVal is None:
        outName = fileName
    else:
        outName = outNameVal

    # Create the output file.
    with open('%s.bed' % outName, 'w', encoding='utf-8') as output:
        output.write('\n'.join(outList))


def main():
    """Converts a .bed file to a .fastq file, taking the filenames as
    command line arguments."""
    
    userInput = argparse.ArgumentParser(description=\
        '%s version %s. Requires a .fastq file containing chr, start, stop '
        'information in the sequence ID field for each entry in the format '
        '@chr:start-stop. Returns a .bed file.' % (scriptName, Version))
        
    requiredNamed = userInput.add_argument_group('required arguments')
    requiredNamed.add_argument('-f', '--file', action='store', required=True,
                               help='The .fastq file to convert to .bed')
    userInput.add_argument('-s', '--salt', action='store', default=390,
                           type=float, # Changed to float to be safe, though int works
                           help='The mM Na+ concentration, default is 390')
    userInput.add_argument('-F', '--formamide', action='store', default=50,
                           type=float,
                           help='The percent formamide being used, default is 50')
    userInput.add_argument('-o', '--output', action='store', default=None,
                           type=str,
                           help='Specify the name prefix of the output file')

    # Import user-specified command line values
    args = userInput.parse_args()
    inputFile = args.file
    saltConc = args.salt
    formConc = args.formamide
    outNameVal = args.output

    convertFastqToBed(inputFile, saltConc, formConc, outNameVal)


if __name__ == '__main__':
    main()
