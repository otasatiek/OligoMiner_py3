#!/usr/bin/env python
# --------------------------------------------------------------------------
# OligoMiner
# probeRC.py
#
# Modified for Python 3 and modern Biopython compatibility.
# --------------------------------------------------------------------------

# Specific script name.
scriptName = 'probeRC'

# Specify script version.
Version = '1.7.1-py3-fix'

# Import module for handling input arguments.
import argparse

# Import Biopython modules.
from Bio.Seq import Seq
# from Bio.Alphabet import IUPAC # Removed: Deprecated in Biopython 1.78+

def createRCs(inputFile, outNameVal):
    """Creates a .bed file with the reverse complements of the given set of
    sequences."""

    # Determine the stem of the input filename.
    fileName = str(inputFile).split('.')[0]

    # Open input file for reading.
    with open(inputFile, 'r') as f:
        file_read = [line.strip() for line in f]

    # Create list to hold output.
    outList = []

    # Parse out probe info, flip sequence to RC, and write to output list.
    for i in range(0, len(file_read), 1):
        # Skip empty lines if any
        if not file_read[i]:
            continue
            
        parts = file_read[i].split('\t')
        if len(parts) < 5:
            continue # Skip lines that don't have enough columns

        chrom = parts[0]
        start = parts[1]
        stop = parts[2]
        probeSeq = parts[3]
        
        # Modified: No longer using IUPAC.unambiguous_dna
        RevSeq = Seq(probeSeq).reverse_complement()
        
        Tm = parts[4]
        outList.append('%s\t%s\t%s\t%s\t%s' % (chrom, start, stop, RevSeq, Tm))

    # Determine the name of the output file.
    if outNameVal is None:
        outName = '%s_RC' % fileName
    else:
        outName = outNameVal

    # Create the output file.
    output = open('%s.bed' % outName, 'w')

    # Write the output file
    output.write('\n'.join(outList))
    output.close()


def main():
    """Produces a .bed file with the reverse complements of the given
    chromosome sequences."""

    # Allow user to input parameters on command line.
    userInput = argparse.ArgumentParser(description=\
        '%s version %s. Requires a .bed file with first four columns in the '
        'format chromosome <tab> start <tab> stop <tab> sequence <tab> Tm such '
	    ' as the .bed files produced by outputClean. Returns a .bed file that is '
        'identical to the input file except that the probe sequences have been '
        'replaced with their reverse complements.' % (scriptName, Version))
    requiredNamed = userInput.add_argument_group('required arguments')
    requiredNamed.add_argument('-f', '--file', action='store', required=True,
                               help='The .bed file containing the probe '
                                    'sequences to take the reverse complements '
                                    'of')
    userInput.add_argument('-o', '--output', action='store', default=None,
                           type=str, help='Specify the name prefix of the '
                                          'output file')

    # Import user-specified command line values
    args = userInput.parse_args()
    inputFile = args.file
    outNameVal = args.output

    createRCs(inputFile, outNameVal)

if __name__ == '__main__':
    main()