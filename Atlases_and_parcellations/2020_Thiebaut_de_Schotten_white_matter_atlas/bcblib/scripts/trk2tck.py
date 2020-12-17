#!/usr/bin/env python3
# -*-coding:Utf-8 -*

""" trk2tck will convert trk files (from TrackVis) into tck files (MRtrix)
This script has been written by Marc-Alexandre Côté, you can find his git
repository on https://gist.github.com/MarcCote
"""

import os
import argparse
import nibabel as nib

__author__ = "Marc-Alexandre Côté"

def build_argparser():
    DESCRIPTION = "Convert tractograms (TRK -> TCK)."
    p = argparse.ArgumentParser(description=DESCRIPTION)
    p.add_argument('tractograms', metavar='bundle', nargs="+",
                   help='list of tractograms.')
    p.add_argument('-f', '--force', action="store_true",
                   help='overwrite existing output files.')
    return p


def main():
    print("This scripts has been written by Marc-Alexandre Côté")
    print("You can find his repository on github at: ")
    print("https://gist.github.com/MarcCote")
    parser = build_argparser()
    args = parser.parse_args()
    for tractogram in args.tractograms:
        if (nib.streamlines.detect_format(tractogram) is not
            nib.streamlines.TrkFile):
            print("Skipping non TRK file: '{}'".format(tractogram))
            continue

        output_filename = tractogram[:-4] + '.tck'
        if os.path.isfile(output_filename) and not args.force:
            print("Skipping existing file: '{}'. Use -f to overwrite.".format(
                output_filename))
            continue

        trk = nib.streamlines.load(tractogram)
        nib.streamlines.save(trk.tractogram, output_filename)

if __name__ == '__main__':
    main()
