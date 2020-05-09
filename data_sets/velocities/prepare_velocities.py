import os
import itertools
from uafgi import make,glaciers
import sys
sys.path.append('.')
import data


def main():

    # Create a blank makefile
    makefile = make.Makefile()

    # Add a rule to it, to create 
    filter_attrs=dict(source='TSX', grid='W69.10N')
    rule = glaciers.merge(makefile, 'data', 'outputs',
        os.path.join('outputs', '{source}_{grid}_2008_2020.nc'),
        ('vx','vy'),
        #max_files=3,
        filter_attrs=filter_attrs,
        blacklist=data.blacklist).rule

    # Build the outputs of that rule
    make.build(makefile, rule.outputs)

    # Remove intermediate files
    make.cleanup(makefile, rule.outputs)


main()
