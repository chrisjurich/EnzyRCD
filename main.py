import argparse
from collections import namedtuple

ValidatedArgs=namedtuple('ValidatedArgs', 'mutations')



def parse_args():
    parser = argparse.ArgumentParser()
   

    parser.add_argument(
        'structure',
        type=str,
        help='File to the .pdb file containing the ligand/structure complex. Must exist and be in .pdb format'
    )

    parser.add_argument(
        '--ph',
        type=float,
        default=7.0,
        help='Desired pH for protonation. Default is 7.0 and must be on range [0,14]'

    )
    parser.add_argument(
        '--mutation',
        type=str,
        help='The desired Mutations to apply in EnzyHTP format. If multiple, split by commas'
    )
    parser.add_argument(
        '--conformer_engine',
        type=str,
        default='BCL',
        help='Engine for generating ligand conformers. Default is \'BCL\'. Allowed values are \'BCL\' and \'MOE\''
    )
    parser.add_argument(
        '--n_conformers',
        type=int,
        default=100,
        help='Number of conformers to generate. Default is 100. Must be positive'
    )

    #TODO(CJ): add the ability to add constraints for the structure

    args = parser.parse_args()
    print(args)
    
    
    #TODO(CJ): add argument checks

    
    validated = ValidatedArgs(
        mutations=None

    )


    return validated






def main( params ):

    pass
    #1. Mutate if necessary

    #2. protonate ligands/enzyme

    #3. relax

    #4. make conformers

    #5. Generate input files
    #   + what do I do here
    #   + param files for ligands
    #   + encode the xml for rosettascripts

    #6. Run it




if __name__ == '__main__':
    main( parse_args() ) 
