import os
import argparse
import shutil
from pymol import cmd
import enzy_htp as eh
from typing import List
from pathlib import Path
from enzy_htp.core import file_system as fs
from collections import namedtuple, defaultdict

import pymol2 #TODO(CJ): add this to the main enzy_htp part
import pymol.invocation
pymol.invocation.parse_args(['pymol', '-q']) # optional, for quiet flag
pymol2.SingletonPyMOL().start()

ValidatedArgs=namedtuple('ValidatedArgs', 'structure mutations pH conformer_engine n_conformers ligand_1 ligand_2 ligand_1_name ligand_2_name constraints')
ValidatedArgs.__doc__="""namedtuple() that holds the arguments from the commandline parser after being treated and validated."""


def check_residue_id( res_id : str ) -> None:
    """Helper function that checks if a supplied string is a valid residue id in the PDB format. Supplied
    code must: 1) be three characters in length and 2) be all uppercase. Strips whitespace.

    Args:
        res_id: The residue id to check as a str().

    Returns:
        Nothing.

    """ 
    res_id = ''.join(res_id.split())

    if len(res_id) != 3:
        eh.core._LOGGER.error(f"The supplied residue id '{res_id}' is invalid. Must be 3 characters long. Exiting...")
        exit( 1 )

    if not res_id.isupper():
        eh.core._LOGGER.error(f"The supplied residue id '{res_id}' is invalid. Must be uppercase. Exiting...")
        exit( 1 )

def parse_args() -> ValidatedArgs:
    """Commandline argument parser. Performs basic checks on the inputs and slightly rearranges them.

    Returns:
        A namedtuple with validated and treated arguments from the commandline.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'structure',
        type=str,
        help='File to the .pdb file containing the ligand/structure complex. Must exist and be in .pdb format.'
    )
    parser.add_argument(
        '--pH',
        type=float,
        default=7.0,
        help='Desired pH for protonation. Default is 7.0 and must be on range [0,14].'

    )
    parser.add_argument(
        '--mutations',
        type=str,
        help='The desired Mutations to apply in EnzyHTP format. If multiple, split by commas.'
    )
    parser.add_argument(
        '--conformer_engine',
        type=str,
        default='BCL',
        help='Engine for generating ligand conformers. Default is \'BCL\'. Allowed values are \'BCL\' and \'MOE\'.'
    )
    parser.add_argument(
        '--n_conformers',
        type=int,
        default=100,
        help='Number of conformers to generate. Default is 100. Must be positive.'
    )
    parser.add_argument(
        '--ligand_1',
        type=str,
        help='The path to ligand_1. Note this is REQUIRED.'
    )
    parser.add_argument(
        '--ligand_1_name',
        type=str,
        help='The three letter code for the residue. Note this REQUIRED.'
    )
    parser.add_argument(
        '--ligand_2',
        type=str,
        help='The path to ligand_2. Note this is NOT REQUIRED.'
    )
    parser.add_argument(
        '--ligand_2_name',
        type=str,
        help='The three letter code for the residue. Note this is NOT REQUIRED unless --ligand_2 is supplied.'
    )
    parser.add_argument( #TODO(CJ): try to make it so that this can also be specified
        '--constraints',
        type=str,
        help='The path to the constraints file. Not required but must exist if supplied.'
    )


    args = parser.parse_args()
    
    fs.check_file_exists( args.structure )

    if not Path(args.structure).suffix == '.pdb':
        eh.core._LOGGER.error(f'The supplied file "{args.structure}" is NOT in a .pdb format. Exiting...')
        exit( 1 )

    if args.pH < 0.0 or args.pH > 14.0:
        eh.core._LOGGER.error(f'The supplied pH {args.pH:0.2f} is not in the range [0.00, 14,00]. Exiting...')
        exit( 1 )


    mutations:List[str] = list()
    
    if args.mutations:
        assert False, "Need to implement this -CJ"
        #TODO(CJ): parser out the mutations and check that they are valid 
    
    if args.conformer_engine not in ["BCL", "MOE"]:
        eh.core._LOGGER.error(f'The supplied conformer engine "{args.conformer_engine}" is not supported. Allowed Options are BCL or MOE. Exiting...')
        exit( 1 )

    if args.n_conformers < 0:
        eh.core._LOGGER.error(f'--n_conformers must be a postiive value. Supplied value {args.n_conformers} is invalid. Exiting...')
        exit( 1 )

    if not args.ligand_1:
        eh.core._LOGGER.error("--ligand_1 is a required argument.")
        exit( 1 )
    
    if not args.ligand_1_name:
        eh.core._LOGGER.error("--ligand_1_name is a required argument.")
        exit( 1 )


    if args.ligand_2:
        if not args.ligand_2_name:
            eh.core._LOGGER.error("--ligand_2_name is a required argument. When --ligand_2 is supplied.")
            exit( 1 )

    ligand_1, ligand_2 = str(), str()
    ligand_1_name, ligand_2_name = str(), str()

    fs.check_file_exists( args.ligand_1 )

    ligand_1 = args.ligand_1

    check_residue_id( args.ligand_1_name )
    ligand_1_name = args.ligand_1_name

    if args.ligand_2:
        fs.check_file_exists( args.ligand_2 )
        ligand_2 = args.ligand_2
        check_residue_id( args.ligand_2_name )
        ligand_2_name = args.ligand_2_name


    constraints = str()

    if args.constraints:
        constraints = args.constraints
        fs.check_file_exists( constraints )

    validated = ValidatedArgs(
        structure=args.structure,
        pH=args.pH,
        mutations=mutations,
        conformer_engine=args.conformer_engine,
        n_conformers=args.n_conformers,
        ligand_1=ligand_1,
        ligand_2=ligand_1,
        ligand_1_name=ligand_1_name,
        ligand_2_name=ligand_2_name,
        constraints=constraints
    )


    return validated


def parameterize_ligand( fname ):
    #TODO(CJ): should probably deduce residue name from the file contents
    fpath = Path(fname)
    prefix=fpath.parent / fpath.stem.replace('_conformers','') #TODO(CJ): check that it replaces the one at the end
    
    cmd.delete('all')
    cmd.load(fpath)
    cmd.split_states('all')
    
    print(fname)
    
    name = fname.stem[:3]
    param_file = f"{name}.params"
    fs.safe_rm( param_file )

    cmd_str=f"$ROSETTA3/source/scripts/python/public/molfile_to_params.py {fname} --name={name} >/dev/null"
    print(cmd_str)
    result = os.system(cmd_str)
    if result:
        #TODO(CJ): add the error stuff here
        pass
    
    lig_pdbs = sorted(Path('.').glob(f'{name}_????.pdb'))
    shutil.copy(lig_pdbs[0],fpath.parent / f"{name}_base.pdb")
    conformers_lines:List[str] = list()

    for lp in lig_pdbs:
        for ll in fs.lines_from_file(lp):
            if not ll.startswith('ATOM') and not ll.startswith('HETATM'):
                continue
            conformers_lines.append( ll )

            fs.safe_rm( lp )

        conformers_lines.append('TER')
    
    conformers_lines.append('END')

    conf_file = fpath.parent/f"{name}_conformers.pdb"
    fs.write_lines( conf_file, conformers_lines )
    
    pfile_lines = fs.lines_from_file( param_file )
    pfile_lines.append(f"PDB_ROTATMERS {conf_file}\n")

    fs.write_lines( param_file, pfile_lines )

    shutil.move(param_file, fpath.parent / param_file )

    return param_file, conf_file 

def make_options_file( params : ValidatedArgs, lig_params_1:str, lig_params_2:str=None, work_dir:str=None ) -> str:
    """Function that creates the options file for a RosettaScripts run for the RosettaLigand protocol.

    Args:
        

    Returns:
        The path to the options file as a str().

    """
    if not work_dir:
        work_dir = str(Path(params.structure).parent)

    content:List[str] = [
        "-in:file",
        f"  -s '{params.structure}'",
        f"  -extra_res_fa '{lig_params_1}'",
    ]

    if lig_params_2:
        content.append(
            f"  -extra_res_fa '{lig_params_2}'",
        )
    
    content.extend("-run:preserve_header")
    content.append("-packing")
    content.append("    -ex1")
    content.append("    -ex2aro")
    content.append("    -ex2 ")
    content.append("    -no_optH false")
    content.append("    -flip_HNQ true")
    content.append("    -ignore_ligand_chi true")

    if constraints:
        content.extend(["-enzdes",f"    -cstfile '{params.constraints}'")

    content.extend(
        [
            "-parser",
            "   -protocol ligand_dock.xml", #TODO(CJ): fix this
            "-out",
           f"   -file:scorefile score.sc", #TODO(CJ): fix this
            "   -level 700",
            "   -nstruct 50", #TODO(CJ): change this
            "   -overwrite",
           f"   -path",
           f"       -all work_dir/complexes", #TODO(CJ): fix this

        ]
    )
   
    fname = Path(work_dir)/"options.txt"

    fs.write_lines( fname, content )

    return fname


def protonate_ligand( molfile : str ) -> str:
    """ """
    local = eh.interface.pymol.convert( molfile, new_ext='.mol2')
    local = eh.interface.pymol.de_protonate( local )
    protonated = eh.interface.moe.protonate( local )
    return eh.interface.pymol.convert( protonated, new_ext='.sdf') 

def main( params : ValidatedArgs ):
    """ Main routine """
    #TODO(CJ): give option to give the conformer files ahead of time 
    #1. Mutate if necessary
    if params.mutations:
        eh.mutate_pdb(
                params.structure,
                mutations=params.mutation,
                engine='rosetta'
            )

    #2. protonate ligands/enzyme
    #2.a protonated the structure
    parser = eh.structure.PDBParser()
    structure = parser.get_structure( params.structure )
    structure = eh.preparation.protonate_stru( structure )
    #2.b protonating the ligands 
    ligand_files:List[str] = list()

    ligand_1, ligand_2 = str(), str()

    ligand_1 = protonate_ligand( params.ligand_1 )
    ligand_1_conformers = eh.interface.bcl.generate_conformers( ligand_1 )
    ligand_1_parameters = eh.interface.rosetta.parameterize_ligand( ligand_1, params.ligand_1_name )

    if params.ligand_2:
        ligand_2 = protonate_ligand( params.ligand_2 )
        ligand_2_conformers = eh.interface.bcl.generate_conformers( ligand_2 )
        ligand_2_parameters = eh.interface.rosetta.parameterize_ligand( ligand_2, params.ligand_2_name )


    options_file:str = make_options_file( params, lig_params_1=ligand_1_parameters, lig_2_params=ligand_2_parameters, work_dir=path(params.structure).parent ) 

    
    #eh.interface.rosetta.run_rosetta_scripts(

if __name__ == '__main__':
    main( parse_args() ) 

