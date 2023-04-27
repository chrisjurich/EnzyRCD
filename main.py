import os
import time
import shutil
import argparse
import pandas as pd
import enzy_htp as eh
from rdkit import Chem
from pathlib import Path
from pymol import cmd, stored
from typing import List, Tuple
from enzy_htp.core import file_system as fs
from collections import namedtuple, defaultdict

import pymol2 #TODO(CJ): add this to the main enzy_htp part
import pymol.invocation
pymol.invocation.parse_args(['pymol', '-q']) # optional, for quiet flag
pymol2.SingletonPyMOL().start()

ValidatedArgs=namedtuple(
    'ValidatedArgs', 
    'structure mutations pH conformer_engine n_conformers ligand_1 ligand_2 ligand_1_name ligand_2_name constraints n_structures'
)
ValidatedArgs.__doc__="""namedtuple() that holds the arguments from the commandline parser after being treated and validated."""
#TODO(CJ): add documentation for the attributes 

def make_df( sele:str ) -> pd.DataFrame:
    stored.holder = []
    
    cmd.iterate_state(-1, sele, 'stored.holder.append( (name, elem, x, y, z, ID, chain, resn ) )')
    
    df =  pd.DataFrame(
        columns='aname elem x y z ID chain resn'.split(),
        data=stored.holder
    )
    df.sort_values(by='ID', inplace=True)
    return df.reset_index(drop=True)


def check_mutations( params ) -> List[str]:
    """Function that checks if the supplied mutations are valid within the EnzyHTP framework. Takes
    in the parameters namedtuple() from the commandline parser and returns a list() of validated 
    mutation codes. Assumes that codes are ',' delimited. Will error and exit if any codes are invalid.
    Args:
        params: The namedtuple() from the ArgumentParser().
    Returns:
        The list() of validated EnzyHTP codes in correct format as str().
    """
    mutations:List[str] = list()
    
    if not params.mutations:
        return mutations
        
    error = False
    invalid_mutations:List[str] = list()
    
    for raw_mut in params.mutations.split(','):
        if not eh.mutation.valid_mutation( raw_mut ):
            invalid_mutations.append( raw_mut )
            error = True
        mutations.append( raw_mut )
    
    if error:
        eh.core._LOGGER.error(f"The following supplied mutations are invalid: {','.join(invalid_mutations)}. Exiting...")
        exit( 1 )
    
    return mutations


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
    parser.add_argument(#TODO(CJ): finish this and add checking for valid values
        '--n_struct',
        type=int,
        help='Number of structures to make in the RosettaScripts portion.',
        default=50
    )

    args = parser.parse_args()
    
    fs.check_file_exists( args.structure )

    if not Path(args.structure).suffix == '.pdb':
        eh.core._LOGGER.error(f'The supplied file "{args.structure}" is NOT in a .pdb format. Exiting...')
        exit( 1 )

    if args.pH < 0.0 or args.pH > 14.0:
        eh.core._LOGGER.error(f'The supplied pH {args.pH:0.2f} is not in the range [0.00, 14,00]. Exiting...')
        exit( 1 )

    if args.n_struct <= 0.0:
        eh.core._LOGGER.error(f'The supplied nstruct is {argsn.n_struct}. Must be greater than 0. Exiting...')
        exit( 1 )

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

    mutations:List[str] = check_mutations( args )

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

    return ValidatedArgs(
        structure=args.structure,
        pH=args.pH,
        mutations=mutations,
        conformer_engine=args.conformer_engine,
        n_conformers=args.n_conformers,
        ligand_1=ligand_1,
        ligand_2=ligand_2,
        ligand_1_name=ligand_1_name,
        ligand_2_name=ligand_2_name,
        constraints=constraints,
        n_structures=args.n_struct
    )

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
        params:
        lig_params_1:
        lig_params_2:
        work_dir:
        
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
    
    content.extend(
        [
            "-run:preserve_header",
            "-packing",
            "    -ex1",
            "    -ex2aro",
            "    -ex2 ",
            "    -no_optH false",
            "    -flip_HNQ true",
            "    -ignore_ligand_chi true",
        ]
    )

    if params.constraints:
        content.extend(
            [
                "-enzdes",
                f"    -cstfile '{params.constraints}'"
            ]
        )

    content.extend(
        [
            "-parser",
            "   -protocol xmls/lig_dock2.xml", #TODO(CJ): fix this
            "-out",
           f"   -file:scorefile 'score.sc'", 
            "   -level 200",
           f"   -nstruct {params.n_structures}", #TODO(CJ): change this
            "   -overwrite",
            "   -path",
           f"       -all '{work_dir}/complexes'",

        ]
    )
   
    fname = Path(work_dir)/"options.txt"

    fs.write_lines( fname, content )

    scorefile:str = f"{work_dir}/complexes/score.sc"
    return (fname, scorefile)

def fix_conformers( template: str, conformers:str) -> str:
    """ """
    cmd.delete('all')
    cmd.load(template)
    template_df:pd.DataFrame = make_df('all')
    cmd.delete('all')

    cmd.load(conformers)
    
    original = cmd.get_object_list()
    assert len(original) == 1
    original = original[0]
    content:List[str]=list()
    cmd.split_states('all')
    for oidx, oo in enumerate(cmd.get_object_list()):
        if oo == original:
            continue
        df = make_df( oo )
        assert len(df) == len(template_df), f"{len(df)} {len(template_df)}"
        for (tidx, trow), (idx, row) in zip(template_df.iterrows(), df.iterrows()):
            assert trow.elem == row.elem
            cmd.alter(f"{oo} and ID {row.ID}", f"name='{trow.aname}'")
            cmd.alter(f"{oo} and ID {row.ID}", f"chain='{trow.chain}'")
            cmd.alter(f"{oo} and ID {row.ID}", f"resn='{trow.resn}'")

        temp_fname = f"state_{oidx}.pdb"
        cmd.save(temp_fname, oo )

        for ll in fs.lines_from_file( temp_fname ):
            if ll.startswith('HETATM') or ll.startswith('ATOM'):
                content.append( ll )
        
        content.append('TER')
        fs.safe_rm( temp_fname )

    content.pop()
    content.append('END')

    outfile = Path(conformers).with_suffix('.pdb')
    fs.write_lines(outfile, content)

    return str( outfile )


def protonate_ligand( molfile : str ) -> str:
    """ """
    local = eh.interface.pymol.convert( molfile, new_ext='.mol2')
    local = eh.interface.pymol.de_protonate( local )
    protonated = eh.interface.moe.protonate( local )
    return eh.interface.pymol.convert( protonated, new_ext='.sdf') 

def log_elapsed_time( elapsed : int ) -> None:
    """ """ 
    days, hours, minutes, seconds = 0, 0, 0, 0
    
    minutes_denom = 60
    hours_denom = minutes_denom*60
    days_denom = hours_denom*24
    
    if elapsed >= days_denom:
        days = int(elapsed / days_denom)
        elapsed //= days_denom

    if elapsed >= hours_denom:
        hours = int(elapsed / hours_denom)
        elapsed //= hours_denom

    if elapsed >= minutes_denom:
        minutes = int(elapsed / minutes_denom)
        elapsed //= minutes_denom

    seconds = int(elapsed)

    eh.core._LOGGER.info(f"Elapsed time: {days} days {hours} hours {minutes} minutes {seconds} seconds")


def kekulize( in_file:str )->str:   
    """ """ 
    mol2_infile = Path(in_file).with_suffix('.mol2')
    eh.interface.pymol.convert( in_file, new_ext='.mol2')
    mol = Chem.MolFromMol2File(str(mol2_infile), removeHs=False)
    mol_file = Path(in_file).with_suffix('.mol')
    Chem.MolToMolFile(mol, str(mol_file),  kekulize=True) 
    eh.interface.pymol.convert( mol_file, new_ext='.sdf')
    #w = Chem.SDWriter(in_file )
    #w.SetKekulize(True)
    #w.write(mol)
    #w.close()
    #return in_file

def main( params : ValidatedArgs ) -> None:
    """ Main routine """
    start = time.time()
    #TODO(CJ): give option to give the conformer files ahead of time 
    #1. Mutate if necessary
    if params.mutations:
        eh.mutate_pdb(
                params.structure,
                mutations=params.mutation,
                engine='rosetta'
            )
        #TODO(CJ): I think i need to fast relax here?                

    parser = eh.structure.PDBParser()
    structure = parser.get_structure( params.structure )
    structure = eh.preparation.protonate_stru( structure )
    
    ligand_files:List[str] = list()

    ligand_1, ligand_2 = str(), str()

    ligand_1:str = protonate_ligand( params.ligand_1 )
    ligand_1_parameters:str = eh.interface.rosetta.parameterize_ligand( ligand_1, params.ligand_1_name )
    kekulize( ligand_1 )
    ligand_1_conformers:str = eh.interface.bcl.generate_conformers( ligand_1 )
    ligand_1_conformers:str = fix_conformers( ligand_1_parameters[1], ligand_1_conformers)

    if params.ligand_2:
        ligand_2:str = protonate_ligand( params.ligand_2 )
        ligand_2_parameters:str = eh.interface.rosetta.parameterize_ligand( ligand_2, params.ligand_2_name )
        kekulize( ligand_2 )
        ligand_2_conformers:str = eh.interface.bcl.generate_conformers( ligand_2 )
        ligand_2_conformers:str = fix_conformers( ligand_2_parameters[1], ligand_2_conformers)


    (options_file, score_file) = make_options_file(
            params,
            lig_params_1=ligand_1_parameters[0],
            lig_params_2=ligand_2_parameters[0],
            work_dir=Path(params.structure).parent
        ) 
    
    fs.safe_rm( score_file )    

    eh.interface.rosetta.run_rosetta_scripts( [f"@{options_file}"] )
    df:pd.DataFrame = eh.interface.rosetta.parse_score_file( score_file )

    #show top 5
    df.sort_values(by='total_score',inplace=True)
    df.reset_index(drop=True, inplace=True)
    idx = min(20, len(df))    
    eh.core._LOGGER.info(f"EnzyRCD run complete. Top {idx} structures are:")
    eh.core._LOGGER.info(f"No. \t    REU\t Structure")
    
    for i, row in df.iterrows():
        if i == idx:
            break
        eh.core._LOGGER.info(f"{i+1: 3}\t{row.total_score:6.2f}\t{row.description}")

    elapsed = time.time() - start
    log_elapsed_time( elapsed )

if __name__ == '__main__':
    main( parse_args() ) 
