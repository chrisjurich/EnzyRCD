from pathlib import Path
from typing import Any, Union, List


import pandas as pd
import enzy_htp as eh
from enzy_htp.core import file_system as fs

import utils as uu

class RosettaCst:
    """TODO"""
    def __init__(self, **kwargs):
        """"""
        self.rname_1 = None 
        self.rnum_1 = None
        self.ratoms_1 = None
        self.rchain_1 = None 
        self.rname_2 = None
        self.rnum_2 = None
        self.ratoms_2 = None
        self.rchain_2 = None
        self.constraints = None

        ALLOWED:List[str]='rname_1 rnum_1 ratoms_1 rchain_1 rname_2 rnum_2 ratoms_2 rchain_2 constraints'.split()
        
        for k,v in kwargs.items():
            if k not in ALLOWED:
                eh.core._LOGGER.error(f"The supplied keyword '{k}' is not valid for the RosettaCst() ctor. Exiting...")
                exit( 1 )
            setattr(self, k, v)

            ALLOWED.remove( k )

        if ALLOWED:
            eh.core._LOGGER.error(f"Missing the following keyword arguments for RosettaCst() ctor: {', '.join(ALLOWED)}. Exiting...")
            exit( 1 )


    def set(self, key:str, value: Any) -> None:
        """
        """
        if key not in self.__dict__:
            assert False

        self.__dict__[key] = value


    def satisfied(self, file:str) -> bool:
        """ 
        """
        pass




def validate_cst(start: Union[str, Path, pd.DataFrame],
                 cst: RosettaCst) -> None:
    """Function that checks if all the atoms for the supplied enzymatic constraint exist in the 
    supplied .pdb file. Logs all invalid atoms for the specific constraint and exits afterwards.

    Args:
        start: Either a str()/Path() to a .pdb file or a pandas DataFrame from one.
        cst: A RosettaCst() namedtuple().

    Returns:
        Nothing.
    """

    bad_atom: bool = False

    if type(start) == pd.DataFrame:
        df: pd.DataFrame = start
    else:
        fs.check_file_exists(str(start))
        df = eh.interface.pymol.collect(
            start, "resi name chain resn".split()
        )

    for ratom in cst.ratoms_1:

        row: pd.Series = df[(df['name'] == ratom) & (df.resn == cst.rname_1) &
                            (df.resi == str(cst.rnum_1)) &
                            (df.chain == cst.rchain_1)]

        if len(row) == 0:
            eh.core._LOGGER.error(
                f'Atom {cst.rchain_1}.{cst.rnum_1}.{cst.rname_1}.{ratom} does not exist'
            )
            bad_atom = True

    for ratom in cst.ratoms_2:

        row: pd.Series = df[(df['name'] == ratom) & (df.resn == cst.rname_2) &
                            (df.resi == str(cst.rnum_2)) &
                            (df.chain == cst.rchain_2)]

        if len(row) == 0:
            eh.core._LOGGER.error(
                f'Atom {cst.rchain_2}.{cst.rnum_2}.{cst.rname_2}.{ratom} does not exist'
            )
            bad_atom = True

    if bad_atom:
        eh.core._LOGGER.error("Errors in supplied constraints. Exiting...")
        exit(1)


def parse_rosetta_cst(raw: str) -> List[RosettaCst]:
    """Takes raw constraints and converts them into RosettaCst() namedtuple()'s,  
    
    Args:
        raw:

    Returns:
        A list of RosettaCst() namedtuple()'s.

    """
    raw = raw[1:-1]
    var = {'constraints': []}
    ALLOWED_CSTS: Set[str] = set(
        "distanceAB angle_A angle_B torsion_A torsion_B torsion_AB".split())

    tokens: List[str] = list(filter(len, raw.split(')(')))

    if len(tokens) < 3:
        eh.core._LOGGER.info(
            f"There must be at least 3 blocks in an individual constraint. There are only {len(tokens)} in '{raw}'. Exiting..."
        )
        exit(1)

    for tidx, tk in enumerate(tokens):

        spl: List[str] = tk.split(',')
        if tidx < 2:
            var[f"rchain_{tidx+1}"] = spl[0]
            var[f"rnum_{tidx+1}"] = int(spl[1])
            var[f"rname_{tidx+1}"] = spl[2]
            var[f"ratoms_{tidx+1}"] = spl[3:]
        else:
            if spl[0] not in ALLOWED_CSTS:
                eh.core._LOGGER.error(
                    f"The supplied constraint type {spl[0]} is not supported. Allowed are: {', '.join(sorted(list(ALLOWED_CSTS)))}. Exiting..."
                )
                exit(1)

            temp = [spl[0]]
            for tt in spl[1:]:
                if tt.find('.') == -1:
                    temp.append(int(tt))
                else:
                    temp.append(float(tt))

            var['constraints'].append(temp)
    
    return RosettaCst(rname_1=var['rname_1'],
                      rnum_1=var['rnum_1'],
                      ratoms_1=var['ratoms_1'],
                      rchain_1=var['rchain_1'],
                      rname_2=var['rname_2'],
                      rnum_2=var['rnum_2'],
                      ratoms_2=var['ratoms_2'],
                      rchain_2=var['rchain_2'],
                      constraints=var['constraints'])



def parse_constraints(raw:str) -> List[RosettaCst]:
    """

    """
    result: List[RosettaCst] = list()
    c_file = Path(raw)

    if uu.likely_a_file(raw):
        fs.check_file_exists(c_file)

    if c_file.exists():
        content: str = fs.content_from_file(c_file)
        content = ''.join(content.split())
        if not content:
            eh.core._LOGGER.warning(
                f"The supplied file '{raw}' contained no constraints. Continuing..."
            )
            return result
    else:
        content: str = raw 
        if not content:
            eh.core._LOGGER.warning(
                "The supplied argument for --constraints is empty. Continuing..."
            )
            return result

    for chunk in content.split('),('):
        if chunk[0] != '(':
            chunk = '(' + chunk
        if chunk[-1] != ')':
            chunk += ')'

        result.append(parse_rosetta_cst(chunk))

    return result



