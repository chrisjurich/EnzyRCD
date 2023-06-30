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
        """Sets an attribute specified by the 'key' to the given 'value'. No checks are performed
        on the supplied value.

        Args:
            key: The str() name of the attribute to change.
            value: What you want to set the attribute to.

        Returns:
            Nothing.

        Raises:
            KeyError() if the supplied attribute does not exist. 
        """
        if key not in self.__dict__:
            raise KeyError(f"{key} is not in the constraint")

        self.__dict__[key] = value


    def evaluate(self, file:str) -> List[float]:
        """ 

        Args:
            file: Name of the file to analyze.

        Returns:
            
        """
        eh.interface.pymol.execute([
            ('delete','all'),
            ('load', file)
        ])

        sele1, sele2, sele3, sele4 = None, None, None, None

        differences:List[float] = list() 

        for cst in self.constraints:
            #print(cst)
            #TODO(CJ): do a tuple expansion here
            cst_type:str = cst[0]
            target = cst[1]
            tolerance = cst[2]
            #TODO(CJ): need to check if the angle is weird for this

            if cst_type.startswith('distance'):
                sele1 = f"chain {self.rchain_1} and resi {self.rnum_1} and name {self.ratoms_1[0]}"
                sele2 = f"chain {self.rchain_2} and resi {self.rnum_2} and name {self.ratoms_2[0]}"
                dist:float = eh.interface.pymol.execute( [('distance', None, sele1, sele2)] )[0]
                differences.append( abs(dist - target ) / tolerance) 
            elif cst_type.startswith('angle'):
                if cst_type == 'angle_A':
                    sele1 = f"chain {self.rchain_1} and resi {self.rnum_1} and name {self.ratoms_1[1]}"
                    sele2 = f"chain {self.rchain_1} and resi {self.rnum_1} and name {self.ratoms_1[0]}"
                    sele3 = f"chain {self.rchain_2} and resi {self.rnum_2} and name {self.ratoms_2[0]}"
                elif cst_type == 'angle_B':
                    sele1 = f"chain {self.rchain_1} and resi {self.rnum_1} and name {self.ratoms_1[0]}"
                    sele2 = f"chain {self.rchain_2} and resi {self.rnum_2} and name {self.ratoms_2[0]}"
                    sele3 = f"chain {self.rchain_2} and resi {self.rnum_2} and name {self.ratoms_2[1]}"
                angle:float = eh.interface.pymol.execute([('angle', None, sele1, sele2, sele3)])[0]
                differences.append( abs(dist-target)/tolerance )

            elif cst_type.startswith('dihedral'):
                if cst_type == 'torsionA':
                    sele1 = f"chain {self.rchain_1} and resi {self.rnum_1} and name {self.ratoms_1[2]}"
                    sele2 = f"chain {self.rchain_1} and resi {self.rnum_1} and name {self.ratoms_1[1]}"
                    sele3 = f"chain {self.rchain_1} and resi {self.rnum_1} and name {self.ratoms_1[0]}"
                    sele4 = f"chain {self.rchain_2} and resi {self.rnum_2} and name {self.ratoms_2[0]}"
                elif cst_type == 'torsionAB':
                    sele1 = f"chain {self.rchain_1} and resi {self.rnum_1} and name {self.ratoms_1[1]}"
                    sele2 = f"chain {self.rchain_1} and resi {self.rnum_1} and name {self.ratoms_1[0]}"
                    sele3 = f"chain {self.rchain_2} and resi {self.rnum_2} and name {self.ratoms_2[0]}"
                    sele4 = f"chain {self.rchain_2} and resi {self.rnum_2} and name {self.ratoms_2[1]}"
                elif cst_type == 'torsionB':
                    sele1 = f"chain {self.rchain_1} and resi {self.rnum_1} and name {self.ratoms_1[0]}"
                    sele2 = f"chain {self.rchain_2} and resi {self.rnum_2} and name {self.ratoms_2[0]}"
                    sele3 = f"chain {self.rchain_2} and resi {self.rnum_2} and name {self.ratoms_2[1]}"
                    sele4 = f"chain {self.rchain_2} and resi {self.rnum_2} and name {self.ratoms_2[2]}"
                assert False
                args.append(
                    ('dihedral', sele1, sele2, sele3, sele4 )
                )

        return differences

    def contains(self, chain:str, res_num:int) -> bool:
        """Does the constraint contain the residue at <chain>.<res_num>?"""
        if chain == self.rchain_1 and res_num == self.rnum_1:
            return True

        if chain == self.rchain_2 and res_num == self.rnum_2:
            return True

        return False



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


def read_cst(raw: str) -> List[RosettaCst]:
    """Takes a raw constraint and converts it into a RosettaCst object. This is a modified version of the information
    required to create an enzyme constraint in Rosetta, with the exception that some parameters for distance and angle
    constraints can be ommitted. Each constraint is converted into a single constraint block in the Rosetta constraint
    system. Consider the below example:

    Residue 1: Chain A, ResNum 1, ResName L01, Atoms A1 A2 A3
    Residue 2: Chain B, ResNum 2, ResName L02, Atoms A4 A5 A6

    FormattedCst:

    (A,1,L01,A1,A2,A3)(B,2,L02,A4,A5,A6)(distanceAB,2.00)(angle_A,180.00)(angle_B)
   


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
            cst_type:str=spl[0]
            if cst_type not in ALLOWED_CSTS:
                eh.core._LOGGER.error(
                    f"The supplied constraint type {cst_type} is not supported. Allowed are: {', '.join(sorted(list(ALLOWED_CSTS)))}. Exiting..."
                )
                exit(1)
            
            temp = [cst_type]
            
            for tt in spl[1:]:
                if tt.find('.') == -1:
                    temp.append(int(tt))
                else:
                    temp.append(float(tt))
            
            t_len:int=len(temp)

            if t_len < 5:

                #TODO(CJ): add in some discussion of using a database here
                # and mention that we are using default parameters
                if cst_type == 'distanceAB':
                    temp.extend(
                        [2.00,0.25,1000.00,0][t_len-1:]
                    )
                elif cst_type in 'angle_A angle_B'.split():
                    temp.extend(
                        [180.0,5.0,1000.0,360.0,1][t_len-1:]
                    )

                else:
                    raise TypeError()

            var['constraints'].append(temp)
    
    return RosettaCst(**var)



def parse_constraints(raw:str) -> List[RosettaCst]:
    """Method that parses a raw str() into a list() of RosettaCst objects(). Allows for the input
    to be either an actual str() or the path of a file containing constraints. On errors or bad input,
    primarily returns an empty list(). 

    Args:
        raw: The constraints raw str() or path to a file containing constraints.

    Returns:
        A list() of validated RosettaCst().

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
                "The supplied constraints str is empty and is likely not a file. Continuing..."
            )
            return result

    for chunk in content.split('),('):
        if chunk[0] != '(':
            chunk = '(' + chunk
        if chunk[-1] != ')':
            chunk += ')'

        result.append(read_cst(chunk))

    return result



