def check_residue_id(res_id: str) -> None:
    """Helper function that checks if a supplied string is a valid residue id in the PDB format. Supplied
    code must: 1) be three characters in length and 2) be all uppercase. Strips whitespace.
    Args:
        res_id: The residue id to check as a str().
    Returns:
        Nothing.
    """
    res_id = ''.join(res_id.split())

    if len(res_id) != 3:
        eh.core._LOGGER.error(
            f"The supplied residue id '{res_id}' is invalid. Must be 3 characters long. Exiting..."
        )
        exit(1)

    if not res_id.isupper():
        eh.core._LOGGER.error(
            f"The supplied residue id '{res_id}' is invalid. Must be uppercase. Exiting..."
        )
        exit(1)



def redirect_stdout():
    #print "Redirecting stdout"
    sys.stdout.flush()  # <--- important when redirecting to files
    newstdout = os.dup(1)
    devnull = os.open(os.devnull, os.O_WRONLY)
    os.dup2(devnull, 1)
    os.close(devnull)
    sys.stdout = os.fdopen(newstdout, 'w')



def likely_a_file(raw: str) -> bool:
    """Applies heuristics to guess if the supplied str() is a filepath or not. Checks if:
        + has './' 
        + has '//' 
        + has a suffix as defined by pathlib

    Args:
        raw: The candidate str() to check.
    
    Returns:
        Whether the supplied str() is likely a file based on heuristics.
    """
    if raw.find('./') != -1:
        return True

    if raw.find('//') != -1:
        return True

    temp = Path(raw)
    if len(temp.suffix):
        return True

    return False



