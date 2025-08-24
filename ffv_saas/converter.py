"""
ffv.converter
Converts molecules between MDL Mol V3000 and Sybyl Mol2
"""

from openbabel import pybel

def convert(string=None, file_in=None, file_out=None):

    if string:
        print("using string")
        molecule = pybel.readstring(FORMAT, string)
    elif file_in:
        print("using file")
        molecule = pybel.readfile(FORMAT, file_in)
    else:
        raise RuntimeError("string or file_in must be specified")

    return molecule.write("mol2", filename=file_out)
