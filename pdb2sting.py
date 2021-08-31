from pypka import Titration
from fixtermini import convert_pdb

input_file = "tests/4lzt.pdb"
output_file = "out.pdb"  # created by PypKa
delphi_file = "out_delphi.pdb"  # created by fixtermini
pH = "7.0"

params = {
    "structure": input_file,
    "ncpus": -1,  # -1 to use all available
    "pH": pH,
    "epsin": 15,
    "ionicstr": 0.1,
    "pbc_dimensions": 0,
    "ser_thr_titration": False,
    "save_pdb": "delphi_in.pdb",
    "structure_output": (output_file, pH, "AMBER"),
}

if __name__ == "__main__":
    # Run PypKa to assign protonation states
    tit = Titration(params)

    # Run fixtermini to fix nomenclature
    convert_pdb(output_file, delphi_file)
