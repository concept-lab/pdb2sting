NTR_ATOM_TYPE = "H1"
CTR_ATOM_TYPE = "OXT"

RES_TABLE = {
    "ALA": ["AL0", "AL1"],
    "ARG": ["AR0", "AR1"],
    "ASN": ["AS4", "AS5"],
    "ASP": ["AS0", "AS1"],
    "ASH": ["AS2", "AS3"],
    "CYS": ["CY0", "CY1"],
    "CYX": ["CY2", "CY3"],
    "CYM": ["CY4", "CY5"],
    "GLY": ["GL0", "GL1"],
    "GLN": ["GL2", "GL3"],
    "GLU": ["GL4", "GL5"],
    "GLH": ["GL6", "GL7"],
    "HID": ["HI0", "HI1"],
    "HIE": ["HI2", "HI3"],
    "HIP": ["HI4", "HI5"],
    "ILE": ["IL0", "IL1"],
    "LEU": ["LE0", "LE1"],
    "LYS": ["LY0", "LY1"],
    "LYN": ["LY2", "LY3"],
    "MET": ["ME0", "ME1"],
    "PHE": ["PH0", "PH1"],
    "PRO": ["PR0", "PR1"],
    "SER": ["SE0", "SE1"],
    "THR": ["TH0", "TH1"],
    "TRP": ["TR0", "TR1"],
    "TYR": ["TY0", "TY1"],
    "VAL": ["VA0", "VA1"],
}

NUCLEIC_ACIDS = {
    "DTH": "DT",
    "DCY": "DC",
    "DAD": "DA",
    "DGU": "DG",
    "RUR": "RU",
    "RCY": "RC",
    "RAD": "RA",
    "RGU": "RG",
}

NUCLEIC3_ATOM_TYPE = "H3T"
NUCLEIC5_ATOM_TYPE = "H5T"

QUICK_CHANGE = {"DTH": ["C5M", "C7  "]}
QUICK_DEL = {"RAD": "H2'", "RGU": "H2'"}


def read_file(fin):
    with open(fin) as f:
        for line in f:
            yield line


def identify_termini(fin):
    termini_is = []
    termini_type = []
    cur_res = None
    res_is = []
    add_trigger = False
    for i, line in enumerate(read_file(fin)):
        if line.startswith("ATOM "):
            atom = line[12:16].strip()
            res = line[22:26]
            resname = line[17:20]

            if res != cur_res:
                if add_trigger:
                    termini_is += res_is
                    termini_type += [res_type for i in range(len(res_is))]
                cur_res = res
                res_is = [i]
                add_trigger = False
            else:
                res_is.append(i)

            if atom == NTR_ATOM_TYPE and resname not in NUCLEIC_ACIDS.keys():
                add_trigger = True
                res_type = "NTR"
            elif atom == CTR_ATOM_TYPE and resname not in NUCLEIC_ACIDS.keys():
                add_trigger = True
                res_type = "CTR"
            elif atom == NUCLEIC3_ATOM_TYPE:
                add_trigger = True
                res_type = "NU3"
            elif atom == NUCLEIC5_ATOM_TYPE:
                add_trigger = True
                res_type = "NU5"

    if add_trigger:
        termini_is += res_is
        termini_type += [res_type for i in range(len(res_is))]

    return termini_is, termini_type


def convert_pdb(fin, fout):
    termini_is, termini_type = identify_termini(fin)
    new_pdb = ""
    for i, line in enumerate(read_file(fin)):
        if line.startswith("ATOM "):
            res = line[17:20]
            atom = line[12:16].strip()

            if res in NUCLEIC_ACIDS:
                line = f"{line[:17]}{NUCLEIC_ACIDS[res]} {line[20:]}"

            if i in termini_is:
                ter_i = termini_is.index(i)
                if res in RES_TABLE:
                    ter_type = 0 if termini_type[ter_i] == "NTR" else 1
                    line = f"{line[:17]}{RES_TABLE[res][ter_type]}{line[20:]}"
                elif res in NUCLEIC_ACIDS and termini_type[ter_i] not in ("NTR", "CTR"):
                    ter_type = termini_type[ter_i][-1]
                    line = f"{line[:17]}{NUCLEIC_ACIDS[res]}{ter_type}{line[20:]}"

            if res in QUICK_CHANGE.keys() and atom == QUICK_CHANGE[res][0]:
                line = f"{line[:12]}{QUICK_CHANGE[res][1]}{line[16:]}"

            if res in QUICK_DEL.keys():

                if atom == QUICK_DEL[res]:
                    line = ""

        new_pdb += line
    with open(fout, "w") as f:
        f.write(new_pdb)


if __name__ == "__main__":
    fname = "prot_ph7.pdb"
    fout = "prot_ph7_delphi.pdb"
    convert_pdb(fname, fout)