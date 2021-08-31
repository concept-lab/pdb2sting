from delphi4py.delphi4py import DelPhi4py

delphi_file = "prot_ph7_delphi.pdb"  # created by fixtermini

delphimol = DelPhi4py(
    "parameters/amber.crg",
    "parameters/amber.siz",
    delphi_file,
    121,
    2,
    "single",
    conc=0.1,
    ibctyp=4,
    res2=0.01,
    nlit=500,
    pbx=False,
    pby=False,
    outputfile="LOG_",
    perfil=80,
)

print("DelPhi Read Completed")
print(delphimol.igrid)

delphimol.runDelPhi(
    nonit=0,
    nlit=500,
    outputfile="LOG_run",
)
print("DelPhi Run Completed")
print("corrected reaction field energy:", delphimol.getSolvation())  # float
print("returned potential size:", len(delphimol.getSitePotential()))  # array