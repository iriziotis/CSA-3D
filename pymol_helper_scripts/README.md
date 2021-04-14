# PyMol CSA-3D helper scripts
### These scripts implement some useful PyMol functions to better visualize catalytic sites in PyMol

- `color_catres.py`: This is to color all catalytic residues consistently, according to their index\
e.g. All 1st catalytic residues of the homologous family are colored red, all 2nd green etc (or\
whatever colour palette the user provides).\
Usage:\
	`color_catres selection [,palette]`\
Args:\
	`selection`: Selection of objects to color\
	`palette`: Colour palette to use. Colours are equally space in the spectrum according to the\
		   number of catalytic residues. {default: 'rainbow'}\

- `select_catres.py`: This is to select all groups of homologous catalytic residues. Makes a selection\
for each catalytic residue, containing all residues from the homologous family that have the same index.\  
Usage:\
    `select_catres selection`\
Args:\
    `selection`: Selection of objects to make catalytic residue selection to
    