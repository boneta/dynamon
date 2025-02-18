#!/usr/bin/env python3

"""
  PyMOL plug-in with DYNAMON functionalities
  ==========================================

  Add extended capabilities to the PyMOL molecular viewer
  to manage files in fDynamo/DYNAMON formats

  - Write a selection of atoms to a file:
        write_sele filename [, selection_name [, selection [, resolution ]]]

  - Write directly QM atoms or NOFIX residues from selection:
        write_qm  filename [, selection ]
        write_nofix  filename [, selection ]

  - Load/read extended file formats:  load  filename [...]
        .ff    -  fDynamo's force field file
                  re-bond ligands accordingly and display properties:
                  atomtypes as 'text_type' and charges as 'partial_charges'
        .crd   -  fDynamo's coordinates file
        .dynn  -  DYNAMON options/selection file


  by Sergio Boneta Martinez
  GPLv3 2025

"""

__version__ = '0.1.0'


##  DEPENDENCIES  #####################################################

import os
from copy import deepcopy
from glob import glob
from textwrap import dedent

from pymol import cmd, importing, plugins, CmdException
from chempy import atomic_number, Atom
from chempy.models import Indexed
from chempy.protein_residues import normal as aa_dict


##  DYNAMON CLASS  ####################################################

class DynnConfigSele():
    """
    DESCRIPTION

        DYNAMON configuration class - Selection oriented

    ATTRIBUTES

        selection = dict
            'sele_name' = dict
                'segi' = dict
                    'resi' = list
    """

    def __init__(self):
        self.opt_raw = []
        self.selection = dict()

    def read_selection(self, filename):
        """
        DESCRIPTION

            Read selection blocks from .dynn file

        ARGUMENTS

            filename = string: file to read selection
        """

        # check file exists
        if not os.path.isfile(filename): return

        # read file as list of strings
        with open(filename, 'r') as f:
            dynn_file = f.readlines()
            # remove empty lines, comments and space-split
            dynn_file = [line.split() for line in dynn_file if line.strip() and not line.startswith(("!","#"))]

        n = 0   # line number
        while n < len(dynn_file):
            if dynn_file[n][0].upper() == 'SELECTION':
                sele_name = dynn_file[n][1].upper()
                n += 1
                # read whole selection to dict of dict of list
                sele = dict()
                segi = ""
                resi = ""
                while dynn_file[n][0].upper() != 'SELECTION' and n < len(dynn_file):
                    subsect = dynn_file[n][0].upper()
                    select  = dynn_file[n][1].upper()
                    if subsect == "S":
                        segi = select
                        resi = ""
                        sele[select] = dict()
                    elif subsect == "R":
                        resi = int(select)
                        sele[segi][int(select)] = []
                    elif subsect == "A":
                        sele[segi][resi].append(select)
                    n += 1
                else:
                    n += 1
                self.selection[sele_name] = sele
            else:
                self.opt_raw.append("  ".join(dynn_file[n]))
                n += 1

    def write_selection(self, filename, resolution='atom'):
        """
        DESCRIPTION

            Write selection blocks to .dynn file

        ARGUMENTS

            filename = string: file to create/append and write selection

            resolution = 'atom'/'residue'/'subsystem' : minimum entity size to treat no-whole when writting {default: 'atom'}
        """

        if not self.selection:
            raise CmdException("No selection to write", "pyDYNAMON")

        # build strings for each selection
        sele_str = []
        for sele_name, sele in self.selection.items():
            s = f"\nSELECTION {sele_name}\n"
            for segi, resis in sele.items():
                s += " "*4+f"S {segi}\n"
                if resolution == 'subsystem': continue
                for resi, atoms in resis.items():
                    s += " "*8+f"R {resi}\n"
                    if resolution == 'residue': continue
                    for name in atoms:
                        s += " "*12+f"A {name}\n"
            s += "SELECTION\n"
            sele_str.append(s)

        # write selections to file
        with open(filename, 'w') as f:
            f.write("\n"+"\n".join(self.opt_raw)+"\n")
            f.write("".join(sele_str))


##  FUNCTIONS  ########################################################

def write_sele(filename, selection_name='', selection='sele', resolution='atom'):
    """
    DESCRIPTION

        Append selection to file and overwrite section if existing

    USAGE

        write_sele filename [, selection_name [, selection [, resolution ]]]

    ARGUMENTS

        filename = string: file path to create/append and write selection

        selection_name = string: name of selection to write (i.e. 'QM'/'NOFIX') {default: taken from selection argument}

        selection = string: name of a PyMOL selection object {default: 'sele'}

        resolution = 'atom'/'residue'/'subsystem' : minimum entity size to treat no-whole when writting {default: 'atom'}

    EXAMPLES

        write_sele protein.dynn
        write_sele protein.dynn, QM, QM

    SEE ALSO

        write_qm, write_nofix
    """

    selection_name = selection_name.upper() or selection.upper()

    # check selection exists
    if not selection in cmd.get_names('selections'):
        raise CmdException(f"Selection '{selection}' not found", "pyDYNAMON")

    # get selection with DynnConfigSele structure
    natoms = 0
    sele = dict()
    obj_list = cmd.get_object_list(selection)
    for obj in obj_list:
        model = cmd.get_model(selection+" and "+obj)
        natoms += model.nAtom
        for a in model.atom:
            segi = str(a.segi)
            resi = int(a.resi)
            name = str(a.name)
            sele.setdefault(segi, {})
            if resolution == 'subsystem': continue
            sele[segi].setdefault(resi, [])
            if resolution == 'residue': continue
            sele[segi][resi].append(name)

    # read file to overwrite a section if already exists
    dynn = DynnConfigSele()
    dynn.read_selection(filename)

    # assign to object and write
    dynn.selection[selection_name] = sele
    dynn.write_selection(filename, resolution='atom')

    print(f" pyDYNAMON: {selection_name} with {natoms} written to \"{os.path.abspath(filename)}\"")

def load_dynn(filename):
    """
    DESCRIPTION

        Load a DYNAMON selection file and read QM/NOFIX to sele objects

    ARGUMENTS

        filename = string: file path
    """

    # check file exists
    if not os.path.isfile(filename):
        raise CmdException(f"File '{filename}' not found.", "pyDYNAMON")

    print(f" pyDYNAMON: reading \"{filename}\"")
    dynn = DynnConfigSele()
    dynn.read_selection(filename)

    if not dynn.selection:
        raise CmdException("No selections to read in file. Aborting.", "pyDYNAMON")

    # build selection algebra
    for sele_name, sele in dynn.selection.items():
        sele1_list = []
        for segi, resdict in sele.items():
            sele1 = f"segi {segi}"
            if resdict:
                sele1 += " & ("
                sele2_list = []
                for resi, namelist in resdict.items():
                    sele2 = f"resi {resi}"
                    if namelist:
                        sele2 += " & name "+"+".join(namelist)
                    sele2_list.append(sele2)
                sele1 += " | ".join(sele2_list)+" )"
            sele1_list.append(sele1)
        sele_final = " | ".join(sele1_list)

        # selection command
        cmd.select(sele_name, sele_final, enable=0, quiet=1)
        natoms = cmd.count_atoms(sele_name)
        print(f" pyDYNAMON: selection \"{sele_name}\" defined with {natoms} atoms.")

def load_ff(filename):
    """
    DESCRIPTION

        Load a fDynamo force field file (.ff)

        Re-bond ligands accordingly
        Alter properties of all atoms:
         - 'text_type' : atomtypes
         - 'partial_charges' : charges

    ARGUMENTS

        filename = string: file path
    """

    print(f" pyDYNAMON: reading \"{filename}\"")

    # residue names excluded from re-bonding
    rebond_excluded = []
    rebond_excluded.extend(aa_dict.keys())
    rebond_excluded.extend(["SOL", "NA", "CL"])

    # read file as list of strings
    with open(filename, "rt") as f:
        ff_file = f.readlines()
        # remove comment lines and blank lines
        ff_file = [line.strip() for line in ff_file if line.strip() and not line.startswith("!")]

    # read residue atoms and bonds
    topology = dict()
    current_sect = None
    for line in ff_file:
        line = line.split("!")
        keyword = line[0].split()[0].lower()
        content = line[0].split()
        if keyword == 'end':
            current_sect = None
        elif keyword == 'residues':
            current_sect = keyword
        elif current_sect == 'residues' or keyword == 'residue':
            if keyword == 'residue':  # new residue
                current_sect = 'residues'
                resname = content[1]
                topology[resname] = dict()
                continue
            top = topology[resname]
            if not topology[resname]: # first line (numbers)
                top['natoms'] = int(content[0])
                top['nbonds'] = int(content[1])
                top['nimpropers'] = int(content[2])
                top['atoms'] = dict()
                top['bonds'] = []
                top['impropers'] = []
            elif len(top['atoms']) < top['natoms']:
                param = { 'atomtype' : str(content[1]).upper(),
                          'charge' : float(content[2]) }
                top['atoms'][str(content[0]).upper()] = param
            elif len(top['bonds']) < top['nbonds']:
                bonds = [(b.split()[0].upper(), b.split()[1].upper())
                          for b in " ".join(content).split(";") if b.strip()]
                top['bonds'].extend(bonds)
            elif len(top['impropers']) < top['nimpropers']:
                impropers = [(i.split()[0].upper(), i.split()[1].upper(),
                              i.split()[2].upper(), i.split()[3].upper())
                              for i in " ".join(content).split(";") if i.strip()]
                top['impropers'].extend(impropers)

    for resname, top in topology.items():
        # atomtypes and partial charges
        for name, param in top['atoms'].items():
            cmd.alter(f"resn {resname} & name {name}", f"partial_charge={param['charge']}")
            cmd.alter(f"resn {resname} & name {name}", f"text_type=\'{param['atomtype']}\'")
        # unbond all
        if resname in rebond_excluded: continue
        cmd.unbond(f"resn {resname}", f"resn {resname}")
        # re-bond by individual residues
        resn_model = cmd.get_model(f"resn {resname}")
        for atom_at_residues in resn_model.get_residues():
            ndx_first = resn_model.atom[atom_at_residues[0]].index
            ndx_last = resn_model.atom[atom_at_residues[1]-1].index
            for bond in top['bonds']:
                cmd.bond(f"resn {resname} & name {bond[0]} & index {ndx_first}-{ndx_last}",
                         f"resn {resname} & name {bond[1]} & index {ndx_first}-{ndx_last}")

def load_crd(filename, object=''):
    """
    DESCRIPTION

        Load a fDynamo coordinates file (.crd)

    ARGUMENTS

        filename = string: file path

        object = string: name of the new object {default: file basename}
    """

    if not object:
        object = "".join(os.path.basename(filename).rpartition('.')[:-2])

    # read file as list of strings
    with open(filename, "rt") as f:
        crd_file = f.readlines()
        # remove comment lines, trailing comments and space split
        crd_file = [line.split("!")[0].split() for line in crd_file
                    if line.strip() and not line.startswith("!")]

    # reversed dictionary of atomic numbers to elements
    atomic_number_inv = { n:elem for elem, n in atomic_number.items() }

    # create new model and append crd's atoms
    model = Indexed()
    a = Atom()
    a.hetatm  = False
    for line in crd_file:
        if line[0].lower() == "subsystem":
            a.segi = str(line[2])
        elif line[0].lower() == "residue":
            a.resi_number  = int(line[1])
            a.resn = str(line[2])
        elif len(line) != 6:
            continue
        else:
            a.name   = str(line[1])
            a.symbol = atomic_number_inv[int(line[2])]
            a.coord  = ( float(line[3]), float(line[4]), float(line[5]) )
            model.add_atom(deepcopy(a))
    model.update_index()
    cmd.load_model(model, object)
    cmd.rebond(object)
    cmd.dss(object)

    print(f" pyDYNAMON: \"{filename}\" loaded as \"{object}\"")

def load_ext(filename, object='', state=0, format='', finish=1,
             discrete=-1, quiet=1, multiplex=None, zoom=-1, partial=0,
             mimic=1, object_props=None, atom_props=None, *, _self=cmd):
        """
        REMARK

            This is a wrapper to load function with extended funtionality
            for fDynamo's files.

            .ff
                fDynamo's force field file
            .crd
                fDynamo's coordinates file
            .dynn
                DYNAMON options/selection file
        """

        if not format:
            format = os.path.basename(filename).rpartition('.')[-1]

        # fDynamo's force field
        if format.lower() == "ff":
            load_ff(filename)

        # DYNAMON dynn options/selection
        elif format.lower() == "dynn":
            load_dynn(filename)

        # fDynamo's crd coordinates
        elif format.lower() == "crd":
            load_crd(filename, object)

        # original load function
        else:
            importing.load(filename, object, state, format, finish,
                           discrete, quiet, multiplex, zoom, partial)

load_ext.__doc__ = importing.load.__doc__ + dedent(load_ext.__doc__)


##  WRAPPERS  #########################################################

def write_qm(filename, selection="sele"):
    """
    DESCRIPTION

        Write QM atom selection to a file

    USAGE

        write_qm filename [, selection]

    ARGUMENTS

        filename = string: file path to create/append and write selection

        selection = string: name of a PyMOL selection object {default: 'sele'}

    EXAMPLES

        write_qm protein.dynn
        write_qm protein.dynn, QM

    SEE ALSO

        write_sele, write_nofix
    """
    write_sele(filename, "QM", selection, resolution='atom')

def write_nofix(filename, selection="sele"):
    """
    DESCRIPTION

        Write NOFIX residue selection to a file

    USAGE

        write_nofix filename [, selection]

    ARGUMENTS

        filename = string: file path to create/append and write selection

        selection = string: name of a PyMOL selection object {default: 'sele'}

    EXAMPLES

        write_nofix protein.dynn
        write_nofix protein.dynn, NOFIX

    SEE ALSO

        write_sele, write_qm
    """
    write_sele(filename, "NOFIX", selection, resolution='residue')


##  PYMOL FUNCTIONS  ##################################################
# auto-completion
cmd.auto_arg[0]['write_sele'] = [lambda: cmd.Shortcut(glob('*.dynn')), 'filename', ', ']
cmd.auto_arg[1]['write_sele'] = [cmd.selection_sc, 'selection_name', ', ']
cmd.auto_arg[2]['write_sele'] = [cmd.selection_sc, 'selection', ', ']
cmd.auto_arg[3]['write_sele'] = [lambda: cmd.Shortcut(['atom', 'residue', 'subsystem']), 'resolution', '']
cmd.auto_arg[0]['write_qm'] = [lambda: cmd.Shortcut(glob('*.dynn')), 'filename', ', ']
cmd.auto_arg[1]['write_qm'] = [cmd.selection_sc, 'selection', '']
cmd.auto_arg[0]['write_nofix'] = [lambda: cmd.Shortcut(glob('*.dynn')), 'filename', ', ']
cmd.auto_arg[1]['write_nofix'] = [cmd.selection_sc, 'selection', '']
# add functions to PyMOL run environment
cmd.extend("write_sele", write_sele)
cmd.extend("write_qm", write_qm)
cmd.extend("write_nofix", write_nofix)
cmd.load = load_ext
cmd.extend("load", cmd.load)


##  GUI  ##############################################################
def __init_plugin__(app=None):
    """Add an entry to the PyMOL 'Plugin' menu"""
    plugins.addmenuitemqt("pyDYNAMON", lambda: print(f"  pyDYNAMON: {__version__}"))
