"""
python-babel - A pythonic interface to openbabel SWIG bindings

Global variables:
  ob - the underlying SWIG bindings for Open Babel
  informats - a dictionary of supported input formats
  outformats - a dictionary of supported output formats
  fps - a list of supported fingerprint types
  forcefields - a list of supported forcefields
"""

import sys
import os.path
import xml.etree.ElementTree as ET
import openbabel as ob
import numpy as np
import pubchempy as pubchem

_obfuncs = _obconsts = ob
_obconv = ob.OBConversion()
_elements = ob.OBElementTable()


input_formats = [format.split('--')[0] for format in
                 _obconv.GetSupportedInputFormat()]
output_formats = [format.split('--')[0] for format in
                  _obconv.GetSupportedOutputFormat()]


def readfile(format, filename, options=None):
    """ Read a molecule from a file.

        format: a string specifying the format to interpret the
                file. e.g. "smi", "mol", "xyz" etc
        filename: the path to the file
        options: A dict object of the openbabel options analogous
                 to the obabel command line options, with key being
                 the option and the value (if it exists)
                 e.g. -xa would be {'a': None}
    """
    if not options:
        options = {}

    obconversion = ob.OBConversion()
    supported_format = obconversion.SetInFormat(format)

    if not supported_format:
        raise ValueError("{:s} is not a recognised Open \
                         Babel format".format(format))

    if not os.path.isfile(filename):
        raise IOError("No such file: '{:s}'".format(filename))

    for key, val in options.items():
        if val is None:
            obconversion.AddOption(key, obconversion.INOPTIONS)
        else:
            obconversion.AddOption(key, obconversion.INOPTIONS, str(val))

    mol = ob.OBMol()

    success = obconversion.ReadFile(mol, filename)
    if not success:
        raise IOError("Failed to read {:s}. Check the \
                       file for errors.".format(filename))

    return Molecule(mol)


def readstring(format, string, options=None):
    """ Read a molecule from a string.

        format: a string specifying the format to interpret the
                string. e.g. "smi", "mol", "xyz" etc
        string: the string specifying the molecule
        options: A dict object of the openbabel options analogous
                 to the obabel command line options, with key being
                 the option and the value (if it exists)
                 e.g. -xa would be {'a': None}
    """
    if not options:
        options = {}

    mol = ob.OBMol()
    obconversion = ob.OBConversion()

    supported_format = obconversion.SetInFormat(format)
    if not supported_format:
        raise ValueError("{:s} is not a recognised Open \
                          Babel format".format(format))
    for key, val in options.items():
        if val is None:
            obconversion.AddOption(key, obconversion.INOPTIONS)
        else:
            obconversion.AddOption(key, obconversion.INOPTIONS, str(val))

    success = obconversion.ReadString(mol, string)
    if not success:
        raise IOError("Failed to convert '{}' to format '{}'. Check \
                       the string for errors".format(string, format))

    return Molecule(mol)


class Molecule(object):
    def __init__(self, OBMol):
        self.OBMol = OBMol

    @property
    def atoms(self):
        return [Atom(a) for a in ob.OBMolAtomIter(self.OBMol)]

    @property
    def bonds(self):
        return [Bond(b) for b in ob.OBMolBondIter(self.OBMol)]

    @property
    def charge(self):
        return self.OBMol.GetTotalCharge()

    @property
    def conformers(self):
        return self.OBMol.GetConformers()

    @property
    def data(self):
        return MoleculeData(self.OBMol)

    @property
    def dim(self):
        return self.OBMol.GetDimension()

    @property
    def energy(self):
        return self.OBMol.GetEnergy()

    @property
    def exactmass(self):
        return self.OBMol.GetExactMass()

    @property
    def formula(self):
        return self.OBMol.GetFormula()

    @property
    def molwt(self):
        return self.OBMol.GetMolWt()

    @property
    def spin(self):
        return self.OBMol.GetTotalSpinMultiplicity()

    @property
    def sssr(self):
        return self.OBMol.GetSSSR()

    def _gettitle(self):
        return self.OBMol.GetTitle()

    def _settitle(self, val):
        self.OBMol.SetTitle(val)
    title = property(_gettitle, _settitle)

    @property
    def clone(self):
        return Molecule(ob.OBMol(self.OBMol))

    @property
    def _exchange(self):
        if self.OBMol.HasNonZeroCoords():
            return (1, self.write("mol"))
        else:
            return (0, self.write("can").split()[0])

    def add_bond(self, begin_id, end_id, order):
        self.OBMol.AddBond(begin_id, end_id, order)

    def __iter__(self):
        """Iterate over the Atoms of the Molecule.

        This allows constructions such as the following:
           for atom in mymol:
               print atom
        """
        return iter(self.atoms)

    def _repr_svg_(self):
        """For IPython notebook, renders 2D pybabel.Molecule SVGs."""

        # Open babel returns a nested svg, which IPython unpacks and treats as
        # two SVGs, messing with the display location. This parses out the
        # inner svg before handing over to IPython.
        namespace = "http://www.w3.org/2000/svg"
        ET.register_namespace("", namespace)
        obsvg = self.clone.write("svg", opt={'a': None})
        tree = ET.fromstring(obsvg)
        svg = tree.find("{{{ns}}}g/{{{ns}}}svg".format(ns=namespace))
        return ET.tostring(svg).decode("utf-8")

    def _repr_html_(self):
        """For IPython notebook, renders 3D pybabel.Molecule webGL objects."""
        return None

    def write(self, format="smi", filename=None, overwrite=False, opt=None):
        """Write the molecule to a file or return a string.

        Optional parameters:
           format -- see the informats variable for a list of available
                     output formats (default is "smi")
           filename -- default is None
           overwite -- if the output file already exists, should it
                       be overwritten? (default is False)
           opt -- a dictionary of format specific options
                  For format options with no parameters, specify the
                  value as None.

        If a filename is specified, the result is written to a file.
        Otherwise, a string is returned containing the result.

        To write multiple molecules to the same file you should use
        the Outputfile class.
        """
        if opt is None:
            opt = {}
        obconversion = ob.OBConversion()
        formatok = obconversion.SetOutFormat(format)
        if not formatok:
            raise ValueError("%s is not a recognised Open Babel format" %
                             format)
        if filename and filename.split('.')[-1] == 'gz':
            obconversion.AddOption('z', self.obConversion.GENOPTIONS)
        for k, v in opt.items():
            if v is None:
                obconversion.AddOption(k, obconversion.OUTOPTIONS)
            else:
                obconversion.AddOption(k, obconversion.OUTOPTIONS, str(v))

        if filename:
            if not overwrite and os.path.isfile(filename):
                raise IOError(("%s already exists. Use 'overwrite=True' to "
                               "overwrite it.") % filename)
            obconversion.WriteFile(self.OBMol, filename)
            obconversion.CloseOutFile()
        else:
            return obconversion.WriteString(self.OBMol)

    def addh(self):
        """Add hydrogens."""
        self.OBMol.AddHydrogens()

    def removeh(self):
        """Remove hydrogens."""
        self.OBMol.DeleteHydrogens()

    def convertdbonds(self):
        """Convert Dative Bonds."""
        self.OBMol.ConvertDativeBonds()

    def find_symmetry_unique_bonds(self):
        unique_bonds = []
        bonds = self.bonds[:]
        while len(bonds) > 0:
            a = bonds.pop()
            if not any(a.symmetry_equivalent(b) for b in unique_bonds):
                unique_bonds.append(a)
        return unique_bonds

    def __str__(self):
        return self.write()


class Atom(object):
    def __init__(self, OBAtom):
        self.OBAtom = OBAtom

    @property
    def coords(self):
        return (self.OBAtom.GetX(), self.OBAtom.GetY(), self.OBAtom.GetZ())

    @property
    def atomicmass(self):
        return self.OBAtom.GetAtomicMass()

    @property
    def atomicnum(self):
        return self.OBAtom.GetAtomicNum()

    @property
    def cidx(self):
        return self.OBAtom.GetCIdx()

    @property
    def coordidx(self):
        return self.OBAtom.GetCoordinateIdx()

    @property
    def exactmass(self):
        return self.OBAtom.GetExactMass()

    @property
    def formalcharge(self):
        return self.OBAtom.GetFormalCharge()

    @property
    def heavyvalence(self):
        return self.OBAtom.GetHvyValence()

    @property
    def heterovalence(self):
        return self.OBAtom.GetHeteroValence()

    @property
    def hyb(self):
        return self.OBAtom.GetHyb()

    @property
    def idx(self):
        return self.OBAtom.GetIdx()

    @property
    def implicitvalence(self):
        return self.OBAtom.GetImplicitValence()

    @property
    def isotope(self):
        return self.OBAtom.GetIsotope()

    @property
    def partialcharge(self):
        return self.OBAtom.GetPartialCharge()

    @property
    def spin(self):
        return self.OBAtom.GetSpinMultiplicity()

    @property
    def type(self):
        return self.OBAtom.GetType()

    @property
    def valence(self):
        return self.OBAtom.GetValence()

    @property
    def vector(self):
        return self.OBAtom.GetVector()

    @property
    def symbol(self):
        return _elements.GetSymbol(self.atomicnum)

    def distance_to(self, other):
        return np.linalg.norm(self.coords - other.coords)

    def distance_to_origin(self):
        return np.linalg.norm(self.coords)

    def __str__(self):
        c = self.coords
        return "Atom: {0}({4}) ({1:.2f} {2:.2f} {3:.2f})"\
                .format(self.symbol, c[0], c[1], c[2], self.idx)
    def __repr__(self):
        return self.__str__()


class Smarts(object):
    """A Smarts Pattern Matcher

    Required parameters:
       smartspattern

    Methods:
       findall(molecule)

    Example:
    >>> mol = readstring("smi","CCN(CC)CC") # triethylamine
    >>> smarts = Smarts("[#6][#6]") # Matches an ethyl group
    >>> print smarts.findall(mol)
    [(1, 2), (4, 5), (6, 7)]

    The numbers returned are the indices (starting from 1) of the atoms
    that match the SMARTS pattern. In this case, there are three matches
    for each of the three ethyl groups in the molecule.
    """

    def __init__(self, smartspattern):
        """Initialise with a SMARTS pattern."""
        self.obsmarts = ob.OBSmartsPattern()
        success = self.obsmarts.Init(smartspattern)
        if not success:
            raise IOError("Invalid SMARTS pattern")

    def findall(self, molecule):
        """Find all matches of the SMARTS pattern to a particular molecule.

        Required parameters:
           molecule
        """
        self.obsmarts.Match(molecule.OBMol)
        vector = self.obsmarts.GetUMapList()
        return list(vector)


class MoleculeData(object):
    """Store molecule data in a dictionary-type object

    Required parameters:
      obmol -- an Open Babel OBMol

    Methods and accessor methods are like those of a dictionary except
    that the data is retrieved on-the-fly from the underlying OBMol.

    Example:
    >>> mol = readfile("sdf", 'head.sdf').next() # Python 2
    >>> # mol = next(readfile("sdf", 'head.sdf')) # Python 3
    >>> data = mol.data
    >>> print data
    {'Comment': 'CORINA 2.61 0041  25.10.2001', 'NSC': '1'}
    >>> print len(data), data.keys(), data.has_key("NSC")
    2 ['Comment', 'NSC'] True
    >>> print data['Comment']
    CORINA 2.61 0041  25.10.2001
    >>> data['Comment'] = 'This is a new comment'
    >>> for k,v in data.items():
    ...    print k, "-->", v
    Comment --> This is a new comment
    NSC --> 1
    >>> del data['NSC']
    >>> print len(data), data.keys(), data.has_key("NSC")
    1 ['Comment'] False
    """

    def __init__(self, obmol):
        self._mol = obmol

    def _data(self):
        data = self._mol.GetData()
        answer = [x for x in data if
                  x.GetDataType() == _obconsts.PairData or
                  x.GetDataType() == _obconsts.CommentData]
        return answer

    def _testforkey(self, key):
        if key not in self:
            raise KeyError("'%s'" % key)

    def keys(self):
        return [x.GetAttribute() for x in self._data()]

    def values(self):
        return [x.GetValue() for x in self._data()]

    def items(self):
        return iter(zip(self.keys(), self.values()))

    def __iter__(self):
        return iter(self.keys())

    def iteritems(self):  # Can remove for Python 3
        return self.items()

    def __len__(self):
        return len(self._data())

    def __contains__(self, key):
        return self._mol.HasData(key)

    def __delitem__(self, key):
        self._testforkey(key)
        self._mol.DeleteData(self._mol.GetData(key))

    def clear(self):
        for key in self:
            del self[key]

    def has_key(self, key):
        return key in self

    def update(self, dictionary):
        for k, v in dictionary.items():
            self[k] = v

    def __getitem__(self, key):
        self._testforkey(key)
        answer = self._mol.GetData(key)
        return answer.GetValue()

    def __setitem__(self, key, value):
        if key in self:
            pairdata = self._mol.GetData(key).Downcast[ob.OBPairData]()
            pairdata.SetValue(str(value))
        else:
            pairdata = ob.OBPairData()
            pairdata.SetAttribute(key)
            pairdata.SetValue(str(value))
            self._mol.CloneData(pairdata)

    def __repr__(self):
        return dict(self.items()).__repr__()


class Bond(object):

    def __init__(self, obbond):
        self.obbond = obbond

    @property
    def order(self):
        return self.obbond.GetBondOrder()

    @property
    def order_char(self):
        if self.order == 1:
            return '-'
        elif self.order == 2:
            return '='
        elif self.order == 3:
            return '#'
        else:
            return '?'

    @property
    def begin_atom(self):
        return Atom(self.obbond.GetBeginAtom())

    @property
    def end_atom(self):
        return Atom(self.obbond.GetEndAtom())

    @property
    def atoms(self):
        return sorted([self.begin_atom, self.end_atom], key=lambda x: x.symbol)

    @property
    def length(self):
        return self.obbond.GetLength()

    @property
    def aromatic(self):
        return self.obbond.IsAromatic()

    @property
    def in_ring(self):
        return self.obbond.IsInRing()

    @property
    def equilibrium_length(self):
        return self.obbond.GetEquibLength()

    @property
    def id(self):
        return self.obbond.GetIdx()

    @property
    def smiles(self):
        return self.write('smi').split('\t')[0]

    def same_kind(self, other):
        return (set(a.symbol for a in self.atoms) ==
                set(a.symbol for a in other.atoms))

    def same_equilibrium_length(self, other):
        return (self.equilibrium_length == other.equilibrium_length)

    def same_id(self, other):
        return self.id == other.id

    def same_atoms(self, other):
        return (set(a.idx for a in self.atoms) ==
                set(a.idx for a in other.atoms))

    def is_equivalent_to(self, other):
        return (self.same_equilibrium_length(other) and
                self.same_kind(other) and
                self.same_id(other))

    def symmetry_equivalent(self, other):
        return (self.same_kind(other) and
                self.same_length(other))

    def same_length(self, other, eps=0.00005):
        return abs(self.length - other.length) < eps


    def __str__(self):
        return "{0}{1}{2}".format(self.atoms[0].symbol,
                                  self.order_char,
                                  self.atoms[1].symbol)

    def __repr__(self):
        return self.__str__() + ": {}".format(self.length)

    def __sub__(self, other):
        return self.length - other.length

    def __format(self):
        return self.__str__()
