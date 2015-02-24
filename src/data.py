"""
    A module to store any global information
    such as the Van Der Waal's radii for various
    elements etc.
"""
import logging
import progressbar as pb
from timeit import default_timer

FORMAT = "%(levelname)s%(message)s"
logging.basicConfig(level=logging.INFO, format=FORMAT)
logging.addLevelName(logging.ERROR, 'error: ')
logging.addLevelName(logging.WARNING, 'warning: ')
logging.addLevelName(logging.INFO, '')
logging.addLevelName(logging.DEBUG, 'debug: ')
logging.addLevelName(logging.CRITICAL, 'CRITICAL: ')

logger = logging.getLogger("sarlacc")


class Timer(object):
    """ A context manager timer class, to measure wall clock time"""
    def __init__(self):
        self.timer = default_timer

    def __enter__(self):
        self.start = self.timer()
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.end = self.timer()
        self.elapsed_s = self.elapsed()
        self.elapsed_s = self.elapsed_s * 1000.0

    def elapsed(self):
        return self.timer() - self.start


def log_traceback(e):
    logger.exception(e)

# Dict containing Van Der Waal's radii in angstroms
vdw_radii = {
    'Ag': 1.72,
    'Ar': 1.88,
    'As': 1.85,
    'Au': 1.66,
    'Br': 1.85,
    'C': 1.70,
    'Cd': 1.58,
    'Cl': 1.75,
    'Cu': 1.40,
    'D': 1.20,  # Deuterium added to be treated as Hydrogen in interactions
    'F': 1.47,
    'Ga': 1.87,
    'H': 1.20,
    'He': 1.40,
    'Hg': 1.55,
    'I': 1.98,
    'In': 1.93,
    'K': 2.75,
    'Kr': 2.02,
    'Li': 1.82,
    'Mg': 1.73,
    'N': 1.55,
    'Na': 2.27,
    'Ne': 1.54,
    'Ni': 1.63,
    'O': 1.52,
    'P': 1.80,
    'Pb': 2.02,
    'Pd': 1.63,
    'Pt': 1.72,
    'S': 1.80,
    'Se': 1.90,
    'Si': 2.10,
    'Sn': 2.17,
    'Te': 2.06,
    'Tl': 1.96,
    'U': 1.86,
    'Xe': 2.16,
    'Zn': 1.39
}

elements = {'Ac', 'Ag', 'Al', 'Am', 'Ar', 'As', 'At', 'Au', 'B', 'Ba', 'Be',
            'Bh', 'Bi', 'Bk', 'Br', 'C', 'Ca', 'Cd', 'Ce', 'Cf', 'Cl', 'Cm',
            'Cn', 'Co', 'Cr', 'Cs', 'Cu', 'Db', 'Ds', 'Dy', 'Er', 'Es', 'Eu',
            'F', 'Fe', 'Fl', 'Fm', 'Fr', 'Ga', 'Gd', 'Ge', 'H', 'He', 'Hf',
            'Hg', 'Ho', 'Hs', 'I', 'In', 'Ir', 'K', 'Kr', 'La', 'Li', 'Lr',
            'Lu', 'Lv', 'Md', 'Mg', 'Mn', 'Mo', 'Mt', 'N', 'Na', 'Nb', 'Nd',
            'Ne', 'Ni', 'No', 'Np', 'O', 'Os', 'P', 'Pa', 'Pb', 'Pd', 'Pm',
            'Po', 'Pr', 'Pt', 'Pu', 'Ra', 'Rb', 'Re', 'Rf', 'Rg', 'Rh', 'Rn',
            'Ru', 'S', 'Sb', 'Sc', 'Se', 'Sg', 'Si', 'Sm', 'Sn', 'Sr', 'Ta',
            'Tb', 'Tc', 'Te', 'Th', 'Ti', 'Tl', 'Tm', 'U', 'Uuo', 'Uup', 'Uus',
            'Uut', 'V', 'W', 'Xe', 'Y', 'Yb', 'Zn', 'Zr'}
elnames = {
    'Ac': 'Actinium',
    'Ag': 'Silver',
    'Al': 'Aluminium',
    'Am': 'Americium',
    'Ar': 'Argon',
    'As': 'Arsenic',
    'At': 'Astatine',
    'Au': 'Gold',
    'B': 'Boron',
    'Ba': 'Barium',
    'Be': 'Beryllium',
    'Bh': 'Bohrium',
    'Bi': 'Bismuth',
    'Bk': 'Berkelium',
    'Br': 'Bromine',
    'C': 'Carbon',
    'Ca': 'Calcium',
    'Cd': 'Cadmium',
    'Ce': 'Cerium',
    'Cf': 'Californium',
    'Cl': 'Chlorine',
    'Cm': 'Curium',
    'Cn': 'Copernicium',
    'Co': 'Cobalt',
    'Cr': 'Chromium',
    'Cs': 'Caesium',
    'Cu': 'Copper',
    'Db': 'Dubnium',
    'Ds': 'Darmstadtium',
    'Dy': 'Dysprosium',
    'Er': 'Erbium',
    'Es': 'Einsteinium',
    'Eu': 'Europium',
    'F': 'Fluorine',
    'Fe': 'Iron',
    'Fl': 'Flerovium',
    'Fm': 'Fermium',
    'Fr': 'Francium',
    'Ga': 'Gallium',
    'Gd': 'Gadolinium',
    'Ge': 'Germanium',
    'H': 'Hygrogen',
    'He': 'Helium',
    'Hf': 'Hafnium',
    'Hg': 'Mercury',
    'Ho': 'Holmium',
    'Hs': 'Hassium',
    'I': 'Iodine',
    'In': 'Indium',
    'Ir': 'Iridium',
    'K': 'Potassium',
    'Kr': 'Krypton',
    'La': 'Lanthanum',
    'Li': 'Lithium',
    'Lr': 'Lawrencium',
    'Lu': 'Lutetium',
    'Lv': 'Livermorium',
    'Md': 'Mendelevium',
    'Mg': 'Magnesium',
    'Mn': 'Manganese',
    'Mo': 'Molybdenum',
    'Mt': 'Meitnerium',
    'N': 'Nitrogen',
    'Na': 'Sodium',
    'Nb': 'Niobium',
    'Nd': 'Neodymium',
    'Ne': 'Neon',
    'Ni': 'Nickel',
    'No': 'Nobelium',
    'Np': 'Neptunium',
    'O': 'Oxygen',
    'Os': 'Osmium',
    'P': 'Phosphorus',
    'Pa': 'Protactinium',
    'Pb': 'Lead',
    'Pd': 'Palladium',
    'Pm': 'Promethium',
    'Po': 'Polonium',
    'Pr': 'Praseodymium',
    'Pt': 'Platinum',
    'Pu': 'Plutonium',
    'Ra': 'Radium',
    'Rb': 'Rubidium',
    'Re': 'Rhenium',
    'Rf': 'Rutherfordium',
    'Rg': 'Roentgenium',
    'Rh': 'Rhodium',
    'Rn': 'Radon',
    'Ru': 'Ruthenium',
    'S': 'Sulfur',
    'Sb': 'Antimony',
    'Sc': 'Scandium',
    'Se': 'Selenium',
    'Sg': 'Seaborgium',
    'Si': 'Silicon',
    'Sm': 'Samarium',
    'Sn': 'Tin',
    'Sr': 'Strontium',
    'Ta': 'Tantalum',
    'Tb': 'Terbium',
    'Tc': 'Technetium',
    'Te': 'Tellurium',
    'Th': 'Thorium',
    'Ti': 'Titanium',
    'Tl': 'Thallium',
    'Tm': 'Thulium',
    'U': 'Uranium',
    'Uuo': 'Ununoctium',
    'Uup': 'Ununpentium',
    'Uus': 'Ununseptium',
    'Uut': 'Ununtrium',
    'V': 'Vanadium',
    'W': 'Tungsten',
    'Xe': 'Xenon',
    'Y': 'Yttrium',
    'Yb': 'Ytterbium',
    'Zn': 'Zinc',
    'Zr': 'Zirconium'
}


def getWidgets(msg, color='white'):
    return [msg, pb.Percentage(), ' ',
            pb.Bar(marker=chr(0x2500), left='',
            right=''), ' ', pb.ETA(), ' ']


def log(s):
    logger.info(s)

def log_error(s):
    logger.error(s)
