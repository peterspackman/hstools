from pymatgen.core.structure import Structure as PMGStructure

COVALENT_RADII = [
    0.23, 1.50, 1.28, 0.96, 0.83, 0.68, 0.68, 0.68, 0.64, 1.50, 
    1.66, 1.41, 1.21, 1.20, 1.05, 1.02, 0.99, 1.51, 2.03, 1.76, 
    1.70, 1.60, 1.53, 1.39, 1.61, 1.52, 1.26, 1.24, 1.32, 1.22, 
    1.22, 1.17, 1.21, 1.22, 1.21, 1.50, 2.20, 1.95, 1.90, 1.75, 
    1.64, 1.54, 1.47, 1.46, 1.45, 1.39, 1.45, 1.44, 1.42, 1.39, 
    1.39, 1.47, 1.40, 1.50, 2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 
    1.99, 1.98, 1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.90, 1.87, 
    1.87, 1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.50, 1.32, 
    1.45, 1.46, 1.48, 1.40, 1.21, 1.50, 2.60, 2.21, 2.15, 2.06, 
    2.00, 1.96, 1.90, 1.87, 1.80, 1.69, 1.54, 1.83, 1.50, 1.50, 
    1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50  
]

ATOM_BOND_RANGE_FACTOR = 0.211671

def bonded(n1, n2, distance):
    threshold = (
            COVALENT_RADII[n1 - 1] + 
            COVALENT_RADII[n2 - 1])
    bond_min = max(threshold - ATOM_BOND_RANGE_FACTOR, 0)
    bond_max = threshold + ATOM_BOND_RANGE_FACTOR
    bonded = bond_min < distance and distance < bond_max
    return bonded


def merge_common_sets(sets):
    merged = True
    while merged:
        merged = False
        results = []
        while sets:
            common, rest = sets[0], sets[1:]
            sets = []
            for x in rest:
                if x.isdisjoint(common):
                    sets.append(x)
                else:
                    merged = True
                    common |= x
            results.append(common)
        sets = results
    return sets

class Structure(PMGStructure):

    def get_bonded_sites(self, site_number, max_radius=max(COVALENT_RADII)):
        bonded_site_indices = set()
        site = self.sites[site_number]
        element = site.specie.number
        r = COVALENT_RADII[element - 1]
        neighbors = self.get_neighbors(site, r + max_radius, include_index=True)

        for neighbor, distance, index in neighbors:
            if bonded(element, neighbor.specie.number, distance):
                bonded_site_indices.add(index)
        return bonded_site_indices


# NOT WORKING YET
#   def expand_connections(self, sites):
#       # coorce to set
#       if isinstance(sites, int):
#           sites = {sites}
#       if not isinstance(sites, set):
#           sites = set(sites)

#       elements = set(x.specie.number for x in self.sites)
#       max_radius = max(COVALENT_RADII[e -1] for e in elements)
#       unexplored = set(range(len(sites)))
#       explored = set()
#       fragments = []

#       while unexplored:
#           site = unexplored.pop()
#           if site in explored:
#               continue
#           explored.add(site)
#           bonded = self.get_bonded_sites(site, max_radius=max_radius)
#           unexplored |= {x for x in bonded if x not in explored}

#       return explored
       


    def get_connected_fragments(self):
        """Get all bonded fragments within this structure"""
        molecules = []

        elements = set(x.specie.number for x in self.sites)
        max_radius = max(COVALENT_RADII[e -1] for e in elements)
        unexplored = set(range(len(self.sites)))

        while unexplored:
            site = unexplored.pop()
            bonded = self.get_bonded_sites(site, max_radius=max_radius)
            bonded.add(site)
            molecules.append(bonded)

        molecules = merge_common_sets(molecules)
        return molecules


    def get_symmetry_unique(self):
        pass


