import re
import xml.etree.ElementTree as ET
from pathlib import Path
from collections import defaultdict

from pint.errors import DefinitionSyntaxError
from openff.units import unit
from openff.toolkit.typing.engines.smirnoff import ForceField

minimal = False
offxml_src = Path('/home/timber/Documents/Python/openff-workspace/polymer_examples/xml examples')

if minimal:
    ff = ForceField(str(offxml_src/'minimal xml.offxml')) # causes error depending on choice of Angle and ProperTorsion entries including
else:
# define paths here
    offxml_paths = [path for path in offxml_src.iterdir() if path.suffix in ('.xml', '.offxml')]
    # print(offxml_paths)

    # xml parsing code begins
    for xml_path in offxml_paths:
        # xml_path = Path('xml examples/base_library_charges.offxml')
        # xml_path = Path('xml examples/openff_unconstrained-2.0.0.offxml')
        print('\n', xml_path.stem)
        tree = ET.parse(xml_path)
        root = tree.getroot()

        error_types = defaultdict(list)
        all_smirks = {
            'Problem' : defaultdict(list),
            'Working' : defaultdict(list),
        }

        # sort entries based on which type of error they raise
        for entry_type in ('Angle', 'Proper'):
            for entry in root.iter(entry_type):
                try:
                    u = unit.Quantity(entry.get('smirks')) 
                except Exception as e:
                    error_types[type(e)].append(entry.get('id')) 
                    err_label = 'Problem' if type(e) == DefinitionSyntaxError else 'Working'
                    all_smirks[err_label][entry_type].append(entry.get('smirks'))
        error_counts = {err : len(instances) for err, instances in error_types.items()}

        print(error_counts)
        # print(all_smirks)

        # check if all atoms in smirks are of wild type
        ATOM_TYPE = re.compile('\[(.*?)[:|;]')
        for working_status, smirks_subset in all_smirks.items():
            print('\t', working_status)
            atom_types = []
            for entry_type, coll in smirks_subset.items():
                # print('\t', entry_type)
                for entry in coll:
                    atom_types.append( set(re.findall(ATOM_TYPE, entry)) )
                    # print('\t\t', set(re.findall(ATOM_TYPE, entry)) )#, entry )

            print('\t\tAll atoms wild type:' , all(typeset == {'*'} for typeset in atom_types) and atom_types != [] )