from pathlib import Path
from openff.toolkit.typing.engines.smirnoff import ForceField

offxml_src = Path('/home/timber/Documents/Python/openff-workspace/polymer_examples/xml examples/base_library_charges.offxml')
ff = ForceField(str(offxml_src))