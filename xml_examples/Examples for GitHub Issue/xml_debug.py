from pathlib import Path
from openff.toolkit.typing.engines.smirnoff import ForceField

offxml_src = Path.cwd()
ff1 = ForceField(str(offxml_src/'minimal working xml.offxml')) # will load without a problem
ff2 = ForceField(str(offxml_src/'minimal erring xml.offxml'))  # will raise DefinitionSyntaxError and fail
