#!/usr/bin/env python3
# encoding: utf-8

"""
Copied from https://github.com/xiaoruiDong/easy_rmg_model/blob/master/easy_rmg_model/template_writer/input/arkane_input.py
"""

import json

from arkane.encorr.corr import assign_frequency_scale_factor
from arkane.input import process_model_chemistry
from arkane.modelchem import LOT

from .base import BaseTemplateWriter

class ArkaneSpecies(BaseTemplateWriter):

    default_settings = {
        'model_chemistry': 'cbs-qb3',
        'freq_scale_factor': None,
        'use_bond_corrections': True,
        'atom_dict': {},
        'bond_dict': {},
        'rotors_dict': {},
        'multiplicity': 1,
        'charge': 0,
        'linear': False,
        'external_symmetry': 1,
        'optical_isomers': 1,
        'freq': './freq/output.out',
        'composite': './composite/output.out',
        'sp': None,
        'template_file': None,
        'save_path': './input.py'
    }

    default_template = """#!/usr/bin/env python3
# encoding: utf-8

modelChemistry = "{{ model_chemistry }}"
frequencyScaleFactor = {{ freq_scale_factor }}
useHinderedRotors = {{ rotors_dict | length > 0 }}
useBondCorrections = {{ use_bond_corrections }}

atoms = {{ atom_dict | safe }}

bonds = {{ bond_dict | safe }}

spinMultiplicity = {{ multiplicity }}

linear = {{ linear }}

externalSymmetry = {{ external_symmetry }}

opticalIsomers = {{ optical_isomers }}

energy = {
    '{{ model_chemistry }}': Log('{{ energy | safe }}'),
}

geometry = Log('{{ freq | safe }}')

frequencies = Log('{{ freq | safe }}')

{% if rotors_dict | length > 0 -%}
rotors = [
{%- for rotor in rotors_dict.values() %}
{%- if rotor['success'] %}
    HinderedRotor(scanLog=Log('{{ rotor['scan_path'] | safe }}'),
                  pivots={{ rotor['pivots'] }},
                  top={{ rotor['top'] }},
                  symmetry={{ rotor['symmetry'] }},
                  fit='fourier'),
{%- endif %}
{%- endfor %}
]
{%- endif %}
"""

    @property
    def model_chemistry(self):
        return self._model_chemistry

    @model_chemistry.setter
    def model_chemistry(self, value):
        if isinstance(value, str):
            self._model_chemistry = process_model_chemistry(value)
        elif isinstance(value, LOT):
            self._model_chemsitry = value

    @property
    def freq_scale_factor(self):
        if not self._freq_scale_factor:
            self._freq_scale_factor = assign_frequency_scale_factor(
                self._model_chemistry)
        return self._freq_scale_factor

    @freq_scale_factor.setter
    def freq_scale_factor(self, value):
        if value == None or \
           (isinstance(value, (int, float)) and
                value > 0 and value < 5):
            self._freq_scale_factor = value
        else:
            raise ValueError(f'Not valid frequency scale factor, got: {value}')

    @property
    def use_bond_corrections(self):
        if self.bond_dict == '{}':
            # Add a warning / or automatically generate one
            # False, otherwise will fail the job
            return False
        else:
            return self._use_bond_corrections

    @use_bond_corrections.setter
    def use_bond_corrections(self, value):
        if isinstance(value, bool):
            self._use_bond_corrections = value

    @property
    def atom_dict(self):
        return json.dumps(self._atom_dict, sort_keys=True,
                          indent=4, separators=(',', ': '))

    @atom_dict.setter
    def atom_dict(self, value):
        if isinstance(value, dict):
            self._atom_dict = value
        elif value == None:
            self._atom_dict = {}
        else:
            raise ValueError(f'Invalid atom dictionary, got {value}')

    @property
    def bond_dict(self):
        return json.dumps(self._bond_dict, sort_keys=True,
                          indent=4, separators=(',', ': '))

    @bond_dict.setter
    def bond_dict(self, value):
        if isinstance(value, dict):
            self._bond_dict = value
        elif value == None:
            self._bond_dict = {}
        else:
            raise ValueError(f'Invalid bond dictionary, got {value}')

    @property
    def energy(self):
        # TODO: a method check, check if model chemistry is the
        # same as indicated
        if hasattr(self, '_energy') and self._energy:
            return self._energy
        if hasattr(self, 'composite')\
           and self.composite:
            return self.composite
        elif hasattr(self, 'sp')\
                and self.sp:
            return self.sp
        elif hasattr(self, 'opt')\
                and self.opt:
            return self.opt
        elif hasattr(self, 'freq')\
                and self.freq:
            return self.freq

    @energy.setter
    def energy(self, value):
        if isinstance(value, str) or value == None:
            self._energy = value

    @property
    def freq(self):
        if self._freq:
            return self._freq
        elif hasattr(self, 'composite') and self.composite:
            return self.composite

    @freq.setter
    def freq(self, value):
        if isinstance(value, str) or value == None:
            self._freq = value
        else:
            raise ValueError(f'Not valid frequency, got {value}')

    def to_dict(self):
        return {'model_chemistry': self. model_chemistry.to_model_chem(),
                'freq_scale_factor': self.freq_scale_factor,
                'use_bond_corrections': self.use_bond_corrections,
                'atom_dict': self.atom_dict,
                'bond_dict': self.bond_dict,
                'rotors_dict': self.rotors_dict,
                'multiplicity': self.multiplicity,
                'charge': self.charge,
                'linear': self.linear,
                'external_symmetry': self.external_symmetry,
                'optical_isomers': self.optical_isomers,
                'energy': self.energy,
                'freq': self.freq,
                'save_path': self.save_path,
                }


class ArkaneThermo(BaseTemplateWriter):

    default_settings = {
        'model_chemistry': 'cbs-qb3',
        'freq_scale_factor': None,
        'use_bond_corrections': True,
        'use_hindered_rotors': True,
        'species_label': '',
        'species_file': 'species.py',
        'species_smiles': '',
        'thermo_type': 'NASA',
        'calc_statmech': True,
        'template_file': None,
        'save_path': './input.py'
    }

    default_template = """#!/usr/bin/env python3
# encoding: utf-8

modelChemistry = "{{ model_chemistry }}"
frequencyScaleFactor = {{ freq_scale_factor }}
useHinderedRotors = {{ use_hindered_rotors }}
useBondCorrections = {{ use_bond_corrections }}


species("{{ species_label }}",
        "{{ species_file }}",
        structure=SMILES("{{ species_smiles }}"),
        )

{%- if calc_statmech %}
statmech("{{ species_label }}")
{%- endif %}

thermo("{{ species_label }}",
       "{{ thermo_type }}")
"""

    @property
    def model_chemistry(self):
        return self._model_chemistry

    @model_chemistry.setter
    def model_chemistry(self, value):
        if isinstance(value, str):
            self._model_chemistry = process_model_chemistry(value)
        elif isinstance(value, LOT):
            self._model_chemsitry = value

    @property
    def freq_scale_factor(self):
        if not self._freq_scale_factor:
            self._freq_scale_factor = assign_frequency_scale_factor(
                self._model_chemistry)
        return self._freq_scale_factor

    @freq_scale_factor.setter
    def freq_scale_factor(self, value):
        if value == None or \
           (isinstance(value, (int, float)) and
                value > 0 and value < 5):
            self._freq_scale_factor = value
        else:
            raise ValueError(f'Not valid frequency scale factor, got: {value}')

    @property
    def species_smiles(self):
        return self._smiles

    @species_smiles.setter
    def species_smiles(self, value):
        self._smiles = value
        if hasattr(self, 'species_label') and getattr(self, 'species_label'):
            pass
        else:
            setattr(self, 'species_label', value)

    def to_dict(self):
        return {'model_chemistry': self. model_chemistry.to_model_chem(),
                'freq_scale_factor': self.freq_scale_factor,
                'use_bond_corrections': self.use_bond_corrections,
                'freq_scale_factor': self.freq_scale_factor,
                'use_hindered_rotors': self.use_hindered_rotors,
                'species_label': self.species_label,
                'species_file': self.species_file,
                'species_smiles': self.species_smiles,
                'thermo_type': self.thermo_type,
                'calc_statmech': self.calc_statmech,
                'save_path': self.save_path,
                }

class ArkaneKinetics(BaseTemplateWriter):

    default_settings = {
        'model_chemistry': 'cbs-qb3',
        'freq_scale_factor': None,
        'use_bond_corrections': False,
        'use_hindered_rotors': True,
        'reactants': {'': ''},
        'products': {'': ''},
        'TS': 'TS.py',
        'tunneling': 'Eckart',
        'template_file': None,
        'save_path': './input.py'
    }

    default_template = """#!/usr/bin/env python3
# encoding: utf-8

title = "{{ reaction_label | safe }}"

description = \
'''
{{ reaction_label | safe }}
'''
modelChemistry = "{{ model_chemistry }}"
frequencyScaleFactor = {{ freq_scale_factor }}
useHinderedRotors = {{ use_hindered_rotors }}
useBondCorrections = {{ use_bond_corrections }}
{% for smi, path in reactants %}
species("{{ smi | safe}}",
        "{{ path | safe }}",
        structure=SMILES("{{ smi | safe }}"),
        )
{%- endfor %}
{% for smi, path in products %}
species("{{ smi | safe }}",
        "{{ path | safe }}",
        structure=SMILES("{{ smi | safe }}"),
        )
{%- endfor %}

transitionState('TS', '{{ TS | safe }}')

reaction(
    label='{{ reaction_label | safe }}',
    reactants={{ reactant_list | safe }},
    products={{ product_list | safe }},
    transitionState = 'TS',
    tunneling='{{ tunneling }}',
)

kinetics(
    label='{{ reaction_label | safe }}',
    Tmin = (300,'K'), Tmax = (2000,'K'), Tcount = 200,
)

"""

    @property
    def model_chemistry(self):
        return self._model_chemistry

    @model_chemistry.setter
    def model_chemistry(self, value):
        if isinstance(value, str):
            self._model_chemistry = process_model_chemistry(value)
        elif isinstance(value, LOT):
            self._model_chemsitry = value

    @property
    def freq_scale_factor(self):
        if not self._freq_scale_factor:
            self._freq_scale_factor = assign_frequency_scale_factor(
                self._model_chemistry)
        return self._freq_scale_factor

    @freq_scale_factor.setter
    def freq_scale_factor(self, value):
        if value == None or \
           (isinstance(value, (int, float)) and
                value > 0 and value < 5):
            self._freq_scale_factor = value
        else:
            raise ValueError(f'Not valid frequency scale factor, got: {value}')

    @property
    def reaction_label(self):
        label = " + ".join(self.reactant_list)
        label += " <=> "
        label += " + ".join(self.product_list)
        return label

    @property
    def reactant_list(self):
        return [reactant[0] for reactant in self.reactants]

    @property
    def product_list(self):
        return [product[0] for product in self.products]

    def to_dict(self):
        return {'model_chemistry': self. model_chemistry.to_model_chem(),
                'freq_scale_factor': self.freq_scale_factor,
                'use_bond_corrections': self.use_bond_corrections,
                'freq_scale_factor': self.freq_scale_factor,
                'use_hindered_rotors': self.use_hindered_rotors,
                'reaction_label': self.reaction_label,
                'reactants': set(self.reactants),
                'products': set(self.products),
                'TS': self.TS,
                'reactant_list': self.reactant_list,
                'product_list': self.product_list,
                'tunneling': self.tunneling,
                'save_path': self.save_path,
                }