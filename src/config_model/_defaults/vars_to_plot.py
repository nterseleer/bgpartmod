"""Fallback variable sets. See src/config_model/_defaults/__init__.py.

Only what `src/utils/plotting.py` needs to work without a user-supplied module.
"""
from src.utils.plotted_variables_sets import PlottedVariablesSet

nutvars = PlottedVariablesSet(
    name="nutrients",
    variables=['NO3_concentration', 'NH4_concentration',
               'DIP_concentration', 'DSi_concentration'],
    ncols=2,
    description="Nutrient concentrations"
)

stoichioPhy = PlottedVariablesSet(
    name="phytoplankton_stoichiometry",
    variables=['Phy_QN', 'Phy_QP', 'Phy_QSi', 'thetaC'],
    ncols=2,
    description="Phytoplankton stoichiometric ratios"
)

phyvars = PlottedVariablesSet(
    name="phytoplankton_pools",
    variables=['Phy_C', 'Phy_Chl', 'Phy_N', 'Phy_P', 'Phy_Si'],
    ncols=3,
    description="Phytoplankton elemental pools"
)

# Light-attenuation contributions, grouped by family. Component names follow
# base_config.Kerimoglu2022; add the mineral pools if the flocculation module is used.
kd_contrib_groups = {
    'kd_contrib_Detritus':     ['kd_contrib_DetL', 'kd_contrib_DetS'],
    'kd_contrib_Heterotrophs': ['kd_contrib_BacA', 'kd_contrib_BacF',
                                'kd_contrib_HF', 'kd_contrib_Cil'],
    'kd_contrib_Flocs':        ['kd_contrib_Microflocs', 'kd_contrib_Micro_in_Macro'],
}


def build_kd_contributions_list(grouped_heterotrophs=True, grouped_detritus=True,
                                grouped_flocs=True):
    """Return the ordered k_d-contribution variable list, with each family either
    merged into its single grouped series or expanded to its per-component members.
    Order is preserved: Phy, detritus, heterotrophs, mineral flocs."""
    lst = ['kd_contrib_Phy']
    lst += (['kd_contrib_Detritus'] if grouped_detritus
            else kd_contrib_groups['kd_contrib_Detritus'])
    lst += (['kd_contrib_Heterotrophs'] if grouped_heterotrophs
            else kd_contrib_groups['kd_contrib_Heterotrophs'])
    lst += (['kd_contrib_Flocs'] if grouped_flocs
            else kd_contrib_groups['kd_contrib_Flocs'])
    return lst
