"""Minimal example: run the reference biogeochemical configuration and plot the results.

Run from the repository root:  python src/main_example.py
"""
import matplotlib.pyplot as plt

from src.config_model import base_config
from src.core import model, phys
from src.utils import plotting as plotres


def run_reference():
    """Kerimoglu et al. (2022) biogeochemistry, 30 days, synthetic light and temperature."""
    setup = phys.Setup(tmax=30., dt=1e-2, dt2=1e-3)
    return model.Model(base_config.Kerimoglu2022, setup=setup, name='reference_run')


def run_with_modified_parameter():
    """Same configuration with one parameter changed, to compare against the reference."""
    from src.utils import config_tools as cfg

    config = cfg.deep_update(base_config.Kerimoglu2022,
                             {'Phy': {'parameters': {'mu_max': 3.5}}})
    setup = phys.Setup(tmax=30., dt=1e-2, dt2=1e-3)
    return model.Model(config, setup=setup, name='lower_mu_max')


if __name__ == '__main__':
    simulations = [run_reference(), run_with_modified_parameter()]
    # observations=None: no observation overlay (plot_results otherwise loads a default dataset).
    plotres.plot_results(simulations, ['Phy_C', 'Phy_Chl', 'NO3_concentration', 'DIP_concentration'],
                         observations=None)
    plt.show()
