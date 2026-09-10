"""Minimal example: calibrate a few parameters against observations, then inspect the result.

Run from the repository root:  python src/optim_main_example.py

Needs an observation file in Observations/ (see src/utils/observations_example.py).
"""
import matplotlib.pyplot as plt

from src.config_model import base_config
from src.core import phys
from src.utils import optimization as optim
from src.utils import observations_example as observations
from src.utils import plotting as plotres

# (parameter, lower bound, upper bound). The parameter path is '<component>+<parameter>'.
PARAMETERS = [
    ('Phy+mu_max', 2.0, 5.5),
    ('Phy+mortrate', 0.02, 0.15),
    ('Phy+lysrate', 0.05, 0.30),
]

CALIBRATED_VARS = ['Phy_C', 'Phy_Chl', 'DIN_concentration', 'DIP_concentration', 'DSi_concentration']


def run_optimization():
    """Differential evolution over PARAMETERS. Keep the population small for a first try."""
    setup = phys.Setup(tmax=30., dt=1e-2, dt2=1e-3)
    return optim.Optimization.run_new(
        modkwargs={'setup': setup, 'verbose': False, 'do_diagnostics': True},
        dconf=base_config.Kerimoglu2022,
        obs=observations.Obs(),
        optimized_parameters=[p[0] for p in PARAMETERS],
        bounds=([p[1] for p in PARAMETERS], [p[2] for p in PARAMETERS]),
        population_size=20,
        num_cpus=4,
        num_generations=50,
        calibrated_vars=CALIBRATED_VARS,
    )


def process_optimization(name):
    """Reload a finished optimization and plot its best model."""
    opt = optim.Optimization.load_existing(name)
    summary = opt.process_results()
    print(f"Best score: {summary['best_score']}")
    print(f"Best parameters: {summary['best_parameters']}")

    plotres.plot_results(opt.get_best_model(), CALIBRATED_VARS,
                         observations=opt.obs, calibrated_vars=opt.calibrated_vars)
    plotres.plot_optimization_evolution(opt.df, name=opt.name)
    return opt


if __name__ == '__main__':
    optimization = run_optimization()
    process_optimization(optimization.name)
    plt.show()
