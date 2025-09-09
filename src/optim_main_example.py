import time
import matplotlib.pyplot as plt
import warnings
import numpy as np

from utils import optimization as optim
from Config_model import config
from core import phys
from utils import evaluation
from utils import observations
from utils import plotting as plotres


def configure_warning_behavior(silence_warnings=False, warnings_as_errors=False):
    """
    Configure warning behavior for model runs.

    Args:
        silence_warnings: If True, suppress all runtime warnings
        warnings_as_errors: If True, convert warnings to errors
    """
    if silence_warnings and warnings_as_errors:
        raise ValueError("Cannot both silence warnings and treat them as errors")

    if silence_warnings:
        # Suppress all warnings
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
    elif warnings_as_errors:
        # Treat warnings as errors
        warnings.filterwarnings("error", category=RuntimeWarning)
        warnings.filterwarnings("error", category=np.VisibleDeprecationWarning)
    else:
        # Default behavior: show warnings
        warnings.filterwarnings("default")



def process_optimization(opt_name, savefig=False):
    # Process results
    opt = optim.Optimization.load_existing(opt_name)
    summary = opt.process_results()
    print(f"Best score: {summary['best_score']}")
    print(f"Best parameters: {summary['best_parameters']}")
    print(f"Convergence info: {summary['convergence']}")
    print(f"Calibrated variables: {opt.calibrated_vars}")

    # Get best model
    best_model = opt.get_best_model()
    lnl = evaluation.calculate_likelihood(best_model.df,
                                          opt.obs,
                                          calibrated_vars=opt.calibrated_vars,
                                          verbose=True, plot=True)

    plotres.plot_results(best_model, plotres.phy_nuts_TEP_flocs,
                         observations=opt.obs,
                         calibrated_vars=opt.calibrated_vars,
                         ncols=4,
                         save=savefig, filename=opt.name)
    # Compare scores
    # opt.compare_model_vs_optim(obs = Observations.Obs(station='MOW1-St330-TEP-D50'))

    plotres.plot_optimization_evolution(opt.df, name=opt.name, savefig=savefig)

    # Plot parameter summary
    plotres.plot_optimization_summary(
        opt.df,
        opt.config['optimized_parameters'],
        name=opt.name,
        savefig=savefig
    )

    # Create parameter table
    # fig = plotres.create_parameter_table(
    #     opt.df,
    #     parameters=[(p, min[i], max[i])
    #                 for i, p in enumerate(opt.Config_model['optimized_parameters'])]
    # )

    return opt





def run_optim_202503(silence_warnings=False, warnings_as_errors=False):
    # Configure warning behavior
    configure_warning_behavior(silence_warnings, warnings_as_errors)

    dconf = config.MOW1
    setup = phys.Setup(**phys.DEFAULT_SETUPS['onur22'], PARfromfile=True, tmax=80, dt=1e-2, dt2=1e-3)
    params_to_optimize = [
        ('Phy+mu_max', 2., 5.5),
        ('Phy+mortrate', 0.02, 0.15), # 0.05
        ('Phy+lysrate', 0.05, 0.3),  # 0.1
    ]

    calibrated_vars = [
        'Phy_C',
        'Phy_Chl',
        'DIN_concentration',
        'DIP_concentration',
        'DSi_concentration',
    ]

    start_time = time.time()
    opt = optim.Optimization.run_new(
        modkwargs={'setup': setup, 'verbose': False, 'do_diagnostics': True, 'full_diagnostics': False},
        dconf=dconf,
        obs=observations.Obs(),
        optimized_parameters=[p[0] for p in params_to_optimize],
        bounds=([p[1] for p in params_to_optimize],
                [p[2] for p in params_to_optimize]),
        population_size=75,
        num_cpus=25,
        num_generations=2000,
        calibrated_vars=calibrated_vars
    )

    print(f"Optimization {opt.name} completed in {time.time() - start_time:.2f} seconds")
    return opt


if __name__ == "__main__":
    # Run new optimization
    # run_optim_202503(warnings_as_errors=True)

    # Process existing optimization
    opt = process_optimization('OPT060', savefig=True)
    plt.show()
