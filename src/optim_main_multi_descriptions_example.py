import copy
import warnings

from matplotlib import pyplot as plt

from src.config_model.phys_setup import YOURI_MEDIUM_BASE
from src.config_model.strains_config import Strains, get_strain_config
from src.core import phys
from src.optimizations.description import Description
from src.optimizations.optimization import OptimizationConfig, Optimization
from src.optimizations.optimization_metadata import OptimizationMetadata as MetaData
from src.optimizations.plot_description import plot_timeseries_comparison
from src.utils import observations, evaluation, plotting


def run_test_opti(silence_warnings=False, warnings_as_errors=False):
    if silence_warnings:
        warnings.filterwarnings('ignore')
    elif warnings_as_errors:
        warnings.filterwarnings('error')

    params_to_optimize: list[tuple[str, float, float]] = [
        ('Phy+mu_max', 0.1, 3.),
        ('Phy+alpha', 2e-6, 3e-4),
    ]

    calibrated_vars = ['Phy_Chl', 'TEPC_C', "N_part_tot"]

    opt_config = OptimizationConfig(
        calibrated_vars=tuple(calibrated_vars),
        optimized_parameters=tuple(p[0] for p in params_to_optimize),
        bounds=([p[1] for p in params_to_optimize], [p[2] for p in params_to_optimize]),
        solver_kwargs={"diffScale": 0.5},
        population_size=10,
        num_cpus=3,
        num_generations=3,
        badlnl=-10000)

    setup = phys.Setup(**YOURI_MEDIUM_BASE)
    selected_strains: list[Strains] = [Strains.ODONTELLA, Strains.LAUDERIA]

    descriptions: list[Description] = []

    for strain in selected_strains:
        full_name = f"{strain.value}_medium"
        obs = observations.Obs(name=full_name, station=full_name)
        obs.format_auria_obs()

        config = get_strain_config(strain)

        descriptions.append(Description(name=full_name,
                                        setup=copy.deepcopy(setup),
                                        config=copy.deepcopy(config),
                                        observation=copy.deepcopy(obs)))

    diagnostic_vars = ['C_tot', 'N_tot', 'P_tot', 'Si_tot', 'Chl_tot',
                       'Cphy_tot', 'POC', 'PON', 'POP', 'DOC', "C_part_tot", "N_part_tot"]
    modkwargs = {'verbose': False, 'do_diagnostics': True, 'full_diagnostics': False,
                 "aggregate_vars": diagnostic_vars}

    optimization = Optimization.run_multi_config(modkwargs=modkwargs,
                                                 optimization_config=opt_config,
                                                 descriptions=descriptions)

    return optimization


def load_opt_test(optimization_name):
    metadata = MetaData.load(optimization_name)
    optimization = Optimization.build_from_metadata(metadata)
    return optimization


def optimization_details_plots(optimization: Optimization, savefig: bool = False):
    plotting.plot_optimization_evolution(optimization.optimization_results, name=optimization.name, savefig=savefig)
    plotting.plot_optimization_summary(
        optimization.optimization_results,
        optimization.optimization_config.optimized_parameters,
        name=optimization.name,
        savefig=savefig
    )

    for description in optimization.descriptions:
        evaluation.calculate_likelihood(description.best_model.df,
                                        description.observation,  # Use first obs for likelihood calculation
                                        calibrated_vars=optimization.optimization_config.calibrated_vars,
                                        verbose=True,
                                        plot=True,
                                        save_plots=savefig,
                                        file_name=f"likelihood_comparison_description_{description.name}")


if __name__ == '__main__':
    new_opt_test = run_test_opti()
    loaded_opt = load_opt_test("OPT151")
    loaded_opt.process_results()
    optimization_details_plots(loaded_opt)
    plot_timeseries_comparison(loaded_opt, loaded_opt.optimization_config.calibrated_vars)
    plt.show()
