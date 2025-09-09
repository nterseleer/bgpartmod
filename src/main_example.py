import matplotlib.pyplot as plt

from Config_model import config
from utils import plotting as plotres
from core import phys
from utils import simulation_manager as sim_manager
from utils import observations


def sensitivity_analysis():
    conf = config.MOW1 | config.Flocs
    setup = phys.Setup(**phys.DEFAULT_SETUPS['onur22'], PARfromfile=True, tmax=15, dt=1e-2, dt2=1e-3)

    sim0 = sim_manager.run_or_load_simulation(
        conf, setup,
        name="base_config_test",
        user_notes="Testing base configuration with flocs"
    )
    sim1 = sim_manager.run_sensitivity(sim0,
                                       {'Microflocs+K_glue': [5, 15], 'Microflocs+alpha_FF_ref': 0.1, },
                                       name="alphas",
                                       user_notes="Checking alphas for Micro",
                                       save=False,
                                       setup_changes={'tmax': 10}
                                       )

    obs = observations.Obs()
    simus = [sim0, sim1, ]
    plotres.plot_results(simus, plotres.phy_nuts, climatology=True, observations=obs)
    plotres.plot_results(simus, plotres.flocsrelatedvars, climatology=True, observations=obs)
    plotres.plot_results(simus, plotres.flocsvar, climatology=True, observations=obs)
    plotres.plot_results(simus, plotres.flocsvar2, climatology=True, observations=obs)

def process_optim():
    import utils.optimization as optim
    obs = observations.Obs(station='MOW1_202503')
    opt = optim.Optimization.load_existing('OPT060')
    sim0 = opt.get_best_model()

    simus = [sim0, ]
    plotres.plot_results(simus, plotres.phy_nuts, climatology=True, observations=obs)


if __name__ == "__main__":
    # sensitivity_analysis()
    process_optim()

    plt.show()