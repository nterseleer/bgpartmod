# Marine Biogeochemical Model

A modular Python framework for simulating marine ecosystem dynamics, with a focus on phytoplankton, heterotrophs, and particulate matter aggregation processes.

## Overview

This marine biogeochemical model implements an ecosystem with phytoplankton, multiple types of heterotrophs (bacteria, ciliates), dissolved organic matter, detritus, and flocculation processes. It is designed to simulate the interactions between biological components and their environment in aquatic ecosystems, with particular emphasis on:

- Phytoplankton growth and nutrient uptake
- Bacterial activity and organic matter decomposition
- Dissolved and particulate organic matter dynamics
- Flocculation and aggregation processes
- Nutrient cycling (N, P, Si, C)

The model framework features a modular design that allows for components to be enabled or disabled, and parameter values to be easily calibrated through optimization tools.

## Model Structure

The model is organized into components, with each component representing a specific group of organisms or chemical species:

### Core Components

- **Phytoplankton (`Phyto`)**: Models primary producers with multi-nutrient (N, P, Si) limitation, variable stoichiometry, and photosynthetic dynamics
- **Heterotrophs (`Heterotrophs`)**: Models bacteria, flagellates, and ciliates with different feeding preferences
- **Dissolved Organic Matter (`DOM`)**: Models various pools of dissolved organic carbon and nutrients
- **Detritus (`Detritus`)**: Models particulate organic matter of different size classes
- **Dissolved Inorganic Matter (`DIM`)**: Models nutrients like NH4, NO3, DIP, DSi
- **Flocculation (`Flocs`)**: Models aggregation of particles into larger flocs with size-dependent properties

### Framework Features

- **Modular design**: Components can be added or removed via configuration
- **Flexible coupling**: Components can be coupled in various ways to represent different food webs
- **Multiple time scales**: Support for components operating at different timescales
- **Parameter optimization**: Built-in capabilities for model calibration with observations
- **Comprehensive diagnostics**: Extensive diagnostic outputs for model analysis

## Installation

### Prerequisites

- Python 3.8+
- NumPy, Pandas, Matplotlib, SciPy
- dill, pickle, json

### Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/username/marine-biogeochemical-model.git
   cd marine-biogeochemical-model
   ```

2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Create necessary directories:
   ```bash
   mkdir -p Figs Observations Simulations/Model_runs Simulations/References_simulations Simulations/Optimizations
   ```

## Getting Started

### Running a Basic Simulation

```python
from config import config
from utils import phys
from utils import simulation_manager as sim_manager

# Load a base configuration and setup
conf = config.MOW1
setup = phys.Setup(
    PARfromfile=True,  # Load light data from file
    tmax=20,           # Run for 20 days
    dt=1e-2,           # Main timestep (days)
    dt2=1e-3,          # Secondary timestep for fast processes
    T=10.0             # Temperature (Celsius)
)

# Run the simulation
sim = sim_manager.run_or_load_simulation(
    conf, setup,
    name="my_first_simulation",
    user_notes="Testing base configuration"
)

# Plot results
from utils import plotting
import matplotlib.pyplot as plt
plotting.plot_results(sim, ['Phy_Chl', 'NO3_concentration', 'NH4_concentration'])
plt.show()
```

### Running a Sensitivity Analysis

```python
# Starting from a base simulation, test parameter sensitivity
sim_sensitivity = sim_manager.run_sensitivity(
    sim,
    {'Phy+mu_max': [3.0, 4.0, 5.0]},  # Test different maximum growth rates
    name="growth_rate_sensitivity",
    user_notes="Testing phytoplankton growth rate sensitivity"
)

# Plot comparison
plotting.plot_results(sim_sensitivity, ['Phy_Chl', 'Phy_C'])
plt.show()
```

### Running Model Optimization

```python
from utils import optimization
from utils import observations

# Define parameters to optimize and their bounds
params_to_optimize = [
    ('Phy+mu_max', 2.0, 5.5),
    ('Phy+mortrate', 0.02, 0.15)
]

# Run optimization
opt = optimization.Optimization.run_new(
    dconf=config.MOW1,
    modkwargs={'setup': setup},
    obs=observations.Obs(),  # Load Observations
    optimized_parameters=[p[0] for p in params_to_optimize],
    bounds=([p[1] for p in params_to_optimize], 
            [p[2] for p in params_to_optimize]),
    population_size=50,
    num_generations=100
)
```

## Model Configuration

The model uses dictionary-based configuration with nested parameters for each component:

```python
model_config = {
    'formulation': 'Onur22',  # Model formulation (implementation)
    'Phy': {
        'class': phyto.Phyto,
        'parameters': {
            'mu_max': 5.2,    # Maximum growth rate (d-1)
            'QN_max': 0.15,   # Maximum N:C ratio (molN:molC)
            # ... more parameters
        },
        'coupling': {
            'coupled_NH4': 'NH4',
            'coupled_NO3': 'NO3',
            # ... coupling to other components
        },
        'initialization': {
            'C': 20.0,        # Initial carbon concentration
            'N': 3.0,         # Initial nitrogen concentration
            # ... initial conditions
        },
        'aggregate': {
            'C_tot': 'C',     # How this component contributes to totals
            # ... aggregation declarations
        }
    },
    # ... other components
}
```

## Model Output

The model produces a pandas DataFrame with all state variables and diagnostic variables. Output can be saved in various formats (pickle, feather, CSV) and includes:

- State variables over time for all components
- Diagnostic variables (e.g., growth rates, nutrient limitation factors)
- Aggregate variables (e.g., total C, N, P in the system)

## Acknowledgements

This model is based on the work of Kerimoglu et al. (2022) and other sources, with modifications and extensions.

## License

[Insert license information here]

## File Structure

The codebase is organized into several key directories and modules:

```
marine-biogeochemical-model/
├── src/                    # Source code directory
│   ├── components/         # Model components
│   │   ├── heterotrophs.py   # Heterotrophic organisms (bacteria, ciliates)
│   │   ├── phytoplankton.py  # Phytoplankton component
│   │   ├── dim.py            # Dissolved inorganic matter
│   │   ├── dom.py            # Dissolved organic matter
│   │   ├── detritus.py       # Particulate detritus
│   │   └── flocs.py          # Flocculation processes
│   ├── core/               # Core model functionality
│   │   ├── base.py           # Base classes for model components
│   │   └── model.py          # Main model class
│   ├── config/             # Configuration files
│   │   ├── base_config.py    # Base model configurations
│   │   └── varinfos.py       # Variable information and units
│   ├── utils/              # Utility modules
│   │   ├── functions.py      # Utility functions
│   │   ├── phys.py           # Physical setup configuration
│   │   ├── plotting.py       # Plotting utilities
│   │   ├── evaluation.py     # Model evaluation methods
│   │   ├── desolver.py       # Differential evolution solver
│   │   ├── optimization.py   # Optimization framework
│   │   └── simulation_manager.py  # Simulation management
│   ├── main_example.py     # Example script for running simulations
│   └── optim_main_example.py  # Example script for optimizations
├── Figs/                   # Directory for generated figures
├── Observations/           # Directory for observation data
└── simulations/            # Directory for simulation output
    ├── Simulations/        # Regular simulations
    ├── References/         # Reference simulations
    └── Optimizations/      # Optimization results
```

### Key Files and Their Roles

#### Components
- **heterotrophs.py**: Implements bacteria, flagellates, and ciliates with feeding preferences
- **phytoplankton.py**: Implements photosynthetic organisms with variable stoichiometry
- **dim.py**: Handles dissolved inorganic nutrients (NH4, NO3, DIP, DSi)
- **dom.py**: Manages dissolved organic matter pools
- **detritus.py**: Manages particulate organic matter and aggregation
- **flocs.py**: Implements flocculation processes for particles

#### Core
- **base.py**: Contains base classes for state variables and organisms
- **model.py**: Main model class that coordinates components and runs simulations

#### Utilities
- **functions.py**: Utility functions for data handling and calculations
- **phys.py**: Physical setup configurations (light, temperature, etc.)
- **plotting.py**: Comprehensive plotting utilities for model output
- **evaluation.py**: Model evaluation against observations
- **optimization.py**: Parameter optimization framework
- **simulation_manager.py**: Manages simulation runs, storage, and retrieval

## Contributing

[Insert contribution guidelines here]