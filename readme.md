# BGC Model

Biogeochemical model for marine ecosystems with configurable components and optimization capabilities.

## Project Structure

```
bgc-model/
├── main.py               # Main script for running simulations
├── optim_main.py        # Main script for running optimizations
├── simulations/         # All simulation outputs
│   ├── Simulations/    # Individual simulation folders
│   ├── References/     # Reference simulation folders
│   ├── Optimizations/  # Optimization results folders
│   │   ├── OPT000/    # Individual optimization results
│   │   └── ...
│   ├── Simulations_log.xlsx    # Log file for simulations
│   └── Optimizations_log.xlsx  # Log file for optimizations
├── Figs/               # Symbolic link to project figures
├── observations/       # Observation data
│   ├── processed/     # Processed observation files
│   └── raw/          # Raw observation data
├── src/               # Source code
│   ├── bgc/
│   │   ├── core/     # Core model functionality
│   │   │   ├── model.py    # Main model implementation
│   │   │   └── base.py     # Base classes
│   │   ├── components/     # Model components
│   │   │   ├── phytoplankton.py
│   │   │   ├── heterotrophs.py
│   │   │   ├── detritus.py
│   │   │   ├── dom.py
│   │   │   ├── dim.py
│   │   │   └── flocs.py
│   │   └── utils/          # Utility modules
│   │       ├── functions.py       # General utilities
│   │       ├── plotting.py        # Plotting functions
│   │       ├── config.py          # Model configurations
│   │       ├── varinfos.py        # Variable information
│   │       ├── phys.py            # Physical setup
│   │       ├── observations.py    # Observation handling
│   │       ├── simulation_manager.py  # Simulation management
│   │       ├── desolver.py        # DE solver
│   │       └── optimization.py    # Optimization framework
├── tests/             # Test files
│   ├── test_model.py
│   └── test_components/
├── examples/          # Example scripts
│   └── basic_simulation.py
└── docs/             # Documentation

```

## Directory Descriptions

- **main.py**: Primary script for running model simulations
- **optim_main.py**: Script for running model optimizations
- **simulations/**: Contains all simulation outputs and logs
  - Organized into Simulations, References, and Optimizations
  - Includes log files tracking all runs
- **observations/**: Contains observation data for model validation
- **src/bgc/**: Main source code directory
  - **core/**: Core model functionality
  - **components/**: Individual model components
  - **utils/**: Utility functions and configurations
- **tests/**: Test files
- **examples/**: Example scripts demonstrating usage
- **docs/**: Documentation files
- **Figs/**: Symbolic link to figure storage

## Notes

- The `Figs/` directory is typically a symbolic link to an external figure storage location
- Simulation results are automatically organized in the `simulations/` directory
- Observations are kept separate from simulation results
- Core functionality is separated from components for modularity
