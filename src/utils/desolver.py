import concurrent.futures
import numpy as np
import pandas as pd
import time
import os
from datetime import datetime
from typing import Any, Callable, List, Optional, Tuple, Union


class DESolver:
    """Differential Evolution solver with parallel processing support."""

    def __init__(self,
                 job: Any,
                 populationSize: int,
                 maxGenerations: int,
                 minInitialValue: Union[List[float], np.ndarray],
                 maxInitialValue: Union[List[float], np.ndarray],
                 diffScale: float = 0.5,
                 crossoverProb: float = 0.9,
                 initialpopulation: Optional[np.ndarray] = None,
                 num_cpus: Optional[int] = 30,
                 n_unused_cpu: int = 2):
        """
        Initialize the Differential Evolution solver.

        Args:
            job: Object with evaluate_fitness method
            populationSize: Size of population per generation
            maxGenerations: Maximum number of generations to run
            minInitialValue: Minimum parameter values
            maxInitialValue: Maximum parameter values
            diffScale: Differential weight (F)
            crossoverProb: Crossover probability (CR)
            initialpopulation: Optional initial population
            num_cpus: Number of CPUs to use (default: max - n_unused_cpu)
            n_unused_cpu: Number of CPUs to leave unused
        """
        # Store job object
        self.job = job

        # Parameter constraints and population settings
        self.minInitialValue = np.asarray(minInitialValue)
        self.maxInitialValue = np.asarray(maxInitialValue)
        assert len(self.minInitialValue) == len(self.maxInitialValue), 'Min/max parameter lengths do not match'

        self.parameterCount = len(self.minInitialValue)
        self.populationSize = populationSize
        self.maxGenerations = maxGenerations

        # DE algorithm parameters
        self.scale = diffScale
        self.crossOverProbability = crossoverProb
        self.strictbounds = True

        # Parallel processing settings
        if num_cpus is None:
            num_cpus = max(1, os.cpu_count() - n_unused_cpu)
        self.num_cpus = num_cpus
        print('Number of used cpus: ', self.num_cpus)

        # Initial population handling
        self.initialpopulation = initialpopulation
        if self.initialpopulation is not None:
            assert self.initialpopulation.ndim == 2, 'Initial population must be a 2D array'
            assert self.initialpopulation.shape[1] == self.parameterCount, \
                f'Initial population must have {self.parameterCount} columns'

        # Initialize random state
        self.randomstate = np.random.RandomState()

        # Adjust population size to match CPU count if needed
        jobsperworker = int(np.round(self.populationSize / float(self.num_cpus)))
        if self.populationSize != jobsperworker * self.num_cpus:
            self.populationSize = jobsperworker * self.num_cpus
            print(f'Adjusted population size to {self.populationSize} to match {self.num_cpus} CPUs')

    def evaluate_trial(self, trial: np.ndarray) -> Tuple[np.ndarray, float]:
        """
        Evaluate a trial solution.

        Args:
            trial: Array of parameter values to evaluate

        Returns:
            Tuple of (parameter array, likelihood)
        """
        # Get the actual log likelihood (will be negative)
        likelihood = self.job.evaluate_model(trial)

        # Return the actual likelihood - DE will handle negation internally
        return trial, likelihood

    def Solve(self) -> bool:
        """Run the differential evolution optimization."""
        # Initialize population
        if self.initialpopulation is not None:
            self.population = self.initialpopulation[-self.populationSize:, :]
        else:
            self.population = self.randomstate.uniform(
                self.minInitialValue,
                self.maxInitialValue,
                size=(self.populationSize, self.parameterCount)
            )

        cost = np.full(self.population.shape[0], np.inf)
        ibest = None
        file_exists = False

        try:
            with concurrent.futures.ProcessPoolExecutor(max_workers=self.num_cpus) as executor:
                for generation in range(self.maxGenerations):
                    print(f'Processing generation {generation}')
                    results = []

                    # Generate trials
                    trials = [
                        self.generateNew(itarget, ibest, F=self.scale, CR=self.crossOverProbability)
                        for itarget in range(self.population.shape[0])
                    ]

                    # Evaluate trials in parallel
                    future_to_trial = {
                        executor.submit(self.evaluate_trial, trial): (trial, i)
                        for i, trial in enumerate(trials)
                    }

                    # Process results as they complete
                    for future in concurrent.futures.as_completed(future_to_trial):
                        original_trial, itarget = future_to_trial[future]
                        try:
                            trial, likelihood = future.result()

                            # Store parameters and actual likelihood
                            param_dict = {
                                name: value
                                for name, value in zip(self.job.config['optimized_parameters'], trial)
                            }
                            param_dict['cost'] = likelihood  # Store actual likelihood
                            results.append(param_dict)

                            # For DE optimization, minimize negative likelihood
                            neg_likelihood = -likelihood

                            if neg_likelihood < cost[itarget]:  # Using < because we're minimizing the negative
                                self.population[itarget, :] = trial
                                cost[itarget] = neg_likelihood
                                if ibest is None or neg_likelihood < cost[ibest]:
                                    ibest = itarget

                        except Exception as e:
                            print(f'Trial evaluation failed: {e}')

                    # Create results dataframe
                    newdf = pd.DataFrame(results)
                    newdf['generation'] = generation + 1

                    # Save results using actual likelihoods
                    if generation == 0:
                        resdf = newdf
                    else:
                        newwinning = newdf['cost'] > olddf['cost']  # Higher likelihood is better
                        resdf = newdf.where(newwinning, olddf)
                        resdf['generation'] = generation + 1

                    # Add timestamp column
                    resdf['timestamp'] = datetime.now().strftime('%Y/%m/%d %H:%M:%S')
                    
                    # Append to results file
                    resdf.to_csv(
                        self.job.files['results'],
                        mode='a',
                        header=not file_exists,
                        index=False
                    )
                    file_exists = True
                    olddf = resdf

                    # if generation % 10 == 0:
                    best_likelihood = -cost[ibest] if ibest is not None else -np.inf
                    print(f'Generation {generation}: Best likelihood = {best_likelihood}')

        except Exception as e:
            print(f'Optimization failed: {e}')
            return False

        return True

    def generateNew(self, itarget: int, ibest: Optional[int],
                    CR: float = 0.9, F: float = 0.5,
                    ndiffvector: int = 1,
                    randomancestor: bool = True) -> np.ndarray:
        """
        Generate new trial vector using differential evolution algorithm.

        Args:
            itarget: Index of target vector
            ibest: Index of best vector (or None)
            CR: Crossover probability
            F: Differential weight
            ndiffvector: Number of difference vectors to use
            randomancestor: Whether to use random ancestor instead of best
        """
        # Select base vector (ancestor)
        available = list(range(self.population.shape[0]))
        available.remove(itarget)

        if ibest is not None and not randomancestor:
            available.remove(ibest)
            ancestor = self.population[ibest]
        else:
            idx = self.randomstate.choice(available)
            available.remove(idx)
            ancestor = self.population[idx]

        # Create difference vector(s)
        difference = np.zeros_like(ancestor)
        for _ in range(ndiffvector):
            idx1, idx2 = self.randomstate.choice(available, size=2, replace=False)
            difference += self.population[idx1] - self.population[idx2]
            # Remove used indices to prevent reuse
            available.remove(idx1)
            available.remove(idx2)

        # Create mutant vector
        mutant = ancestor + F * difference

        # Crossover
        cross = self.randomstate.random(self.parameterCount) <= CR
        cross[self.randomstate.randint(self.parameterCount)] = True
        trial = np.where(cross, mutant, self.population[itarget])

        # Enforce bounds
        if self.strictbounds:
            while np.any((trial < self.minInitialValue) | (trial > self.maxInitialValue)):
                trial = self.minInitialValue + np.abs(self.minInitialValue - trial)
                trial = self.maxInitialValue - np.abs(self.maxInitialValue - trial)

        return trial