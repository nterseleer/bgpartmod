from __future__ import annotations

import copy
import os
from dataclasses import asdict
from datetime import datetime
from pathlib import Path

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from src.optimizations.optimization import Optimization

from src.optimizations.description import Description
from src.optimizations.optimization_metastructure import OptimizationMetaStructure as MetaStructure
from src.utils.functions import load_json, deep_update, save_to_json


class OptimizationMetadata:

    def __init__(self):
        self._is_optimization_processed: bool = False
        self._is_optimization_run: bool = False
        self._is_override_enable: bool = False
        self._is_metadata_updated: bool = False
        self._optimization_name: str = ""
        self._meta_structure: MetaStructure | None = None
        self._descriptions: list[dict] = []
        self._configs: list[dict] = []
        self._modkwargs: str = ""
        self._results_solver: str = ""
        self._observations: list[Any] = []
        self._optimization_config: dict[str, Any] | None = None
        self._setups: list[dict] = []
        self._user_note: str = ""
        self._creation_date: datetime | None = None
        self._log_structure: dict = {}

    def _refactor_meta_descriptions_as_fct_process_status(self):
        if not self._is_optimization_processed:
            for description in self._descriptions:
                description.pop("winning_config", None)
                description.pop("best_model", None)

    def _from_optimizations_to_meta_descriptions(self,
                                                 optimization: Optimization):  # todo think about better implementation to catch "update case"
        self._meta_structure.correct_optimization_description_file_names(optimization.descriptions)
        self._meta_structure.generate_files_path_from_optimization_descriptions(optimization.descriptions,
                                                                                optimization.is_optimization_processed)
        path_description_dicts: list[dict[str, str]] = self._meta_structure.generate_deep_update_path_dict_list()

        self._descriptions.clear()
        for opt_description, path_description_dict in zip(optimization.descriptions, path_description_dicts):
            opt_description: Description
            path_description_dict: dict[str, str]

            description_summary = opt_description.description_summary()
            description_summary = deep_update(description_summary, path_description_dict)
            self._descriptions.append(description_summary)

    def _update_log_structure(self):
        self._log_structure["is_optimization_run"] = self._is_optimization_run
        self._log_structure["is_optimization_processed"] = self._is_optimization_processed
        self._log_structure["descriptions"] = self._descriptions

    def _from_optimizations_to_log_structure(self):

        self._log_structure["name"] = self._optimization_name
        self._log_structure["creation_date"] = datetime.strftime(self._creation_date,'%Y-%m-%d %H:%M:%S')
        self._log_structure["user_note"] = self._user_note
        self._update_log_structure()
        self._log_structure["modkwargs_path"] = self._modkwargs
        self._log_structure["optim_solver_config"] = self._optimization_config


    @classmethod
    def from_optimization(cls, optimization: Optimization, is_optimization_override_enable: bool = False):

        meta_data = cls()
        meta_data._optimization_name = optimization.name
        meta_data._meta_structure = MetaStructure(optimization.name)
        meta_data._creation_date = optimization.creation_date

        meta_data._is_override_enable = is_optimization_override_enable
        meta_data._is_optimization_processed = optimization.is_optimization_processed
        meta_data._is_optimization_run = optimization.is_optimization_run
        meta_data._user_note = optimization.user_note
        meta_data._is_metadata_updated = not os.path.exists(
            meta_data._meta_structure.optimization_path)  # todo check meta_data_integrity

        meta_data._from_optimizations_to_meta_descriptions(optimization)
        meta_data._modkwargs = meta_data._meta_structure.modkwargs_file
        meta_data._optimization_config = asdict(optimization.optimization_config)  # type: ignore[arg-type]
        meta_data._from_optimizations_to_log_structure()
        meta_data._decompose_meta_descriptions()
        meta_data._results_solver = meta_data.meta_structure.results_file
        meta_data._refactor_meta_descriptions_as_fct_process_status()

        return meta_data

    @classmethod
    def load(cls, optimization_name, is_optimization_override_enable: bool = False):

        meta_structure: MetaStructure = MetaStructure(optimization_name)

        if not os.path.exists(meta_structure.optimization_path):
            raise (FileExistsError(
                f"Optimization  \"{optimization_name}\" does not exist, please change correct Optimization name"))

        meta_data = cls()
        meta_data._optimization_name = optimization_name
        meta_data._meta_structure = meta_structure
        meta_data._is_new_meta_data = False
        meta_data._log_structure = load_json(meta_structure.log_file)[0]

        meta_data._descriptions = meta_data._log_structure["descriptions"]
        meta_data._decompose_meta_descriptions()
        meta_data._modkwargs = meta_structure.modkwargs_file
        meta_data._results_solver = meta_structure.results_file
        meta_data._is_optimization_processed = meta_data._log_structure["is_optimization_processed"]
        meta_data._is_optimization_run = meta_data._log_structure["is_optimization_run"]
        meta_data._optimization_config = meta_data._log_structure["optim_solver_config"]
        meta_data._user_note = meta_data._log_structure["user_note"]
        meta_data._creation_date = datetime.strptime(meta_data._log_structure["creation_date"],'%Y-%m-%d %H:%M:%S')
        meta_data._meta_structure.generate_files_path_from_meta_descriptions(meta_data._descriptions,
                                                                             meta_data._is_optimization_processed,
                                                                             meta_data._is_optimization_run)

        return meta_data

    def _decompose_meta_descriptions(self) -> None:
        for description in self._descriptions:
            self._configs.append(description["config"])
            self._setups.append(description["setup"])
            self._observations.append(description["observation"])

    def save(self):
        if self._is_metadata_updated:
            self._update_log_structure()
            save_to_json(self._log_structure, self.meta_structure.log_file)
            self._is_metadata_updated = False

    @property
    def optimization_name(self) -> str:
        return self._optimization_name

    @property
    def descriptions(self) -> list[dict[str, Any]]:
        return copy.deepcopy(self._descriptions)

    @property
    def modkwargs(self) -> str:
        return self._modkwargs

    @property
    def optimization_config(self) -> dict[str, Any]:
        return copy.deepcopy(self._optimization_config)

    @property
    def user_note(self) -> str:
        return self._user_note

    @property
    def meta_structure(self) -> MetaStructure | None:
        return copy.deepcopy(self._meta_structure)

    @property
    def creation_date(self) -> datetime | None:
        return copy.copy(self._creation_date)

    @property
    def is_optimization_run(self) -> bool:
        return self._is_optimization_run

    @property
    def is_optimization_processed(self) -> bool:
        return self._is_optimization_processed

    @is_optimization_run.setter
    def is_optimization_run(self, value) -> None:
        self._is_optimization_run = value
        self._is_metadata_updated = True

    def update_optimization_process_status(self, optimization: Optimization) -> None:

        if not self._is_optimization_run:
            raise ValueError("Unable to define optimization as processed when it is not considered to be run.")
        self._is_optimization_processed = optimization.is_optimization_processed
        self._is_metadata_updated = True
        self._from_optimizations_to_meta_descriptions(optimization)
        self._refactor_meta_descriptions_as_fct_process_status()
        self._update_log_structure()
