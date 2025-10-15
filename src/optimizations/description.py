from __future__ import annotations

import warnings

from typing import Optional, Any, TYPE_CHECKING

from src.core.model import Model
from src.core import phys
from src.utils.functions import load_pkl, load_json, write_dict_to_file, save_to_pkl, load_feather
from src.utils.observations import Obs

if TYPE_CHECKING:
    from src.optimizations.optimization_metastructure import OptimizationMetaStructure as MetaStructure


class Description:

    def __init__(self, name: str, setup: phys.Setup, config: dict, observation: Obs,
                 winning_config: dict = None, best_model: Model = None):
        self.name: str = name
        self.setup: phys.Setup = setup
        self.config: dict = config
        self.observation: Obs = observation
        self.winning_config: dict | None = winning_config
        self.best_model: Model | None = best_model

    @classmethod
    def from_meta_description(cls, meta_description: dict[str, Any]):
        name = meta_description["name"]
        setup = load_pkl(meta_description["setup"]["path"])
        config = load_pkl(meta_description["config"]["path"])
        try:
            observation = load_pkl(meta_description["observation"]["path"])
        except FileNotFoundError:
            warnings.warn("observation file.pkl did not exist, try to load origin file.feather")
            observation = load_feather(meta_description["observation"]["origin"])
        winning_config = None if meta_description.get("winning_config") is None else load_pkl(
            meta_description["winning_config"]["path"])
        best_model = None if meta_description.get("best_model") is None else load_pkl(
            meta_description["best_model"]["path"])

        return cls(name, setup, config, observation, winning_config, best_model)

    def description_summary(self) -> dict[str, Any]:
        dict_summary = {
            "name": self.name,
            "config": {
                "name": self.config.get("name"),
                "formulation": self.config.get("formulation"),
            },
            "setup": {
                "name": self.setup_name,
            },
            "observation": {
                "name": self.observation.name,
                "origin": self.observation.data_file
            },
            "winning_config": {
                "name": self.winning_config_name,
            },
            "best_model": {
                "name": self.best_model_name
            }

        }
        return dict_summary

    @staticmethod
    def _check_if_name_is_none(attribute_to_check: Any) -> None:
        if attribute_to_check is None:
            raise ValueError("The attribute must be initialized before returning its name.")
        if isinstance(attribute_to_check, dict) and attribute_to_check.get("name") is None:
            raise ValueError(
                "The name's attribute is not set, please correct optimization description file names (MetaStructure)")

    def save(self, meta_structure: MetaStructure):
        self._check_if_name_is_none(self.config_name)
        write_dict_to_file(self.config, meta_structure.config_file(self.config_name))
        # Save setup as pickle to preserve phys.Setup object
        if self.setup:
            self._check_if_name_is_none(self.setup_name)
            write_dict_to_file(self.setup, meta_structure.setup_file(self.setup_name))

        self.observation.save(meta_structure.observation_file(self.observation_name), is_save_summary=True)

        if isinstance(self.winning_config, dict) and len(self.winning_config) > 0:
            self._check_if_name_is_none(self.winning_config_name)
            write_dict_to_file(self.winning_config,meta_structure.winning_config_file(self.winning_config_name))
        if isinstance(self.best_model, Model):
            self._check_if_name_is_none(self.best_model_name)
            save_to_pkl(self.best_model, meta_structure.best_model_file(self.best_model_name))

    @property
    def config_name(self):
        return self.config.get("name")

    @config_name.setter
    def config_name(self, name: str):
        self.config = {"name": name} | self.config

    @property
    def setup_name(self):
        return self.setup.name if self.setup else None

    @setup_name.setter
    def setup_name(self, name: str):
        self.setup.name = name

    @property
    def observation_name(self):
        return self.observation.name

    @observation_name.setter
    def observation_name(self, name: str):
        self.observation.name = name

    @property
    def winning_config_name(self):
        if isinstance(self.winning_config, dict):
            return self.winning_config.get("name")
        return None

    @winning_config_name.setter
    def winning_config_name(self, name: str):
        if self.winning_config is None:
            self.winning_config = {}
        if isinstance(self.winning_config, dict) and len(self.winning_config) > 0:
            self.winning_config = {"name": name} | self.winning_config

    @property
    def best_model_name(self) -> Optional[str]:
        if isinstance(self.best_model, Model):
            return self.best_model.name
        return None

    @best_model_name.setter
    def best_model_name(self, name: str):
        # Can only set name if best_model exists
        if isinstance(self.best_model, Model):
            self.best_model.name = name
