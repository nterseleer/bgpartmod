import copy
import os
from collections import Counter
from itertools import zip_longest

from typing import Any

from src.config_system import path_config as pc
from src.optimizations.optimization_description import OptimizationDescription as Description


class OptimizationMetaStructure:

    def __init__(self, optimization_name: str):
        self.optimization_path: str = os.path.join(pc.OPTIMIZATIONS_ROOT_DIR, optimization_name)

        self.readme_file: str = os.path.join(self.optimization_path, pc.OPTIMIZATIONS_README_SUB_FILE)
        self.summary_file: str = os.path.join(self.optimization_path, pc.OPTIMIZATIONS_SUMMARY_SUB_FILE)
        self.results_file: str = os.path.join(self.optimization_path, pc.OPTIMIZATIONS_RESULTS_SUB_FILE)
        self.log_file: str = os.path.join(self.optimization_path, pc.OPTIMIZATIONS_LOG_SUB_FILE)
        self.calibration_file: str = os.path.join(self.optimization_path, pc.OPTIMIZATIONS_CALIBRATION_SUB_FILE)

        self._optimization_src_path: str = os.path.join(self.optimization_path, pc.OPTIMIZATIONS_SRC_SUB_DIR)
        self.observation_dir: str = os.path.join(self._optimization_src_path, pc.OPTIMIZATIONS_OBS_SUB_DIR)
        self.config_dir: str = os.path.join(self._optimization_src_path, pc.OPTIMIZATIONS_CONFIG_SUB_DIR)
        self.setup_dir: str = os.path.join(self._optimization_src_path, pc.OPTIMIZATIONS_SETUP_SUB_DIR)
        self.modkwargs_file: str = os.path.join(self._optimization_src_path, pc.OPTIMIZATIONS_MOD_KWARGS_SUB_FILE)
        self._config_files: dict = {}
        self._observation_files: dict = {}
        self._setup_files: dict = {}

        self._optimization_out_path: str = os.path.join(self.optimization_path, pc.OPTIMIZATIONS_OUT_SUB_DIR)
        self.winning_config_dir: str = os.path.join(self._optimization_out_path,
                                                    pc.OPTIMIZATIONS_WINNING_CONFIG_SUB_DIR)
        self.best_model_dir: str = os.path.join(self._optimization_out_path, pc.OPTIMIZATIONS_BEST_MODEL_SUB_DIR)
        self._winning_config_files: dict = {}
        self._best_model_files: dict = {}

    def generate_files_path_from_meta_descriptions(self, meta_descriptions: list[dict[str, Any]],
                                                   _is_optim_processed: bool,
                                                   _is_optim_run: bool) -> None:

        for meta_desc in meta_descriptions:
            self._config_files[meta_desc["config"]["name"]] = meta_desc["config"]["file"]
            self._observation_files[meta_desc["observation"]["name"]] = meta_desc["observation"]["file"]
            self._setup_files[meta_desc["setup"]["name"]] = meta_desc["setup"]["file"]

            if _is_optim_run:
                self._winning_config_files[meta_desc["winning_config"]["name"]] = meta_desc["winning_config"]["file"]
            if _is_optim_processed:
                self._best_model_files[meta_desc["best_model"]["name"]] = meta_desc["best_model"]["file"]

        self._generate_missing_directory(_is_optim_processed)

    @staticmethod
    def _ensure_valid_name(current_name: str | None, alternative_name: str) -> str:
        if isinstance(current_name, str) and len(current_name.strip()) > 0:
            return current_name

        return alternative_name

    @staticmethod
    def correct_optimization_description_file_names(opt_descriptions: list[Description]):
        counter_description_names = Counter([opt.name for opt in opt_descriptions])
        duplicates: list[str] = [name for name in counter_description_names if counter_description_names[name] > 1]
        if len(duplicates) > 0:
            raise ValueError(f"Duplicate optimization names: {duplicates}")
        for i, opt_description in enumerate(opt_descriptions):
            opt_description.config_name = OptimizationMetaStructure._ensure_valid_name(
                opt_description.config_name, f"CONFIG_{i}")
            opt_description.setup_name = OptimizationMetaStructure._ensure_valid_name(
                opt_description.setup_name, f"SETUP_{i}")
            opt_description.observation_name = OptimizationMetaStructure._ensure_valid_name(
                opt_description.observation_name, f"OBSERVATION_{i}")
            opt_description.winning_config_name = OptimizationMetaStructure._ensure_valid_name(
                opt_description.winning_config_name, f"WINNING_CONFIG_{i}")
            opt_description.best_model_name = OptimizationMetaStructure._ensure_valid_name(
                opt_description.best_model_name, f"BEST_MODEL_{i}")

    def _generate_missing_directory(self, is_optimization_processed: bool):

        src_dirs = (self.observation_dir, self.config_dir, self.setup_dir)
        for src_dir in src_dirs:
            os.makedirs(src_dir, exist_ok=True)

        if is_optimization_processed:
            out_dirs = (self.winning_config_dir, self.best_model_dir)
            for out_dir in out_dirs:
                os.makedirs(out_dir, exist_ok=True)

    @property
    def config_files(self) -> dict:
        return copy.copy(self._config_files)

    def config_file(self, name) -> str:
        return self._config_files.get(name)

    @property
    def setup_files(self) -> dict:
        return copy.copy(self._setup_files)

    def setup_file(self, name) -> str:
        return self._setup_files.get(name)

    @property
    def observation_files(self) -> dict:
        return copy.copy(self._observation_files)

    def observation_file(self, name) -> str:
        return self._observation_files.get(name)

    @property
    def winning_config_files(self) -> dict:
        return copy.copy(self._winning_config_files)

    def winning_config_file(self, name) -> str:
        return self._winning_config_files.get(name)

    @property
    def best_model_files(self) -> dict:
        return copy.copy(self._best_model_files)

    def best_model_file(self, name) -> str:
        return self._best_model_files.get(name)

    def generate_files_path_from_optimization_descriptions(self, opt_descriptions: list[Description]):

        for opt_description in opt_descriptions:
            self._config_files[opt_description.config_name] = os.path.join(self.config_dir, opt_description.config_name)
            self._setup_files[opt_description.setup_name] = os.path.join(self.setup_dir, opt_description.setup_name)
            self._observation_files[opt_description.observation_name] = opt_description.observation.datadir
            if opt_description.winning_config_name:
                self._winning_config_files[opt_description.winning_config_name] = os.path.join(
                    self.winning_config_dir, opt_description.winning_config_name)
            if opt_description.best_model_name:
                self._best_model_files[opt_description.best_model_name] = os.path.join(
                    self.best_model_dir, opt_description.best_model_name)

    def generate_deep_update_path_dict_list(self) -> list[dict[str, str]]:
        path_dict_list = []
        opt_description_path_zip = zip_longest(
            self.config_files.values(),
            self.setup_files.values(),
            self.observation_files.values(),
            self.winning_config_files.values(),
            self.best_model_files.values(),
            fillvalue=None
        )
        for conf_p, setup_p, obs_p, win_p, bm_p in opt_description_path_zip:
            path_dict: dict[str, str] = {
                "config": {"path": conf_p},
                "setup": {"path": setup_p},
                "observation": {"path": obs_p},
                "winning_config": {"path": win_p},
                "best_model": {"path": bm_p}
            }
            path_dict_list.append(path_dict)

        return path_dict_list
