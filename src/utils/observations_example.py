import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import os
from typing import Dict

from src.utils import functions as fns
from src.Config_system import path_config as path_cfg


class Obs:
    def __init__(self,
                 name='ObservationData',
                 station='MOW1',
                 datadir=None,
                 read_df=True
                 ):
        """
        Initialize observation data handler.

        Args:
            name: Name identifier for the dataset
            station: Station identifier
            datadir: Optional custom data directory path. If None, uses project's Observations directory
            read_df: Whether to immediately read the data file
        """
        self.name = name
        self.station = station
        self.datadir = datadir if datadir is not None else path_cfg.OBSERVATION_DIR

        # Ensure the data directory exists
        if not os.path.exists(self.datadir):
            raise FileNotFoundError(
                f"Data directory not found: {self.datadir}\n"
                f"Please ensure the Observations directory exists at the project root"
            )

        if read_df:
            self.df = pd.read_feather(os.path.join(self.datadir, f"{station}.feather"))

    def create_summary(self) -> Dict:
        """
        Create metadata summary of observation data.

        Returns:
            Dictionary containing observation metadata and statistics
        """
        if not hasattr(self, 'df'):
            raise AttributeError("No data loaded to summarize")

        summary = {
            "station": self.station,
            "name": self.name,
            "timespan": {
                "start": str(self.df.index.min()),
                "end": str(self.df.index.max())
            },
            "variables": {}
        }

        for column in self.df.columns:
            col_data = self.df[column]
            summary["variables"][column] = {
                "count": int(len(col_data)),  # Convert numpy types to native Python
                "valid": int(col_data.count()),
                "mean": fns.serialize_for_json(col_data.mean()),
                "std": fns.serialize_for_json(col_data.std()),
                "min": fns.serialize_for_json(col_data.min()),
                "max": fns.serialize_for_json(col_data.max())
            }

        return summary


    def read_obs(self, fname, readdic={}):
        df = pd.read_csv(self.datadir + fname, **readdic)
        df['date'] = pd.to_datetime(df['year'].astype(str) + ' ' + df['julian_day'].astype(str), format='%Y %j.%f')
        df.set_index('date', inplace=True)
        return df

    def convert_cols_from_obs_to_mod(self, df, dconv):
        oldnamesindata = df.columns
        # Translate between name in dataset and name used in model
        usedVars = []
        dtranslation = {}
        for nameindata, dico in dconv.items():
            if nameindata in df.columns:
                if isinstance(dico['correspondingVars'], list):
                    newnames = dico['correspondingVars']
                else:
                    newnames = (dico['correspondingVars'],)
                for newname in newnames:
                    usedVars.append(newname)
                    dtranslation[newname] = nameindata
                    df[newname] = df[nameindata]
                if not nameindata in newnames:
                    df.drop(columns=nameindata, inplace=True)

        # Keep transformed and untransformed data
        for col in df.columns:
            if col in usedVars:
                df['m' + col] = df[col]
                df[col] = df[col] * dconv.get(dtranslation[col], {}).get('trsfrm', 1)

        # Compute additional variables from others
        for var in list(set(dconv.keys() - set(oldnamesindata))):
            df[var] = fns.eval_expr(dconv[var]['oprt'], df, None)
        # Keep only relevant variables
        df = df[usedVars]
        df.plot()
        return df

    def compute_climatology(self, df):
        df = df.assign(julian_day=df.index.dayofyear)
        df = df.groupby('julian_day').mean()
        return df

    def write_climatology(self, df,):
        df.to_feather(self.datadir + self.station + '.feather')

    def plot_df(self, columns = None):
        columns = columns or self.df.columns
        for col in columns:
            plt.figure()
            plt.scatter(self.df.index, self.df[col])
            plt.plot(self.df.index, self.df[col])
            plt.title(col + '_' + self.station)
            plt.xlabel(self.df.index.name)
            plt.tight_layout()
        # plt.show()


if __name__ == "__main__":
    myobs = Obs()
