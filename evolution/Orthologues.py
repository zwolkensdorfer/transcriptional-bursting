import contextlib
import os
import shutil
import urllib.parse
import urllib.request
from pathlib import Path
import sys

import pandas as pd
import scanpy

from rp2.environment import make_semver, get_package_version
from rp2.paths import get_data_path
from rp2 import hagai_2018, load_one_to_one_mouse_orthologues, load_biomart_gene_symbols_df
from evolution import Load_counts

class ort_analysis:
    def __init__(self, species, gene_ids, orthologues_df=None):
        self._condition_columns = ["replicate", "treatment", "time_point"]
        self._species = species
        self._gene_ids = gene_ids
        self._orthologues_df = orthologues_df

    @property
    def condition_columns(self):
        return self._condition_columns

    @property
    def index_columns(self):
        return ["gene"] + self._condition_columns

    @property
    def species(self):
        return self._species

    @property
    def gene_ids(self):
        return self._gene_ids

    def create_mouse_orthologue_analysis(self, species):
        if species == self._species:
            return self

        if self._orthologues_df is None:
            self._orthologues_df = load_one_to_one_mouse_orthologues().reset_index()

        gene_mask = self._orthologues_df[f"{self._species}_gene"].isin(self._gene_ids)
        orthologue_gene_ids = self._orthologues_df.loc[gene_mask, f"{species}_gene"].to_list()
        return ort_analysis(species, orthologue_gene_ids, self._orthologues_df)

    def create_rat_orthologue_analysis(self, species):
        if species == self._species:
            return self

        if self._orthologues_df is None:
            self._orthologues_df = load_one_to_one_rat_orthologues().reset_index()

        gene_mask = self._orthologues_df[f"{self._species}_gene"].isin(self._gene_ids)
        orthologue_gene_ids = self._orthologues_df.loc[gene_mask, f"{species}_gene"].to_list()
        return ort_analysis(species, orthologue_gene_ids, self._orthologues_df)



def create_default_mouse_analysis():

    mouse_lps_responsive_genes = hagai_2018.load_lps_responsive_genes()
    mouse_analysis_IDs = sorted(mouse_lps_responsive_genes)

    return ort_analysis("mouse", mouse_analysis_IDs)


def create_default_rat_analysis():

    rat_lps_responsive_genes = Load_counts.load_lps_responsive_genes("rat")
    rat_analysis_IDs = sorted(rat_lps_responsive_genes)

    return ort_analysis("rat", rat_analysis_IDs)

def create_default_pig_analysis():

    pig_lps_responsive_genes = Load_counts.load_lps_responsive_genes("pig")
    pig_analysis_IDs = sorted(pig_lps_responsive_genes)

    return ort_analysis("pig", pig_analysis_IDs)

def create_default_rabbit_analysis():

    rabbit_lps_responsive_genes = Load_counts.load_lps_responsive_genes("rabbit")
    rabbit_analysis_IDs = sorted(rabbit_lps_responsive_genes)

    return ort_analysis("rabbit", rabbit_analysis_IDs)
