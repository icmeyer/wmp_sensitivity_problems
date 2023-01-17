import pandas as pd
import numpy as np
import glob
import os
import openmc

from scalepy.isotopes.name_conversions import ElemA_to_latex
from mpullse.data import write_df

def make_capture_rate_table(sp_path, table_dir, benchmark_name):
    sp = openmc.StatePoint(sp_path)

    total_abs = sp.get_tally(name="total absorption")
    df = total_abs.get_pandas_dataframe()
    total_mean = df['mean'].values[0]
    total_sd = df['std. dev.'].values[0]

    dfs = []
    for tallyid in sp._tallies:
        tally = sp._tallies[tallyid]
        if tally._name == 'total absorption':
            continue
        if tally.scores == ['absorption']:
            dfs.append(tally.get_pandas_dataframe())
    df = pd.concat(dfs, ignore_index=True)
    df = df.sort_values(by='mean', ascending=False)
    df['relative_rate'] = df['mean']/total_mean * 100

    table_name = 'relative_absorption'
    table_path = table_dir + table_name + '.txt'
    write_df(df, table_path)

    df = df[['nuclide', 'relative_rate']]

    header_names = ['Nuclide', r'Relative Absorption [\%]']

    formatters = [ElemA_to_latex,
                  '{:.3f}'.format]

    column_format = '|c|c|'

    caption = 'The relative absorption rate for the {:s} benchmark problem'.format(benchmark_name)
    latex_str = df.to_latex(header = header_names,
                            formatters = formatters,
                            label='tab:{:s}_rel_abs'.format(benchmark_name),
                            escape=False,
                            index=False,
                            column_format=column_format,
                            caption=caption)

    latex_path = table_dir + table_name + '.tex'
    latex_str = latex_str.replace("\\\n", "\\ \hline\n")
    latex_str = latex_str.replace(r"\toprule", r"\hline")
    latex_str = latex_str.replace(r"\midrule", "")
    latex_str = latex_str.replace(r"\bottomrule", r"")
    with open(latex_path, 'w') as f:
        f.write(latex_str)



if __name__=='__main__':
    make_capture_rate_table('/home/icmeyer/research/wmp_sensitivity_problems/problems/kritz-2/sensitivity_run/statepoint.100.h5', '/home/icmeyer/research/wmp_sensitivity_problems/problems/kritz-2/sensitivity_run/', 'KRITZ-2')





