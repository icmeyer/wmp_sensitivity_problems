import pandas

from scalepy.isotopes.name_conversions import ElemA_to_latex
from mpullse.data import write_df

df_fresh = pandas.read_csv('/home/icmeyer/research/wmp_sensitivity_problems/problems/partial_bwr/fresh/sensitivity_run/relative_absorption.csv')
df_4MWd = pandas.read_csv('/home/icmeyer/research/wmp_sensitivity_problems/problems/partial_bwr/4MWd/sensitivity_run/relative_absorption.csv')

print(df_fresh)
print(df_4MWd)
df = df_fresh.merge(df_4MWd, on='nuclide')
df = df.rename(columns = {'relative_rate_x': 'fresh_relative_rate', 'relative_rate_y': '4MWd_relative_rate'})
df = df.sort_values('fresh_relative_rate', ascending=False)

print(df)

df = df[['nuclide', 'fresh_relative_rate', '4MWd_relative_rate']]

header_names = ['Nuclide', r'Fresh Relative Absorption [\%]', r'4MWd Relative Absorption [\%]']

formatters = [ElemA_to_latex,
              '{:.3f}'.format,
              '{:.3f}'.format]

column_format = '|c|c|c|'

caption = 'The relative absorption rate for the partial BWR assembly benchmark problem'
latex_str = df.to_latex(header = header_names,
                        formatters = formatters,
                        label='tab:pbwr_rel_abs',
                        escape=False,
                        index=False,
                        column_format=column_format,
                        caption=caption)

latex_path = 'pbwr_capture_rates.tex'
latex_str = latex_str.replace("\\\n", "\\ \hline\n")
latex_str = latex_str.replace(r"\toprule", r"\hline")
latex_str = latex_str.replace(r"\midrule", "")
latex_str = latex_str.replace(r"\bottomrule", r"")
with open(latex_path, 'w') as f:
    f.write(latex_str)


