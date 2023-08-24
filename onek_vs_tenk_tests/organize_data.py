import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import matplotlib.patches as mpatches


profiler_dump = Path("polymer_examples/onek_vs_tenk_tests/profile_data")

times_df = pd.DataFrame(columns=["run_name", "n_atoms", "n_chains", "time_searching", "time_else", "total_time"])
cum_df_1 = pd.DataFrame(columns=["function"])
cum_df_10 = pd.DataFrame(columns=["function"])
for file in profiler_dump.iterdir():
    file_str = file.stem
    data = file_str.split("E_")[1]
    n_atoms = int(data.split("x")[0])
    chains_trip = data.split("x")[1]
    n_chains = int(chains_trip.split("_t")[0])
    sample = chains_trip.split("_t")[1]
    print(f"n atoms: {n_atoms}   n_chains: {n_chains}    sample: {sample}")

    with open(file, "r") as f:
        contents = f.readline()
        words = contents.split()
        total_seconds = float(words[-2])

    names = ["ncalls",
             "tottime",
             "percall1",
             "cumtime",
             "percall2",
             "function_name"]
    # widths = [9, 9, 9, 9, 9, 100]
    # data = pd.read_fwf(file, widths=widths, names=names, skiprows=5)
    
    # manually parse junk pstats output
    data = pd.DataFrame(columns = names)
    i = 0
    with open(file, "r") as file:
        [file.readline() for i in range(5)]
        line = file.readline()
        
        while line:
            if line == '\n':
                line = file.readline()
                continue
            line_list = line.split()
            if len(line_list) > 6:
                function_name = " ".join(line_list[5:])
            else:
                function_name = line_list[5]
            
            data.loc[len(data), :] = [*line_list[:5], function_name]

            line = file.readline()

    time_searching = float(data[data["function_name"] == "molecule.py:3490(chemical_environment_matches)"].iloc[0,3])
    times_df.loc[len(times_df), :] = [file_str, n_atoms, n_chains, time_searching, total_seconds-time_searching, total_seconds]

    if n_chains == 1:
        if len(cum_df_1) == 0:
            cum_df_1 = data[["cumtime", "function_name"]]
            cum_df_1.set_index("function_name", inplace=True)
            cum_df_1["cumtime"] = cum_df_1["cumtime"].astype(float)
            cum_df_1.columns = [n_atoms]
        else:
            new_df = data[["cumtime", "function_name"]]
            new_df.set_index("function_name", inplace=True)
            new_df["cumtime"] = new_df["cumtime"].astype(float)
            new_df.columns = [n_atoms]
            cum_df_1 = pd.concat([cum_df_1, new_df], axis=1, join="inner")
    elif n_chains == 10:
        if len(cum_df_10) == 0:
            cum_df_10 = data[["cumtime", "function_name"]]
            cum_df_10.set_index("function_name", inplace=True)
            cum_df_10["cumtime"] = cum_df_10["cumtime"].astype(float)
            cum_df_10.columns = [n_atoms]
        else:
            new_df = data[["cumtime", "function_name"]]
            new_df.set_index("function_name", inplace=True)
            new_df["cumtime"] = new_df["cumtime"].astype(float)
            new_df.columns = [n_atoms]
            cum_df_10 = pd.concat([cum_df_10, new_df], axis=1, join="inner")




# ------------------------------------------------------
# code to generate large line plot of function times:
# ------------------------------------------------------
# fig = plt.figure(figsize=(10, 100))

# names_df = pd.Series()
# for function, data in cum_df_1.iterrows():

#     data = data.groupby(by=data.index).mean()
#     if data[3000] < 10:
#         continue
#     # if data[300] < 10:
#     #     continue

#     data.sort_index(inplace=True)
#     plt.plot(data.index, data.values / total_times.values, label=function)
#     plt.text(data.index.max() + 30, data.values[-1] / total_times.values[-1], function, fontsize = 1)
#     names_df.loc[data.values[-1]] = function
# names_df.sort_index(inplace=True, ascending=False)
# names_df.to_csv("names.txt")
# plt.yscale("log")
# plt.legend()
# plt.savefig("figure.png", dpi=600)

# -----------------------------------------------------
# Code to generate bar graph of function times based 
# on simplified function run times
# -----------------------------------------------------
total_times_series1 = cum_df_1.max(axis=0).sort_index()
total_times1 = total_times_series1.groupby(by=total_times_series1.index).mean()

total_times_series10 = cum_df_10.max(axis=0).sort_index()
total_times10 = total_times_series10.groupby(by=total_times_series10.index).mean()

bottom1 = np.array([0] * 9)
bottom10 = np.array([0] * 9)
functions = {"to_openmm": "interchange.py:359(to_openmm)",
             "validate_topology": "interchange.py:170(validate_topology)",
             "bonds": "_create.py:133(_bonds)",
             "angles": "_create.py:182(_angles)",
             "propers": "_create.py:196(_propers)",
             "impropers": "_create.py:216(_impropers)",
             "vdw": "_create.py:230(_vdw)",
             "electrostatics": "_create.py:244(_electrostatics)"
             }
fig = plt.figure(dpi=600)
mpl.rcParams['hatch.linewidth'] = 0.17
cmap = mpl.cm.tab10
multiple = 10
bar_width = 0.36 * multiple
x_pos = np.arange(9).astype(float) * multiple
i = 0
legend_handles = []
for label, function in functions.items():
    cum_time_1 = cum_df_1.loc[function, :]
    cum_time_10 = cum_df_10.loc[function, :]

    cum_time_1 = cum_time_1.groupby(by=cum_time_1.index).mean().sort_index()
    cum_time_10 = cum_time_10.groupby(by=cum_time_10.index).mean().sort_index()


    plt.bar(x_pos+0.5*bar_width, cum_time_1, width=bar_width, linewidth=0.25, bottom=bottom1,  color=cmap.colors[i], hatch="///", edgecolor="k")
    plt.bar(x_pos-0.5*bar_width, cum_time_10, width=bar_width, linewidth=0.25, bottom=bottom10,            color=cmap.colors[i], hatch="oo", edgecolor="k")

    bottom1 = bottom1 + cum_time_1.to_numpy()
    bottom10 = bottom10 + cum_time_10.to_numpy()
    legend_handles.append(mpatches.Patch(color=cmap.colors[i], label=label))
    i += 1

# final bar graph for all other 
plt.bar(x_pos+0.5*bar_width, total_times1 - bottom1, width=bar_width, linewidth=0.25, bottom=bottom1, color=cmap.colors[i], hatch="///", edgecolor="k")
plt.bar(x_pos-0.5*bar_width, total_times10 - bottom10, width=bar_width, linewidth=0.25, bottom=bottom10,              color=cmap.colors[i], hatch="oo", edgecolor="k")
legend_handles.append(mpatches.Patch(color=cmap.colors[i], label="other"))

legend_handles.append(mpatches.Patch(hatch="///", label="10*N"))
legend_handles.append(mpatches.Patch(hatch="oo", label="10xN"))

plt.legend(handles=legend_handles)

plt.xticks(x_pos, labels = total_times1.index)
plt.xlabel("Total System Size (DOP)")
plt.ylabel("Loading Time (sec)")
plt.tight_layout()
plt.savefig("condensed_function_splits.png", dpi=600)

functions_of_interest = []

# times_df.to_csv("polymer_examples/onek_vs_tenk_tests/compiler_condensed_results_variable_N.txt")