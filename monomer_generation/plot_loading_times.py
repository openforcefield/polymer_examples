import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.font_manager as fm
import pandas as pd
import os
from pathlib import Path
import numpy as np

current_dir = Path(__file__).parent.resolve()
os.chdir(current_dir)

labels = {"Polymer": ["PEO_PLGA",
                    "atactic_styrene",
                    "atactic_styrene-s9",
                    "bisphenolA",
                    "naturalrubber",
                    "naturalrubber-s49",
                    "paam_modified",
                    "paam_modified-s64",
                    "PAMAM",
                    "peg_modified",
                    "peg_modified-s49",
                    "polyethylene",
                    "polyethylene-s9",
                    "polyethylmethacrylate",
                    "polyethylmethacrylate-s81",
                    "polymethylketone",
                    "polyphenyleneI",
                    "polyphenyleneII",
                    "PolyphenyleneIII",
                    "polyphenylenesulfone",
                    "polyphenylenesulfone-s16",
                    "polythiophene",
                    "polyvinylchloride",
                    "polyvinylchloride-s81",
                    "syntactic_styrene",
                    "syntactic_styrene-s49",
                    "pnipam_modified",
                    "pnipam_modified-s49",
                    "PET",
                    "vulcanizedrubber",],
          "Sugar": ["cellulose",
                    "chitin",
                    "dextran",
                    "messy_sugar",
                    "trimannose",],
          "Nucleic Acid": ["2q1r",
                        "122d",
                        "130d",
                        "133d",
                        "144d",],
          "Protein": ["6cww",
                    "7qt2",
                    "7wcc",
                    "6mtg",
                    "7xjf",
                    "7fse",
                    "7pvu",
                    "8ovp",
                    "8gt9",
                    "8fy3",
                    "8f0x",
                    "8e8i",
                    "8d1b",
                    "8bhw",
                    "8ciq",
                    "7yb4",
                    "1lyd"],
          "Peptide": ["bip23267-sup-0002-appendixs1",
                    "ja7b02319_si_002",
                    "ja7b02319_si_003",
                    "Nspe5_gaff-phi_top6pops",]
}

# Load Data
data = pd.read_csv("polymer_energies.txt", sep=", ", header=0, engine="python")
# for name in data["name"]:
#     print(f"\"{name}\",")
# print(data)

fig, (ax2, ax1) = plt.subplots(1, 2, figsize=(10, 5))

# Left Plot
# font_path = fm.findfont(fm.FontProperties(family='Times New Roman'))
# font_prop = fm.FontProperties(fname=font_path, size=12)
# plt.rcParams['font.family'] = 'Times New Roman'
# plt.rcParams['font.size'] = 12

cmap = mpl.cm.tab20
colors = [cmap.colors[i] for i in [0, 2, 4, 6, 8]]
for data_type in ["time_to_parameterize"]:
    i = 0
    for label, mols in labels.items():
        cat = data[data["name"].isin(mols)]
        ax1.scatter(cat["num_atoms"], cat[data_type], label=label, alpha=0.7, color=colors[i])
        i += 1

ax1.set_yscale("log")
ax1.set_xscale("log")

ax1.set_title('(b)', y=-0.25)
ax1.set_xlabel("Number of Atoms")
ax1.set_ylabel("Time to Parameterize (sec)")
ax1.legend()

# Right Plot
# data = data.sort_values(by="num_atoms", inplace=True)
data = data[["num_atoms", "time_to_load", "time_to_parameterize", "time_to_energy_minimize"]].to_numpy()

bins = [0, 100, 1000, 10000, 100000]
binned_data_list = []
for i in range(1, len(bins)):
    bin_start = bins[i-1]
    bin_end = bins[i]
    ids = np.where(np.logical_and(data[:,0] < bin_end, data[:,0] >= bin_start))
    bin_data = data[ids]
    binned_data_list.append(np.average(bin_data, axis=0))

data = np.array(binned_data_list)

categories = ["System Loading", "Parameterizing"]
x_indices = np.arange(1, len(bins))
times = data[:,1:].T
bar_width = 0.33
cmap = mpl.cm.tab20c

# b1 = ax2.bar(x_indices+bar_width, times[2], edgecolor='k', width=bar_width, label=categories[2], color=cmap.colors[3])
b2 = ax2.bar(x_indices+0.5*bar_width,           times[1], edgecolor='k', width=bar_width, label=categories[1], color=cmap.colors[15])
b3 = ax2.bar(x_indices-0.5*bar_width, times[0], edgecolor='k', width=bar_width, label=categories[0], color=cmap.colors[11])
ax2.set_yscale("log")

# for bars in [b1, b2, b3]:
#     bar_color = bars[0].get_facecolor()
#     for bar in bars:
#         ax2.text(
#             bar.get_x() + bar.get_width() / 2,
#             bar.get_height() + 0.3,
#             round(bar.get_height(), 1),
#             horizontalalignment='center',
#             color=bar_color,
#             weight='bold'
#         )


ax2.set_title('(a)', y=-0.25)
ax2.set_xlabel('System Size (Num. Atoms)')
ax2.set_ylabel('Average Time (sec)')
ax2.set_xticks(x_indices)
ax2.set_xticklabels(["Small \n($<10^2$)", "Medium \n($10^2$-$10^3$)", "Large\n($<10^3$-$10^4$)", "Very Large\n($>10^4$)"])
ax2.legend()

fig.tight_layout()
fig.savefig("loading_times.png", dpi=300)