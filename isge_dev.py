# --------------------------------------------------------------------------- #
#                        in-silico gel electrophoresis                        #
# --------------------------------------------------------------------------- #

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# goal: render sequences lengths as dna bands in a gel

# for 1kb ladder vertical displacement = -
# [7, 12, 19, 24, 33, 51, 81, 109, 133, 153, 163, 175, 189, 204, 220, 241, 266, 294, 327]

# or + and plot inverted

v_disp_sizekb_dict = {7: 10.0, 12: 8.0, 19: 6.0, 24: 5.0, 33: 4.0, 51: 3.0,
                      81: 2.0, 109: 1.5, 133: 1.2, 153: 1.0, 163: 0.9,
                      175: 0.8, 189: 0.7, 204: 0.6, 220: 0.5, 241: 0.4,
                      266: 0.3, 294: 0.2, 327: 0.1}

# max v displacement = 350

size_to_frac_dict = {}
for key, value in v_disp_sizekb_dict.items():
    size_to_frac_dict[value * 1000] = key / 350

df = pd.DataFrame({'DNA_size': size_to_frac_dict.keys(),
                   'vertical_pos': size_to_frac_dict.values()})

sns.lineplot(data=df, x='DNA_size', y='vertical_pos')
plt.title("DNA Size vs. Vertical Displacement of 1kb+ ladder")
plt.xlabel("DNA Size (bp)")
plt.ylabel("Normalized Vertical Position")

plt.show()

# order dictionary for correct interpolation:
ladder_1kb = dict(sorted(size_to_frac_dict.items()))
new_dict = {}
for k, v in ladder_1kb.items():
    new_dict[k] = round(number=v, ndigits=3)

ladder_1kb = new_dict


def gel_migrate(x, ladder: str = "1kb"):

    """
    Convert DNA lengths to migration distances based on the specified ladder.

    Args:
        x (float or list of floats): DNA size(s) in base pairs.
        ladder (str): Which ladder to use, ladders are internal dictionaries
            of DNA fragment length keys to migration distance proportion
            values.

    Returns:

        float or list of floats: Normalized vertical position(s).

    """
    if ladder == '1kb':
        ladder = {100.0: 0.934, 200.0: 0.84, 300.0: 0.76, 400.0: 0.689,
                  500.0: 0.629, 600.0: 0.583, 700.0: 0.54, 800.0: 0.5,
                  900.0: 0.466, 1000.0: 0.437, 1200.0: 0.38, 1500.0: 0.311,
                  2000.0: 0.231, 3000.0: 0.146, 4000.0: 0.094, 5000.0: 0.069,
                  6000.0: 0.054, 8000.0: 0.034, 10000.0: 0.02}
    # space for other ladder types...
    else:
        raise ValueError(f'ladder type {ladder} unrecognized!. Expected '
                         f'ladders: "1kb".')

    # restrict domain conditional on ladder min and max:
    domain = (min(ladder.keys()), max(ladder.keys()))

    # ensure x is an array for consistent handling:
    x = np.array(x, dtype=np.float64)
    # clip values to the specified domain:
    x = np.clip(x, domain[0], domain[1])

    # interpolate values using the ladder dictionary:
    sizes = np.array(list(ladder.keys()))
    positions = np.array(list(ladder.values()))
    interpolated_positions = np.interp(x, sizes, positions)

    return interpolated_positions


def gel_plot(x: dict[str, list], title: str, ladder: str = '1kb',
             row_len: int = 12):
    """
    Plot simulated gel electrophoresis with multiple rows of lanes.

    Args:
        x (dict): A dictionary of sample names and their corresponding DNA
        sizes (list of floats).
        title (str): A descriptor to display as plot title.
        ladder (str): Type of ladder to use (currently only "1kb" recognized).
            Defaults to "1kb".
        row_len (int): Number of samples per row (excluding flanking ladders).

    Returns:
        None
    """
    if ladder == '1kb':
        ladder_dict = {100.0: 0.934, 200.0: 0.84, 300.0: 0.76, 400.0: 0.689,
                       500.0: 0.629, 600.0: 0.583, 700.0: 0.54, 800.0: 0.5,
                       900.0: 0.466, 1000.0: 0.437, 1200.0: 0.38,
                       1500.0: 0.311, 2000.0: 0.231, 3000.0: 0.146,
                       4000.0: 0.094, 5000.0: 0.069, 6000.0: 0.054,
                       8000.0: 0.034, 10000.0: 0.02}
    else:
        raise ValueError(f'Ladder type {ladder} unrecognized!')

    # calculate the number of rows needed, rounding up:
    total_samples = len(x)
    rows = (total_samples + row_len - 1) // row_len

    # fill out row if needed:
    sample_names = list(x.keys())
    padding = row_len - (
        len(sample_names) % row_len
        ) if len(sample_names) % row_len != 0 else 0
    sample_names += ["Empty"] * padding

    # convert ladder and sample DNA sizes to migration positions:
    ladder_positions = gel_migrate(list(ladder_dict.keys()), ladder=ladder)
    samp_dict = {name: gel_migrate(sizes, ladder=ladder)
                 for name, sizes in x.items()}

    # empty lanes no bands:
    samp_dict.update({"Empty": []})

    # create the figure with subplots for each row:
    fig, axes = plt.subplots(rows, 1, figsize=(15, rows * 5),
                             facecolor='dimgrey', sharex=True)

    # ensure axes iterable for single row:
    if rows == 1:
        axes = [axes]

    for row_idx, ax in enumerate(axes):
        start_idx = row_idx * row_len
        end_idx = start_idx + row_len
        current_row_samples = sample_names[start_idx:end_idx]

        # add flanking ladders to row:
        lane_positions = ['Ladder'] + current_row_samples + ['Ladder']
        lane_migrations = (
            [ladder_positions] +
            [samp_dict.get(name, []) for name in current_row_samples] +
            [ladder_positions]
        )

        # plot current row:
        for col_idx, (lane, migrations) in enumerate(zip(lane_positions,
                                                         lane_migrations)):
            for pos in migrations:
                ax.plot([col_idx - 0.4, col_idx + 0.4],
                        [pos, pos], color="white", linewidth=5)
                if lane == "Ladder":
                    for band_size, band_pos in zip(ladder_dict.keys(),
                                                   ladder_positions):
                        # match position with a small tolerance:
                        if np.isclose(pos, band_pos, atol=1e-3):
                            ax.text(col_idx, pos, f"{int(band_size)}",
                                    ha="center", va="center", color="red",
                                    fontsize=8)

            # add lane labels above plot (below because inverted):
            ax.text(col_idx, -0.05, lane, ha="center", va="bottom",
                    color="white")

        # set appearance of current row:
        ax.set_ylim(1, 0)
        ax.set_xlim(-0.5, len(lane_positions) - 0.5)
        ax.axis("off")

    # add title for entire figure:
    fig.suptitle(f"{title} insilico gel", color="white", fontsize=16)
    plt.tight_layout()
    return plt


samples = {
    "Sample 1": [200, 500, 1500],
    "Sample 2": [100, 400, 5000],
    "Sample 3": [300, 1000],
    "Sample 4": [600, 800, 1200],
    "Sample 5": [3000],
    "Sample 6": [700, 8000],
    "Sample 7": [900],
    "Sample 8": [4000],
    "Sample 9": [500],
    "Sample 10": [1500, 2000],
    "Sample 11": [600],
    "Sample 12": [10000],
    "Sample 13": [800],
    "Sample 14": [300, 400],
}

x = gel_plot(samples, title='title', ladder="1kb", row_len=12)
x.show()

samples = {'A': [], 'B': [100, 50, 500, 5000]}

x = gel_plot(samples, title="two example", ladder="1kb", row_len=12)
x.show()

samples = {
    "Sample 1": [200, 500, 1500],
    "Sample 2": [100, 400, 5000],
    "Sample 3": [300, 1000],
    "Sample 4": [600, 800, 1200],
    "Sample 5": [3000],
    "Sample 6": [700, 8000],
    "Sample 7": [900],
    "Sample 8": [4000],
    "Sample 9": [500],
    "Sample 10": [1500, 2000],
    "Sample 11": [600],
    "Sample 12": [10000],
    "Sample 13": [800],
    "Sample 14": [300, 400],
    "Sample 15": [500],
    "Sample 16": [500],
    "Sample 1aj": [500],
    "Sample 1vc": [500],
    "Sample 1fg": [500],
    "Sample 2f": [500],
    "Sample 1a": [500],
    "Sample 1b": [500],
    "Sample 1c": [500],
    "Sample 1d": [500],
    "Sample 1e": [500],
    "Sample 3f": [500],
}

x = gel_plot(samples, title='26 by 24', ladder="1kb", row_len=24)
x.show()
x = gel_plot(samples, title='26 by 12', ladder="1kb", row_len=12)
x.show()
