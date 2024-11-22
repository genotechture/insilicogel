import matplotlib.pyplot as plt
import numpy as np


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
    fig.suptitle(f"{title} in-silico gel", color="white", fontsize=16)
    plt.tight_layout()
    return plt
