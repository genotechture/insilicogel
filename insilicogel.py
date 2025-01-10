import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
import plotly.subplots as sp


def gel_migrate(x, ladder: str = "1kb+"):
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
    if ladder == '1kb+':
        ladder_dict = {100.0: 0.934, 200.0: 0.84, 300.0: 0.76, 400.0: 0.689,
                       500.0: 0.629, 600.0: 0.583, 700.0: 0.54, 800.0: 0.5,
                       900.0: 0.466, 1000.0: 0.437, 1200.0: 0.38,
                       1500.0: 0.311, 2000.0: 0.231, 3000.0: 0.146,
                       4000.0: 0.094, 5000.0: 0.069, 6000.0: 0.054,
                       8000.0: 0.034, 10000.0: 0.02}
    # elif add other ladder types:
    else:
        raise ValueError(f'ladder type {ladder} unrecognized!.')

    # restrict domain conditional on ladder min and max:
    domain = (min(ladder_dict.keys()), max(ladder_dict.keys()))

    # ensure x is an array for consistent handling:
    x = np.array(x, dtype=np.float64)
    # clip values to the specified domain:
    x = np.clip(x, domain[0], domain[1])

    # interpolate values using the ladder dictionary:
    sizes = np.array(list(ladder_dict.keys()))
    positions = np.array(list(ladder_dict.values()))
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


def gel_plotly(samples: dict[str, list], title: str = "In-silico gel",
               ladder: str = '1kb+', visible_samples: int = 20):
    """
    Creates an interactive gel image using plotly. Reference ladder is assigned
    to a separate subplot for easy comparison. DNA bands are hover responsive,
    their color is associated with 
    Samples are scrollable along x-axis allowing scalability.

    Args:
        samples (dict[str, list]): Samples are stored as key: value pairs,
            where each value is a list of all the DNA band sizes associated
            with that sample experimental stage. I.e. since PCR band, many
            restriction digest bands.
        title (str, optional): The title of the plot.
            Defaults to "In-silico gel".
        ladder (str, optional): Which ladder to use, informing the y-axis range
            and interpolation of band sizes. Defaults to '1kb+'.
        visible_samples (int, optional): How many samples to show initially,
            prior to horizontal scrolling and plotly selections/adjustments.
            Defaults to 20.

    Raises:
        ValueError: unrecognized ladder.

    Returns:
        json: A JSON object of the plot data for rendering in HTML.
    """

    # define ladder positions:
    if ladder == '1kb+':
        ladder_sizes = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
                        1200, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000]
        ladder_positions = gel_migrate(x=ladder_sizes, ladder=ladder)

    # elif space for more ladders:

    else:
        raise ValueError(f'ladder type:{ladder} unrecognized!. '
                         'Expected ladders: "1kb+"')

    # assign hover color based on migration distance using ladder range:
    min_size, max_size = ladder_positions.min(), ladder_positions.max()
    cmap = cm.get_cmap('viridis')
    norm = mcolors.Normalize(vmin=min_size, vmax=max_size)

    # create subplots:
    fig = sp.make_subplots(
        rows=1, cols=2,
        column_widths=[0.1, 0.9],
        horizontal_spacing=0.001,
        shared_yaxes=True  # Share the y-axis to align bands
    )

    # add ladder to the first subplot:
    lane_width = 0.4
    x_pos = 0
    for band, size in zip(ladder_positions, ladder_sizes):
        fig.add_trace(
            go.Scatter(
                x=[x_pos - lane_width / 2, x_pos + lane_width / 2],
                y=[band, band],
                mode='lines',
                line=dict(color='white', width=3),
                text=[f'Ladder {ladder}: {size} bp'] * 2,
                hoverinfo='none',
                name=ladder,
                showlegend=False
            ),
            row=1, col=1
        )

        # assign band a color:
        color = mcolors.to_hex(cmap(norm(band)))

        # invisible ladder trace for informative color hover behavior:
        hover_x = np.linspace(x_pos - lane_width / 2,
                              x_pos + lane_width / 2,
                              num=5)
        hover_y = [band] * len(hover_x)
        fig.add_trace(
            go.Scatter(
                x=hover_x,
                y=hover_y,
                mode='markers',
                marker=dict(
                    size=1,
                    color=color,
                    opacity=0.0
                ),
                text=[f'{ladder}: {size} bp'] * len(hover_x),
                hoverinfo='text',
                showlegend=False
            ),
            row=1, col=1
        )

    # add sample bands to the second subplot:
    lane_width = 0.8
    lanes = list(samples.keys())
    for i, lane in enumerate(lanes):
        x_pos = i
        sizes = samples[lane]
        bands = gel_migrate(samples[lane], ladder=ladder)

        for band, size in zip(bands, sizes):
            if size >= 1:
                fig.add_trace(
                    go.Scatter(
                        x=[x_pos - lane_width / 2, x_pos + lane_width / 2],
                        y=[band, band],
                        mode='lines',
                        line=dict(color='white', width=3),
                        text=[f'{size} bp'] * 2,
                        hoverinfo='none',
                        showlegend=False
                    ),
                    row=1, col=2
                )

                # assign band a color:
                color = mcolors.to_hex(cmap(norm(band)))

                # invisible ladder trace for informative color hover behavior:
                hover_x = np.linspace(x_pos - lane_width / 2,
                                    x_pos + lane_width / 2,
                                    num=5)
                hover_y = [band] * len(hover_x)
                fig.add_trace(
                    go.Scatter(
                        x=hover_x,
                        y=hover_y,
                        mode='markers',
                        marker=dict(
                            size=1,
                            color=color,
                            opacity=0.0
                        ),
                        text=[f'{lane}: {size} bp'] * len(hover_x),
                        hoverinfo='text',
                        showlegend=False
                    ), row=1, col=2
                )

    # figure layout:
    fig.update_xaxes(
        tickangle=45,  # Rotate x-axis labels by 45 degrees
        range=[-0.5, visible_samples + 0.5],
        rangeslider=dict(visible=True),
        row=1,
        col=2
    )

    fig.update_yaxes(
        autorange='reversed',
        showgrid=False,
        zeroline=False,
        range=[0, 1],
        row=1, col=2
    )

    fig.update_layout(
        margin=dict(t=150),
        title=dict(text=title, x=0.5, y=0.95, xanchor='left', yanchor='top'),
        plot_bgcolor='dimgrey',
        xaxis1=dict(
            tickmode='array',
            tickvals=[0],
            ticktext=[f'{ladder} ladder'],
            showgrid=False,
            zeroline=False,
            range=[-0.5, 0.5],
            side='top',
            ),
        xaxis2=dict(
            tickmode='array',
            tickvals=list(range(len(lanes))),
            ticktext=lanes,
            showgrid=False,
            zeroline=False,
            side='top',
        ),
        yaxis=dict(
            title='Migration distance',
            ticks='',
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            autorange='reversed'
        ),
        height=800,
        width=1200,
        showlegend=False
    )

    return pio.to_json(fig)