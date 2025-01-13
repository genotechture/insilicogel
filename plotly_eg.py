from insilicogel import gel_plotly
import plotly.io as pio

# example data:
example_dna_bands = {'Sample A': [500, 1000, 5000],
                     'Sample B': [1500, 5000],
                     'Sample C': [6500],
                     'Sample D': [0],  # no band
                     'Sample E': [15000],  # size capped at ladder max (10 kbp)
                     'Sample F': [50]}  # size capped at ladder min (100 bp)

# generate JSON string of plotly figure:
gel_json = gel_plotly(samples=example_dna_bands, title='Example gel',
                      ladder='1kb+', visible_samples=20)

# render on website, but since this is github we will save it as an HTML file:
gel_fig = pio.from_json(gel_json)
pio.write_html(gel_fig, file='plotly_eg.html', auto_open=True)

# Note: band migration interpolation is confined to ladder range
