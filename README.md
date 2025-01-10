# In-silico gel electrophoresis visualization

- Take a dict of {samples : [DNA lengths list]}
- Migrate according to a reference ladder's DNA size: vertical displacement
  proportion (1kb+ only)
- assign samples to howevermany rows of length `row_len`
- Visualize hypothetical plot
- Could go downstream of PCR or restriction digestion simulation

## Updates:

### 2025-01-10:

- Tested some visualization alternatives and settled on Plotly.
- See [here](./plotly_eg.html) 
- Data generation detailed in `plotly_eg.py`.

## Example static image

![example multiple rows](26x24.png)

## Requirements:

- numpy
- matplotlib
- plotly

## shortcomings:

- band intensity is singular
- unusable standalone without some sequence parsing know-how
