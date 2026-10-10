import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# My colors
nice_colors = ['peru', 'maroon', 'darkgreen',
               'indigo', 'plum']  # autumn at IAR

# My plot config
nice_fonts = {
    # Use LaTeX to write all text
    "text.usetex": True,
    "font.family": 'serif',
    # Use 10pt font in plots, to match 10pt font in the document
    "axes.labelsize": 18,
    "font.size": 20,
    "axes.linewidth": 1.5,
    "axes.titlesize": 16,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 14,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "xtick.major.size": 5.5,  # major tick size in points
    "xtick.minor.size": 4,    # minor tick size in points
    "ytick.major.size": 5.5,  # major tick size in points
    "ytick.minor.size": 3.5,  # minor tick size in points
    "xtick.major.width": 1.4,  # major tick width in points
    "xtick.minor.width": 1.3,  # minor tick width in points
    "ytick.major.width": 1.4,  # major tick width in points
    "ytick.minor.width": 1.3,  # minor tick width in points
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
}

matplotlib.rcParams.update(nice_fonts)

def color_from_cmap(cmap,ncol=10, vmin=0.0, vmax=1.0):
    cmap = matplotlib.cm.get_cmap(cmap)
    x    = np.linspace(vmin, vmax, ncol)
    norm = matplotlib.colors.BoundaryNorm(x, cmap.N)#, extend='both')
    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm),
    return cmap(norm(x)), sm[0]
