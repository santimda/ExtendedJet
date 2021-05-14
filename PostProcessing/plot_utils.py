import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# My colors
iar_aut = ['peru','maroon','darkgreen','indigo'] # otoño en el IAR

# My plot config
def luc_std_plot(ax=None):
    luc_config = {
            # Use LaTeX to write all text
            "text.usetex": True,
            "font.family": "serif",
            # Use 10pt font in plots, to match 10pt font in document
            "axes.labelsize": 25,
            "font.size": 30,
            "axes.linewidth": 1.2,
            # Make the legend/label fonts a little smaller
            "legend.fontsize": 20,
            "xtick.labelsize": 20,
            "ytick.labelsize": 20,
            "axes.linewidth" : 1.5,
            "axes.grid": True,
            "xtick.major.size":    15,     # major tick size in points
            "xtick.minor.size":    15,       # minor tick size in points
            "xtick.major.width":   1.8,     # major tick width in points
            "xtick.minor.width":   1.6,     # minor tick width in points
            "grid.color":     'grey',  # grid color
            "grid.linestyle": '--',       # solid
            "grid.linewidth": 0.5,   # in points
            "grid.alpha":     0.5,    # transparency, between 0.0 and 1.0
    }

    # Update the fonts
    matplotlib.rcParams.update(luc_config)
    if ax is not None:
       ax.tick_params(axis="y", direction="in", length=10, width=1, color="black")
       ax.tick_params(which='minor', axis="y", direction="in", length=5, width=1, color="black")
       ax.tick_params(axis="x", direction="in", length=10, width=1, color="black")
       ax.tick_params(which='minor', axis="x", direction="in", length=5, width=1, color="black")

def color_from_cmap(cmap,ncol=10, vmin=0.0, vmax=1.0):
    cmap = matplotlib.cm.get_cmap(cmap)
    x    = np.linspace(vmin, vmax, ncol)
    norm = matplotlib.colors.BoundaryNorm(x, cmap.N)#, extend='both')
    #norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm),
    return cmap(norm(x)), sm[0]
