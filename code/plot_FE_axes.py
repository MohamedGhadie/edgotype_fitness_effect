#----------------------------------------------------------------------------------------
# Plot fitness effect axes with no data points.
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from plot_tools import bar_plot

def main():
        
    # figure directory
    figDir = Path('../figures') / 'combined'
    
    # create output directories if not existing
    if not figDir.exists():
        os.makedirs(figDir)
    
    bar_plot ([0, 0, 0],
              xlabels = ('Effectively\nneutral', 'Mildly\ndeleterious', 'Strongly\ndetrimental'),
              ylabel = 'Fraction (%)',
              msize = 0,
              ewidth = 0,
              edgecolor = 'white',
              fontsize = 18,
              xlim = None,
              ylim = [0, 100],
              xticks = None,
              yMinorTicks = 4,
              show = False,
              figdir = figDir,
              figname = 'mut_fitness_effect_axis')

if __name__ == "__main__":
    main()