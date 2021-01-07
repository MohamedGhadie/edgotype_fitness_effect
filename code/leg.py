
import os
from pathlib import Path
from plot_tools import multi_bar_plot

def main():
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # figure directory
    figDir = Path('../figures')
    
    # create output directories if not existing
    if not figDir.exists():
        os.makedirs(figDir)
    
    multi_bar_plot ([[1,1,1], [1,1,1]],
                    #xlabel = 'RSA',
                    #ylabel = 'Difference in fraction\nfrom all residues (%)',
                    #xlabels = xticklabels,
                    colors = ['magenta', 'turquoise'],
                    hatches = ['//', '..'],
                    barwidth = 0.4,
                    fontsize = 22,
                    ylim = [0, 3],
                    #xticks = xticks,
                    #opacity = None,
                    leg = ('Pathogenic mutations', 'Non-pathogenic mutations'),
                    overlap = False,
                    show = False,
                    figdir = figDir,
                    figname = 'legend')

if __name__ == "__main__":
    main()