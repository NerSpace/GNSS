
import sys, os
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FuncFormatter
import numpy as np
import conda
CondaFileDir = conda.__file__
CondaDir = CondaFileDir.split('lib')[0]
ProjLib = os.path.join(os.path.join(CondaDir, 'share'), 'proj')
os.environ["PROJ_LIB"] = ProjLib
from mpl_toolkits.basemap import Basemap

import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)

def createFigure(PlotConf):
    try:
        fig, ax = plt.subplots(1, 1, figsize = PlotConf["FigSize"])
    
    except:
        fig, ax = plt.subplots(1, 1)

    return fig, ax

def saveFigure(fig, Path):
    Dir = os.path.dirname(Path)
    try:
        os.makedirs(Dir)
    except: pass
    fig.savefig(Path, dpi=150., bbox_inches='tight')

def prepareAxis(PlotConf, ax):
    def format_func(value, tick_number):
        return f'{value:.3f}'


    for key in PlotConf:
        if key == "Title":
            ax.set_title(PlotConf["Title"])

        for axis in ["x", "y"]:
            if axis == "x":
                if key == axis + "Label":
                    ax.set_xlabel(PlotConf[axis + "Label"])

                if key == axis + "Ticks":
                    ax.set_xticks(PlotConf[axis + "Ticks"])

                if key == axis + "TicksLabels":
                    ax.set_xticklabels(PlotConf[axis + "TicksLabels"])
                
                if key == axis + "Lim":
                    ax.set_xlim(PlotConf[axis + "Lim"])

            if axis == "y":
                if key == axis + "Label":
                    ax.set_ylabel(PlotConf[axis + "Label"])

                if key == axis + "Ticks":
                    ax.set_yticks(PlotConf[axis + "Ticks"])

                if key == axis + "TicksLabels":
                    if "FontSize" in PlotConf:
                        ax.set_yticklabels(PlotConf[axis + "TicksLabels"], fontsize = PlotConf["FontSize"])
                    else:
                        ax.set_yticklabels(PlotConf[axis + "TicksLabels"])
                if key == axis + "Lim":
                    ax.set_ylim(PlotConf[axis + "Lim"])

                if "NumberFormat" in PlotConf and PlotConf["NumberFormat"] == 1:
                    ax.yaxis.set_major_formatter(FuncFormatter(format_func))

        if key == "Grid" and PlotConf[key] == True:
            ax.grid(linestyle='--', linewidth=0.5, which='both')

def prepareColorBar(PlotConf, ax, Values):
    try:
        Min = PlotConf["ColorBarMin"]
    except:
        Mins = []
        for v in Values.values():
            Mins.append(min(v))
        Min = min(Mins)
    try:
        Max = PlotConf["ColorBarMax"]
    except:
        Maxs = []
        for v in Values.values():
            Maxs.append(max(v))
        Max = max(Maxs)

    divider = make_axes_locatable(ax)
    color_ax = divider.append_axes("right", size="3%", pad="2%")
    cmap = mpl.cm.get_cmap(PlotConf["ColorBar"])

    if "DiscreteColorBar" in PlotConf:
        bounds = np.linspace(Min, Max, PlotConf["DiscreteColorBar"]+1)
        normalize = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=256)
        cbar = mpl.colorbar.ColorbarBase(color_ax, 
                                         cmap=cmap, 
                                         norm=normalize,
                                         label=PlotConf["ColorBarLabel"],
                                         boundaries=bounds,
                                         ticks=bounds)
    else:
        normalize = mpl.cm.colors.Normalize(vmin=Min, vmax=Max)
        cbar = mpl.colorbar.ColorbarBase(color_ax, 
        cmap=cmap,
        norm=normalize,
        label=PlotConf["ColorBarLabel"])

    if "zTicks" in PlotConf:
        cbar.set_ticks(np.arange(Min, Max+1, 1))

    return normalize, cmap

def drawMap(PlotConf, ax,):
    Map = Basemap(projection = 'cyl',
    llcrnrlat  = PlotConf["LatMin"]-0,
    urcrnrlat  = PlotConf["LatMax"]+0,
    llcrnrlon  = PlotConf["LonMin"]-0,
    urcrnrlon  = PlotConf["LonMax"]+0,
    lat_ts     = 10,
    resolution = 'l',
    ax         = ax)

    # Draw map meridians
    Map.drawmeridians(
    np.arange(PlotConf["LonMin"],PlotConf["LonMax"]+1,PlotConf["LonStep"]),
    labels = [0,0,0,1],
    fontsize = 6,
    linewidth=0.2)
        
    # Draw map parallels
    Map.drawparallels(
    np.arange(PlotConf["LatMin"],PlotConf["LatMax"]+1,PlotConf["LatStep"]),
    labels = [1,0,0,0],
    fontsize = 6,
    linewidth=0.2)

    # Draw coastlines
    Map.drawcoastlines(linewidth=0.5)

    # Draw countries
    Map.drawcountries(linewidth=0.25)

def generateLinesPlot(PlotConf):
    fig, ax = createFigure(PlotConf)

    prepareAxis(PlotConf, ax)

    for key in PlotConf:
        if key == "LineWidth":
            LineWidth = PlotConf["LineWidth"]
        if key == "ColorBar":
            normalize, cmap = prepareColorBar(PlotConf, ax, PlotConf["zData"])
        if key == "Map" and PlotConf[key] == True:
            drawMap(PlotConf, ax)

    for Label in PlotConf["yData"].keys():
        LineColor = PlotConf.get("LineColor", {}).get(Label, None) 
        LegendLabel = PlotConf.get("LegendLabels", {}).get(Label, Label)

        if "ColorBar" in PlotConf:
            ax.scatter(PlotConf["xData"][Label], PlotConf["yData"][Label], 
            marker = PlotConf["Marker"],
            s = PlotConf["MarkerSize"],
            linewidth = LineWidth,
            c = cmap(normalize(np.array(PlotConf["zData"][Label]))))

        else:
            if "LineColor" in PlotConf: 
                ax.plot(PlotConf["xData"][Label], PlotConf["yData"][Label],
                        PlotConf["Marker"],
                        linewidth=LineWidth,
                        color=LineColor,
                        label=LegendLabel)
            else: 
                ax.plot(PlotConf["xData"][Label], PlotConf["yData"][Label],
                        PlotConf["Marker"],
                        linewidth=LineWidth,
                        label=LegendLabel)

    if "Legend" in PlotConf and PlotConf["Legend"] == True:
        ax.legend() 

    # For adding the labels to the markers
    if "TickName" in PlotConf:
        for i in range(len(PlotConf["xData"][0])):
            color = cmap(normalize(np.array(PlotConf["zData"][0][i])))
            if i == 0:
                ax.annotate(PlotConf["TickName"][i], (float(PlotConf["xData"][0][i]), float(PlotConf["yData"][0][i])),
                            textcoords="offset points", xytext=(0, -10), ha='center', fontsize = 6, color = color)
            else:
                if (PlotConf["xData"][0][i]-PlotConf["xData"][0][i-1]) < 20 and PlotConf["yData"][0][i] == PlotConf["yData"][0][i-1] and PlotConf["zData"][0][i] == PlotConf["zData"][0][i-1]:
                    ax.annotate(" ", (float(PlotConf["xData"][0][i]), float(PlotConf["yData"][0][i])),
                            textcoords="offset points", xytext=(0, 5), ha='center', fontsize = 6, color = color)
                else:
                    if i % 2 == 0:
                        ax.annotate(PlotConf["TickName"][i], (float(PlotConf["xData"][0][i]), float(PlotConf["yData"][0][i])),
                            textcoords="offset points", xytext=(0, 5), ha='center', fontsize = 6, color = color)
                    else:
                        ax.annotate(PlotConf["TickName"][i], (float(PlotConf["xData"][0][i]), float(PlotConf["yData"][0][i])),
                            textcoords="offset points", xytext=(0, -10), ha='center', fontsize = 6, color = color)
    
    saveFigure(fig, PlotConf["Path"])

# Function for the plot with the status 0 in gray
def generateLinesPlotStatus(PlotConf):
    fig, ax = createFigure(PlotConf)
    prepareAxis(PlotConf, ax)

    if "LineWidth" in PlotConf:
        LineWidth = PlotConf["LineWidth"]

    if "status_1" in PlotConf["yData"]:
        normalize, cmap = prepareColorBar(PlotConf, ax, PlotConf["zData"]["status_1"])

    for Label in PlotConf["yData"].keys():
        LegendLabel = PlotConf.get("LegendLabels", {}).get(Label, Label)

        if Label == "status_1":
            ax.scatter(PlotConf["xData"][Label], PlotConf["yData"][Label], 
                       c=cmap(normalize(np.array(PlotConf["zData"][Label]))), 
                       marker=PlotConf["Marker"], linewidth=LineWidth, label=LegendLabel, zorder = PlotConf["Zorder"][Label], s = PlotConf["MarkerSize"])
            
        else:
            ax.scatter(PlotConf["xData"][Label], PlotConf["yData"][Label], 
                       color='gray', marker=PlotConf["Marker"], linewidth=LineWidth, label=LegendLabel, zorder = PlotConf["Zorder"][Label], s = PlotConf["MarkerSize"], alpha = 0.2)

    if "Legend" in PlotConf and PlotConf["Legend"]:
        ax.legend()

    saveFigure(fig, PlotConf["Path"])

# Plot with two y-axis

def generateDualYAxisPlot(PlotConf):
    fig, ax1 = createFigure(PlotConf)

    prepareAxis(PlotConf, ax1)

    LineWidth = PlotConf.get("LineWidth", 1.0)
    Marker = PlotConf["Marker"]

    # Plot left Y-axis data
    for Label in PlotConf["yDataLeft"].keys():
        LineColor = PlotConf.get("LineColor", {}).get(Label, None)
        LegendLabel = PlotConf.get("LegendLabels", {}).get(Label, Label)
        
        ax1.plot(PlotConf["xData"][Label], PlotConf["yDataLeft"][Label],
                 Marker,
                 linewidth=LineWidth,
                 color=LineColor,
                 label=LegendLabel)
        
    ax1.set_ylabel(PlotConf["yLabelLeft"])
    ax1.tick_params(axis='y')

    # Create a second y-axis
    ax2 = ax1.twinx()

    # Plot right Y-axis data
    for Label in PlotConf["yDataRight"].keys():
        LineColor = PlotConf.get("LineColor", {}).get(Label, None)
        LegendLabel = PlotConf.get("LegendLabels", {}).get(Label, Label)
        
        ax2.plot(PlotConf["xData"][Label], PlotConf["yDataRight"][Label],
                 Marker,
                 linewidth=LineWidth,
                 color=LineColor,
                 label=LegendLabel)
        
    ax2.set_ylabel(PlotConf["yLabelRight"])
    ax2.tick_params(axis='y')
    ax2.set_ylim(PlotConf["yLimR"])
    ax2.set_yticks(PlotConf["yTicksR"])


    # Add legend
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(handles1 + handles2, labels1 + labels2, loc='upper left')

    saveFigure(fig, PlotConf["Path"])


# Define the kind of plot used

def generatePlot(PlotConf):
    if(PlotConf["Type"] == "Lines"):
        generateLinesPlot(PlotConf)
    if(PlotConf["Type"] == "LinesStatus"):
        generateLinesPlotStatus(PlotConf)
    if(PlotConf["Type"] == "DualYAxis"):
        generateDualYAxisPlot(PlotConf)