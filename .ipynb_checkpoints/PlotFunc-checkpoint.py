import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cm
import scipy
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import matplotlib.ticker as ticker
from brokenaxes import brokenaxes



def CanvasStyle(ax, x_min=0, y_min=0):
    ax.patch.set_facecolor('white')
    ax.grid(False)
    #plt.tight_layout()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.axhline(linewidth=3, y=y_min, color='k')
    ax.axvline(linewidth=3, x=x_min, color='k')
    ax.yaxis.set_tick_params(width=2)
    ax.xaxis.set_tick_params(width=2)
    return ax

# label p value and significance.
def significance(arrs, ax, columns=[], mode='box'):
    
    def stars(p):
        if p < 0.0001:
            return "****"
        elif (p < 0.001):
            return "***"
        elif (p < 0.01):
            return "**"
        elif (p < 0.05):
            return "*"
        else:
            return "n.s."
    
    def ttest(arr1, arr2):
        # Calculate t-test and p value.
        # Use scipy.stats.ttest_ind.
        t, p = scipy.stats.ttest_ind(arr1, arr2)
        s = stars(p)
        print("ttest_ind:            t = %g  p = %g" % (t, p))
        return s
    
    trans = ax.get_xaxis_transform()
    label_min, label_max = ax.get_ylim()
    props = {'connectionstyle':"bar,fraction={0}".format(0.1),
             'arrowstyle':'-',
             'linewidth':2,
             'ec':'#000000'}
    rank = np.array([abs(v[0]-v[1]) for v in columns])
    cols = np.array([columns[ind] for ind in np.argsort(rank)])
    overlap_record = np.zeros(len(arrs))
    if mode == 'bar':
        x_position = np.arange(len(arrs))/2+0.25
        y_max = []
        norm = []
        for i in range(len(arrs)):
            y = (np.mean(arrs[i])+np.std(arrs[i]))
            norm.append(max(arrs[i]))
            y_max.append(y)
        y_max = 1#max(y_max)/max(norm)
        standard = 0.025*(len(arrs)-1)
        
    else:
        x_position = np.arange(len(arrs))+0.75+0.25
        y_max = []
        DataMax = []
        for i in range(len(arrs)):
            DataMax.append(np.percentile(arrs[i], 98))
            Q1 = np.percentile(arrs[i], 75)
            Q3 = np.percentile(arrs[i], 25)
            max_bound = Q1+1.5*(Q1-Q3)
            y = np.sort(arrs[i], axis=None)[np.nonzero(np.sort(arrs[i], axis=None)<=max_bound)[0][-1]]
            y_max.append(y)
        if mode == "box":
            y_max = max(y_max)
            standard = 0.05*(len(arrs)-1)
        else:
            y_max = max(y_max)
            standard = 0.05*2
        
    for col in cols:
        # Calculate t-test and p value.
        s = ttest(arrs[col[0]], arrs[col[1]])
        # Update props.
        props['connectionstyle'] = "bar,fraction={0}".format(
            standard/(abs(x_position[col[0]]-x_position[col[1]])))
        # Process the label.
        passby = np.arange(col[0], col[1]+1)
        adj = (y_max-label_min)/(label_max-label_min)
        adj = adj+0.04*(1+np.max(overlap_record[passby]))
        overlap_record[passby] +=1
        ax.annotate("", xy=(x_position[col[0]], adj),
                    xycoords=trans,
                    xytext=(x_position[col[1]], adj),
                    textcoords=trans,
                    arrowprops=props)
        if s=="n.s.":
            text_adj = y_max+(label_max-label_min)*0.04*(
            1+np.max(overlap_record[passby]))+0.005
        else:
            text_adj = y_max+(label_max-label_min)*0.04*(
                1+np.max(overlap_record[passby]))-(label_max-label_min)*0.01
        ax.text((x_position[col[0]]+x_position[col[1]])/2,
                text_adj,
                s,
                horizontalalignment='center',
                verticalalignment='center',
                weight='bold',)
#               backgroundcolor='white')
    return ax

class PlotMethods():
    
    def __init__(self, path, columns=[]):
        plt.rcParams['font.weight'] = 'normal'#'bold'
        plt.rcParams['font.size'] = 12
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['axes.labelsize'] = 14
        plt.rcParams['axes.facecolor'] = 'white'
        plt.rcParams['axes.grid'] = False
        plt.style.use('seaborn-ticks')
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42
        self.save_format = ['png', 'svg', 'pdf']
        self.columns = columns
        self.lw = 3
        self.errbarattr = {'lw':3,
                           'capthick':3,
                           'capsize':10,
                           'ecolor':'black'
                          }
        self.ylabelDict = {'PULScore':"Colocalization Ratio",
                           'LRratio':"Ratio of PopZ pole to opposite pole",
                           'LRMax':"Ratio of PopZ pole to opposite pole",
                           'PolarityScore':"Polarity score",
                           'PeakScore':"Peak score",
                           'Pearson':"Pearson Coefficient",
                           'Bowman':"Localization index",
                           'TotalIntensity':"Normalized Flourescence Intensity",
                           'PolarIntensity':"Flourescent intensity (A.U.)",
                           'CellLen':"Cell length",
                           'Others':"Flourescent intensity (A.U.)"
                          }
        self.path = path
        
    def ProfileStatPlot(self, avg, std, names):
        # Draw the average values with standard deviation of profiles.
        # Open a empty figure for new plot.
        x_scale = np.linspace(0, 1, len(avg))
        fig = plt.figure(figsize=(8, 3))
        ax = fig.add_subplot(111)
        # Draw the plot.
        ax.plot(x_scale, avg, 'g', label='{}'.format(names), linewidth=3)
        # Draw the error bar with a filled region.
        ax.fill_between(x_scale,
                        avg-std,
                        avg+std,
                        alpha=0.2,
                        facecolor='g',
                        edgecolor='g',
                        linewidth=0.0)
        # Change some arributes of the plot.
        ax.set_xlim([0, 1])
        ax.set_ylim(0, np.max(std+avg))
        ax.set_xlabel('Relative Position of Cell')
        ax.set_ylabel('Normalized Fluorescence Intensity')
        ax.legend()
        ax = CanvasStyle(ax)
        for figformat in self.save_format:
            fig.savefig(r'{0}/{1}_{2}.{3}'.format(
                self.path, names, 'ProfileStatPlot', figformat), format='{}'.format(figformat))
#         fig.savefig('{0}\{1}_{2}.png'.format(self.path, names, 'ProfileStatPlot'))
        plt.show()
        
    def ProfileGroupBar(self, avg_arr, std_arr, names_arr):
        # Grouped bar plot against relative positions.
        # Open a new figure.
        # Specifically for comparing constitutive DivIVA to NSCCDivIVA.
        fig = plt.figure(figsize=(8, 3))
        ax = fig.add_subplot(111)
        index = np.arange(len(avg_arr[0]))
        bar_width = 0.35

        rects1 = ax.bar(index/np.max(index), avg_arr[0], bar_width,
                         color='orange',
                         label=names_arr[0])

        rects2 = ax.bar((index + bar_width)/np.max(index + bar_width), avg_arr[1], bar_width,
                         color='skyblue',
                         label=names_arr[1])

        ax.set_xlabel('Relative Position of Cell')
        ax.set_ylabel('Fluorescence Intensity')
        plt.legend()

        ax = CanvasStyle(ax)
        ax.set_xlim(np.min(index)/np.max(index), 1)
        fig.savefig('{0}\{1}_{2}.png'.format(self.path, names_arr[0], 'ProfileGroupBar'))
        plt.show()
        
    
    def ProfilePlot(self, profiles, label, mean_line=[], norm=1):
        n = len(profiles)
        col = plt.cm.Greys(np.linspace(0.95,0.45,n))
        fig = plt.figure(figsize=(8, 3))
        ax = fig.add_subplot(111)
        y_upper = 0
#         norm = 13.791527876117886
        for line, c in zip(profiles, col): # 4.0420503680338
            x_scale = np.linspace(0, 1, len(line))#2.270227918688883 #13.3802942173593
            ax.plot(x_scale, [v/norm for v in line], color=c, linewidth=2)
            y_upper = np.max(line) if y_upper<np.max(line) else y_upper
        if len(mean_line)>0:
            x_scale = np.linspace(0, 1, len(mean_line))
            ax.plot(x_scale, [v/norm for v in mean_line], color='r', linewidth=4, label='{}'.format(label))
        print('info', y_upper)
        ax.set_xlabel('Relative Position of Cell')
        ax.set_ylabel('Fluorescence Intensity (A.U.)')
        if np.any(mean_line)==0:
            ax.set_ylim([1, -1])
        else:
            ax.set_ylim([0, 1])#y_upper
        ax.set_xlim([0, 1])
        ax.legend()
        ax = CanvasStyle(ax)
        ax.grid(False)
        for figformat in self.save_format:
            fig.savefig(r'{0}/{1}_{2}.{3}'.format(
                self.path, label, 'SingleTraces', figformat), format='{}'.format(figformat))
        plt.show()
        
    def TwoProfilePlot(self, profile1_arr, profile2_arr, names):
        fig = plt.figure(figsize=(8, 3))
        ax = fig.add_subplot(111)
        x_scale = np.linspace(0, 1, len(profile1_arr[0]))
        max1 = (profile1_arr[0]).max()
        max2 = (profile2_arr[0]).max()
        maxmax = np.max([max1, max2])
        line1 = ax.plot(x_scale, profile1_arr[0]/maxmax,
                        color='tomato', linewidth=3, label='{}'.format(names[0]))
        ax.fill_between(x_scale, (profile1_arr[0]-profile1_arr[1])/maxmax, (profile1_arr[0]+profile1_arr[1])/maxmax, alpha=0.2, facecolor='tomato', edgecolor='tomato', linewidth=0.0)

        x_scale = np.linspace(0, 1, len(profile2_arr[0]))
        line2 = ax.plot(x_scale, profile2_arr[0]/maxmax,
                        color='grey', linewidth=3, label='{}'.format(names[1]))
        ax.fill_between(x_scale, (profile2_arr[0]-profile2_arr[1])/maxmax, (profile2_arr[0]+profile2_arr[1])/maxmax, alpha=0.2, facecolor='grey', edgecolor='grey', linewidth=0.0)
        
        ax = CanvasStyle(ax)
        ax.set_xlabel('Relative Position of Cell')
        ax.set_ylabel('Normalized Intensity')
        plt.legend()
        a = np.max((profile1_arr[0]+profile1_arr[1])/maxmax)
        b = np.max((profile2_arr[0]+profile2_arr[1])/maxmax)
        ax.set_ylim([0, np.max([a, b])])
        ax.set_xlim([0, 1])
        plt.tight_layout()
        for figformat in self.save_format:
            fig.savefig(r'{0}/{1}.{2}'.format(
                self.path, 'TwoProfileComparison', figformat), format='{}'.format(figformat))
        plt.show()
            
            
    def BarPlot(self, avg_arr, std_arr, names_arr, FigName='Others', arrs=False):
        n = len(avg_arr)
        col = plt.cm.Greys(np.linspace(0.95,0.25,n))
        width = 0.25  # the width of the bars
        ticks = np.arange(len(names_arr))/2+width
        fig = plt.figure(figsize=(4, 5))
        ax = fig.add_subplot(111)
        std_neg = np.zeros(len(std_arr))
        stds = np.array([std_neg, std_arr])
        bars = ax.bar(np.arange(len(avg_arr))/2+width,
                      avg_arr,
                      width,
                      color=col,
                      yerr=stds,
                      edgecolor="black",
                      lw=self.lw,
                      error_kw=self.errbarattr)
        if FigName == 'LRratio' or FigName == 'LRMax':
            ax.plot(np.arange(n), np.ones(n), color='grey', linestyle='--', linewidth=2)
        # add some text for labels, title and axes ticks
        ax.set_xlim([0, np.arange(n)[-1]/2+2*width])
        ax.set_ylim([0, 1.1*np.max(np.array(avg_arr)+np.array(std_arr))])
        ax.set_ylabel("{}".format(self.ylabelDict[FigName]))
        ax.set_xticks(ticks)
        ax.set_xticklabels((names_arr))
        ax = CanvasStyle(ax)
        ax = significance(arrs,
                     ax,
                     columns=self.columns,
                     mode='bar')
        
        for figformat in self.save_format:
            fig.savefig(r'{0}/{1}_{2}.{3}'.format(
                self.path, 'Promoter_Strength', 'BoxPlot', figformat), format='{}'.format(figformat))
        plt.show()
        
        
    def BoxPlot(self, arrs, names_arr, FigName='Others'):
        # Parameters
        width = 0.75
        x_position = np.arange(len(names_arr))+width+0.25
        x_flank = np.arange(len(names_arr)+1)+width-0.25
        # Open a figure
        figbox = plt.figure(figsize=(5, 5))
        ax = figbox.add_subplot(111)
        # cmap
        colors = plt.cm.bone(np.linspace(0.45, 0.85, len(names_arr)))
        # props
        capprops = dict(linestyle='-', linewidth=2.5, color='k')
        whiskerprops = dict(linestyle='-', linewidth=2.5, color='k')
        medianprops = dict(linestyle='-', linewidth=4, color='k')
        boxprops = dict(linestyle='-', linewidth=2.5, color='k')
        # rectangular box plot.
        
        global_max = []
        for arr in arrs:
            arr = np.array(arr)
            Q1 = np.percentile(arr, 75)
            Q3 = np.percentile(arr, 25)
            max_bound = Q1+1.5*(Q1-Q3)
            min_bound = Q3-1.5*(Q1-Q3)
            index1 = np.nonzero(arr<=max_bound)[0]
            tmp_arr = arr[index1]
            global_max.append(np.max(tmp_arr))
        global_max = np.max(global_max)
        arrs = np.array([np.log10(arr/global_max) for arr in arrs])
        bplot = ax.boxplot([arr for arr in arrs],
                           widths = 0.4,
                           medianprops=medianprops,
                           boxprops=boxprops,
                           whiskerprops=whiskerprops,
                           capprops=capprops,
                           showfliers=False,
                           vert=True,   # vertical box aligmnent
                           patch_artist=True)   # fill with color
        
        # fill with colors
        for patch, color in zip(bplot['boxes'], colors):
            patch.set_facecolor(color)
            
        # Add additional line
        if FigName=='LRratio':
            ax.plot(x_flank,
                    np.ones(len(x_flank)),
                    color='grey',
                    linestyle='--',
                    linewidth=2)
        
        # add some text for labels, title and axes ticks
        ax.set_ylabel("{}".format(self.ylabelDict[FigName]))
        ax.set_xticks(x_flank)
        ax.xaxis.set_major_locator(ticker.FixedLocator(x_flank))
        ax.xaxis.set_minor_locator(ticker.FixedLocator(x_position))
        ax.set_yticklabels(['$10^{{{0}}}$'.format(v) for v in ax.get_yticks()])
        ax.xaxis.set_major_formatter(ticker.NullFormatter())
        ax.xaxis.set_minor_formatter(ticker.FixedFormatter(names_arr))
        #ax.text(0.0, 0.1, fontsize=15, transform=ax.transAxes)
        for tick in ax.xaxis.get_minor_ticks():
            tick.tick1line.set_markersize(0)
            tick.tick2line.set_markersize(0)
            tick.label1.set_horizontalalignment('center')
        y_min, _ = ax.get_ylim()
        x_min, _ = ax.get_xlim()
        ax = CanvasStyle(ax, x_min=x_min, y_min=y_min)
        ax = significance(arrs,
                     ax,
                     columns=self.columns,
                     mode='box')
        plt.tight_layout()
        for figformat in self.save_format:
            figbox.savefig(r'{0}/{1}_{2}.{3}'.format(
                self.path, 'Promoter_Strength', 'BoxPlot', figformat), format='{}'.format(figformat))
        plt.show()

        
    def ViolinPlot(self, arrs, names_arr, FigName='Others', norm=False):
        # Parameters
        width = 0.75
        x_position = np.arange(len(names_arr))+width+0.25
        x_flank = np.arange(len(names_arr)+1)+width-0.25
        # Open a figure
        figbox = plt.figure(figsize=(3+1*(len(names_arr)-2), 5))
        ax = figbox.add_subplot(111)
        # cmap
        colors = plt.cm.bone(np.linspace(0.45, 0.85, len(names_arr)))
        # props
        capprops = dict(linestyle='-', linewidth=2.5, color='k')
        whiskerprops = dict(linestyle='-', linewidth=2.5, color='k')
        medianprops = dict(linestyle='-', linewidth=4, color='k')
        boxprops = dict(linestyle='-', linewidth=2.5, color='k')
        # rectangular box plot
        
        global_max = []
        for arr in arrs:
            arr = np.array(arr)
            Q1 = np.percentile(arr, 75)
            Q3 = np.percentile(arr, 25)
            max_bound = Q1+1.5*(Q1-Q3)
            min_bound = Q3-1.5*(Q1-Q3)
            index1 = np.nonzero(arr<=max_bound)[0]
            tmp_arr = arr[index1]
            global_max.append(np.max(tmp_arr))
        global_max = np.max(global_max)
        #np.array([np.log10(arr) for arr in arrs])/global_max
        if norm==True:
            data = [np.array(arr)/global_max for arr in arrs]
        else:
            data = [np.array(arr) for arr in arrs]
        bp = bplot = ax.boxplot(data,
                           widths = 0.1, # 0.4
                           medianprops=medianprops,
                           boxprops=boxprops,
                           whiskerprops=whiskerprops,
                           capprops=capprops,
                           showfliers=False,
                           vert=True,   # vertical box aligmnent
                           patch_artist=True)   # fill with color
        
        caps = bp['caps']
        for cap in caps:
            cap.set(xdata=cap.get_xdata() + (-0.1,+0.1))
        
        # fill with colors
        for patch, color in zip(bplot['boxes'], colors):
            patch.set_facecolor('grey')
        
        arr_show = []
        for arr in arrs:
            arr = np.array(arr)
            Q1 = np.percentile(arr, 75)
            Q3 = np.percentile(arr, 25)
            max_bound = Q1+1.5*(Q1-Q3)
            min_bound = Q3-1.5*(Q1-Q3)
            index1 = np.nonzero(arr<=max_bound)[0]
            tmp_arr = arr[index1]
            index2 = np.nonzero(tmp_arr>=min_bound)[0]
            arr_show.append(tmp_arr[index2])#/global_max)
        vplot = ax.violinplot(arr_show,
                              showmeans=False,
                              showmedians=False,
                              showextrema=False,
                              widths = 0.75)
        
        for prop in vplot['bodies']:
            prop.set_facecolor(colors)
            prop.set_edgecolor('white')
            prop.set_alpha(0.5)
        
        # Add additional line
        if FigName=='LRratio' or FigName == 'LRMax':
            ax.plot(x_flank,
                    np.ones(len(x_flank)),
                    color='grey',
                    linestyle='--',
                    linewidth=2)

        # add some text for labels, title and axes ticks
        ax.set_ylabel("{}".format(self.ylabelDict[FigName]))
        ax.set_xticks(x_flank)
        ax.xaxis.set_major_locator(ticker.FixedLocator(x_flank))
        ax.xaxis.set_minor_locator(ticker.FixedLocator(x_position))
        #ax.set_yticklabels(['$10^{{{0}}}$'.format(v) for v in ax.get_yticks()])
        ax.xaxis.set_major_formatter(ticker.NullFormatter())
        ax.xaxis.set_minor_formatter(ticker.FixedFormatter(names_arr))
        #ax.text(0.0, 0.1, fontsize=15, transform=ax.transAxes)
        for tick in ax.xaxis.get_minor_ticks():
            tick.tick1line.set_markersize(0)
            tick.tick2line.set_markersize(0)
            tick.label1.set_horizontalalignment('center')
        
        if FigName=='TotalIntensity':
            ax.set_ylim([0,1])
        y_min, _ = ax.get_ylim()
        x_min, _ = ax.get_xlim()
        ax = CanvasStyle(ax, x_min=x_min, y_min=y_min)
        #/global_max
        significance([np.array(arr) for arr in arrs],
                     ax,
                     columns=self.columns,
                     mode='violin')
        
        ax.grid(False)
        plt.tight_layout()
        for figformat in self.save_format:
            figbox.savefig(r'{0}/{1}_{2}.{3}'.format(
                self.path, FigName, 'ViolinPlot', figformat), format='{}'.format(figformat))
        plt.show()

    def BrokenBoxPlot(self, arrs, names_arr):
        ## plot
        figbox = plt.figure(figsize=(8, 5))
        bax = brokenaxes(xlims=((.25, .34), (.35, .65), (.66, .75)), hspace=.05)
        ## cmap
        values = np.arange(len(names_arr))
        cool = plt.get_cmap('bone') #gist_rainbow
        cNorm  = colors.Normalize(vmin=0, vmax=values[-1]*1.2)
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cool)
        # BOX plot
        capprops = dict(linestyle='-', linewidth=2.5, color='k')
        whiskerprops = dict(linestyle='--', linewidth=2.5, color='k')
        medianprops = dict(linestyle='-', linewidth=4, color='tomato')
        boxprops = dict(linestyle='-', linewidth=2.5, color='k')
        #ax_box
        box = bax.boxplot(arrs,
                          widths = 0.4,
                        medianprops=medianprops,
                                boxprops=boxprops,
                                whiskerprops=whiskerprops,
                                capprops=capprops,
                                showfliers=False,
                                vert=False,
                         patch_artist=False)   # vertical box aligmnent


        # Dots plot
        for i in range(len(names_arr)):
            colorVal = scalarMap.to_rgba(values[i])
            y = arrs[i]
            x = np.random.normal(1+i, 0.04, size=len(y))
            if i == 7:
                marker = '.'
                point = colorVal
                #point = 'r.'
            else:
                marker = '.'
                point = colorVal
                #point = 'k.'
                #ax_box
            dots = bax.plot(y, x, color = point, marker = marker, linestyle = 'None', alpha=0.8)

        
        # add a line to label the position of midcell.
        width = 0.75
        #ax_box
        bax.plot(0.5*np.ones(len(names_arr)*4),
                    np.arange(len(names_arr)*4)*0.5,
                    color = 'black',
                    linestyle = '--',
                    alpha=0.8)
        
        # Add some text for labels, title and axes ticks
        bax.axs[0].set_yticks(np.arange(len(names_arr))+width+0.25)
        bax.axs[1].set_yticks(np.arange(len(names_arr))+width+0.25)
        bax.axs[2].set_yticks(np.arange(len(names_arr))+width+0.25)
        NAME = ['', '']
        for N in names_arr:
            NAME.append(N)
        bax.axs[0].set_yticklabels((NAME))
        bax.set_ylim([0, 6])
        bax.axhline(linewidth=3, y=0, color='k')
        bax.axvline(linewidth=3, x=0.25, color='k')
        bax.axs[0].tick_params(axis = 'x', colors = 'w')
        bax.axs[1].tick_params(axis = 'x', colors = 'w')
        bax.axs[2].tick_params(axis = 'x', colors = 'w')
        bax.axs[0].yaxis.set_tick_params(width=2)
        for figformat in self.save_format:
            figbox.savefig(r'{0}/{1}_{2}.{3}'.format(
                self.path, 'Comparison', 'MassCenterPlot', figformat), format='{}'.format(figformat),)
        plt.show()
        
        
    def HORBoxPlot(self, arrs, names_arr):
        ## plot
        figbox = plt.figure(figsize=(8, 5))
        ax_box = figbox.add_subplot(111)
        ## cmap
        values = np.arange(len(names_arr))
        cool = plt.get_cmap('bone') #gist_rainbow
        cNorm  = colors.Normalize(vmin=0, vmax=values[-1]*1.2)
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cool)
        # BOX plot
        capprops = dict(linestyle='-', linewidth=2.5, color='k')
        whiskerprops = dict(linestyle='--', linewidth=2.5, color='k')
        medianprops = dict(linestyle='-', linewidth=4, color='tomato')
        boxprops = dict(linestyle='-', linewidth=2.5, color='k')
        #ax_box
        box = ax_box.boxplot(arrs,
                          widths = 0.4,
                        medianprops=medianprops,
                                boxprops=boxprops,
                                whiskerprops=whiskerprops,
                                capprops=capprops,
                                showfliers=False,
                                vert=False,
                         patch_artist=True)   # vertical box aligmnent


        for ind, patch, cap, wk, medians in zip(range(len(names_arr)),
                                                box['boxes'],
                                                box['caps'],
                                                box['whiskers'],
                                                box['medians']):
            colorVal = scalarMap.to_rgba(values[ind])
            coloredit=list(colorVal)
            coloredit[3] = 0.6
            wk.set(color = coloredit)
            coloredit[3] = 0.4
            patch.set_facecolor(coloredit)

        # Dots plot
        for i in range(len(names_arr)):
            colorVal = scalarMap.to_rgba(values[i])
            y = arrs[i]
            x = np.random.normal(1+i, 0.04, size=len(y))
            if i == 7:
                marker = '.'
                point = colorVal
            else:
                marker = '.'
                point = colorVal

            dots = ax_box.plot(y, x, color = point, marker = marker, linestyle = 'None', alpha=0.8)

        # add a line to label the position of midcell.
        width = 0.75
        #ax_box
        ax_box.plot(0.5*np.ones(len(names_arr)*2),
                    np.arange(len(names_arr)*2),
                    color = 'black',
                    linestyle = '--',
                    alpha=0.8)
        # Add some text for labels, title and axes ticks
        ax_box.set_xlabel('Relative Position')
        ax_box.set_xticks([0.25, 0.5, 0.75])
        ax_box.set_xticklabels(['PopZ pole', 'Midcell', 'Opposite pole'])
        ax_box.set_yticks(np.arange(len(names_arr))+width+0.25)

        NAME = [0]
        for N in names_arr:
            NAME.append(N)
        ax_box.set_yticklabels((names_arr))
        ax_box.set_ylim([0,6])
        ax_box.set_xlim([0.25,0.75])
        y_min, _ = ax_box.get_ylim()
        x_min, _ = ax_box.get_xlim()
        ax_box = CanvasStyle(ax_box, x_min=x_min, y_min=y_min)
        plt.tight_layout()
        for figformat in self.save_format:
            figbox.savefig(r'{0}/{1}_{2}.{3}'.format(
                self.path, 'Comparison', 'OriMassCenterPlot', figformat), format='{}'.format(figformat))
        plt.show()
        
        
    def kymograph(self, data, max_divide_t, bound):
        # Open a new plot fig.
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)
        # Show data by Heatmap.
        data = data[:10, :]
        heatmap = ax.imshow(data[:10, :], aspect='auto',
                            cmap='viridis', interpolation='gaussian', vmin=0, vmax=1)
        # Label the cell division time point.
        ax.plot(np.arange(1000), np.ones(1000)*(max_divide_t-1), color='white', linestyle='--', linewidth=3)
        # Label the boundary of the cells.
        ax.plot(bound[0][0], bound[1], color='white', linestyle='-', linewidth=3)
        ax.plot(bound[0][1], bound[1], color='white', linestyle='-', linewidth=3)
        # Limit the boundary of the axis.
        r_len, c_len = data.shape
        ax.set_xlim([0, c_len])
        ax.set_ylim([r_len-2., 0])
        # Labeling
        ax.set_xlabel('Relative cell length')
        ax.set_ylabel('Time (minutes)')
        ax.set_xticks(np.around(np.linspace(0, len(data[9, :]), 11), 2))
        ax.set_xticklabels(np.around(np.linspace(0, 1, 11), 2))
        labels = [item.get_text() for item in ax.get_yticklabels()]
        labels = np.linspace(0, 54, 10)
        ax.set_yticklabels(labels)
        plt.colorbar(heatmap)
        # Standardize the style of plotting. 
        y_min, _ = ax.get_ylim()
        x_min, _ = ax.get_xlim()
        ax = CanvasStyle(ax, x_min=x_min, y_min=y_min)
        for figformat in self.save_format:
            fig.savefig(r'{0}/{1}.{2}'.format(
                self.path, 'kymograph', figformat), format='{}'.format(figformat))
        plt.show()
        
        
    def divergent(self, mother, daughter):
        # Open a new plot figure.
        fig = plt.figure(figsize=(8, 5))
        ax = fig.add_subplot(111)
        mother_mean = np.mean(mother, axis=0)[:-1]
        mother_std = np.std(mother, axis=0)[:-1]
        daughter_mean = np.mean(daughter, axis=0)[:-1]
        daughter_std = np.std(daughter, axis=0)[:-1]
        print('len bug', daughter_mean[5:], len(np.arange(6, 12)))
        ax.plot(np.arange(0, 10), mother_mean[0:10], color='tomato', linestyle='-',
                linewidth=3, label='Mother cell')
        ax.fill_between(np.arange(10),
                        mother_mean-mother_std,
                        mother_mean+mother_std,
                        alpha=0.2,
                        facecolor='tomato',
                        edgecolor='tomato',
                        linewidth=0.0)

        ax.plot(np.arange(4, 10), daughter_mean[4:10], color='steelblue',
                linestyle='-', linewidth=3, label='Daughter cell')
        ax.fill_between(np.arange(4, 10),
                        daughter_mean[4:]-daughter_std[4:],
                        daughter_mean[4:]+daughter_std[4:],
                        alpha=0.2,
                        facecolor='steelblue',
                        edgecolor='steelblue',
                        linewidth=0.0)

        ax.set_ylabel('Normalized fluorescence Intensity')
        ax.set_xlabel('Time (minutes)')
        y_min, y_max = ax.get_ylim()
        x_min, _ = ax.get_xlim()
        y_indicator = np.arange(0, int(y_max+2))
        ax.plot(np.ones(len(y_indicator))*4,
                y_indicator,
                color='grey',
                linestyle='--',
                linewidth=3)
        ax.set_ylim([y_min, y_max])
        ax.set_xlim([0, 10])
        plt.legend()
        ax = CanvasStyle(ax, x_min=x_min, y_min=y_min)
        labels = np.linspace(0, 54, 10)
        ax.set_xticks(np.arange(10))
        ax.set_xticklabels(labels)
        for figformat in self.save_format:
            fig.savefig(r'{0}/{1}.{2}'.format(
                self.path, 'divergent', figformat), format='{}'.format(figformat))
        plt.show()