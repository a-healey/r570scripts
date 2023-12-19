#get the stuff you need
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import pandas as pd
from matplotlib import collections  as mc
import matplotlib as mpl
#setting teh fonttype to 42 is useful as it uses a font in the pdf that is editable in Adobe Illustrator
mpl.rcParams['pdf.fonttype'] = 42
import seaborn as sns
plt.style.use("dark_background")
from sys import argv
#======================================================
def real_main():
    script, coords, refOfInterest, queryOfInterest, minLen = argv
    plotDict = {}
    fileName = coords.split(".")[0]
    keySet = set()
    minLen = int(minLen)
    fig = plt.figure(1, figsize=(10, 8))
    with open(coords) as f:
        for _ in range(5):
            next(f)
        for line in f:
            s = line.split(None)
            start1 = int(s[0])
            end1 = int(s[1])
            start2 = int(s[3])
            end2 = int(s[4])
            len1 = int(s[6])
            len2 = int(s[7])
            if len1 < minLen or len2 < minLen: continue
            ID = float(s[9])
            ref = s[14]
            qry = s[15]
            if ref ==  refOfInterest and qry == queryOfInterest:
                key1 = "_".join((str(start1),str(start2),str(end1),str(end2)))
                key2 = "_".join((str(start2),str(start1),str(end2),str(end1)))
                plotDict[key1] = [start1,start2,ID,key1]
                plotDict[key2] = [end1,end2,ID,key1]
                keySet.add(key1)
        #####
    #####
#####
    #plotter makes a scatter plot of alignment starts and stops,
    #then draws a dotted line between them
    #color palette used is based on percent identity from the mummmer alignment
    plotDF = pd.DataFrame.from_dict(plotDict, orient ='index')
    #df.sort_index(inplace=True)
    plotDF= plotDF.reset_index()
    plotDF.head()
    plotDF.columns=['index','ref', 'qry', 'ID','key']
    plotDF.head()
    print(len(keySet))
    len(plotDF['ID'].unique())
    palette = sns.color_palette("viridis",len(plotDF['ID'].unique()))
    points = plt.scatter(data=plotDF, x = 'ref', y= 'qry',
                         c='ID', s=20, cmap="viridis",)
    plt.colorbar(points,label='Percent Identity')
    sns.lineplot(data=plotDF, 
                 x="ref", 
                 y="qry", 
                 hue="ID",
                 units="key",
                 estimator=None,lw=1, linestyle='--',
                palette=palette, legend = False)
                #####
            #####
        #####
    #####
    plt.title("%s" % (fileName))
    plt.xlabel("Reference:%s" % (refOfInterest))
    plt.ylabel("Query:%s" % (queryOfInterest))
    #plt.show()
    plt.savefig('%s.%s_%s.heatmapPlot.L_%d_bp.pdf' % (fileName,refOfInterest,queryOfInterest, minLen))
#======================================================
if ( __name__ == '__main__' ):
    real_main()