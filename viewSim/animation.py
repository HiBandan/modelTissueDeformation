"""
Written by bandan chakrabortty @IBDM
"""
import os
import userModule
import numpy as np
from pylab import Line2D 
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon as PolyPatch

# input-path
dataPath =  './../simulationData'
if os.path.exists(dataPath):
    # output-patrh
    moviePath = dataPath + '/' + 'plots' 
    userModule.recreateDirectory(moviePath)
    #*********************************************************#
    #             reference-embryo-shape-information          #
    #*********************************************************#
    fileName = open(dataPath + '/refFrame', 'r')
    embryoInfo,cellVerticesInfo = fileName.readlines()
    fileName.close()
    # embryo-information
    embryoGeometryInfo = [item for item in embryoInfo.rstrip('\n').split(',')]
    epitheliumType = embryoGeometryInfo[0]
    semi_a,semi_b,radius,vitellineWidth,vitellineReferenceLine,epithelialReferenceLine,epithelialHeight,L = [float(cor) for cor in embryoGeometryInfo[1:]]
    # vertex-information
    coordinatesAllVertexType = cellVerticesInfo.rstrip('\n').split('#')
    coordinatesPerVertexType = [np.array_split(np.array([float(coordinate) for coordinate in indexCoordinate.split(',')]),3) for indexCoordinate in coordinatesAllVertexType[1].split('@')[:-1]]
    vitellineCoordinates_ref = [vtype[0] for vtype in  coordinatesPerVertexType]
    vitellineCoordinates_ref = np.reshape(vitellineCoordinates_ref,(-1, 2)) 
    apicalCoordinates_ref = [vtype[1] for vtype in  coordinatesPerVertexType]
    apicalCoordinates_ref = np.reshape(apicalCoordinates_ref,(-1, 2)) 
    vertexNumber = len(apicalCoordinates_ref[:,0])
    #*********************************************************#
    #           time-series-embryo-shape-information          #
    #*********************************************************#
    for frameNumber,frameName in enumerate(['initialFrame','finalFrame']):
        fig, (ax,ax_y) = plt.subplots(2,1,figsize = (4, 2))
        fileName = open(dataPath + '/' + frameName, 'r')
        embryoInfo,cellVerticesInfo = fileName.readlines()
        fileName.close()
        ###########################
        # free-vertex-information #
        ###########################
        coordinatesAllVertexType = cellVerticesInfo.rstrip('\n').split('#') 
        coordinatesFixed = [int(node) for node in coordinatesAllVertexType[0].split('@')[:-1]]
        coordinatesPerVertexType = [np.array_split(np.array([float(coordinate) for coordinate in indexCoordinate.split(',')]),3) for indexCoordinate in coordinatesAllVertexType[1].split('@')[:-1]]
        # vitelline-vertices 
        vitellineCoordinates = [vtype[0] for vtype in  coordinatesPerVertexType]
        vitellineCoordinates = np.reshape(vitellineCoordinates,(-1, 2)) 
        xv = vitellineCoordinates[:,0]
        yv = vitellineCoordinates[:,1]
        # apical-vertices 
        apicalCoordinates = [vtype[1] for vtype in  coordinatesPerVertexType]
        apicalCoordinates = np.reshape(apicalCoordinates,(-1, 2)) 
        xa = apicalCoordinates[:,0]
        ya = apicalCoordinates[:,1]
        # basal-vertices 
        basalCoordinates = [vtype[2] for vtype in  coordinatesPerVertexType]
        basalCoordinates = np.reshape(basalCoordinates,(-1, 2))
        xb = basalCoordinates[:,0]
        yb = basalCoordinates[:,1]
        ####################
        # edge-information #
        ####################
        # apical/basal-edges 
        edgesApicalBasal = [[int(node) for node in Nodes.split(',')] for Nodes in coordinatesAllVertexType[2].split('@')[:-1]]
        edgesApicalBasal = np.reshape(edgesApicalBasal,(-1, 2))
        for edge in edgesApicalBasal:
            ax.add_line(Line2D([xa[edge[0]], xa[edge[1]]], [ya[edge[0]], ya[edge[1]]], color='k',lw = 0.2))
            ax.add_line(Line2D([xb[edge[0]], xb[edge[1]]], [yb[edge[0]], yb[edge[1]]], color='k',lw = 0.2))
            ax_y.add_line(Line2D([xa[edge[0]], xa[edge[1]]], [ya[edge[0]], ya[edge[1]]], color='k',lw = 0.2))
        # lateral-edges 
        lateralEdgeCoordinates = [None]*(2*vertexNumber)
        lateralEdgeCoordinates[::2] = apicalCoordinates
        lateralEdgeCoordinates[1::2] = basalCoordinates
        lateralEdgeCoordinates = np.reshape(lateralEdgeCoordinates,(-1, 2))
        xl = lateralEdgeCoordinates[:,0]
        yl = lateralEdgeCoordinates[:,1]
        edgesLateral = [(i*2,i*2+1) for i in range(len(apicalCoordinates))]
        edgesLateral = np.reshape(edgesLateral, (-1,2))
        for edge in edgesLateral:
            ax.add_line(Line2D([xl[edge[0]], xl[edge[1]]], [yl[edge[0]], yl[edge[1]]], color='k',lw = 0.2))
        ######################
        # myosin-information #
        ######################
        edgesMyo = [[int(edge) for edge in Edges.split(',')] for Edges in coordinatesAllVertexType[3].split('@')[:-1]]
        if edgesMyo:
            edgesMyo = np.reshape(edgesMyo,(-1, 2))
            for edge in edgesMyo:
                cellsAsPolygons   = []   
                va_1 = [xa[edge[0]],ya[edge[0]]]
                va_2 = [xa[edge[1]],ya[edge[1]]]
                vb_1 = [xb[edge[0]],yb[edge[0]]]
                vb_2 = [xb[edge[1]],yb[edge[1]]]
                vm_1 = [0.5*(e1 + e2) for e1, e2 in zip(va_1, vb_1)] 
                vm_2 = [0.5*(e1 + e2) for e1, e2 in zip(va_2, vb_2)]
                apicalVerts = [va_1,va_2]
                midVerts = [vm_1,vm_2]
                cellsAsPolygons = [[apicalVerts[i],apicalVerts[i+1],midVerts[i+1],midVerts[i]] for i in range(len(apicalVerts)-1)] 
                collection = PatchCollection([PolyPatch(cellsAsPolygons[i],closed=True) for i in range(len(cellsAsPolygons))])
                collection.set_color(['g' for i in range(len(cellsAsPolygons))])
                ax.add_collection(collection) 
        ############################
        # fixed-vertex-information #
        ############################
        xa_fixed = [xa[node] for node in coordinatesFixed]
        ya_fixed = [ya[node] for node in coordinatesFixed]
        ax.plot(xa_fixed, ya_fixed, 'ro', ms = 1.0)
        ax_y.plot(xa_fixed, ya_fixed, 'ro', ms = 1.0)
        refTrackNodes = None
        ###################
        # flat-epithelium #
        ###################
        if(epitheliumType == 'flat'): 
            refTrackNodes = [0]
            vitellineReferenceLine_x = np.linspace(-0.1*semi_a+epithelialReferenceLine,0.1*semi_a+epithelialReferenceLine+1.0*L,50)
            vitellineMembraneOuter_y  = [vitellineReferenceLine]*len(vitellineReferenceLine_x)
            vitellineMembraneInner_y = [vout - vitellineWidth for vout in vitellineMembraneOuter_y]
            # reference-line-for-x-displacement
            ax.plot([-L/2,-L/2],[-epithelialHeight,vitellineReferenceLine-vitellineWidth],  c = 'k', ms = 2.0, ls = '--', lw = 0.3)
            ax_y.plot([-L/2,-L/2],[1.0,vitellineReferenceLine-vitellineWidth],  c = 'k', ms = 2.0, ls = '--', lw = 0.3)
            # reference-line-for-tissue-centre
            ax.plot(vitellineReferenceLine_x, vitellineMembraneOuter_y, color = 'm', ls = '--', lw = 0.3)
            ax.plot(vitellineReferenceLine_x, vitellineMembraneInner_y, color = 'm', ls = '--', lw = 0.3)
            ax_y.plot(vitellineReferenceLine_x, vitellineMembraneOuter_y, c = 'm', ms = 2.0, ls = '--', lw = 0.5)
            ax_y.plot(vitellineReferenceLine_x, vitellineMembraneInner_y, c = 'm', ms = 2.0, ls = '--', lw = 0.5)
        #####################
        # curved-epithelium #
        #####################
        else:
            refTrackNodes = [int(vertexNumber/4)]
            refTrackCoordinates = [[apicalCoordinates_ref[:,0][node],apicalCoordinates_ref[:,1][node]] for node in refTrackNodes]
            xv = vitellineCoordinates_ref[:,0]
            yv = vitellineCoordinates_ref[:,1]
            vit_x = np.append(xv,xv[0])
            vit_y = np.append(yv,yv[0])
            ax.add_patch(Ellipse((0, 0),width=2*semi_a,height=2*semi_b,fill=False,color= 'm',ls='--',lw=0.3))
            ax_y.add_patch(Ellipse((0, 0),width=2*semi_a,height=2*semi_b,fill=False,color= 'm',ls='--',lw=0.3))
            ax_y.axis('equal')
        # axes
        ax.axis('equal')  
        ax.axis('off')
        ax_y.axis('off')
        # save-figure
        figName = '/F_' + str(frameNumber) + '.png' 
        fig.savefig(moviePath + figName,figsize = (10, 3), dpi = 500, transparent=False)
        plt.close(fig)     
else:
    print ("no input data found !")
                            
