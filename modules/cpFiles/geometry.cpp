#include "drosoSymBreak.h"
#include "algebra.h"

void System:: initializeParameter()
{
    // load-parameters
    Switches sw(vary_PARAMETERS_File);
    // update-yolk-pressure
    p.pressure = laplasPressure;
    // update-myosin
    int myoEdgeIndex(0);
    myosin_Edges.clear();
    for(int l = 1; l <= sw.myoCellNumber; l++)
    {
        myoEdgeIndex = (sw.refVertex - l + numNode)%numNode;
        // increase-apical-myosin
        geometry->CELL[myoEdgeIndex].alphA = 1.0;
        // decrease-basal-myosin
        geometry->CELL[myoEdgeIndex].alphB = 1.0;
        // coloring-myosin-induced-cells
        myosin_Edges.push_back((myoEdgeIndex+1+numNode)%numNode);
    }
    myosin_Edges.push_back(myoEdgeIndex);
    // update-fixed-vertex
    p.fixedVertex = -1;
    // update-movie-switch
    p.movieFrameInterval = 0;
    p.movieFrameNumber = 0;

    return;
}

void System:: updateParameter()
{
    // load-parameters
    Switches sw(vary_PARAMETERS_File);
    // update-yolk-pressure
    p.pressure += sw.pressure;
    // update-myosin
    int myoEdgeIndex(0);
    myosin_Edges.clear();
    double alpha_m_plus = 1.0 + 0.5*sw.delAlpha;
    for(int l = 1; l <= sw.myoCellNumber; l++)
    {
        myoEdgeIndex = (sw.refVertex - l + numNode)%numNode;
        // increase-apical-myosin
        geometry->CELL[myoEdgeIndex].alphA = alpha_m_plus;
        // decrease-basal-myosin
        geometry->CELL[myoEdgeIndex].alphB = 2.0 - alpha_m_plus;
        // coloring-myosin-induced-cells
        if(alpha_m_plus > 1.0)
        myosin_Edges.push_back((myoEdgeIndex+1+numNode)%numNode);
    }
    myosin_Edges.push_back(myoEdgeIndex);
    // update-fixed-vertex
    if(sw.apicalAttachment > 0)
    p.fixedVertex = sw.refVertex;
    // update-movie-switch
    p.movieFrameInterval = sw.movieFrameInterval;
    p.movieFrameNumber = sw.movieFrameNumber;

    return;
}

void System:: initializeGeometry()
{
    // reset-vertices/edges/cells
    Vv.clear();
    Va.clear();
    Vb.clear();
    lateral_Edg.clear();
    geometry->CELL.clear();
    // read-input-parameters
    p.initialize(constant_PARAMETERS_File);
    // type-of-epithelium
    epitheliumType = p.typeEpithelium;
    // number-of-vertices
    numNode = p.vertex_Number;
    // initialize-vertices/lateral-edges
    for(int i = 0; i < numNode; i++)
    {
        Vv.push_back(Nodes(0.0,0.0,1.0,1.0));
        Va.push_back(Nodes(0.0,0.0,0.0,1.0));
        Vb.push_back(Nodes(0.0,0.0,-1.0,1.0));
        lateral_Edg.push_back(lateralEdge(0.0,0.0));
    }
    // initialize-cells
    for(int i = 0; i < numNode; i++)
    {
        geometry->CELL.push_back(crossSectionRegion(geometry));
    }
    // update-default-lateral-tension
    default_alphaL  = p.lateralTension;
    ///////////////////////
    /// flat-epithelium ///
    ///////////////////////
    if(p.typeEpithelium== "flat")
    {
        // epithelial-width(initial)
        epithelialWidth = sqrt(p.cellRefArea/p.cellAspectRatio);
        // epithelial-height(initial)
        epithelialHeight = p.cellAspectRatio*epithelialWidth;
        // system-size
        L = (numNode*p.cellRefArea)/epithelialHeight;
        // midline-radius
        midline_radius  = L/(2.0*M_PI);
        initialTension = p.default_alphaA + p.default_alphaB - default_alphaL*numNode*(epithelialHeight/L);
        flatEpithelium_refX = -L/2.0;//0.0;
        flatEpithelium_refY = epithelialHeight;
        // vertices(CLOCK-WISE)-and-lateral-edges
        double xvalues(flatEpithelium_refX);
        for(int i = 0; i < numNode; i++)
        {
            // vertices
            Vv[i] = Nodes(xvalues,epithelialHeight,1.0,1.0);
            Va[i] = Nodes(xvalues,epithelialHeight,0.0,1.0);
            Vb[i] = Nodes(xvalues,0.0,-1.0,1.0);
            xvalues += epithelialWidth;
            // lateral-edge
            lateral_Edg[i] = lateralEdge(0.0,default_alphaL);
            // assign-default-value-of-apical/basal-tension-to-cell-edges
            geometry->CELL[i].alphA = p.default_alphaA;;
            geometry->CELL[i].alphB = p.default_alphaB;
        }
        // yolk-height
        H = - epithelialHeight;
        // cut-off-length-for-cell-edges
        lambda = (L*p.cutOffEdgeLength)/numNode;
        // equilibrate-epithelium
        updateGeometry();
        run("noUseFile",false);
        // update-vitelline-coordinates
        for(int i = 0; i < numNode; i++)
        {
            Vv[i].C(0,0) = Va[i].C(0,0);
            Vv[i].C(1,0) = Va[i].C(1,0);
        }
        // update-system-size
        L = 0.0;
        for(int i = 0; i < numNode-1; i++)
        {
            int j1 = (i+numNode)%numNode; // i
            int j2 = (i+1+numNode)%numNode; // i+1

            double dx =  0.5*((Va[j1].C(0,0) + Vb[j1].C(0,0))-(Va[j2].C(0,0) + Vb[j2].C(0,0)));
            double dy =  0.5*((Va[j1].C(1,0) + Vb[j1].C(1,0))-(Va[j2].C(1,0) + Vb[j2].C(1,0)));

            L+= sqrt(dx*dx + dy*dy);
        }
        L += epithelialWidth;
        // update-cut-off-length-for-cell-edges
        lambda = (L*p.cutOffEdgeLength)/numNode;
        epithelialHeight = Va[numNode/2].C(1,0) - Vb[numNode/2].C(1,0);
        initialTension = p.default_alphaA + p.default_alphaB - default_alphaL*numNode*(epithelialHeight/L);
        updateGeometry();
    }

    ///////////////////////////
    /// curved-epithelium ///
    ///////////////////////////
    else
    {
        // epithelial-width(initial)
        epithelialWidth = sqrt(p.cellRefArea/p.cellAspectRatio);
        // epithelial-height(initial)
        epithelialHeight = p.cellAspectRatio*epithelialWidth;
        // system-size
        L = (numNode*p.cellRefArea)/epithelialHeight;
        //  midline-radius: circle-pass-thgrough-the-midline-of-epithelium
        midline_radius  = L/(2.0*M_PI);
        // update-yolk-pressure
        initialTension = p.default_alphaA + p.default_alphaB - default_alphaL*numNode*(epithelialHeight/L);
        // vitelline-radius: taking-mid-line-of-apical-to-basal-surface-as-circle-reference
        vitelline_radius = midline_radius + (epithelialHeight/2.0);
        // consider-vitelline-radius: defalut-radius-go-thgrough-the-midline-of-epithelium
        double K_emb = M_PI*(3.0*(1.0 + p.embAspectRatio) - sqrt((3.0 + p.embAspectRatio)*(1.0 + 3.0*p.embAspectRatio)));
        double midline_semiMinorAxis = (numNode*p.cellRefArea)/(epithelialHeight*K_emb);
        double midline_semiMajorAxis = midline_semiMinorAxis*p.embAspectRatio;

        double axisCorrection = (epithelialHeight/2.0);
        semiMajorAxis = midline_semiMajorAxis + axisCorrection;
        semiMinorAxis = midline_semiMinorAxis + axisCorrection;

        // vertices(CLOCK-WISE)-and-lateral-edges
        vector<Vector2d> V;
        vector<Vector2d> normalDirector;
        Vector2d vertexUpdate;
        double angleStep = 2.0*M_PI/numNode;
        for(int i = 0; i < numNode; i++)
        {
            double angle = 2.0*M_PI  - i*angleStep;
            V.push_back(Vector2d(semiMajorAxis*cos(angle),semiMinorAxis*sin(angle)));
        }
        // normals-on-vitelline: apical/basal-vertices
        for(int i = 0; i < numNode; i++)
        {
            // avoid: angle = PI/2, 3PI/
            V[i] = (V[i] + V[(i+1+numNode)%numNode])/2.0;
            Vector2d normalToCurvedLine(semiMinorAxis*V[i](0,0)/semiMajorAxis, semiMajorAxis*V[i](1,0)/semiMinorAxis);
            normalDirector.push_back(Vector2d(normalToCurvedLine));
            normalDirector[i]/=normalDirector[i].norm();
            vertexUpdate = V[i];
            Vv[i] = Nodes(vertexUpdate(0,0),vertexUpdate(1,0),1.0,1.0);
            vertexUpdate = V[i];
            Va[i] = Nodes(vertexUpdate(0,0),vertexUpdate(1,0),0.0,1.0);
            vertexUpdate = V[i] - epithelialHeight*normalDirector[i];
            Vb[i] = Nodes(vertexUpdate(0,0),vertexUpdate(1,0),-1.0,1.0);
            // lateral-edge
            lateral_Edg[i] = lateralEdge(0.0,default_alphaL);
            // assign-default-value-of-apical/basal-tension-to-cell-edges
            geometry->CELL[i].alphA = p.default_alphaA;
            geometry->CELL[i].alphB = p.default_alphaB;
        }
        // cut-off-length-for-cell-edges
        lambda = (L*p.cutOffEdgeLength)/numNode;
        // update-system-size
        double yolkArea_ref = 0.0;
        for(int i = 0; i < numNode; i++)
        {
            int j1 = (i+numNode)%numNode; // i
            int j2 = (i-1+numNode)%numNode; // i-1
            int j3 = (i+1+numNode)%numNode; // i+1

            yolkArea_ref+= Vb[j1].C(0,0)*(Vb[j2].C(1,0)-Vb[j3].C(1,0));
        }
        // update-yolk-reference-area
        geometry->A_Y_ref = 0.5*yolkArea_ref;
        if(geometry->A_Y_ref < 0)
        {
            cout << "initial yolk-area is NEGATIVE !!! -> EXIT " << endl;
            exit(0);
        }
        // equilibrate-epithelium
        double laplas_L = M_PI*(3*(semiMajorAxis+semiMinorAxis) - sqrt((3*semiMajorAxis+semiMinorAxis)*(semiMajorAxis+3*semiMinorAxis)));
        double laplas_Tension = p.default_alphaA + p.default_alphaB - default_alphaL*numNode*(epithelialHeight/laplas_L);
        laplasPressure =  laplas_Tension/semiMajorAxis;
        p.pressure = laplasPressure;
        updateGeometry();
        run("noUseFile",false);
        // update-vitelline-coordinates
        for(int i = 0; i < numNode; i++)
        {
            Vv[i].C(0,0) = Va[i].C(0,0);
            Vv[i].C(1,0) = Va[i].C(1,0);
        }
        // update-system-size
        double length_Epithelium = 0.0;
        for(int i = 0; i < numNode; i++)
        {
            int j1 = (i+numNode)%numNode; // i
            int j2 = (i-1+numNode)%numNode; // i-1
            int j3 = (i+1+numNode)%numNode; // i+1
            double dx =  0.5*((Va[j1].C(0,0) + Vb[j1].C(0,0))-(Va[j3].C(0,0) + Vb[j3].C(0,0)));
            double dy =  0.5*((Va[j1].C(1,0) + Vb[j1].C(1,0))-(Va[j3].C(1,0) + Vb[j3].C(1,0)));
            length_Epithelium+= sqrt(dx*dx + dy*dy);
        }
        // update-cut-off-length-for-cell-edges
        midline_radius  = length_Epithelium/(2.0*M_PI);
        lambda = (length_Epithelium*p.cutOffEdgeLength)/numNode;
        double L_Theo = M_PI*(3.0*(semiMajorAxis+semiMinorAxis) -sqrt((3*semiMajorAxis+semiMinorAxis)*(3*semiMinorAxis+semiMajorAxis)));
        updateGeometry();
    }

    return;
}

/// scale-geometry
double System:: scaleGeometry(double dL, double lateral_Tension)
{

    double L_new = L + dL;

    double epithelialWidth_new = L_new/numNode;

    double epithelialHeight_new = p.cellAspectRatio*epithelialWidth_new;

    // update-vertices(CLOCK-WISE)
    double xvalues(-L_new/2.0);
    for(int i = 0; i < numNode; i++)
    {
        Vv[i].C << xvalues,epithelialHeight_new;
        Va[i].C << xvalues,epithelialHeight_new;
        Vb[i].C << xvalues,0.0;
        xvalues+=epithelialWidth_new;
    }

    // update-system-size
    L = 0.0;
    for(int i = 0; i < numNode-1; i++)
    {
        int j1 = (i+numNode)%numNode; // i
        int j2 = (i+1+numNode)%numNode; // i+1

        double dx =  0.5*((Va[j1].C(0,0) + Vb[j1].C(0,0))-(Va[j2].C(0,0) + Vb[j2].C(0,0)));
        double dy =  0.5*((Va[j1].C(1,0) + Vb[j1].C(1,0))-(Va[j2].C(1,0) + Vb[j2].C(1,0)));

        L+= sqrt(dx*dx + dy*dy);
    }
    L += epithelialWidth_new;

    // update-the-system
    updateGeometry();

    return ((numNode)/L_new)*(p.cellRefArea - (numNode*lateral_Tension/(2.0*p.betaCell*L_new)));
}

