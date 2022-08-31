#include "drosoSymBreak.h"

/// save-movie-frame
void System::initializeOutput(string fileName)
{
    string FileName = outputDir + "/" + fileName;
    movieFile.open(FileName.c_str());

    if (movieFile.is_open())
    {
        movieFile << epitheliumType << "," << semiMajorAxis << "," << semiMinorAxis << "," << midline_radius << "," << vitellineWidth
        << "," << flatEpithelium_refY << "," << flatEpithelium_refX << "," << epithelialHeight << "," << L << endl;
    }
    else
    {
        cout << "ERROR: cannot open file for storing movie frames -> " << "[" << fileName << "]" << endl;
        exit(-1);
    }

    return;
}
void System::saveFrame(void)
{
    ////////////////////
    /// fixed-vertex ///
    ////////////////////
    if (p.fixedVertex > 0)
    {
        for(int i = p.fixedVertex; i < p.fixedVertex + 1; i++)
        movieFile << i <<"@";
    }
    else
    {

    }
    movieFile << "#";

    ////////////////////
    /// all-vertices ///
    ////////////////////
    for(int i = 0; i < numNode; i++)
    {
        movieFile << Vv[i].C(0,0) << "," << Vv[i].C(1,0) << ",";
        movieFile << Va[i].C(0,0) << "," << Va[i].C(1,0) << ",";
        movieFile << Vb[i].C(0,0) << "," << Vb[i].C(1,0) <<"@";
    }

    /////////////////
    /// all-edges ///
    /////////////////
    int numEdges = (epitheliumType == "flat") ? numNode-1 : numNode;
	movieFile << "#";
	for(int i = 0; i < numEdges; i++)
    movieFile << (i+numNode)%numNode << "," << (i+1+numNode)%numNode <<"@";

    ////////////////////
    /// myosin-edges ///
    ////////////////////
    movieFile << "#";
    if(myosin_Edges.size() > 0)
    {
        // myosin-apical-edges
        for(int i = 0; i < myosin_Edges.size()-1; i++)
        movieFile << myosin_Edges[i] << "," << myosin_Edges[i+1] <<"@";
    }

    ///////////////////////////////////////
    /// tissue-tension/detachment-force ///
    ///////////////////////////////////////
    movieFile << "#";
    movieFile << epithelialHeight << "," << F_x << "," << F_y << "," << v_Fx[0]<< "," << v_Fy[0]<< "," << v_Fx[1]<< "," << v_Fy[1];
    movieFile << endl;
}
void System::closeFiles()
{
    movieFile.close();
    return;
}

/// save-default-embryo-information
void System:: saveEmbryoConfiguration(string paramFile)
{
    ofstream pF;
    pF.open(paramFile.c_str());

    // default-system-information
    pF << "Number of vertices:" << numNode << endl;
	pF << "Epithelial width:" << L/numNode << endl;
	pF << "Epithelial height:" << epithelialHeight << endl;
	pF << "Vitelline default width:" << vitellineWidth << endl;
    pF << "Length of the epithelium:" << L << endl;
    pF << "YolkArea " << geometry->A_Y_ref << endl;

    // vitelline-vertex-coordinates
    pF << "vitelline-vertices" << endl;
    for(int i = 0; i < numNode; i++)
    pF << i << " " << Vv[i].C(0,0) << " " << Vv[i].C(1,0) << " "<< Vv[i].EG(0,0) << " " << Vv[i].EG(1,0) << endl;

    // apical-vertex-coordinates
    pF << "apical-vertices" << endl;
    for(int i = 0; i < numNode; i++)
    pF << i << " " << Va[i].C(0,0) << " " << Va[i].C(1,0) << " "<< Va[i].EG(0,0) << " " << Va[i].EG(1,0) << endl;

    // basal-vertex-coordinates
    pF << "basal-vertices" << endl;
    for(int i = 0; i < numNode; i++)
    pF << i << " " << Vb[i].C(0,0) << " " << Vb[i].C(1,0) << " " << Vb[i].EG(0,0) << " " << Vb[i].EG(1,0) << endl;

    // number-of-edges
    int numEdges = (epitheliumType == "flat") ? numNode-1 : numNode;
    pF << "Number of edges:" << numEdges << endl;
    pF.close();
    return ;
}

