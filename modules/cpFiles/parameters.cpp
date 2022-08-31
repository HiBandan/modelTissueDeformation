#include "drosoSymBreak.h"

Parameters::Parameters(System* s):
system(s),fixedVertex(-1),default_alphaA(0.0),default_alphaB(0.0),cutOffEdgeLength(0.0),pressure(0.0),epsilon(0.0),enegryTolerance(1e-6),betaCell(0.0),cellRefArea(0.0),
cellAspectRatio(0.0),embAspectRatio(1.0),betaYolk(0.0),movieFrameInterval(0),movieFrameNumber(0),typeEpithelium("."),vertex_Number(0),lateralTension(0.0)
{
    return;
}

void Parameters::initialize(const char* pf)
{
    // if-no-input-file-exists/if-unable-to-open-the-input-file
    if (!readFromFile(pf, true))
    {
        cout << "Failed to load parameters from file [" << pf << "]. Aborting\n";
        exit(-1);
    }

    return;
}

template<class T> inline void loadParam(bool& recognized, ifstream& file, string& id, const char* tag, T& target)
{
    // try-to-read-the-target-item-only-if-already-not-detected
    if (!recognized)
    {
        // check-if-name-of-the-target-item-is-same-as-defined-in-the-input-file
        if (id == tag)
        {
            // now-read-the-target-item
            file >> target;
            recognized = true;
        }
    }

    return;
}

bool Parameters::readFromFile(const char* pf, bool initialRun)
{
    // if-no-input-file-exists
    if ((pf == NULL) || (*pf == 0))
        return false;

    // try-to-open-the-input-file
    ifstream parFile;
    parFile.open(pf);

    // if-unable-to-open-the-input-file
    if (!parFile.good())
    {
        cerr << "Could not open parameter file [" << pf << "].\n";
        return false;
    }

    // if-input-file-exists-and-able-to-open
    string id;
    bool recognized;
    while (parFile.good() )
    {
        parFile >> id;

        // if-the-current-item-if-not-detectable
        if (!parFile.good())
        break;

        // if-a-commented-line-is-detected-skip-and-move-to-the-next-line
        if (id[0] == '#')
        getline(parFile,id);

        // if-an-item-is-found
        else
        {
            // start-scanning-the-items-assuming-no-item-is-detected-yet
            recognized = false;

            loadParam(recognized, parFile, id, "default_alphaA",default_alphaA);
            loadParam(recognized, parFile, id, "default_alphaB",default_alphaB);


            loadParam(recognized, parFile, id, "cutOffEdgeLength",cutOffEdgeLength);
            loadParam(recognized, parFile, id, "betaCell",betaCell);
            loadParam(recognized, parFile, id, "embAspectRatio",embAspectRatio);

            loadParam(recognized, parFile, id, "cellRefArea",cellRefArea);
            loadParam(recognized, parFile, id, "cellAspectRatio",cellAspectRatio);

            loadParam(recognized, parFile, id, "epsilon",epsilon);
            loadParam(recognized, parFile, id, "enegryTolerance",enegryTolerance);
            loadParam(recognized, parFile, id, "betaYolk",betaYolk);

            loadParam(recognized, parFile, id, "typeEpithelium",typeEpithelium);

            loadParam(recognized, parFile, id, "vertex_Number",vertex_Number);

            loadParam(recognized, parFile, id, "lateralTension",lateralTension);


            // if-undeclared-item-is-detected-in-the-input-file
            if (!recognized)
            {
                cerr << "Parameter not recognized: " << id << "\n";
            }

            // skip-the-last-line-of-the-input-file
            getline(parFile, id);
        }
    }

    return true;
}

void Parameters::writeToFile(string outPutFileName)
{
    ofstream of;
    of.open(outPutFileName);

    of << "embAspectRatio" << " " << embAspectRatio << endl;


    of << "betaCell" << " " << betaCell << endl;
    of << "cellRefArea" << " " << cellRefArea << endl;
    of << "cellAspectRatio" << " " << cellAspectRatio << endl;
    of << "default_alphaA" << " " << default_alphaA << endl;
    of << "default_alphaB" << " " << default_alphaB << endl;


    of << "cutOffEdgeLength" << " " << cutOffEdgeLength << endl;
    of << "epsilon" << " " << epsilon << endl;

    of << "pressure" << " " << pressure << endl;
    of << "betaYolk" << " " << betaYolk << endl;
    of << "enegryTolerance" << " " << enegryTolerance << endl;
    of << "fixedVertex" << " " << fixedVertex << endl;

    of.close();

    return;
}



