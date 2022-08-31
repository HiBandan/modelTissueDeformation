#include "drosoSymBreak.h"

Switches::Switches(char* swFile):
apicalAttachment(0),delAlpha(0.0),
refVertex(0),

movieFrameInterval(0),movieFrameNumber(0),
pressure(0.0),myoCellNumber(0)
{
    // read-input-switches
    initialize(swFile);

    return;
}

void Switches::initialize(const char* sf)
{
    // if-no-input-file-exists/if-unable-to-open-the-input-file
    if (!readFromFile(sf, true))
    {
        cout << "Failed to load parameters from file [" << sf << "]. Aborting\n";
        exit(-1);
    }

    return;
}

/// load-individual-switches
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

/// read-items-from-input-switch-file
bool Switches::readFromFile(const char* sf, bool initialRun)
{
    // if-no-input-file-exists
    if ((sf == NULL) || (*sf == 0))
        return false;

    // try-to-open-the-input-file
    ifstream swFile;
    swFile.open(sf);

    // if-unable-to-open-the-input-file
    if (!swFile.good())
    {
        cerr << "Could not open parameter file [" << sf << "].\n";
        return false;
    }

    // if-input-file-exists-and-able-to-open
    string id;
    bool recognized;
    while (swFile.good() )
    {
        swFile >> id;

        // if-the-current-item-if-not-detectable
        if(!swFile.good())
        break;

        // if-a-commented-line-is-detected-skip-and-move-to-the-next-line
        if (id[0] == '#')
        getline(swFile,id);

        // if-an-item-is-found
        else
        {
            // start-scanning-the-items-assuming-no-item-is-detected-yet
            recognized = false;

            // scan-the-individual-items

            loadParam(recognized, swFile, id, "apicalAttachment",apicalAttachment);
            loadParam(recognized, swFile, id, "refVertex",refVertex);

	        loadParam(recognized, swFile, id, "movieFrameInterval",movieFrameInterval);
	        loadParam(recognized, swFile, id, "movieFrameNumber",movieFrameNumber);


            loadParam(recognized, swFile, id, "pressure",pressure);

            loadParam(recognized, swFile, id,"delAlpha",delAlpha);
            loadParam(recognized, swFile, id,"myoCellNumber",myoCellNumber);


	        // if-undeclared-item-is-detected-in-the-input-file
            if (!recognized)
            {
                cerr << "Parameter not recognized: " << id << "\n";
            }

            // skip-the-last-line-of-the-input-file
            getline(swFile, id);
        }
    }

    return true;
}

