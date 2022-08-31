#include "drosoSymBreak.h"

int main(int argc, char* argv[])
{
	//***************************//
    /// define-reference-system ///
    //***************************//
    System* embryo = new System(argv[1],argv[2],"./simulationData");
    embryo->createOutPutDirectory();
    embryo->saveConfiguration("refFrame");

    //***********************************************//
    /// energy-minimization-with-default-parameters ///
    //***********************************************//
    embryo->initializeParameter();
    embryo->run("noUseFile",false);
    embryo->saveConfiguration("initialFrame");

    //***********************************************//
    /// energy-minimization-with-updated-parameters ///
    //***********************************************//
    // update-parameters
    embryo->updateParameter();
    embryo->run("timeSeriesFrames",false);
    embryo->saveConfiguration("finalFrame");

    //****************//
    /// reset-system ///
    //****************//
    delete embryo;
    embryo = NULL;

	return (0);
}
