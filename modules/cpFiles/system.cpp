#include "drosoSymBreak.h"
#include "MersenneTwister.h"
#include "LmDerivativeNewton.h"
#include "CgMinimizer-inline.h"
#include "CgMinimizerSearchDirectionUpdates.h"

// default-values-for-minimizer: Lmc2(1e-4),LmMinMaxAbsSlope(1e-10),CgMaxItPerDof(1000),CgGradientTolPerDof(1e-6),CgnitialStepSizePerDof(1e-6),

/// constructor/destructor
System::System(char* paramFile_constant, char* paramFile_vary, string outputPath):
p(this),constant_PARAMETERS_File(paramFile_constant),vary_PARAMETERS_File(paramFile_vary),outputDir(outputPath),Lmc2(1e-4),LmMinMaxAbsSlope(1e-10),
CgMaxItPerDof(1000),CgGradientTolPerDof(1e-6),CgnitialStepSizePerDof(1e-6),numNode(0),L(0.0),lambda(0.0),H(0.0),F_x(0.0),semiMajorAxis(0.0),semiMinorAxis(0.0),F_y(0.0),epithelialHeight(0.0),
epithelialWidth(0.0),flatEpithelium_refX(0.0),vitellineWidth(1e-6),flatEpithelium_refY(0.0),midline_radius(0.0),initialTension(0.0),outputFileName("."),geometry(NULL)
{
    // initialize-geometry
    geometry = new Geometry(this);
    initializeGeometry();

    return;
}

System::~System()
{
    return;
}

void System:: createOutPutDirectory(void)
{
    if (boost::filesystem::exists(outputDir))
    {
        boost::filesystem::remove_all(outputDir);
    }
    boost::filesystem::create_directories(outputDir);
}

/// run
void System::run(string allFrameFile, bool checkIntegrity)
{
    initializeOutput(allFrameFile);
    // initialize-the-line-minimizer: newton-method
	LmDerivativeNewton lmDn;
  	lmDn.setC2(Lmc2);
  	lmDn.setMinMaxAbsSlope(LmMinMaxAbsSlope);
    // initialize-the-CG-minimizer: seed-through-the-line-minimizer
	CgMinimizer<LmDerivativeNewton,PolakRibiere> minimizer(lmDn);
	minimizer.setGradientTolerancePerDof(CgGradientTolPerDof);
  	minimizer.setMaxIterationsPerDof(CgMaxItPerDof);
  	minimizer.setInitialStepSizePerDof(CgnitialStepSizePerDof);
  	minimizer.setLogInterval(p.movieFrameInterval);
  	minimizer.setLogNumber(p.movieFrameNumber);
  	minimizer.clearDofs();
    // passing-only-free-apical-vertices-to-the-minimizer
    for(int i = 0; i < p.fixedVertex; i++)
    {
        minimizer.addDof(&Va[i].C(0,0),&Va[i].EG(0,0)); // x
        minimizer.addDof(&Va[i].C(1,0),&Va[i].EG(1,0)); // y
    }

    for(int i = p.fixedVertex + 1; i < numNode; i++)
    {
        minimizer.addDof(&Va[i].C(0,0),&Va[i].EG(0,0)); // x
        minimizer.addDof(&Va[i].C(1,0),&Va[i].EG(1,0)); // y
    }
    // passing-only-free-basal-vertices-to-the-minimizer
    for(int i = 0; i < numNode; i++)
    {
        minimizer.addDof(&Vb[i].C(0,0),&Vb[i].EG(0,0));// x
        minimizer.addDof(&Vb[i].C(1,0),&Vb[i].EG(1,0));// y
    }

    //////////////////////////
    /// start-minimization ///
    //////////////////////////
    bool success = minimizer.minimize
	(
	    // access_iteration-number-on-FLY
        checkIntegrity,
        // pass-as-lambda-function: updateGeometry-and-updateGradient
		[this]()->void
		{
		    this->updateGeometry();
			this->updateGradient();
		},
		// pass-as-lambda-function: saveFrame
		[this]()->void
		{
		    this->saveFrame();
		},
		// pass-as-lambda-function: checkIntegrity
		[this]()->bool
		{
            this->checkEnergyGradient();
		}
	);
	// verify-if-integrity-check-on-energy-gradient-is-sucessful
	if(!success)
	{
        cout << "check for integrity of update gradient failed" << endl;
        //exit(-1);
	}
    // close-all-the-output-files
    closeFiles();

    return;
}

void System:: saveConfiguration(string fileName)
{
    initializeOutput(fileName);
    saveFrame();
    closeFiles();
}

/// update-geometry
void System:: updateGeometry(void)
{
    double dx(0.0), dy(0.0);
    for(int i = 0; i < numNode; i++)
    {
        int j1 = (i+numNode)%numNode; // i
        int j2 = (i-1+numNode)%numNode; // i-1
        int j3 = (i+1+numNode)%numNode; // i+1
        // area-of-cell
        geometry->CELL[i].A = 0.5*(Va[j3].C(0,0)*Va[j1].C(1,0) + Vb[j3].C(0,0)*Va[j3].C(1,0)
                                     + Vb[j1].C(0,0)*Vb[j3].C(1,0) + Va[j1].C(0,0)*Vb[j1].C(1,0)
                                     -Va[j1].C(0,0)*Va[j3].C(1,0) - Va[j3].C(0,0)*Vb[j3].C(1,0)
                                     - Vb[j3].C(0,0)*Vb[j1].C(1,0) - Vb[j1].C(0,0)*Va[j1].C(1,0));
        // length-of-apical-edge
        dx = Va[j3].C(0,0)-Va[j1].C(0,0);
        dy = Va[j3].C(1,0)-Va[j1].C(1,0);
        geometry->CELL[i].l_a = sqrt(dx*dx + dy*dy);
        // length-of-basal-edge
        dx = Vb[j3].C(0,0)-Vb[j1].C(0,0);
        dy = Vb[j3].C(1,0)-Vb[j1].C(1,0);
        geometry->CELL[i].l_b = sqrt(dx*dx + dy*dy);
        // length-of-lateral-edge
        dx = Va[j1].C(0,0)-Vb[j1].C(0,0);
        dy = Va[j1].C(1,0)-Vb[j1].C(1,0);
        lateral_Edg[i].l_l = sqrt(dx*dx + dy*dy);
        // flat-epithelium
        if(epitheliumType == "flat")
        {
            Vv[j1].C(0,0) = Va[j1].C(0,0);
            Va[i].D = DistancePoints(Va[j1].C(0,0),Vv[j1].C(0,0),Va[j1].C(1,0),Vv[j1].C(1,0))*isPointBelowHorizontalLine(Va[j1].C(1,0),Vv[j1].C(1,0));
        }
        // curved-epithelium: circular/elliptic
        else
        {
            Array2d intersectionPoint;
            Array2d targetPoint(Va[j1].C(0,0),Va[j1].C(1,0));
            Va[i].D = DistancePointEllipse(semiMajorAxis,semiMinorAxis,targetPoint,intersectionPoint,1000)*isPointInsideEllipse(Va[j1].C(0,0),Va[j1].C(1,0),semiMajorAxis,semiMinorAxis);
            Vv[j1].C(0,0) = intersectionPoint(0,0);
            Vv[j1].C(1,0) = intersectionPoint(1,0);
        }
    }
    // area-of-yolk
    double yolkArea(0.0);
    for(int i = 0; i < numNode; i++)
    {
        int j1 = (i+numNode)%numNode; // i
        int j2 = (i-1+numNode)%numNode; // i-1
        int j3 = (i+1+numNode)%numNode; // i+1

        yolkArea+= Vb[j1].C(0,0)*(Vb[j2].C(1,0)-Vb[j3].C(1,0));
    }
    yolkArea *=0.5;
    // area-of-vitelline-space
    double viellineSpaceArea(0.0);
    for(int i = 0; i < numNode; i++)
    {
        int j1 = (i+numNode)%numNode; // i
        int j2 = (i-1+numNode)%numNode; // i-1
        int j3 = (i+1+numNode)%numNode; // i+1

        viellineSpaceArea+= Va[j1].C(0,0)*(Va[j2].C(1,0)-Va[j3].C(1,0));
    }
    viellineSpaceArea *=-0.5;

    if(epitheliumType == "flat")
    {
        // corrected-area-of-cell-due-to-apical-vertices
        geometry->CELL[numNode-1].A += 0.5*L*(Va[numNode-1].C(1,0) + Va[0].C(1,0));
        // corrected-area-of-cell-due-to-basal-vertices
        geometry->CELL[numNode-1].A -= 0.5*L*(Vb[numNode-1].C(1,0) + Vb[0].C(1,0));
        // corrected-length-of-apical-edge
        dx = Va[0].C(0,0) + L - Va[numNode-1].C(0,0);
        dy = Va[0].C(1,0) - Va[numNode-1].C(1,0);
        geometry->CELL[numNode-1].l_a = sqrt(dx*dx + dy*dy);
        // corrected-length-of-the-basal-edge
        dx = Vb[0].C(0,0) + L - Vb[numNode-1].C(0,0);
        dy = Vb[0].C(1,0) - Vb[numNode-1].C(1,0);
        geometry->CELL[numNode-1].l_b = sqrt(dx*dx + dy*dy);
        // correction-to-yolk-area
        yolkArea+= 0.5*L*(Vb[numNode-1].C(1,0)+Vb[0].C(1,0));
        yolkArea-= L*H;
        // correction-to-vitelline-area
        viellineSpaceArea-= 0.5*L*(Va[numNode-1].C(1,0)+Va[0].C(1,0));
        viellineSpaceArea+= L*flatEpithelium_refY;
    }
    // corrected-area-of-the-vitelline-space
    geometry->A_V = viellineSpaceArea;
    // corrected-area-of-the-yolk
    geometry->A_Y = yolkArea;
    // weight-factors
    double l_r_inv = 1.0/lambda;
    for(int i = 0; i < numNode; i++)
    {
        int j1 = (i+numNode)%numNode; // i
        int j2 = (i-1+numNode)%numNode; // i-1
        // weight-factor-to-vitelline-distance
        Va[i].w = 1.0;
        // weight-factor-to-the-length-of-apical-edge
        geometry->CELL[i].m_a = 1.0/(exp(geometry->CELL[i].l_a*l_r_inv) - 1.0);
        geometry->CELL[i].g_a = 1.0 - geometry->CELL[i].m_a*(1.0 + geometry->CELL[i].m_a);
        // weight-factor-to-the-length-of-basal-edge
        geometry->CELL[i].m_b = 1.0/(exp(geometry->CELL[i].l_b*l_r_inv) - 1.0);
        geometry->CELL[i].g_b = 1.0 - geometry->CELL[i].m_b*(1.0 + geometry->CELL[i].m_b);
        // weight-factor-to-the-length-of-lateral-edge
        lateral_Edg[i].m_l = 1.0/(exp(lateral_Edg[i].l_l*l_r_inv) - 1.0);
        lateral_Edg[i].g_l = 1.0 - lateral_Edg[i].m_l*(1.0 + lateral_Edg[i].m_l);
    }
}

/// update-gradient
void System::updateGradient(void)
{
    double edgeContribution_apical(0.0),cellContribution_apical(0.0),vitellineContribution_apical(0.0);
    double edgeContribution_basal(0.0),cellContribution_basal(0.0),yolkContribution_basal(0.0);
    double vitelline_EG_x[numNode] = {0.0},vitelline_EG_y[numNode] = {0.0};
    for(int i = 0; i < numNode; i++)
    {
        int j1 = (i+numNode)%numNode; // i
            int j2 = (i-1+numNode)%numNode; // i-1
            int j3 = (i+1+numNode)%numNode; // i+1
            // residual-cell-area
            double A1 = geometry->CELL[j1].A - p.cellRefArea;
            double A2 = geometry->CELL[j2].A - p.cellRefArea;
            double A3 = geometry->A_Y - geometry->A_Y_ref;
            // energy-gradient-at-apical-vertices : w.r.t. -> x
            Va[i].EG(0,0)  = geometry->CELL[j1].g_a*geometry->CELL[j1].alphA*(Va[j1].C(0,0)-Va[j3].C(0,0))/geometry->CELL[j1].l_a
                                     + geometry->CELL[j2].g_a*geometry->CELL[j2].alphA*(Va[j1].C(0,0)-Va[j2].C(0,0))/geometry->CELL[j2].l_a
                                     + lateral_Edg[j1].g_l*lateral_Edg[j1].alphL*(Va[j1].C(0,0)-Vb[j1].C(0,0))/lateral_Edg[j1].l_l
                                     + p.betaCell*(A1*(Vb[j1].C(1,0)-Va[j3].C(1,0))+A2*(Va[j2].C(1,0)-Vb[j1].C(1,0)))
                                     + 2.0*p.epsilon*Va[j1].w*HF(Va[j1].D)*(Va[j1].C(0,0)-Vv[j1].C(0,0));
            vitelline_EG_x[i] = 2.0*p.epsilon*Va[j1].w*HF(Va[j1].D)*(Va[j1].C(0,0)-Vv[j1].C(0,0));
            // energy-gradient-at-apical-vertices : w.r.t. -> y
            Va[i].EG(1,0)  = geometry->CELL[j1].g_a*geometry->CELL[j1].alphA*(Va[j1].C(1,0)-Va[j3].C(1,0))/geometry->CELL[j1].l_a
                                    + geometry->CELL[j2].g_a*geometry->CELL[j2].alphA*(Va[j1].C(1,0)-Va[j2].C(1,0))/geometry->CELL[j2].l_a
                                    + lateral_Edg[j1].g_l*lateral_Edg[j1].alphL*(Va[j1].C(1,0)-Vb[j1].C(1,0))/lateral_Edg[j1].l_l
                                    + p.betaCell*(A1*(Va[j3].C(0,0)-Vb[j1].C(0,0))+A2*(Vb[j1].C(0,0)-Va[j2].C(0,0)))
                                    + 2.0*p.epsilon*Va[j1].w*HF(Va[j1].D)*(Va[j1].C(1,0)-Vv[j1].C(1,0));
            vitelline_EG_y[i] = 2.0*p.epsilon*Va[j1].w*HF(Va[j1].D)*(Va[j1].C(1,0)-Vv[j1].C(1,0));
            // energy-gradient-at-basal-vertices: w.r.t. -> x
            Vb[i].EG(0,0) = geometry->CELL[j1].g_b*geometry->CELL[j1].alphB*(Vb[j1].C(0,0)-Vb[j3].C(0,0))/geometry->CELL[j1].l_b
                                   + geometry->CELL[j2].g_b*geometry->CELL[j2].alphB*(Vb[j1].C(0,0)-Vb[j2].C(0,0))/geometry->CELL[j2].l_b
                                   + lateral_Edg[j1].g_l*lateral_Edg[j1].alphL*(Vb[j1].C(0,0)-Va[j1].C(0,0))/lateral_Edg[j1].l_l
                                   + p.betaCell*(A1*(Vb[j3].C(1,0)-Va[j1].C(1,0))+A2*(Va[j1].C(1,0)-Vb[j2].C(1,0)))
                                   -0.5*p.pressure*(Vb[j2].C(1,0)-Vb[j3].C(1,0)) + p.betaYolk*A3*(Vb[j2].C(1,0)-Vb[j3].C(1,0));
            // energy-gradient-at-basal-vertices: w.r.t. -> y
            Vb[i].EG(1,0) = geometry->CELL[j1].g_b*geometry->CELL[j1].alphB*(Vb[j1].C(1,0)-Vb[j3].C(1,0))/geometry->CELL[j1].l_b
                                   + geometry->CELL[j2].g_b*geometry->CELL[j2].alphB*(Vb[j1].C(1,0)-Vb[j2].C(1,0))/geometry->CELL[j2].l_b
                                   + lateral_Edg[j1].g_l*lateral_Edg[j1].alphL*(Vb[j1].C(1,0)-Va[j1].C(1,0))/lateral_Edg[j1].l_l
                                   + p.betaCell*(A1*(Va[j1].C(0,0)-Vb[j3].C(0,0))+A2*(Vb[j2].C(0,0)-Va[j1].C(0,0)))
                                   -0.5*p.pressure*(Vb[j3].C(0,0)-Vb[j2].C(0,0)) + p.betaYolk*A3*(Vb[j3].C(0,0)-Vb[j2].C(0,0));
    }

    ////////////////////////////////
    /// periodic-boundary-effect ///
    ////////////////////////////////
    if(epitheliumType == "flat")
    {
        double contributionOnBoundary(0.0);
        ///////////////////
        //     i = 0     //
        ///////////////////
        // energy-gradient-at-apical-vertices : w.r.t. -> x
        edgeContribution_apical = geometry->CELL[numNode-1].g_a*L*geometry->CELL[numNode-1].alphA/geometry->CELL[numNode-1].l_a;
        cellContribution_apical = 0.0;
        contributionOnBoundary = edgeContribution_apical + cellContribution_apical;
        Va[0].EG(0,0) += contributionOnBoundary;
        // energy-gradient-at-apical-vertices : w.r.t. -> y
        edgeContribution_apical = 0.0;
        cellContribution_apical = L*p.betaCell*(geometry->CELL[numNode-1].A-p.cellRefArea);
        contributionOnBoundary = edgeContribution_apical + cellContribution_apical;
        Va[0].EG(1,0) += contributionOnBoundary;
        // energy-gradient-at-basal-vertices: w.r.t. -> x
        edgeContribution_basal = geometry->CELL[numNode-1].g_b*L*geometry->CELL[numNode-1].alphB/geometry->CELL[numNode-1].l_b;
        cellContribution_basal = 0.0;
        yolkContribution_basal = 0.0;
        contributionOnBoundary = edgeContribution_basal + cellContribution_basal + yolkContribution_basal;
        Vb[0].EG(0,0) += contributionOnBoundary;
        // energy-gradient-at-basal-vertices: w.r.t. -> y
        edgeContribution_basal = 0.0;
        cellContribution_basal = -L*p.betaCell*(geometry->CELL[numNode-1].A-p.cellRefArea);
        yolkContribution_basal = -0.5*p.pressure*L;
        contributionOnBoundary = edgeContribution_basal + cellContribution_basal + yolkContribution_basal;
        Vb[0].EG(1,0) += contributionOnBoundary;

        /////////////////////
        //     i = n-1     //
        /////////////////////
        // energy-gradient-at-apical-vertices : w.r.t. -> x
        edgeContribution_apical = -geometry->CELL[numNode-1].g_a*L*geometry->CELL[numNode-1].alphA/geometry->CELL[numNode-1].l_a;
        cellContribution_apical = 0.0;
        contributionOnBoundary = edgeContribution_apical + cellContribution_apical;
        Va[numNode-1].EG(0,0) += contributionOnBoundary;
        // energy-gradient-at-apical-vertices : w.r.t. -> y
        edgeContribution_apical = 0.0;
        cellContribution_apical = L*p.betaCell*(geometry->CELL[numNode-1].A-p.cellRefArea);
        contributionOnBoundary = edgeContribution_apical + cellContribution_apical;
        Va[numNode-1].EG(1,0) += contributionOnBoundary;
        // energy-gradient-at-basal-vertices: w.r.t. -> x
        edgeContribution_basal = -geometry->CELL[numNode-1].g_b*L*geometry->CELL[numNode-1].alphB/geometry->CELL[numNode-1].l_b;
        cellContribution_basal = 0.0;
        yolkContribution_basal = 0.0;
        contributionOnBoundary = edgeContribution_basal + cellContribution_basal + yolkContribution_basal;
        Vb[numNode-1].EG(0,0) += contributionOnBoundary;

        // energy-gradient-at-basal-vertices: w.r.t. -> y
        edgeContribution_basal = 0.0;
        cellContribution_basal = -L*p.betaCell*(geometry->CELL[numNode-1].A-p.cellRefArea);
        yolkContribution_basal = -0.5*L*p.pressure;
        contributionOnBoundary = edgeContribution_basal + cellContribution_basal + yolkContribution_basal;
        Vb[numNode-1].EG(1,0) += contributionOnBoundary;
    }

    // force-on-fixed-vertex
    F_x = -Va[p.fixedVertex].EG(0,0);
    F_y = -Va[p.fixedVertex].EG(1,0);

    // force-on-2-vertices-anterior-to-the-fixed-vertex
    for(int i = 0; i < 2; i++)
    {
        v_Fx[i] = -vitelline_EG_x[i +(p.fixedVertex + 1)];
        v_Fy[i] = -vitelline_EG_y[i +(p.fixedVertex + 1)];
    }

    return;
}

/// energy
double System:: energy(void)
{
    double edgeContribution_apical(0.0),edgeContribution_basal(0.0),edgeContribution_lateral(0.0),cellContribuion(0.0),vitellineContribution(0.0),yolkContribution(0.0);
    for(int i = 0; i < numNode; i++)
    {
        // contribution-from-apical-edge
        edgeContribution_apical +=  geometry->CELL[i].alphA*geometry->CELL[i].l_a;
        // contribution-from-basal-edge
        edgeContribution_basal += geometry->CELL[i].alphB*geometry->CELL[i].l_b;
        // contribution-from-lateral-edge
        edgeContribution_lateral += lateral_Edg[i].alphL*lateral_Edg[i].l_l;
        // contribution-from-cell-area
        cellContribuion += p.betaCell*(geometry->CELL[i].A-p.cellRefArea)*(geometry->CELL[i].A-p.cellRefArea);
        // contribution-from-vitelline-membrane
        vitellineContribution += p.epsilon*HF(Va[i].D)*Va[i].w*Va[i].D*Va[i].D;
    }
    // contribution-from-yolk
    yolkContribution = -p.pressure*geometry->A_Y + p.betaYolk*(geometry->A_Y - geometry->A_Y_ref)*(geometry->A_Y - geometry->A_Y_ref);
    // sum-of-all-contributions
    double totalEnergy = edgeContribution_apical + edgeContribution_basal + edgeContribution_lateral + cellContribuion + vitellineContribution + yolkContribution;

    return(totalEnergy);
}

double System:: checkEnergy(void)
{
    double edgeContribution_apical(0.0),edgeContribution_basal(0.0),edgeContribution_lateral(0.0),cellContribuion(0.0),vitellineContribution(0.0),yolkContribution(0.0);
    for(int i = 0; i < numNode; i++)
    {
        // contribution-from-apical-edge
        edgeContribution_apical +=  geometry->CELL[i].alphA*(geometry->CELL[i].l_a + lambda*geometry->CELL[i].m_a);
        // contribution-from-basal-edge
        edgeContribution_basal += geometry->CELL[i].alphB*(geometry->CELL[i].l_b + lambda*geometry->CELL[i].m_b);
        // contribution-from-lateral-edge
        edgeContribution_lateral += lateral_Edg[i].alphL*(lateral_Edg[i].l_l + lambda*lateral_Edg[i].m_l);
        // contribution-from-cell-area
        cellContribuion += p.betaCell*(geometry->CELL[i].A-p.cellRefArea)*(geometry->CELL[i].A-p.cellRefArea);
        // contribution-from-vitelline-membrane
        vitellineContribution += p.epsilon*HF(Va[i].D)*Va[i].w*Va[i].D*Va[i].D;
    }
    // contribution-from-yolk
    yolkContribution = -p.pressure*geometry->A_Y+ p.betaYolk*(geometry->A_Y - geometry->A_Y_ref)*(geometry->A_Y - geometry->A_Y_ref);
    // sum-of-all-contributions
    double totalEnergy = edgeContribution_apical + edgeContribution_basal + edgeContribution_lateral + cellContribuion + vitellineContribution + yolkContribution;

    return(totalEnergy);
}

bool System::checkEnergyGradient(void)
{
    bool valid(true);
    double tolerance(p.enegryTolerance), coordinateShift(CgnitialStepSizePerDof), oldEnergy(checkEnergy()), newEnergy(0.0), diff(0.0);
    for(int i = 0; i < numNode; i++)
    {
        // check: energy-gradient-at-apical-vertices: w.r.t. -> x
        Va[i].C(0,0) = Va[i].C(0,0) + coordinateShift;
        updateGeometry();
        newEnergy = checkEnergy();
        double energyGradient_apical_x  = (newEnergy - oldEnergy)/coordinateShift;
        diff = energyGradient_apical_x - Va[i].EG(0,0);
        Va[i].C(0,0) = Va[i].C(0,0) - coordinateShift;
        valid = diff < tolerance ? true : false;
        if(! valid)
        {
            cout << "vertex-index: " << i << " (" << energyGradient_apical_x << ") -> (" << Va[i].EG(0,0) << ")" << endl;
            break;
        }
        // check: energy-gradient-at-apical-vertices: w.r.t. -> y
        Va[i].C(1,0) = Va[i].C(1,0) + coordinateShift;
        updateGeometry();
        newEnergy = checkEnergy();
        double energyGradient_apical_y  = (newEnergy - oldEnergy)/coordinateShift;
        diff = energyGradient_apical_y - Va[i].EG(1,0);
        Va[i].C(1,0) = Va[i].C(1,0) - coordinateShift;
        valid = diff < tolerance ? true : false;
        if(! valid)
        {
            cout << i << " : (" << energyGradient_apical_y << ") -> (" << Va[i].EG(1,0) << ")" << endl;
            break;
        }
        // check: energy-gradient-at-basal-vertices: w.r.t. -> x
        Vb[i].C(0,0) = Vb[i].C(0,0) + coordinateShift;
        updateGeometry();
        newEnergy = checkEnergy();
        double energyGradient_basal_x  = (newEnergy - oldEnergy)/coordinateShift;
        diff = energyGradient_basal_x - Vb[i].EG(0,0);
        Vb[i].C(0,0) = Vb[i].C(0,0) - coordinateShift;
        valid = diff < tolerance? true : false;
        if(! valid)
        {
            cout << i << " : (" << energyGradient_basal_x << ") -> (" << Vb[i].EG(0,0) << ")" << endl;
            break;
        }
        // check: energy-gradient-at-basal-vertices: w.r.t. -> y
        Vb[i].C(1,0) = Vb[i].C(1,0) + coordinateShift;
        updateGeometry();
        newEnergy = checkEnergy();
        double energyGradient_basal_y  = (newEnergy - oldEnergy)/coordinateShift;
        diff = energyGradient_basal_y - Vb[i].EG(1,0);
        Vb[i].C(1,0) = Vb[i].C(1,0) - coordinateShift;
        valid = diff < tolerance? true : false;
        if(! valid)
        {
            cout << i << " : (" << energyGradient_basal_y << ") -> (" << Vb[i].EG(1,0) << ")" << endl;
            break;
        }
    }

    return valid;
}
