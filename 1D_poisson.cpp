#include "poisson_1D.h"


void phi_setup( MultiBlockLattice2D<T,ADESCRIPTOR>& phiLattice,
                SimulationParams<T> & simParams)
{
    plint nx = simParams.getNx();
    plint ny = simParams.getNy();
    Box2D south_boundary(0,nx-1,0,0);
    Box2D north_boundary(0,nx-1,ny-1,ny-1);

    plint processorLevel = 5;

    initializeAtEquilibrium(    phiLattice,
                                phiLattice.getBoundingBox(),
                                simParams.getPhi_init(),
                                Array<T,2>((T)0.0,(T)0.0));

    OnLatticeAdvectionDiffusionBoundaryCondition2D<T, ADESCRIPTOR>* 
        bc = createLocalAdvectionDiffusionBoundaryCondition2D<T,ADESCRIPTOR>();

    bc->addTemperatureBoundary1N(south_boundary,phiLattice);
    setBoundaryDensity(phiLattice,south_boundary,(T)1.0);

    bc->addTemperatureBoundary1N(north_boundary,phiLattice);
    setBoundaryDensity(phiLattice,north_boundary,(T)1.0);

    integrateProcessingFunctional(  new poisson_1D<T,ADESCRIPTOR>(simParams),
                                    phiLattice.getBoundingBox(),
                                    phiLattice,
                                    processorLevel);


    phiLattice.initialize();

    delete bc;
}


int main(int argc, char* argv[])
//{{{
    {
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    // Instantiate variables
    SimulationParams<T> simParams = assign_params("params.xml");

    IncomprFlowParam<T> parameters( 1.0,
                                    1.0,
                                    simParams.getResolution(),
                                    simParams.getLx(),
                                    simParams.getLy());

    T phiOmega = 1.0/simParams.getTau_phi();
    pcout << "phiOmega = " << phiOmega << endl;
    plint nx = simParams.getNx();
    plint ny = simParams.getNy();

    writeLogFile(parameters,simParams,"output testing square domain");

    // Instantiate Charge Carrier Lattice
    MultiBlockLattice2D<T,ADESCRIPTOR> phiLattice(
                                    nx,
                                    ny,
                                    new ADYNAMICS<T, ADESCRIPTOR>(phiOmega));

    phi_setup(  phiLattice,
                simParams);

    plint convergenceIter = simParams.getConvergenceIter();
    // =========================================================================
    // Time Stepping
    // =========================================================================
    for (plint iT = 0; iT<=simParams.getMaxIter(); iT++)
    {
        if (iT%convergenceIter==0) {
            pcout << "Writing vtk_instance at iT = " << iT << endl;
            writeVTK(   phiLattice,
                        parameters,
                        simParams,
                        iT,
                        "phi Lattice");
        }
        phiLattice.collideAndStream();
    }
    plb_ofstream succesiveProfiles("concentration_final.dat");
    succesiveProfiles << std::setprecision(7)
        << *computeDensity(phiLattice,Box2D(1,1,0,simParams.getResolution())) << endl <<endl;
}
// }}}
