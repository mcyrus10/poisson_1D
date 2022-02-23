# include "output.h"

void population(    Cell<T,ADESCRIPTOR> &cell,
                    T & density)
{
    density = 0.0;
    for (int j =0; j < 5; j++){
        density += cell[j]+ADESCRIPTOR<T>::t[j];
    }
}

template<typename T, template<typename U> class Descriptor>
class poisson_1D : public BoxProcessingFunctional2D_L<T,Descriptor>
{
private:
    SimulationParams<T> simParams;
    T k = 27.79;                // Below Equation 3.2
    T bar_omega_i[5][1] = { {0.0},
                            {(T)1.0/4.0},
                            {(T)1.0/4.0},
                            {(T)1.0/4.0},
                            {(T)1.0/4.0}
                           };
    T dx, dt,conductivity, u, alpha, R;
    plint resolution;


public:
    poisson_1D(   SimulationParams<T> simParams_): 
                    simParams(simParams_)
    {
        alpha = (T)1.0/2.0;                                                 //
        resolution = simParams.getResolution();
        dx = 1.0/(T)resolution;
        // ---------------------------------- Why?
        dt = dx*dx;
        // ---------------------------------- Why?
        conductivity = (1.0/3.0)*(simParams.getTau_phi()-0.5);                    // below equation 2.1
        u = 0.0;
        pcout << "dt = " << dt << endl;
    }
    virtual void process(   Box2D domain,
                            BlockLattice2D<T,Descriptor>& distribution)
    {
        Dot2D offset = distribution.getLocation();
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                // Omega_i' = Delta t * bar-omega_i * R * D
                Cell<T,ADESCRIPTOR>& cell = distribution.get(iX,iY);
                //pcout << cell[0] << endl;
                population(cell,u);
                plint xGlobal = iX+offset.x;
                plint yGlobal = iY+offset.y;
                R = k*k*u;                                                  // Equation 3.2
                //pcout << "iX,iY global = " << xGlobal << "," << yGlobal << "; conductivity = " << conductivity << "; k = " << k << "; R = " << R << "; u = " << u << endl;
                for (int i = 0; i<5; i++){
                    cell[i] = (cell[i]+ADESCRIPTOR<T>::t[i]) - dt*bar_omega_i[i][0]*R*conductivity - ADESCRIPTOR<T>::t[i];   // Equation 2.1
                }
                }
            }
    }

    virtual poisson_1D<T,Descriptor>* clone() const {
        return new poisson_1D(*this);
    }

    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }
    
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
};
