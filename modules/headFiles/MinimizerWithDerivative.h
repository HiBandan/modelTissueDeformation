#ifndef CONJUGATEDGRADIENT_H
#define CONJUGATEDGRADIENT_H

#include <vector>
#include <functional>
#include <fstream>

class MinimizerWithDerivative
{
	public:
        /// constructor
        MinimizerWithDerivative()
        {
            _logInterval = -1;
            _logNumber = -1;
        }

        /// initialize log interval
        void setLogInterval(int logInterval)
        {
            _logInterval = logInterval;
        }

        /// initialize log number
        void setLogNumber(int logNumber)
        {
            _logNumber = logNumber;
        }

        /// clean up degrees of freedoms (dofs)
        void clearDofs()
        {
            _dofs.clear();
            _gradient.clear();
        }

        /// add degrees of freedoms (dofs) via reference through addresses
        void addDof(double &dof, const double &derivative)
        {
            _dofs.push_back(&dof);
            _gradient.push_back(&derivative);
        }

        /// add degrees of freedoms (dofs) via reference through  pointers
        void addDof(double *dof, const double *derivative)
        {
            _dofs.push_back(dof);
            _gradient.push_back(derivative);
        }

        /// minimizer : true = "go on", false = "abort minimization"
        virtual bool minimize( bool checkIntegrity,
                               std::function<void()> updateGradient,
                               std::function<void()> saveFrame = [](){},
                               std::function<bool()> checkGradient = [](){ return false; }
                            ) = 0;
	/// class variables
	protected:
        int _logInterval;
        int _logNumber;
        std::vector<double*> _dofs;
        std::vector<const double*> _gradient;
};

#endif /* CONJUGATEDGRADIENT_H */

