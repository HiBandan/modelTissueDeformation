#ifndef CGMINIMIZER_INLINE_H
#define CGMINIMIZER_INLINE_H

#include <iostream>
#include <math.h>
#include "CgMinimizer.h"
#include <fstream>
#include "Eigen/Dense"
using namespace Eigen;

template <typename LineMinimizer, typename SearchDirectionUpdate>
bool CgMinimizer<LineMinimizer, SearchDirectionUpdate>::minimize(bool checkIntegrity,
                                                                 std::function<void()> updateGradient,
                    						                     std::function<void()> saveFrame,
                                                                 std::function<bool()> checkGradient
                                                                 )
{
    // gradient tolerance for all the nodes
	const double TotalGradientSqCutoff = _dofs.size()*_gradientTolerancePerDof*_gradientTolerancePerDof;

    // maximum interations for all the node
  	const int MaxNumberOfIterations = _dofs.size()*_maxIterationsPerDof;

  	// initialize step-size for all the nodes
  	const double InitialStepSize = _initialStepSizePerDof*sqrt(_dofs.size());

  	// get the gradients at global (system) nodes ( _gradient )
  	updateGradient();

    // check-integrity
    if(checkIntegrity)
    {
        if(!checkGradient())
        {
            _errorCode = Aborted;
            return false;
        }
    }

    // create & initialize local nodes : positions and gradients
  	Eigen::VectorXd position = Eigen::VectorXd::Zero(_dofs.size());
  	Eigen::VectorXd gradient = Eigen::VectorXd::Zero(_dofs.size());

  	for(unsigned int i = 0; i <_dofs.size(); ++i)
	{
    		position(i) = *_dofs[i];
    		gradient(i) = *_gradient[i];
  	}

  	// pass local node gradients to search director (i.e. PolakRibiere)
  	SearchDirectionUpdate searchDirectionUpdate(gradient);

  	// assign next step size
  	double lastStepSizeTimesSlope = InitialStepSize*fabs(searchDirectionUpdate.searchDirectionUnitVector().dot(gradient));

  	// start interation loop for energy minimization
  	int iteration = 0;
  	while(gradient.squaredNorm() > TotalGradientSqCutoff)
	{
        if(iteration > MaxNumberOfIterations)
		{
      			std::cerr << "CgMinimizer::minimize: reached max iterations!" << std::endl;
      			_errorCode = ExceededMaxIterations;
      			return false;
        }

        // line minimization
        Eigen::VectorXd searchDirection = searchDirectionUpdate.searchDirectionUnitVector();
        double slope = searchDirection.dot(gradient);

        if(slope > 0)
		{
      			//std::cout << "CgMinimizer::minimize: positive slope before line minimization, flipping direction!" << std::endl;
      			searchDirection *= -1;
      			slope *= -1;
        }

        double alpha = _lineMinimizer.minimize
		(
            lastStepSizeTimesSlope,

            slope,

            [&](double alpha)->double
			{
                        // update global (system) node coordinates ( _dofs )
              			for(unsigned int i = 0; i<_dofs.size(); ++i)
              			*_dofs[i] = position(i) + alpha*searchDirection(i);

              			// update gradients at the global (system) nodes ( _gradient )
              			updateGradient();

              			double derivative = 0.0;
              			for(unsigned int i = 0; i<_dofs.size(); ++i)
              			derivative += *_gradient[i]*searchDirection(i);

                        // energy derivatives
              			return derivative;
            }
        );

        // after current line search: update of local position coordinates
        lastStepSizeTimesSlope = alpha*fabs(slope);
        position += alpha*searchDirection;

        // after current line search: update global (system) node coordinates ( _dofs )
        for(unsigned int i = 0; i <_dofs.size(); ++i)
        *_dofs[i] = position(i);

        // after current line search: update gradients at the global (system) nodes ( _gradient )
        updateGradient();

        // check-integrity
        if(checkIntegrity)
        {
            if(!checkGradient())
            {
                _errorCode = Aborted;
                return false;
            }
        }

        // next line search: get new search direction at the nodes
        Eigen::VectorXd newGradient = Eigen::VectorXd::Zero(_dofs.size());
        for(unsigned int i = 0; i <_dofs.size(); ++i)
		newGradient(i) = *_gradient[i];
        searchDirectionUpdate.update(alpha, newGradient, gradient);
        // save movie frame (export vertex coordinates)
        if((_logInterval > 0) && (iteration%_logInterval == 0))
        {
            if(iteration/_logInterval < _logNumber)
            saveFrame();

            cout << "intermediate frame number: " << iteration/_logInterval << endl;

            if(iteration/_logInterval > _logNumber)
            exit(0);
        }

   		// next line search: update local gradient with new gradient
        gradient = newGradient;

        // next line search: increment iteration step
        ++iteration;
	}

	// minimization is successful, update error code
  	_errorCode = Success;

  	return true;
}

#endif /* CGMINIMIZER_INLINE_H */

