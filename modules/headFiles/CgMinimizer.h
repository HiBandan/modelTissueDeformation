#ifndef CGMINIMIZER_H
#define CGMINIMIZER_H

#include "Eigen/Dense"
#include "MinimizerWithDerivative.h"


template <typename LineMinimizer, typename SearchDirectionUpdate>
class CgMinimizer : public MinimizerWithDerivative
{
	public:

        /// define a new data types: ErrorCode
  		enum ErrorCode { Success = 0, ExceededMaxIterations, Aborted, Undefined};

        /// constructor
  		CgMinimizer(LineMinimizer &lineMinimizer) : MinimizerWithDerivative(),
                                                    _lineMinimizer(lineMinimizer),
                                                    _gradientTolerancePerDof(1e-12),
                                                    _maxIterationsPerDof(100),
                                                    _initialStepSizePerDof(1e-3),
                                                    _errorCode(Undefined)
        {

        }

        /// specify gradient tolerance per node
  		void setGradientTolerancePerDof(double gradientTolerancePerDof)
		{
			_gradientTolerancePerDof = gradientTolerancePerDof;
		}

		/// specify maximum interations per node
  		void setMaxIterationsPerDof(int maxIterationsPerDof)
		{
			_maxIterationsPerDof = maxIterationsPerDof;
		}

		/// specify starting step-size per node
  		void setInitialStepSizePerDof(double initialStepSizePerDof)
		{
			_initialStepSizePerDof = initialStepSizePerDof;
		}

		/// minimizer (function body in CgMinimizer-inline.h) : true = "go on", false = "abort minimization"
  		virtual bool minimize
		( bool checkIntegrity,
          std::function<void()> updateGradient,
          std::function<void()> saveFrame = [](){},
          std::function<bool()> checkGradient = [](){return false;}
        );

  		ErrorCode errorCode() const
  		{
            return _errorCode;
  		}

        /// class variables
		private:
            LineMinimizer &_lineMinimizer;
            int _maxIterationsPerDof;
            double _initialStepSizePerDof;
            double _gradientTolerancePerDof;
            ErrorCode _errorCode;
};

#endif /* CGMINIMIZER_H */

