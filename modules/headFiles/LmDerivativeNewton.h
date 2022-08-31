#ifndef LMDERIVATIVEBISECTION_H
#define LMDERIVATIVEBISECTION_H

#include <iostream>
#include <functional>
#include <math.h>
#include <sstream>
#include <exception>

class LmDerivativeNewton
{
	public:

        /// constructor
        LmDerivativeNewton(void) : _c2(1e-4),
                                   _minMaxAbsSlope(1e-10),
                                   _maxInitialAlphaFactor(2),
                                   _factorSearchForMaxAlpha(1.5),
                                   _maxIterations(100),
                                   _debug(false),
                                   _iterations(0)
        {

        }

        /// initialize c2 (second, strong wolfe condition)
        void setC2(double c2)
        {
            _c2 = c2;
        }

        /// initialize slope
        void setMinMaxAbsSlope(double minMaxAbsSlope)
        {
            _minMaxAbsSlope = minMaxAbsSlope;
        }

        /// initialize debug
        void debug(void)
        {
            _debug = true;
        }

        /// get number of interations
        long iterations(void) const
        {
            return _iterations;
        }

        double minimize(double lastAlphaTimesLastSlope,
                        double initialSlope,
                        const function<double(double)> &gradient
                        )
        {
	    	/// set abort criterion
            double mas = _c2*fabs(initialSlope);
	    	const double MaxAbsSlope = (mas >_minMaxAbsSlope) ? mas : _minMaxAbsSlope;
	    	const double InitialAlpha = lastAlphaTimesLastSlope/fabs(initialSlope);

	    	const double MaxSecondAlpha = _maxInitialAlphaFactor*InitialAlpha;

            struct Point
            {
                double alpha, slope;

                /// constructor 1
                Point(double a, const function<double(double)> &gradient)
                {
                     alpha = a;
                     slope = gradient(alpha);
                }

                /// constructor 2
      			Point(double a, double s)
      			{
                     alpha = a;
                     slope = s;
      			}

      			/// set the values of alpha/slope
      			void setAlpha(double a, const function<double(double)> &gradient)
      			{
                    alpha = a;
                    slope = gradient(alpha);
                }

      			bool check(const double MaxAbsSlope) const
      			{
                    return fabs(slope) < MaxAbsSlope;
                }

      			bool operator<(const Point &o) const
      			{
                    return alpha < o.alpha;
                }

      			bool operator>(const Point &o) const
      			{
                    return alpha > o.alpha;
                }

      			Point operator-(const Point &o)
      			{
                    Point tmp(*this);
                    tmp.alpha -= o.alpha;
                    tmp.slope -= o.slope;
                    return tmp;
                }

      			void swap(Point &o)
      			{
                    Point temp = *this;
                    *this = o;
                    o = temp;
                }
    		};

            /// minimum-point
    		Point minimum(0, initialSlope);

    		/// middle-point
    		Point middle(InitialAlpha, gradient);

    		if(middle.check(MaxAbsSlope))
    		{
                return middle.alpha;
    		}

    		/// maximum-point
    		double maxAlpha = fabs(middle.alpha*minimum.slope/(middle.slope-minimum.slope)); // linear inter-/extrapolation of slope to get alpha max
    		if(maxAlpha > MaxSecondAlpha)
    		maxAlpha = MaxSecondAlpha;
    		Point maximum(maxAlpha, gradient);
    		if(maximum.check(MaxAbsSlope))
    		{
                return maximum.alpha;
    		}

    		/// ordering: minimum-alpha < middle-alpha < maximum-alpha
    		if(middle > maximum)
    		middle.swap(maximum);

    		/// requirement: at least one slope (middle-slope or maximum-slope) should be > 0
    		if(middle.slope < 0)
            {
      			while(maximum.slope < 0)
                {
        			double maxAlpha = maximum.alpha;
        			minimum = middle;
        			middle = maximum;
       	 			maximum.setAlpha(maxAlpha*_factorSearchForMaxAlpha, gradient);

        			if(maximum.check(MaxAbsSlope))
        			{
                        return maximum.alpha;
                    }
      			};
    		}

    		//////////////////////////////////////
    		/// start loop through iterations ///
    		/////////////////////////////////////
    		_iterations = 0;
    		while(true)
            {
      			if(_iterations > _maxIterations)
                {
        			cout << "LmDerivativeNewton::minimize:  maximum number of iterations!\n";
        			cout << "alpha : " << middle.alpha << " " << _iterations << " " << _maxIterations << std::endl;
        			return middle.alpha;
      			}

      			/// parabolic interpolation of slope (relative to minimum.alpha)
      			Point deltaMid = middle - minimum;
      			Point deltaMax = maximum - minimum;
      			double deltaAlphaSqMid = deltaMid.alpha*deltaMid.alpha;
      			double deltaAlphaSqMax = deltaMax.alpha*deltaMax.alpha;
      			double denominator = deltaMid.alpha*deltaAlphaSqMax - deltaMax.alpha*deltaAlphaSqMid;
      			if(denominator == 0)
                {
        			return middle.alpha;
      			}
      			double c = (deltaMax.slope*deltaMid.alpha - deltaMid.slope*deltaMax.alpha)/denominator;
      			double b = (deltaMid.slope*deltaAlphaSqMax - deltaMax.slope*deltaAlphaSqMid)/denominator;
      			double a = minimum.slope;
      			double offset = 0.5*b/c;
      			double det = offset*offset - a/c;
      			if(det<-1e-14*offset*offset)
                {
                     cout << "LmDerivativeNewton::minimize:  Got negative determinant: " << det << ", quotient: " << det/(offset*offset) << "\n";
      			}
      			double sqrtDet = (det <= 0) ? 0 : sqrt(det);
      			double newAlpha = (c > 0) ? -offset + sqrtDet : -offset - sqrtDet;

      			/// sub-middle-point
      			Point subMiddle(minimum.alpha + newAlpha, gradient);
      			if(subMiddle.check(MaxAbsSlope))
                {
        			return subMiddle.alpha;
      			}

      			if(_debug)
                {
        			cout << "checkMid: " << a + b*middle.alpha + c*middle.alpha*middle.alpha << "\n";
        			cout << "checkMax: " << a + b*maximum.alpha + c*maximum.alpha*maximum.alpha << "\n";
        			cout << "cutoff slope:  " << MaxAbsSlope << "\n";
        			cout << "c:  " << c << "\n";
        			cout << "minimum.alpha*c:  " << minimum.alpha*c << "\n";
        			cout << "b:  " << b << "\n";
        			cout << "minimum.alpha*(b + minimum.alpha*c):  " << minimum.alpha*(b + minimum.alpha*c) << "\n";
        			cout << "minimum.slope:  " << minimum.slope << "\n";
        			cout << "a:  " << a << "\n";
        			cout << "offset:  " << offset << "\n";
        			cout << "offset*offset:  " << offset*offset << "\n";
        			cout << "a/c:  " << a/c << "\n";
        			cout << "sqrtDet:  " << sqrtDet << "\n";
        			cout << "check: " << a + b*newAlpha + c*newAlpha*newAlpha << "\n";
        			cout << "min abs:  " << minimum.alpha << ", slope: " << minimum.slope << "\n";
        			cout << "smid: " << subMiddle.alpha-minimum.alpha << ", slope: " << subMiddle.slope  << ", rel. slope: " << subMiddle.slope-minimum.slope  << "\n";
        			cout << "mid:  " << middle.alpha-minimum.alpha << ", slope: " << middle.slope << ", rel. slope: " << middle.slope-minimum.slope << "\n";
        			cout << "max:  " << maximum.alpha-minimum.alpha << ", slope: " << maximum.slope << ", rel. slope: " << maximum.slope-minimum.slope << "\n";
      			}

      			/// ordering: sub-middle-alpha < middle-alpha
      			if(subMiddle > middle)
      			subMiddle.swap(middle);

      			/// kick out (minimum or maximum):  requirement at least one-slope > 0, which is maximum
      			if(subMiddle.slope > 0)
                {
        			// kick out maximum
        			maximum = middle;
        			middle = subMiddle;
      			}
                else if(middle.slope < 0)
                {
        			// kick out minimum
        			minimum = subMiddle;
      			}
                else
                {
        			/// slope of two above (minimum-point & sub-middle-point) and slope of two below (middle-point & maximum-point): solution ? -> compare slopes
        			if(fabs(minimum.slope) < fabs(maximum.slope))
                    {
         				// kick out maximum: reason minimum is closer to zero
          				maximum = middle;
          				middle = subMiddle;
       				}
                    else
                    {
          				// kick out minimum: reason maximum is closer to zero
          				minimum = subMiddle;
        			}
      			}
      			++_iterations;
    		}

    		return middle.alpha;
  	}

	private:
  		double _c2;  // second, strong wolfe condition
  		double _minMaxAbsSlope;
  		double _maxInitialAlphaFactor;
  		double _factorSearchForMaxAlpha;
  		long _maxIterations;
  		bool _debug;
  		long _iterations;
};

#endif /* LMDERIVATIVEBISECTION_H */

