/***************************************************************************
 *   Copyright (C) 2016 by Саша Миленковић                                 *
 *   sasa.milenkovic.xyz@gmail.com                                         *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *   ( http://www.gnu.org/licenses/gpl-3.0.en.html )                       *
 *									   *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <cmath>
#include "algebra.h"

//---------------------------------------------------------------------------
// solve cubic equation x^3 + a*x^2 + b*x + c
// x - array of size 3
// In case 3 real roots: => x[0], x[1], x[2], return 3
//         2 real roots: x[0], x[1],          return 2
//         1 real root : x[0], x[1] ± i*x[2], return 1
unsigned int solveP3(double *x,double a,double b,double c) {
	double a2 = a*a;
    	double q  = (a2 - 3*b)/9;
	double r  = (a*(2*a2-9*b) + 27*c)/54;
    	double r2 = r*r;
	double q3 = q*q*q;
	double A,B;
    	if(r2<q3)
    	{
    		double t=r/sqrt(q3);
    		if( t<-1) t=-1;
    		if( t> 1) t= 1;
    		t=acos(t);
    		a/=3; q=-2*sqrt(q);
    		x[0]=q*cos(t/3)-a;
    		x[1]=q*cos((t+M_2PI)/3)-a;
    		x[2]=q*cos((t-M_2PI)/3)-a;
    		return 3;
    	}
    	else
    	{
    		A =-pow(fabs(r)+sqrt(r2-q3),1./3);
    		if( r<0 ) A=-A;
    		B = (0==A ? 0 : q/A);

		a/=3;
		x[0] =(A+B)-a;
		x[1] =-0.5*(A+B)-a;
		x[2] = 0.5*sqrt(3.)*(A-B);
		if(fabs(x[2])<eps) { x[2]=x[1]; return 2; }

		return 1;
        }
}

//---------------------------------------------------------------------------
// solve quartic equation x^4 + a*x^3 + b*x^2 + c*x + d
// Attention - this function returns dynamically allocated array. It has to be released afterwards.
DComplex* solve_quartic(double a, double b, double c, double d)
{
	double a3 = -b;
	double b3 =  a*c -4.*d;
	double c3 = -a*a*d - c*c + 4.*b*d;

	// cubic resolvent
	// y^3 − b*y^2 + (ac−4d)*y − a^2*d−c^2+4*b*d = 0

	double x3[3];
	unsigned int iZeroes = solveP3(x3, a3, b3, c3);

	double q1, q2, p1, p2, D, sqD, y;

	y = x3[0];
	// The essence - choosing Y with maximal absolute value.
	if(iZeroes != 1)
	{
		if(fabs(x3[1]) > fabs(y)) y = x3[1];
		if(fabs(x3[2]) > fabs(y)) y = x3[2];
	}

	// h1+h2 = y && h1*h2 = d  <=>  h^2 -y*h + d = 0    (h === q)

	D = y*y - 4*d;
	if(fabs(D) < eps) //in other words - D==0
	{
		q1 = q2 = y * 0.5;
		// g1+g2 = a && g1+g2 = b-y   <=>   g^2 - a*g + b-y = 0    (p === g)
		D = a*a - 4*(b-y);
		if(fabs(D) < eps) //in other words - D==0
			p1 = p2 = a * 0.5;

		else
		{
			sqD = sqrt(D);
			p1 = (a + sqD) * 0.5;
			p2 = (a - sqD) * 0.5;
		}
	}
	else
	{
		sqD = sqrt(D);
		q1 = (y + sqD) * 0.5;
		q2 = (y - sqD) * 0.5;
		// g1+g2 = a && g1*h2 + g2*h1 = c       ( && g === p )  Krammer
		p1 = (a*q1-c)/(q1-q2);
		p2 = (c-a*q2)/(q1-q2);
	}

	DComplex* retval = new DComplex[4];

	// solving quadratic eq. - x^2 + p1*x + q1 = 0
	D = p1*p1 - 4*q1;
	if(D < 0.0)
	{
		retval[0].real( -p1 * 0.5 );
		retval[0].imag( sqrt(-D) * 0.5 );
		retval[1] = std::conj(retval[0]);
	}
	else
	{
		sqD = sqrt(D);
		retval[0].real( (-p1 + sqD) * 0.5 );
		retval[1].real( (-p1 - sqD) * 0.5 );
	}

	// solving quadratic eq. - x^2 + p2*x + q2 = 0
	D = p2*p2 - 4*q2;
	if(D < 0.0)
	{
		retval[2].real( -p2 * 0.5 );
		retval[2].imag( sqrt(-D) * 0.5 );
		retval[3] = std::conj(retval[2]);
	}
	else
	{
		sqD = sqrt(D);
		retval[2].real( (-p2 + sqD) * 0.5 );
		retval[3].real( (-p2 - sqD) * 0.5 );
	}

    return retval;
}

bool intersectionLineSegments(Vector2d v1,Vector2d v2,Vector2d u1,Vector2d u2,Vector2d& interSecPoint)
{
    bool encounter(false);

    Vector2d u;
    u = u2-u1;

    Vector2d v;
    v = v2-v1;

    Vector2d w;
    w = u1-v1;

    double sval(0.0);
    sval = (v(1,0)*w(0,0)-v(0,0)*w(1,0))/(v(0,0)*u(1,0)-v(1,0)*u(0,0));

    double sval2(0.0);
    w = v1-u1;

    sval2 = (u(1,0)*w(0,0)-u(0,0)*w(1,0))/(u(0,0)*v(1,0)-u(1,0)*v(0,0));

    if(((sval>=0.0) && (sval<=1.0)) && ((sval2>=0.0) && (sval2<=1.0)))
    {
        encounter=true;
        interSecPoint << u1(0,0)+sval*(u2(0,0)-u1(0,0)),u1(1,0)+sval*(u2(1,0)-u1(1,0));
    }

    return(encounter);
}

int HF(double D_i)
{
    return (D_i <= 0)? 0 : 1;
}

double RobustLength(double v0, double v1)
{
    double length(sqrt(v0*v0+v1*v1));

    if (fabs(v0) == std::max(fabs(v0),fabs(v1)))
        double length = fabs(v0)*sqrt(1 + (v1/v0)*(v1/v0));

    if (fabs(v1) == std::max(fabs(v0),fabs(v1)))
        double length = fabs(v1)*sqrt(1 + (v0/v1)*(v0/v1));

    return length;
}

double GetRoot(double r0, double z0, double z1, double g, int maxIterations)
{
    double RobustLength(double, double);
    double n0 = r0 * z0 ;
    double s0 = z1 - 1.0;
    double s1 =  g < 0 ? 0 : RobustLength(n0,z1) - 1.0;
    double  s = 0.0;

    for (int i = 0; i < maxIterations; ++i )
    {
        s = (s0 + s1)/2;

        if (s == s0 || s == s1)
        break;

        double ratio0 = n0/(s + r0);
        double ratio1 = z1/(s + 1) ;
        g = ratio0*ratio0 + ratio1*ratio1 - 1;

        if (g > 0)
        s0 = s;

        else if (g < 0)
        s1 = s;

        else
        break;
    }
    return s;
}

double DistancePoints(double x1,double x2,double y1, double y2)
{
    double distance = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
    return distance;
}

double DistancePointEllipse(double semi_a,double semi_b,Array2d targetPoint,Array2d& intersectionPoint,int maxIterations)
{
    // change-in-sign: target-point-not-in-first-quadrant
    Array2d targetPoint_sign(sgn(targetPoint(0,0)),sgn(targetPoint(1,0)));
    targetPoint*=targetPoint_sign;
    // distance-and-intersection-point
    double distance(0.0);
    double v_x = targetPoint(0,0);
    double v_y = targetPoint(1,0);
    if (v_y > 0)
    {
        if (v_x > 0)
        {
            double z0 = v_x / semi_a;
            double z1 = v_y / semi_b;
            double g = z0*z0 + z1*z1 - 1;
            if (g != 0)
            {
                double r0 = (semi_a/semi_b)*(semi_a/semi_b);
                double sbar = GetRoot(r0,z0,z1,g,maxIterations);
                double I_x = r0*v_x/(sbar + r0);
                double I_y = v_y/(sbar + 1);
                distance = sqrt((I_x - v_x)*(I_x - v_x) + (I_y - v_y)*(I_y - v_y));
                intersectionPoint << I_x,I_y;
            }
            else
            {
                double I_x = v_x ;
                double I_y = v_y ;
                distance = 0;
                intersectionPoint << I_x,I_y;
            }
        }
        else // v_x == 0
        {
            double I_x = 0 ;
            double I_y = semi_b ;
            distance = fabs( v_y - semi_b);
            intersectionPoint << I_x,I_y;
        }
    }
    else // v_y == 0
    {
        double numer0 = semi_a*v_x;
        double denom0 = semi_a*semi_a - semi_b*semi_b;
        if (numer0 < denom0)
        {
            double xdsemi_a = numer0/denom0;
            double I_x = semi_a * xdsemi_a ;
            double I_y = semi_b * sqrt(1.0 - xdsemi_a * xdsemi_a);
            distance = sqrt((I_x - v_x)*(I_x - v_x) + I_y*I_y);
            intersectionPoint << I_x,I_y;
        }
        else
        {
            double I_x = semi_a;
            double I_y = 0;
            distance = fabs (v_x - semi_a);
            intersectionPoint << I_x,I_y;
        }
    }
    // change-in-sign: target-point-not-in-first-quadrant
    intersectionPoint*=targetPoint_sign;

    return distance;
}

// Function to check the point
int isPointInsideEllipse(double x, double y, double a, double b)
{
    double p = (pow(x,2)/pow(a,2)) + (pow(y,2)/pow(b,2));
    if (p > 1)
        return 1;
    else
    return -1;
}

int isPointBelowHorizontalLine(double y1,double y2)
{
    double p = y1 - y2;
    if (p > 0)
        return 1;
    else
        return -1;
}

