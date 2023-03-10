#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <sciplot/sciplot.hpp>

using namespace sciplot;
using namespace std;
using namespace boost::numeric::odeint;

std::vector<double> Time;
std::vector<double> Position;

/* we solve the simple ODE x' = 3/(2t^2) + x/(2t)
 * with initial condition x(1) = 0.
 * Analytic solution is x(t) = sqrt(t) - 1/t
 */

void rhs( const double x , double &dxdt , const double t )
{
    dxdt = 3.0/(2.0*t*t) + x/(2.0*t);
}

void write_cout( const double &x , const double t )
{
    cout << t << '\t' << x << endl;
    Time.push_back(t); Position.push_back(x);
}

// state_type = double
typedef runge_kutta_dopri5< double > stepper_type;

int main()
{
    double x = 0.0;
    integrate_adaptive( make_controlled( 1E-12 , 1E-12 , stepper_type() ) ,
                        rhs , x , 1.0 , 10.0 , 0.1 , write_cout );

    Plot2D plot; plot.xlabel("time"); plot.ylabel("position");
    plot.legend().atOutsideBottom().displayHorizontal().displayExpandWidthBy(2);
    plot.drawCurve(Time, Position).label("harmonic");
    Figure fig = {{plot}};
    Canvas canvas = {{fig}};
    canvas.show();
}
