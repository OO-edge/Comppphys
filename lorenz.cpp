#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <boost/array.hpp>
#include <boost/numeric/odeint/integrate/integrate.hpp>
#include <sciplot/sciplot.hpp>

using namespace sciplot;
using namespace std;
using namespace boost::numeric::odeint;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

std::vector<double> Xpos;
std::vector<double> Ypos;
std::vector<double> Zpos;

typedef array< double , 3 > state_type;

void lorenz( const state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = -b * x[2] + x[0] * x[1];
}

void write_lorenz( const state_type &x , const double t )
{
    Xpos.push_back(x[0]); Ypos.push_back(x[1]); Zpos.push_back(x[2]);
}

int main(int argc, char **argv)
{
    state_type x = {{ 10.0 , 1.0 , 1.0 }}; // initial conditions
    integrate( lorenz , x , 0.0 , 40.0 , 0.1 , write_lorenz );

    Plot3D plot; plot.xlabel("x"); plot.ylabel("y"); plot.zlabel("z");
    plot.border().clear();
    plot.border().bottomLeftFront();
    plot.border().bottomRightFront();
    plot.border().leftVertical();

    plot.autoclean(false);

    plot.palette("dark2");

    plot.drawCurve(Xpos, Ypos, Zpos).label("Lorentz attractor").lineColor("black");

    Figure fig = {{plot}}; Canvas canvas = {{fig}}; canvas.size(600, 600);
    canvas.show();
}
