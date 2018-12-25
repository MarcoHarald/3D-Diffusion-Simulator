#include <iostream>
#include <cmath>
using namespace std;

//float xAxis(-x0:x0),yAxis(-y0:y0),zAxis(-z0:z0)
//float dx_temp(-x0:x0,-y0:y0,-z0:z0),dy_temp(-x0:x0,-y0:y0,-z0:z0),dz_temp(-x0:x0,-y0:y0,-z0:z0)
//float tolerance, dt_temp(-x0:x0,-y0:y0,-z0:z0)

//4-dimensional array of REAL with size to be allocated
//float, DIMENSION(:,:,:,:), ALLOCATABLE :: array

//allocatable integers to define dimensions
int x_position,y_position,z_position, t_frame, stepSize;

int xLimit=10, yLimit=10, zLimit=10;
int x,y,z,k, steps, flag;

float bound;
float constDfsn = 1.1, initRad = 2.0, initRad2 = 4.0, stepDist = 1.0, stepTime = 1.0e-1;

//SET ACCURACY & TOLERANCE FOR BOUNDARY LIMITS
float tolerance = pow(10.0,-5.0);

//SET INITIAL GUESS AT BOUNDARY CONDITIONS (BC)
float temp1 = 20.0, temp0 = 100.0;

//SET SIMULATION LENGTH (s)
float t_siml = 10.0;

///DEFINE EXACT VALUE OF COORDINATES ON THREE AXES
///DEFINE EXACT VALUE OF COORDINATES ON THREE AXES

class axis{
        int limit;
    void setLimit(int limit);
    float setAxis() {
        float axis[2*limit];
        int i=-(limit);
        do {
            axis[i]=static_cast<int>(i)*stepDist;
            i++;
        } while(i<limit);
    }
};

void axis::setLimit (int x) {
    limit = x;
}

int main () {

    //set conditions
    axis xAxis;
    xAxis.setLimit (10);

    axis yAxis;
    yAxis.setLimit (10);

    axis yAxis;
    yAxis.setLimit (10);

    return 0;
}