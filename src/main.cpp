#include <string>

#include "../include/mynurbs.h"
#include "../include/GLProgram.h"

#define WINDOW_WIDTH 1600
#define WINDOW_HEIGHT 1200

// declare static members for use in callback functions
int GLProgram::windowWidth = WINDOW_WIDTH;
int GLProgram::windowHeight = WINDOW_HEIGHT;
Camera GLProgram::camera;
bool GLProgram::mousePressed = false;
double GLProgram::prevMouseX, GLProgram::prevMouseY;
glm::mat4 GLProgram::modelMatrix = glm::mat4(1.0f);

typedef tinynurbs::RationalCurve3d RationalCurve3d;
typedef tinynurbs::RationalSurface3d RationalSurface3d;
typedef glm::dvec3 vec3d;

void printSurfaceInfo(const RationalSurface3d &resOfLoftAlgo);

void task11(vector<RationalCurve3d> &curves, vector<RationalSurface3d>& surfaces)
{
    RationalCurve3d crv1;
    crv1.control_points = {vec3d(0, 0, 5), vec3d(0, 3, 5), vec3d(0, 4, 5), vec3d(0, 5, 5), vec3d(0, 7, 5), vec3d(0, 11, 5), vec3d(0, 13, 5), vec3d(0, 15, 5)};
    crv1.knots = {0, 0, 0, 0, 0.2, 0.4, 0.6, 0.8, 1, 1, 1, 1};
    crv1.degree = 3;
    crv1.weights = {1, 1, 1, 1, 1, 1, 1, 1};

    RationalCurve3d crv2;
    crv2.control_points = {vec3d(2, 0, 5), vec3d(2, 3, 4), vec3d(2, 4, 3), vec3d(2, 5, 4), vec3d(2, 7, 6), vec3d(2, 11, 5), vec3d(2, 13, 3), vec3d(2, 15, 5)};
    crv2.knots = {0, 0, 0, 0, 0.2, 0.4, 0.6, 0.8, 1, 1, 1, 1};
    crv2.degree = 3;
    crv2.weights = {1, 1, 1, 1, 1, 1, 1, 1};

    RationalCurve3d crv3;
    crv3.control_points = {vec3d(3, 0, 5), vec3d(3, 3, 3), vec3d(3, 4, 1), vec3d(3, 5, 4), vec3d(3, 7, 7), vec3d(3, 11, 9), vec3d(3, 13, 7), vec3d(3, 15, 5)};
    crv3.knots = {0, 0, 0, 0, 0.2, 0.4, 0.6, 0.8, 1, 1, 1, 1};
    crv3.degree = 3;
    crv3.weights = {1, 1, 1, 1, 1, 1, 1, 1};

    RationalCurve3d crv4;
    crv4.control_points = {vec3d(5, 0, 5), vec3d(5, 3, 4), vec3d(5, 4, 3), vec3d(5, 5, 3), vec3d(5, 7, 2), vec3d(5, 11, 3), vec3d(5, 13, 4), vec3d(5, 15, 5)};
    crv4.knots = {0, 0, 0, 0, 0.2, 0.4, 0.6, 0.8, 1, 1, 1, 1};
    crv4.degree = 3;
    crv4.weights = {1, 1, 1, 1, 1, 1, 1, 1};

    RationalCurve3d crv5;
    crv5.control_points = {vec3d(8, 0, 5), vec3d(8, 3, 5), vec3d(8, 4, 5), vec3d(8, 5, 5), vec3d(8, 7, 5), vec3d(8, 11, 5), vec3d(8, 13, 5), vec3d(8, 15, 5)};
    crv5.knots = {0, 0, 0, 0, 0.2, 0.4, 0.6, 0.8, 1, 1, 1, 1};
    crv5.degree = 3;
    crv5.weights = {1, 1, 1, 1, 1, 1, 1, 1};

    curves.emplace_back(crv1);
    curves.emplace_back(crv2);
    curves.emplace_back(crv3);
    curves.emplace_back(crv4);
    curves.emplace_back(crv5);

    vector<RationalCurve3d> taskCurves;
    taskCurves.emplace_back(crv1);
    taskCurves.emplace_back(crv2);
    taskCurves.emplace_back(crv3);
    taskCurves.emplace_back(crv4);
    taskCurves.emplace_back(crv5);
    
    RationalSurface3d resOfLoftAlgo;
    myNurbs::createLoftSurface(taskCurves, resOfLoftAlgo);

    assert(tinynurbs::surfaceIsValid(resOfLoftAlgo));
    //printSurfaceInfo(resOfLoftAlgo);
    surfaces.emplace_back(resOfLoftAlgo);
}
void task12(vector<RationalCurve3d> &curves, vector<RationalSurface3d>& surfaces)
{
    vector<RationalCurve3d> taskCurves;

    RationalCurve3d crv1;
    crv1.degree = 3;
    crv1.knots = {0, 0, 0, 0, 0.2, 0.4, 0.6, 0.8, 1, 1, 1, 1};
    crv1.control_points = {vec3d(0, 0, 0), vec3d(2, -2, -1), vec3d(4, -4, -2), vec3d(5, -2, -3), vec3d(7, -1, -4), vec3d(8, 0, -5), vec3d(10, 1, -6), vec3d(11, 0, -7)};
    crv1.weights = {1.0, 1.1, 1.0, 0.9, 0.9, 1.0, 1.0, 1.0};
    curves.emplace_back(crv1);
    taskCurves.emplace_back(crv1);

    crv1.control_points = {vec3d(0, 4, 3), vec3d(2, 5, 4), vec3d(5, 1, 3), vec3d(6, 2, 0), vec3d(8, 4, -4), vec3d(9, 3, -6), vec3d(11, 4, -8), vec3d(12, 4, -9)};
    curves.emplace_back(crv1);
    taskCurves.emplace_back(crv1);

    crv1.control_points = {vec3d(0, 8, 6), vec3d(3, 9, 7), vec3d(5, 8, 5), vec3d(6, 8, -1), vec3d(8, 9, -5), vec3d(10, 7, -7), vec3d(12, 9, -10), vec3d(13, 8, -12)};
    crv1.weights = {1.0, 1.2, 1.2, 1.0, 1.0, 1.0, 1.0, 1.0};
    curves.emplace_back(crv1);
    taskCurves.emplace_back(crv1);

    crv1.control_points = {vec3d(1, 12, 2), vec3d(3, 13, 4), vec3d(6, 14, 3), vec3d(7, 14, 0), vec3d(8, 13, -4), vec3d(10, 11, -5), vec3d(12, 14, -8), vec3d(13, 12, -10)};
    curves.emplace_back(crv1);
    taskCurves.emplace_back(crv1);

    crv1.control_points = {vec3d(2, 16, 0), vec3d(3, 17, -1), vec3d(6, 20, -2), vec3d(7, 19, -3), vec3d(8, 18, -4), vec3d(11, 17, -5), vec3d(13, 17, -6), vec3d(14, 16, -7)};
    crv1.weights = {1.0, 1.1, 1.1, 1.1, 1.1, 1.0, 1.0, 1.0};
    curves.emplace_back(crv1);
    taskCurves.emplace_back(crv1);

    RationalSurface3d resOfLoftAlgo;
    myNurbs::createLoftSurface(taskCurves, resOfLoftAlgo);

    assert(tinynurbs::surfaceIsValid(resOfLoftAlgo));
    //printSurfaceInfo(resOfLoftAlgo);
    surfaces.emplace_back(resOfLoftAlgo);
}
void task21(vector<RationalCurve3d> &curves, vector<RationalSurface3d>& surfaces)
{
    vector<RationalCurve3d> taskCurves;

    RationalCurve3d crv1;
    crv1.degree = 3;
    crv1.knots = {0, 0, 0, 0, 1. / 3, 1. / 3, 1. / 3, 2. / 3, 1, 1, 1, 1};
    crv1.control_points = {vec3d(0, 0, 5), vec3d(0, 3, 5), vec3d(0, 4, 5), vec3d(0, 5, 5), vec3d(0, 7, 5), vec3d(0, 11, 5), vec3d(0, 13, 5), vec3d(0, 15, 5)};
    crv1.weights = {1, 1, 1, 1, 1, 1, 1, 1};
    curves.emplace_back(crv1);
    taskCurves.emplace_back(crv1);

    crv1.knots = {0, 0, 0, 0, 1. / 5, 1. / 4, 1. / 3, 1. / 2, 1, 1, 1, 1};
    crv1.control_points = {vec3d(2, 0, 5), vec3d(2, 3, 4), vec3d(2, 4, 3), vec3d(2, 5, 4), vec3d(2, 7, 6), vec3d(2, 11, 5), vec3d(2, 13, 3), vec3d(2, 15, 5)};
    curves.emplace_back(crv1);
    taskCurves.emplace_back(crv1);

    crv1.knots = {0, 0, 0, 0, 1. / 5, 2. / 5, 3. / 5, 4. / 5, 1, 1, 1, 1};
    crv1.control_points = {vec3d(3, 0, 5), vec3d(3, 3, 3), vec3d(3, 4, 1), vec3d(3, 5, 4), vec3d(3, 7, 7), vec3d(3, 11, 9), vec3d(3, 13, 7), vec3d(3, 15, 5)};
    curves.emplace_back(crv1);
    taskCurves.emplace_back(crv1);

    crv1.knots = {0, 0, 0, 0, 1. / 2, 2. / 3, 3. / 4, 4. / 5, 1, 1, 1, 1};
    crv1.control_points = {vec3d(5, 0, 5), vec3d(5, 3, 4), vec3d(5, 4, 3), vec3d(5, 5, 3), vec3d(5, 7, 2), vec3d(5, 11, 3), vec3d(5, 13, 4), vec3d(5, 15, 5)};
    curves.emplace_back(crv1);
    taskCurves.emplace_back(crv1);

    crv1.knots = {0, 0, 0, 0, 1. / 3, 2. / 3, 2. / 3, 2. / 3, 1, 1, 1, 1};
    crv1.control_points = {vec3d(8, 0, 5), vec3d(8, 3, 5), vec3d(8, 4, 5), vec3d(8, 5, 5), vec3d(8, 7, 5), vec3d(8, 11, 5), vec3d(8, 13, 5), vec3d(8, 15, 5)};
    curves.emplace_back(crv1);
    taskCurves.emplace_back(crv1);

    RationalSurface3d resOfLoftAlgo;
    myNurbs::createLoftSurface(taskCurves, resOfLoftAlgo);

    assert(tinynurbs::surfaceIsValid(resOfLoftAlgo));
    //printSurfaceInfo(resOfLoftAlgo);
    surfaces.emplace_back(resOfLoftAlgo);
}
void task22(vector<RationalCurve3d> &curves, vector<RationalSurface3d>& surfaces)
{
    vector<RationalCurve3d> taskCurves;

    RationalCurve3d crv1;
    crv1.degree = 3;
    crv1.knots = {0, 0, 0, 0, 1. / 3, 1. / 3, 1. / 3, 2. / 3, 1, 1, 1, 1};
    crv1.control_points = {vec3d(0, 0, 0), vec3d(2, -2, -1), vec3d(4, -4, -2), vec3d(5, -2, -3), vec3d(7, -1, -4), vec3d(8, 0, -5), vec3d(10, 1, -6), vec3d(11, 0, -7)};
    crv1.weights = {1.0, 0.5, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    curves.emplace_back(crv1);
    taskCurves.emplace_back(crv1);

    crv1.knots = {0, 0, 0, 0, 1. / 5, 1. / 4, 1. / 3, 1. / 2, 1, 1, 1, 1};
    crv1.control_points = {vec3d(0, 2, 3), vec3d(2, 3, 4), vec3d(5, -1, 3), vec3d(6, 0, 0), vec3d(8, 2, -4), vec3d(9, 1, -6), vec3d(11, 2, -8), vec3d(12, 2, -6)};
    crv1.weights = {1.0, 0.6, 0.5, 0.5, 0.6, 0.7, 0.8, 1.0};
    curves.emplace_back(crv1);
    taskCurves.emplace_back(crv1);

    crv1.knots = {0, 0, 0, 0, 1. / 5, 2. / 5, 3. / 5, 4. / 5, 1, 1, 1, 1};
    crv1.control_points = {vec3d(0, 4, 6), vec3d(3, 5, 7), vec3d(5, 4, 5), vec3d(6, 4, -1), vec3d(8, 5, -5), vec3d(10, 3, -7), vec3d(12, 5, -10), vec3d(13, 4, -12)};
    crv1.weights = {1.0, 0.7, 0.6, 0.5, 0.5, 0.6, 0.7, 1.0};
    curves.emplace_back(crv1);
    taskCurves.emplace_back(crv1);

    crv1.knots = {0, 0, 0, 0, 1. / 2, 2. / 3, 3. / 4, 4. / 5, 1, 1, 1, 1};
    crv1.control_points = {vec3d(1, 6, 2), vec3d(3, 7, 4), vec3d(6, 8, 3), vec3d(7, 8, 0), vec3d(8, 7, -4), vec3d(10, 5, -5), vec3d(12, 8, -8), vec3d(13, 6, -10)};
    crv1.weights = {1.0, 0.8, 0.7, 0.6, 0.5, 0.5, 0.6, 1.0};
    curves.emplace_back(crv1);
    taskCurves.emplace_back(crv1);

    crv1.knots = {0, 0, 0, 0, 1. / 3, 2. / 3, 2. / 3, 2. / 3, 1, 1, 1, 1};
    crv1.control_points = {vec3d(2, 8, 0), vec3d(3, 9, -1), vec3d(6, 12, -2), vec3d(7, 11, -3), vec3d(8, 10, -4), vec3d(11, 9, -5), vec3d(13, 9, -6), vec3d(14, 8, -7)};
    crv1.weights = {1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.5, 1.0};
    curves.emplace_back(crv1);
    taskCurves.emplace_back(crv1);

    RationalSurface3d resOfLoftAlgo;
    myNurbs::createLoftSurface(taskCurves, resOfLoftAlgo);

    assert(tinynurbs::surfaceIsValid(resOfLoftAlgo));
    //printSurfaceInfo(resOfLoftAlgo);
    surfaces.emplace_back(resOfLoftAlgo);
}
void task31(vector<RationalCurve3d> &curves, vector<RationalSurface3d>& surfaces)
{
    vector<RationalCurve3d> taskCurves;

    RationalCurve3d crv;
    crv.degree = 1;
    crv.knots = {0, 0, 1, 1};
    crv.control_points = {vec3d(0, 0, 5), vec3d(0, 15, 5)};
    crv.weights = {1, 1};
    curves.emplace_back(crv);
    taskCurves.emplace_back(crv);

    crv.degree = 2;
    crv.knots = {0, 0, 0, 3. / 10, 1. / 2, 7. / 10, 1, 1, 1};
    crv.control_points = {vec3d(2, 0, 5), vec3d(2, 3, 4), vec3d(2, 5, 4), vec3d(2, 7, 6), vec3d(2, 13, 3), vec3d(2, 15, 5)};
    crv.weights = {1, 1, 1, 1, 1, 1};
    curves.emplace_back(crv);
    taskCurves.emplace_back(crv);

    crv.degree = 1;
    crv.knots = {0, 0, 1. / 5, 2. / 5, 3. / 5, 4. / 5, 1, 1};
    crv.control_points = {vec3d(3, 0, 5), vec3d(3, 3, 3), vec3d(3, 5, 4), vec3d(3, 7, 7), vec3d(3, 13, 7), vec3d(3, 15, 5)};
    curves.emplace_back(crv);
    taskCurves.emplace_back(crv);

    crv.degree = 2;
    crv.knots = {0, 0, 0, 3. / 10, 1. / 2, 7. / 10, 1, 1, 1};
    crv.control_points = {vec3d(5, 0, 5), vec3d(5, 3, 4), vec3d(5, 5, 3), vec3d(5, 7, 2), vec3d(5, 13, 4), vec3d(5, 15, 5)};
    curves.emplace_back(crv);
    taskCurves.emplace_back(crv);

    crv.degree = 3;
    crv.knots = {0, 0, 0, 0, 3. / 10, 7. / 10, 1, 1, 1, 1};
    crv.control_points = {vec3d(8, 0, 5), vec3d(8, 3, 5), vec3d(8, 5, 5), vec3d(8, 7, 5), vec3d(8, 13, 5), vec3d(8, 15, 5)};
    curves.emplace_back(crv);
    taskCurves.emplace_back(crv);

    RationalSurface3d resOfLoftAlgo;
    myNurbs::createLoftSurface(taskCurves, resOfLoftAlgo);

    assert(tinynurbs::surfaceIsValid(resOfLoftAlgo));
    //printSurfaceInfo(resOfLoftAlgo);
    surfaces.emplace_back(resOfLoftAlgo);
}
void task32(vector<RationalCurve3d> &curves, vector<RationalSurface3d>& surfaces)
{
    vector<RationalCurve3d> taskCurves;

    RationalCurve3d crv;
    crv.degree = 1;
    crv.knots = {0, 0, 1, 1};
    crv.control_points = {vec3d(0, 0, 0), vec3d(11, 0, -7)};
    crv.weights = {1, 1};
    curves.emplace_back(crv);
    taskCurves.emplace_back(crv);

    crv.degree = 2;
    crv.knots = {0, 0, 0, 3. / 10, 1. / 2, 7. / 10, 1, 1, 1};
    crv.control_points = {vec3d(0, 2, 3), vec3d(2, 3, 4), vec3d(6, 0, 0), vec3d(8, 2, -4), vec3d(11, 2, -8), vec3d(12, 2, -6)};
    crv.weights = {1, 0.9, 0.8, 0.8, 0.9, 1};
    curves.emplace_back(crv);
    taskCurves.emplace_back(crv);

    crv.degree = 1;
    crv.knots = {0, 0, 1. / 5, 2. / 5, 3. / 5, 4. / 5, 1, 1};
    crv.control_points = {vec3d(0, 4, 6), vec3d(3, 5, 7), vec3d(6, 4, -1), vec3d(8, 5, -5), vec3d(12, 5, -10), vec3d(13, 4, -12)};
    crv.weights = {1, 0.8, 0.6, 0.6, 0.8, 1};
    curves.emplace_back(crv);
    taskCurves.emplace_back(crv);

    crv.degree = 2;
    crv.knots = {0, 0, 0, 3. / 10, 1. / 2, 7. / 10, 1, 1, 1};
    crv.control_points = {vec3d(1, 6, 2), vec3d(3, 7, 4), vec3d(7, 8, 0), vec3d(8, 7, -4), vec3d(12, 8, -8), vec3d(13, 6, -10)};
    crv.weights = {1, 0.8, 0.6, 0.6, 0.8, 1};
    curves.emplace_back(crv);
    taskCurves.emplace_back(crv);

    crv.degree = 3;
    crv.knots = {0, 0, 0, 0, 3. / 10, 7. / 10, 1, 1, 1, 1};
    crv.control_points = {vec3d(2, 8, 0), vec3d(3, 9, -1), vec3d(7, 11, -3), vec3d(8, 10, -4), vec3d(13, 9, -6), vec3d(14, 8, -7)};
    crv.weights = {1, 0.9, 0.8, 0.8, 0.9, 1};
    curves.emplace_back(crv);
    taskCurves.emplace_back(crv);

    RationalSurface3d resOfLoftAlgo;
    myNurbs::createLoftSurface(taskCurves, resOfLoftAlgo);

    assert(tinynurbs::surfaceIsValid(resOfLoftAlgo));
    //printSurfaceInfo(resOfLoftAlgo);
    surfaces.emplace_back(resOfLoftAlgo);
}

int main()
{
    vector<RationalCurve3d> curves;
    vector<RationalSurface3d> surfaces;

    task11(curves, surfaces);
    task12(curves, surfaces);
    task21(curves, surfaces);
    task22(curves, surfaces);
    task31(curves, surfaces);
    task32(curves, surfaces);

    //RationalSurface3d resOfLoftAlgo;
    //myNurbs::createLoftSurface(curves, resOfLoftAlgo);
    //
    //assert(tinynurbs::surfaceIsValid(resOfLoftAlgo));
    //printSurfaceInfo(resOfLoftAlgo);

    //surfaces.emplace_back(resOfLoftAlgo);
    
    GLProgram program;
    program.init(curves, surfaces);
    program.setClearColor(0.05f, 0.18f, 0.25f, 1.0f);
    program.run(curves, surfaces);
    program.cleanup();
    return 0;
}

void printSurfaceInfo(const RationalSurface3d& resOfLoftAlgo)
{
   using namespace std;
   cout << "degree u v:" << resOfLoftAlgo.degree_u << ' ' << resOfLoftAlgo.degree_v << endl;
   cout << "knots u:";
   for each (const auto & ele in resOfLoftAlgo.knots_u)
   {
      cout << ' ' << ele;
   }
   cout << endl;
   cout << "knots v:";
   for each (const auto & ele in resOfLoftAlgo.knots_v)
   {
      cout << ' ' << ele;
   }
   cout << endl;
   cout << "ctrpnts.rows cols: " << resOfLoftAlgo.control_points.rows() << ' ' << resOfLoftAlgo.control_points.cols() << endl;
}
