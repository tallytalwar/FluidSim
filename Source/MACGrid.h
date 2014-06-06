#ifndef MACGrid_H_
#define MACGrid_H_

#pragma warning(disable: 4244 4267 4996)
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
using namespace boost::numeric;

#include <windows.h>
#include "GL/gl.h"
#include "vec.h"
#include "GridData.h"

class Camera;
class MACGrid
{
   friend MACGrid;
public:
   MACGrid();
   ~MACGrid();
   MACGrid(const MACGrid& orig);
   MACGrid& operator=(const MACGrid& orig);

   void reset();

   void draw(const Camera& c);
   void updateSources();
   void advectVelocity(double dt);
   void addExternalForces(double dt);
   void project(double dt);
   void advectTemperature(double dt);
   void advectDensity(double dt);

protected:

   // Setup
   void initialize();

   // Simulation
   void computeBouyancy(double dt);
   void computeVorticityConfinement(double dt);

   // Rendering
   struct Cube { vec3 pos; vec4 color; double dist; };
   void drawWireGrid();
   void drawMySolids(const Camera& c);
   void drawSmokeCubes(const Camera& c);
   void drawSmoke(const Camera& c);
   void drawCube(const MACGrid::Cube& c);
   void drawFace(const MACGrid::Cube& c);
   void drawVelocities();
   vec4 getRenderColor(int i, int j, int k);
   vec4 getRenderColor(const vec3& pt);
   void drawZSheets(bool backToFront);
   void drawXSheets(bool backToFront);

   // GridData accessors
   enum Direction { X, Y, Z };
   vec3 getVelocity(const vec3& pt);
   double getVelocityX(const vec3& pt);
   double getVelocityY(const vec3& pt);
   double getVelocityZ(const vec3& pt);
   double getTemperature(const vec3& pt);
   double getDensity(const vec3& pt);
   vec3 getCenter(int i, int j, int k);

   GridDataX mU; // X component of velocity, stored on X faces, size is (dimX+1)*dimY*dimZ
   GridDataY mV; // Y component of velocity, stored on Y faces, size is dimX*(dimY+1)*dimZ
   GridDataZ mW; // W component of velocity, stored on Z faces, size is dimX*dimY*(dimZ+1)
   GridData mP;  // Pressure, stored at grid centers, size is dimX*dimY*dimZ
   GridData mD;  // Density, stored at grid centers, size is dimX*dimY*dimZ
   GridData mT;  // Temperature, stored at grid centers, size is dimX*dimY*dimZ

   GridData mSolid; // TO specify which cell is a solid cell and store it in the MAC grid

public:

   enum RenderMode { CUBES, SHEETS };
   static RenderMode theRenderMode;
   static bool theDisplayVel;

   bool checkBoundary(int i,int j,int k, int direction);

   bool checkSolid(int i, int j, int k);

   int getNeighbor(int i, int j, int k);
};

#endif