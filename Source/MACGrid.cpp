#include "MACGrid.h"
#include "GL/glut.h"
#include "camera.h"
#include "ConjGrad.h"
#include <math.h>
#include <map>
#include <stdio.h>

#include <stdlib.h>
#include <time.h>

// Globals
MACGrid target;
extern int theDim[3];
extern double theCellSize;

// NOTE: x -> cols, z -> rows, y -> stacks
MACGrid::RenderMode MACGrid::theRenderMode = SHEETS;
bool MACGrid::theDisplayVel = false;

#define FOR_EACH_CELL \
   for(int k = 0; k < theDim[MACGrid::Z]; k++)  \
      for(int j = 0; j < theDim[MACGrid::Y]; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_FACE \
   for(int k = 0; k < theDim[MACGrid::Z]+1; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]+1; j++) \
         for(int i = 0; i < theDim[MACGrid::X]+1; i++) 

MACGrid::MACGrid()
{
   initialize();
}

MACGrid::MACGrid(const MACGrid& orig)
{
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;
}

MACGrid& MACGrid::operator=(const MACGrid& orig)
{
   if (&orig == this)
   {
      return *this;
   }
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;   

   return *this;
}

MACGrid::~MACGrid()
{
}

void MACGrid::reset()
{
   srand( time(NULL) );
   mU.initialize();
   mV.initialize();
   mW.initialize();
   mP.initialize();
   mD.initialize();
   mT.initialize(300.0);
   mSolid.initialize(0.0);
}

void MACGrid::initialize()
{
   reset();
}

void MACGrid::updateSources()
{
    // TODO: Set initial values for density, temperature, velocity
	
	vec3 whereto, wheretoSolid;

	double rand_num = (double)rand()/(double)RAND_MAX;

	double angle = rand_num*3.14;
	double cos_ang = cos(angle);
	double sin_ang = sin(angle);

	if(cos_ang < 0.2)
		cos_ang = 0.2;
	if(sin_ang < 0.2)
		sin_ang = 0.2;

	mV(1,2,1) = 1;
	mD(1,2,1) = 1;
	mT(1,2,1) = 350;

	FOR_EACH_CELL
	{
		if(i == 0 || i == theDim[0]-1 || j == 0 || j == theDim[1]-1 || k == 0 || k == theDim[2]-1)
			mSolid(i,j,k) = 1;
	}

	////solids
	////P
	//mSolid(39, 31, 0) = 1;
	//mSolid(39, 32, 0) = 1;
	//mSolid(39, 33, 0) = 1;
	//mSolid(39, 34, 0) = 1;
	//mSolid(39, 35, 0) = 1;
	//mSolid(38, 35, 0) = 1;
	//mSolid(37, 35, 0) = 1;
	//mSolid(36, 35, 0) = 1;
	//mSolid(36, 34, 0) = 1;
	//mSolid(36, 33, 0) = 1;
	//mSolid(37, 33, 0) = 1;
	//mSolid(38, 33, 0) = 1;

	////E
	//mSolid(30, 35, 0) = 1;
	//mSolid(31, 35, 0) = 1;
	//mSolid(32, 35, 0) = 1;
	//mSolid(33, 35, 0) = 1;
	//mSolid(33, 34, 0) = 1;
	//mSolid(33, 33, 0) = 1;
	//mSolid(32, 33, 0) = 1;
	//mSolid(33, 32, 0) = 1;
	//mSolid(33, 31, 0) = 1;
	//mSolid(32, 31, 0) = 1;
	//mSolid(31, 31, 0) = 1;
	//mSolid(30, 31, 0) = 1;

	////N
	//mSolid(26, 31, 0) = 1;
	//mSolid(26, 32, 0) = 1;
	//mSolid(26, 33, 0) = 1;
	//mSolid(26, 34, 0) = 1;
	//mSolid(26, 35, 0) = 1;
	//mSolid(25, 35, 0) = 1;
	//mSolid(24, 34, 0) = 1;
	//mSolid(23, 33, 0) = 1;
	//mSolid(22, 32, 0) = 1;
	//mSolid(21, 31, 0) = 1;
	//mSolid(20, 31, 0) = 1;
	//mSolid(20, 32, 0) = 1;
	//mSolid(20, 33, 0) = 1;
	//mSolid(20, 34, 0) = 1;
	//mSolid(20, 35, 0) = 1;

	////N
	//mSolid(16, 31, 0) = 1;
	//mSolid(16, 32, 0) = 1;
	//mSolid(16, 33, 0) = 1;
	//mSolid(16, 34, 0) = 1;
	//mSolid(16, 35, 0) = 1;
	//mSolid(15, 35, 0) = 1;
	//mSolid(14, 34, 0) = 1;
	//mSolid(13, 33, 0) = 1;
	//mSolid(12, 32, 0) = 1;
	//mSolid(11, 31, 0) = 1;
	//mSolid(10, 31, 0) = 1;
	//mSolid(10, 32, 0) = 1;
	//mSolid(10, 33, 0) = 1;
	//mSolid(10, 34, 0) = 1;
	//mSolid(10, 35, 0) = 1;

	////C
	//mSolid(34, 17, 0) = 1;
	//mSolid(35, 17, 0) = 1;
	//mSolid(36, 17, 0) = 1;
	//mSolid(37, 16, 0) = 1;
	//mSolid(37, 15, 0) = 1;
	//mSolid(37, 14, 0) = 1;
	//mSolid(36, 13, 0) = 1;
	//mSolid(35, 13, 0) = 1;
	//mSolid(34, 13, 0) = 1;

	////G
	//mSolid(27, 17, 0) = 1;
	//mSolid(28, 17, 0) = 1;
	//mSolid(29, 17, 0) = 1;
	//mSolid(30, 16, 0) = 1;
	//mSolid(30, 15, 0) = 1;
	//mSolid(30, 14, 0) = 1;
	//mSolid(29, 13, 0) = 1;
	//mSolid(28, 13, 0) = 1;
	//mSolid(27, 14, 0) = 1;
	//mSolid(27, 15, 0) = 1;
	//mSolid(28, 15, 0) = 1;

	////G
	//mSolid(20, 17, 0) = 1;
	//mSolid(21, 17, 0) = 1;
	//mSolid(22, 17, 0) = 1;
	//mSolid(23, 16, 0) = 1;
	//mSolid(23, 15, 0) = 1;
	//mSolid(23, 14, 0) = 1;
	//mSolid(22, 13, 0) = 1;
	//mSolid(21, 13, 0) = 1;
	//mSolid(20, 14, 0) = 1;
	//mSolid(20, 15, 0) = 1;
	//mSolid(21, 15, 0) = 1;

	////T
	//mSolid(12, 17, 0) = 1;
	//mSolid(13, 17, 0) = 1;
	//mSolid(14, 17, 0) = 1;
	//mSolid(15, 17, 0) = 1;
	//mSolid(16, 17, 0) = 1;
	//mSolid(14, 16, 0) = 1;
	//mSolid(14, 15, 0) = 1;
	//mSolid(14, 14, 0) = 1;
	//mSolid(14, 13, 0) = 1;

	//

	//
	//mD(25,1,0) = 1;
	//mT(25,1,0) = 350;
	//mV(25,1,0) = sin_ang*10;

	//mD(25,50,0) = 1;
	//mT(25,50,0) = 200;
	//mV(25,50,0) = 0.0;

	//mD(7,1,0) = 1;
	//mT(7,1,0) = 350;
	//mV(7,1,0) = sin_ang*10;

	//mD(7,50,0) = 1;
	//mT(7,50,0) = 200;
	//mV(7,50,0) = 0.0;

	//mD(42,1,0) = 1;
	//mT(42,1,0) = 350;
	//mV(42,1,0) = sin_ang*10;

	//mD(42,50,0) = 1;
	//mT(42,50,0) = 200;
	//mV(42,50,0) = 0.0;


	/*mD(whereto[0],whereto[1],whereto[2]) = 1.0;
	mT(whereto[0],whereto[1],whereto[2]) = 500.0;
	
	
	mT(whereto[0]-1,whereto[1],whereto[2]) = 400.0;
	mT(whereto[0]+1,whereto[1],whereto[2]) = 400.0;
	mT(whereto[0],whereto[1]+1,whereto[2]) = 400.0;
	mT(whereto[0]-1,whereto[1]+1,whereto[2]) = 400.0;
	mT(whereto[0]+1,whereto[1]+1,whereto[2]) = 400.0;

	mU(whereto[0],whereto[1],whereto[2]) = cos_ang*5;
	mV(whereto[0],whereto[1],whereto[2]) = sin_ang*5;
	mW(whereto[0],whereto[1],whereto[2]) = 0.0;

	mSolid(theDim[0]/2-8,theDim[1]/2,0) = 1;
	mSolid(theDim[0]/2-7,theDim[1]/2,0) = 1;
	mSolid(theDim[0]/2-6,theDim[1]/2,0) = 1;
	
	mSolid(theDim[0]/2-2, 15, 0) = 1;
	mSolid(theDim[0]/2-1, 15, 0) = 1;
	mSolid(theDim[0]/2, 15, 0) = 1;
	mSolid(theDim[0]/2+1, 15, 0) = 1;
	mSolid(theDim[0]/2+2, 15, 0) = 1;

	mSolid(theDim[0]/2+6,theDim[1]/2,0) = 1;
	mSolid(theDim[0]/2+7,theDim[1]/2,0) = 1;
	mSolid(theDim[0]/2+8,theDim[1]/2,0) = 1;*/

}

void MACGrid::advectVelocity(double dt)
{
    // TODO: Calculate new velocities and store in target
    // Then save the result to our object

	vec3 Xg, face_u_center, face_v_center, face_w_center, U_vec_at_face_u, U_vec_at_face_v, U_vec_at_face_w, Xp_face_u, Xp_face_v, Xp_face_w;
	
	FOR_EACH_FACE
	{
		Xg = getCenter(i,j,k); //get center of the cell
		if(j < theDim[1] && k < theDim[2])
		{
			if(!checkBoundary(i,j,k,X))
			{
				face_u_center.set( (Xg[0]-theCellSize/2), Xg[1], Xg[2] );
				U_vec_at_face_u = getVelocity(face_u_center);
				Xp_face_u = face_u_center - dt*U_vec_at_face_u;
				target.mU(i,j,k) = getVelocityX(Xp_face_u);
			}
			else
				target.mU(i,j,k) = 0.0;
		}
		
		if(i < theDim[0] && k < theDim[2])
		{
			if(!checkBoundary(i,j,k,Y))
			{
				face_v_center.set( Xg[0], (Xg[1]-theCellSize/2), Xg[2] );
				U_vec_at_face_v = getVelocity(face_v_center);
				Xp_face_v = face_v_center - dt*U_vec_at_face_v;
				target.mV(i,j,k) = getVelocityY(Xp_face_v);
			}
			else
				target.mV(i,j,k) = 0.0;
		}

		if(i < theDim[0] && j < theDim[1])
		{
			if(!checkBoundary(i,j,k,Z))
			{
				face_w_center.set( Xg[0], Xg[1], (Xg[2]-theCellSize/2) );
				U_vec_at_face_w = getVelocity(face_w_center);
				Xp_face_w = face_w_center - dt*U_vec_at_face_w;
				target.mW(i,j,k) = getVelocityZ(Xp_face_w);
			}
			else
				target.mW(i,j,k) = 0.0;
		}
	}

	FOR_EACH_CELL
	{
		if(checkSolid(i,j,k))
		{
			target.mU(i,j,k) = 0.0;
			target.mU(i+1,j,k) = 0.0;
		 
			target.mV(i,j,k) = 0.0;
			target.mV(i,j+1,k) = 0.0;

			target.mW(i,j,k) = 0.0;
			target.mW(i,j,k+1) = 0.0;
		}
	}
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
}

void MACGrid::advectTemperature(double dt)
{
    // TODO: Calculate new temp and store in target
    // Then save the result to our object
	FOR_EACH_CELL
	{
		vec3 Xg = getCenter(i,j,k);
		vec3 U_vec_at_Xg = getVelocity(Xg);
		vec3 Xp = Xg - dt*U_vec_at_Xg;
		target.mT(i,j,k) = getTemperature(Xp);
		//cout<<"temperature at "<<i<<','<<j<<','<<k<<':'<<target.mT(i,j,k)<<endl;
	}

	FOR_EACH_CELL
	{
		if(checkSolid(i,j,k))
			target.mT(i,j,k) = 300.0;
	}
	mT = target.mT;
	//cout<<endl;
}

void MACGrid::advectDensity(double dt)
{
    // TODO: Calculate new densitities and store in target
    // Then save the result to our object
	FOR_EACH_CELL
	{
		vec3 Xg = getCenter(i,j,k);
		vec3 U_vec_at_Xg = getVelocity(Xg);
		vec3 Xp = Xg - dt*U_vec_at_Xg;
		target.mD(i,j,k) = getDensity(Xp);
	}

	FOR_EACH_CELL
	{
		if(checkSolid(i,j,k))
			target.mD(i,j,k) = 0.0;
	}
	mD = target.mD;
}

void MACGrid::computeBouyancy(double dt)
{
   // TODO: Calculate bouyancy and store in target
   // and then save the result to our object

	double alpha, beta;
	double temperature, density;
	double tempterature_amb = 300;
	alpha = 1.0;
	beta = 1.0;
	double forcey;

	FOR_EACH_FACE
	{
			if(i < theDim[0] && k < theDim[2])
			{
				if(!checkBoundary(i,j,k,Y))
				{
					vec3 Xg = getCenter(i,j,k);
					Xg[1] = Xg[1] - theCellSize/2;
					temperature = getTemperature(Xg);

					density = getDensity(Xg);
					forcey = -alpha*density + beta*(temperature - tempterature_amb);
					target.mV(i,j,k) = mV(i,j,k) + dt*forcey;
				}
				else
					target.mV(i,j,k) = 0.0;
			}
	}

	FOR_EACH_CELL
	{
		if(checkSolid(i,j,k))
		{
			target.mV(i,j,k) = 0.0;
			target.mV(i,j+1,k) = 0.0;
		}
	}
	mV = target.mV;
}

void MACGrid::computeVorticityConfinement(double dt)
{
   // TODO: Calculate vorticity confinement forces
   // Apply the forces to the current velocity and store the result in target
   // Then save the result to our object

	vec3 vorticity, N_vector;

	vec3 force;

	double etta = 7.0;

	GridData vort_X = target.mU;
	GridData vort_Y = target.mV;
	GridData vort_Z = target.mW;

	GridData grad_vort_X = target.mU;
	GridData grad_vort_Y = target.mV;
	GridData grad_vort_Z = target.mW;

	GridData vort_mag = target.mD;	

	double temp1, temp2, temp3, temp4, temp5, temp6;

	double cell_twice_inv = 1.0/(2.0*theCellSize);

	FOR_EACH_FACE
	{
		if(j < 2)
		{
			temp1 = (mW(i,j+1,k) - 0.0)*cell_twice_inv;
			temp6 = (mU(i,j+1,k) - 0.0)*cell_twice_inv;
		}
		else if((j+2) == theDim[1])
		{
			temp1 = (0.0 - mW(i,j-1,k))*cell_twice_inv;
			temp6 = (0.0 - mU(i,j-1,k))*cell_twice_inv;
		}
		else
		{
			temp1 = (mW(i,j+1,k) - mW(i,j-1,k))*cell_twice_inv;
			temp6 = (mU(i,j+1,k) - mU(i,j-1,k))*cell_twice_inv;
		}

		if(k < 2)
		{
			temp2 = (mV(i,j,k+1) - 0.0)*cell_twice_inv;
			temp3 = (mU(i,j,k+1) - 0.0)*cell_twice_inv;
		}
		else if((k+2) == theDim[2])
		{
			temp2 = (0.0 - mV(i,j,k-1))*cell_twice_inv;
			temp3 = (0.0 - mU(i,j,k-1))*cell_twice_inv;
		}
		else
		{
			temp2 = (mV(i,j,k+1) - mV(i,j,k-1))*cell_twice_inv;
			temp3 = (mU(i,j,k+1) - mU(i,j,k-1))*cell_twice_inv;
		}
		
		if(i < 2)
		{
			temp4 = (mW(i+1,j,k) - 0.0)*cell_twice_inv;
			temp5 = (mV(i+1,j,k) - 0.0)*cell_twice_inv;
		}
		else if( (i+2) == theDim[0])
		{
			temp4 = (0.0 - mW(i-1,j,k))*cell_twice_inv;
			temp5 = (0.0 - mV(i-1,j,k))*cell_twice_inv;
		}
		else
		{
			temp4 = (mW(i+1,j,k) - mW(i-1,j,k))*cell_twice_inv;
			temp5 = (mV(i+1,j,k) - mV(i-1,j,k))*cell_twice_inv;
		}

		vort_X(i,j,k) = (temp1 - temp2);
		vort_Y(i,j,k) = (temp3 - temp4);
		vort_Z(i,j,k) = (temp5 - temp6);
		vorticity.set( vort_X(i,j,k), vort_Y(i,j,k), vort_Z(i,j,k));
		vort_mag(i,j,k) = vorticity.Length();
	}

	FOR_EACH_FACE
	{
		if(i == 0)
			grad_vort_X(i,j,k) = (vort_mag(i+1,j,k) - 0.0)*cell_twice_inv;
		else if(i == theDim[0])
			grad_vort_X(i,j,k) = (0.0 - vort_mag(i-1,j,k))*cell_twice_inv;
		else
			grad_vort_X(i,j,k) = (vort_mag(i+1,j,k) - vort_mag(i-1,j,k))*cell_twice_inv;

		if(j == 0)
			grad_vort_Y(i,j,k) = (vort_mag(i,j+1,k) - 0.0)*cell_twice_inv;
		else if(j == theDim[1])
			grad_vort_Y(i,j,k) = (0.0 - vort_mag(i,j-1,k))*cell_twice_inv;
		else
			grad_vort_Y(i,j,k) = (vort_mag(i,j+1,k) - vort_mag(i,j-1,k))*cell_twice_inv;

		if(k == 0)
			grad_vort_Z(i,j,k) = (vort_mag(i,j,k+1) - 0.0)*cell_twice_inv;
		else if(k == theDim[2])
			grad_vort_Z(i,j,k) = (0.0 - vort_mag(i,j,k-1))*cell_twice_inv;
		else
			grad_vort_Z(i,j,k) = (vort_mag(i,j,k+1) - vort_mag(i,j,k-1))*cell_twice_inv;

		N_vector.set(grad_vort_X(i,j,k), grad_vort_Y(i,j,k), grad_vort_Z(i,j,k));
		double mag = N_vector.Length() + 1e-15;

		vorticity.set(vort_X(i,j,k), vort_Y(i,j,k), vort_Z(i,j,k));
		N_vector = N_vector/mag;

		force = etta*theCellSize*(N_vector.Cross(vorticity));

		if(j < theDim[1] && k < theDim[2])
		{
			if(!checkBoundary(i,j,k,X))
				target.mU(i,j,k) = mU(i,j,k) + dt*force[0];
			else
				target.mU(i,j,k) = 0.0;
		}

		if(i < theDim[0] && k < theDim[2])
		{
			if(!checkBoundary(i,j,k,Y))
				target.mV(i,j,k) = mV(i,j,k) + dt*force[1];
			else
				target.mV(i,j,k) = 0.0;
		}

		if(j < theDim[1] && k < theDim[2])
		{
			if(!checkBoundary(i,j,k,Z))
				target.mW(i,j,k) = mW(i,j,k) + dt*force[2];
			else
				target.mW(i,j,k) = 0.0;
		}
	}

	FOR_EACH_CELL
	{
		if(checkSolid(i,j,k))
		{
			target.mU(i,j,k) = 0.0;
			target.mU(i+1,j,k) = 0.0;

			target.mV(i,j,k) = 0.0;
			target.mV(i,j+1,k) = 0.0;

			target.mW(i,j,k) = 0.0;
			target.mW(i,j,k+1) = 0.0;
		}
	}

   mU = target.mU;
   mV = target.mV;
   mW = target.mW;
}

void MACGrid::addExternalForces(double dt)
{
   computeBouyancy(dt);
   computeVorticityConfinement(dt);
}


bool MACGrid::checkSolid(int i, int j, int k)
{
	if(mSolid(i,j,k) > 0.9)
		return true;
	else 
		return false;
}

int MACGrid::getNeighbor(int i, int j, int k)
{
	int neighbor = 6;
	
	/*if(theDim[0] == 3 || theDim[1] == 3 || theDim[2] == 3)
		neighbor = 4;
	else
		neighbor = 6;
	
	if(theDim[0] > 3)
	{
		if(i == 1)
			neighbor--;
		if(i == (theDim[0] - 2))
			neighbor--;
	}
	if(theDim[1] > 3)
	{
		if(j == 1)
			neighbor--;
		if(j == (theDim[1] - 2))
			neighbor--;
	}
	if(theDim[2] > 3)
	{
		if(k == 1)
			neighbor--;
		if(k == (theDim[2] - 2))
			neighbor--;
	}*/


	//HANDLE SOLID CELL
	if(checkSolid(i,j,k)) //return 0 as the neighbor is the present cell is a solid cell
		return 0;

	if(checkSolid(i-1,j,k)) //if left cell is a solid cell then decrease neighbor by 1
		neighbor--;

	if(checkSolid(i+1,j,k)) //if right cell is a solid cell then decrease neighbor by 1
		neighbor--;

	if(checkSolid(i,j-1,k)) //if down cell is a solid cell then decrease neighbor by 1
		neighbor--;

	if(checkSolid(i,j+1,k)) //if up cell is a solid cell then decrease neighbor by 1
		neighbor--;

	if(checkSolid(i,j,k-1)) //if back cell is a solid cell then decrease neighbor by 1
		neighbor--;

	if(checkSolid(i,j,k+1)) //if front cell is a solid cell then decrease neighbor by 1
		neighbor--;

	return neighbor;
}


void MACGrid::project(double dt)
{
   // TODO: Solve Ax = b for pressure
   // 1. Contruct b
   // 2. Construct A 
   // 3. Solve for p
   // Subtract pressure from our velocity and save in target

   // Then save the result to our object

	//using namespace boost::numeric::ublas;
	
	int size = theDim[0] * theDim[1] * theDim[2];

	boost::numeric::ublas::vector<double> vecB (size);
	boost::numeric::ublas::matrix<double> matA (size, size);
	boost::numeric::ublas::vector<double> press (size);

	double u_left, u_right, v_down, v_up, w_front, w_back;
	double cell_vel_grad;
	double mult_factor;

	double press_u, press_u_left, press_v_down, press_v, press_w, press_w_back;

	matA.clear();
	vecB.clear();
	press.clear();

	int neighbor;

	int mat_index_cell, mat_index_neigh;

	FOR_EACH_CELL
	{
		neighbor = getNeighbor(i,j,k);

		mat_index_cell = theDim[1]*theDim[0]*k + theDim[0]*j + i;
		
		matA(mat_index_cell, mat_index_cell) = neighbor; //solid cell handled in the getNeighbor function

		if(!checkSolid(i,j,k)) //if the present cell is not a solid then only put -1 else leave it
		{
			//	P(i-1,j,k)
			mat_index_neigh = theDim[1]*theDim[0]*k + theDim[0]*j + (i - 1);
			if(!checkSolid(i-1,j,k))
			{
				if(mat_index_neigh >= 0)
				{
					if( (i-1) >= 0)
						matA(mat_index_cell, mat_index_neigh) = -1.0;
				}
			}
				
			//	P(i+1,j,k)
			mat_index_neigh = theDim[1]*theDim[0]*k + theDim[0]*j + (i + 1);
			if(!checkSolid(i+1,j,k))
			{
				if(mat_index_neigh <= (size-1))
				{
					if( (i+1) < theDim[0])
						matA(mat_index_cell, mat_index_neigh) = -1.0;
				}
			}

			//	P(i,j-1,k)
			mat_index_neigh = theDim[1]*theDim[0]*k + theDim[0]*(j-1) + i;
			if(!checkSolid(i,j-1,k))
			{
				if(mat_index_neigh >= 0)
				{
					if( (j-1) >= 0)
						matA(mat_index_cell, mat_index_neigh) = -1.0;
				}
			}

			//	P(i,j+1,k)
			mat_index_neigh = theDim[1]*theDim[0]*k + theDim[0]*(j+1) + i;
			if(!checkSolid(i,j+1,k))
			{
				if(mat_index_neigh <= (size-1))
				{
					if( (j+1) < theDim[1])
						matA(mat_index_cell, mat_index_neigh) = -1.0;
				}
			}

			//	P(i,j,k-1)
			mat_index_neigh = theDim[1]*theDim[0]*(k-1) + theDim[0]*j + i;
			if(!checkSolid(i,j,k-1))
			{
				if(mat_index_neigh >= 0)
				{
					if( (k-1) >= 0)
						matA(mat_index_cell, mat_index_neigh) = -1.0;
				}
			}

			//	P(i,j,k+1)
			mat_index_neigh = theDim[1]*theDim[0]*(k+1) + theDim[0]*j + i;
			if(!checkSolid(i,j,k+1))
			{
				if(mat_index_neigh <= (size-1))
				{
					if( (k+1) < theDim[2])
						matA(mat_index_cell, mat_index_neigh) = -1.0;
				}
			}
		}

		//populate b vector
		if(i == 1)
			u_left = 0.0;
		else
			u_left = mU(i,j,k);
		if( (i+1) == theDim[0]-1)
			u_right = 0.0;
		else
			u_right = mU(i+1,j,k);

		if(j == 1)
			v_down = 0.0;
		else
			v_down = mV(i,j,k);
		if( (j+1) == theDim[1]-1)
			v_up = 0.0;
		else
			v_up = mV(i,j+1,k);
		
		if(k == 1)
			w_front = 0.0;
		else
			w_front = mW(i,j,k);
		if( (k+1) == theDim[2]-1)
			w_back = 0.0;
		else
			w_back = mW(i,j,k+1);

		mult_factor = -1.2*theCellSize/dt;

		cell_vel_grad = mult_factor * ( (u_right - u_left) + (v_up - v_down) + (w_back - w_front));

		if(!checkSolid(i,j,k))
			vecB(mat_index_cell) = cell_vel_grad;
	}

	/*cout<<"A mat values: "<<endl;
	for(int l = 0; l<size; l++)
	{
		for(int m = 0; m<size; m++)
		{
			cout<<matA(l,m)<<' ';
		}
		cout<<endl;
	}
	cout<<endl;
	
	cout<<"B values: "<<endl;
	for(int i = 0 ; i < vecB.size(); i++)
	{
		cout<<vecB(i)<<endl;
	}*/


	cg_solve(matA, vecB, press, 100, 0.1);

	//solve for Pressure
	
	//cout<<"Pressure values: "<<endl;
	FOR_EACH_CELL
	{
		mP(i,j,k) = press(k*theDim[1]*theDim[0] + j*theDim[0] + i);
		if(checkSolid(i,j,k))
			mP(i,j,k) = 0.0;
		//cout<<mP(i,j,k)<<endl;
	}
	
	FOR_EACH_FACE
	{
		press_u = mP(i,j,k);

		press_v = mP(i,j,k);

		press_w = mP(i,j,k);

		if(i > 1)
			press_u_left = mP((i-1),j,k);
		else
			press_u_left = 0.0;

		if(j > 1)
			press_v_down = mP(i,(j-1),k);
		else
			press_v_down = 0.0;
		
		if(k > 1)
			press_w_back = mP(i,j,(k-1));
		else
			press_w_back = 0.0;

		//For X, check for the valid X direction faces
		//cout<<endl;
		if(j < theDim[1]-1 && k < theDim[2]-1)
		{
			if(!checkBoundary(i,j,k,X))
				target.mU(i,j,k) = mU(i,j,k) - dt*(press_u - press_u_left)/1.2;
			else
				target.mU(i,j,k) = 0.0;
		}
		//For Y, check for the valid Y direction faces
		if(i < theDim[0]-1 && k < theDim[2]-1)
		{
			if(!checkBoundary(i,j,k,Y))
				target.mV(i,j,k) = mV(i,j,k) - dt*(press_v - press_v_down)/1.2;
			else
				target.mV(i,j,k) = 0.0;
		}
		//For Z, check for the valid Z direction faces
		if(i < theDim[0]-1 && j < theDim[1]-1)
		{
			if(!checkBoundary(i,j,k,Z))
				target.mW(i,j,k) = mW(i,j,k) - dt*(press_w - press_w_back)/1.2;
			else
				target.mW(i,j,k) = 0.0;
		}
	}
	
	FOR_EACH_CELL
	{
		if(checkSolid(i,j,k))
		{
			target.mU(i,j,k) = 0.0;
			target.mU(i+1,j,k) = 0.0;

			target.mV(i,j,k) = 0.0;
			target.mV(i,j+1,k) = 0.0;

			target.mW(i,j,k) = 0.0;
			target.mW(i,j,k+1) = 0.0;
		}
	}

	mU = target.mU;
    mV = target.mV;
    mW = target.mW;

    /*
	cout<<endl<<endl<<"............PROJECT................"<<endl<<endl;FOR_EACH_CELL
    {
	    cout<<"U vel after for "<<i<<','<<j<<','<<k<<": "<<mU(i,j,k)<<endl;
	    cout<<"V vel after for "<<i<<','<<j<<','<<k<<": "<<mV(i,j,k)<<endl;
	    cout<<"W vel after for "<<i<<','<<j<<','<<k<<": "<<mW(i,j,k)<<endl;
    }*/
    // IMPLEMENT THIS IS A SANITY CHECK: assert (checkDivergence());
}

vec3 MACGrid::getVelocity(const vec3& pt)
{
   vec3 vel;
   vel[0] = getVelocityX(pt); 
   vel[1] = getVelocityY(pt); 
   vel[2] = getVelocityZ(pt); 
   return vel;
}

double MACGrid::getVelocityX(const vec3& pt)
{
   return mU.interpolate(pt);
}

double MACGrid::getVelocityY(const vec3& pt)
{
   return mV.interpolate(pt);
}

double MACGrid::getVelocityZ(const vec3& pt)
{
   return mW.interpolate(pt);
}

double MACGrid::getTemperature(const vec3& pt)
{
   return mT.interpolate(pt);
}

double MACGrid::getDensity(const vec3& pt)
{
   return mD.interpolate(pt);
}

bool MACGrid::checkBoundary(int i, int j, int k, int direction)
{
	//add code for blocks
	
	//for wall
	switch(direction)
	{
	case X:
		if(i == 1 || i == theDim[0]-1)
			return true;
		break;
	case Y:
		if(j == 1 || j == theDim[1]-1)
			return true;
		break;
	case Z:
		if(k == 1 || k == theDim[2]-1)
			return true;
		break;
	}
	
	return false;
}


vec3 MACGrid::getCenter(int i, int j, int k)
{
   double xstart = theCellSize/2.0;
   double ystart = theCellSize/2.0;
   double zstart = theCellSize/2.0;

   double x = xstart + i*theCellSize;
   double y = ystart + j*theCellSize;
   double z = zstart + k*theCellSize;
   return vec3(x, y, z);
}

void MACGrid::drawMySolids(const Camera& c)
{
	std::multimap<double, MACGrid::Cube, std::greater<double>> cubes;
    FOR_EACH_CELL
    {
       if(checkSolid(i,j,k))
	   {
		   MACGrid::Cube cube;
		   cube.color = vec4(0.1,0.1,0.1,1.0);
		   cube.pos = getCenter(i,j,k);
		   cube.dist = DistanceSqr(cube.pos, c.getPosition());
		   cubes.insert(make_pair(cube.dist, cube));
	   }
	} 

    // Draw cubes from back to front
    std::multimap<double, MACGrid::Cube, std::greater<double>>::const_iterator it;
    for (it = cubes.begin(); it != cubes.end(); ++it)
    {
       drawCube(it->second);
    }
}

void MACGrid::draw(const Camera& c)
{   
   drawWireGrid();
   //drawMySolids(c);
   if (theDisplayVel) drawVelocities();   
   if (theRenderMode == CUBES) drawSmokeCubes(c);
   else drawSmoke(c);
}

void MACGrid::drawVelocities()
{
   // draw line at each center
   glColor4f(0.0, 1.0, 0.0, 1.0);
   glBegin(GL_LINES);
      FOR_EACH_CELL
      {
         vec3 pos = getCenter(i,j,k);
         vec3 vel = getVelocity(pos);
         if (vel.Length() > 0.0001)
         {
           vel.Normalize();
           vel *= theCellSize/2.0;
           vel += pos;
           glVertex3dv(pos.n);
           glVertex3dv(vel.n);
         }
      }
   glEnd();
}

vec4 MACGrid::getRenderColor(int i, int j, int k)
{
    double value = mD(i, j, k); 
    return vec4(1.0, 0.9, 1.0, value);
}

vec4 MACGrid::getRenderColor(const vec3& pt)
{
    double value = getDensity(pt); 
    return vec4(1.0, 0.0, 0.0, value);
}

void MACGrid::drawZSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double startk = back - stepsize;
   double endk = 0;
   double stepk = -theCellSize;

   if (!backToFront)
   {
      startk = 0;
      endk = back;   
      stepk = theCellSize;
   }

   for (double k = startk; backToFront? k > endk : k < endk; k += stepk)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double i = 0.0; i <= right; i += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double i = right; i >= 0.0; i -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawXSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double starti = right - stepsize;
   double endi = 0;
   double stepi = -theCellSize;

   if (!backToFront)
   {
      starti = 0;
      endi = right;   
      stepi = theCellSize;
   }

   for (double i = starti; backToFront? i > endi : i < endi; i += stepi)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double k = 0.0; k <= back; k += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double k = back; k >= 0.0; k -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}


void MACGrid::drawSmoke(const Camera& c)
{
   vec3 eyeDir = c.getBackward();
   double zresult = fabs(Dot(eyeDir, vec3(1,0,0)));
   double xresult = fabs(Dot(eyeDir, vec3(0,0,1)));
   double yresult = fabs(Dot(eyeDir, vec3(0,1,0)));

   if (zresult < xresult)
   {      
      drawZSheets(c.getPosition()[2] < 0);
   }
   else 
   {
      drawXSheets(c.getPosition()[0] < 0);
   }
}

void MACGrid::drawSmokeCubes(const Camera& c)
{
   std::multimap<double, MACGrid::Cube, std::greater<double>> cubes;
   FOR_EACH_CELL
   {
      MACGrid::Cube cube;
      cube.color = getRenderColor(i,j,k);
      cube.pos = getCenter(i,j,k);
      cube.dist = DistanceSqr(cube.pos, c.getPosition());
      cubes.insert(make_pair(cube.dist, cube));
   } 

   // Draw cubes from back to front
   std::multimap<double, MACGrid::Cube, std::greater<double>>::const_iterator it;
   for (it = cubes.begin(); it != cubes.end(); ++it)
   {
      drawCube(it->second);
   }
}

void MACGrid::drawWireGrid()
{
   // Display grid in light grey, draw top & bottom

   double xstart = 0.0;
   double ystart = 0.0;
   double zstart = 0.0;
   double xend = theDim[0]*theCellSize;
   double yend = theDim[1]*theCellSize;
   double zend = theDim[2]*theCellSize;

   glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
      glDisable(GL_LIGHTING);
      glColor3f(0.25, 0.25, 0.25);

      glBegin(GL_LINES);
      for (int i = 0; i <= theDim[0]; i++)
      {
         double x = xstart + i*theCellSize;
         glVertex3d(x, ystart, zstart);
         glVertex3d(x, ystart, zend);

         glVertex3d(x, yend, zstart);
         glVertex3d(x, yend, zend);
      }

      for (int i = 0; i <= theDim[2]; i++)
      {
         double z = zstart + i*theCellSize;
         glVertex3d(xstart, ystart, z);
         glVertex3d(xend, ystart, z);

         glVertex3d(xstart, yend, z);
         glVertex3d(xend, yend, z);
      }

      glVertex3d(xstart, ystart, zstart);
      glVertex3d(xstart, yend, zstart);

      glVertex3d(xend, ystart, zstart);
      glVertex3d(xend, yend, zstart);

      glVertex3d(xstart, ystart, zend);
      glVertex3d(xstart, yend, zend);

      glVertex3d(xend, ystart, zend);
      glVertex3d(xend, yend, zend);
      glEnd();
   glPopAttrib();

   glEnd();
}

#define LEN 0.5
void MACGrid::drawFace(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);
      glEnd();
   glPopMatrix();
}

void MACGrid::drawCube(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0, -1.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN, -LEN);         

         glNormal3d( 0.0,  0.0, -0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN,  LEN, -LEN);
         glVertex3d( LEN,  LEN, -LEN);
         glVertex3d( LEN, -LEN, -LEN);

         glNormal3d(-1.0,  0.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d(-LEN,  LEN,  LEN);
         glVertex3d(-LEN,  LEN, -LEN);

         glNormal3d( 0.0, 1.0,  0.0);
         glVertex3d(-LEN, LEN, -LEN);
         glVertex3d(-LEN, LEN,  LEN);
         glVertex3d( LEN, LEN,  LEN);
         glVertex3d( LEN, LEN, -LEN);

         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);

         glNormal3d(1.0,  0.0,  0.0);
         glVertex3d(LEN, -LEN, -LEN);
         glVertex3d(LEN, -LEN,  LEN);
         glVertex3d(LEN,  LEN,  LEN);
         glVertex3d(LEN,  LEN, -LEN);
      glEnd();
   glPopMatrix();
}