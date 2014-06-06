#include "smokeSim.h"
#include <GL/glut.h>
#include <IL/il.h>
#include <IL/ilu.h>
#include <IL/ilut.h>

SmokeSim::SmokeSim() : mFrameNum(0), mRecordEnabled(false)
{
   ilInit();
   iluInit();
   ilEnable(IL_FILE_OVERWRITE);
   ilutRenderer(ILUT_OPENGL);

   reset();
}

SmokeSim::~SmokeSim()
{
}

void SmokeSim::reset()
{
   mGrid.reset();
}

void SmokeSim::setGridDimensions(int x, int y, int z)
{
   extern int theDim[3]; // Naughty globals...
   theDim[0] = x;
   theDim[1] = y;
   theDim[2] = z;
   reset();
}

void SmokeSim::step()
{
   double dt = 0.01;

   // Step0: Gather user forces
   mGrid.updateSources();

   // Step1: Calculate new velocities
   mGrid.advectVelocity(dt);
   mGrid.addExternalForces(dt);
   mGrid.project(dt);

   // Step2: Calculate new temperature
   mGrid.advectTemperature(dt);

   // Step3: Calculate new density 
   mGrid.advectDensity(dt);
}

void SmokeSim::setRecording(bool on)
{
   if (on && ! mRecordEnabled)  // reset counter
   {
      mFrameNum = 0;
   }
   mRecordEnabled = on;
}

bool SmokeSim::isRecording()
{
   return mRecordEnabled;
}

void SmokeSim::draw(const Camera& c)
{
   drawAxes(); 
   mGrid.draw(c);
   if (mRecordEnabled) grabScreen();
}

void SmokeSim::drawAxes()
{
  glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
      glDisable(GL_LIGHTING);

      glLineWidth(2.0); 
      glBegin(GL_LINES);
         glColor3f(1.0, 0.0, 0.0);
         glVertex3f(0.0, 0.0, 0.0);
         glVertex3f(1.0, 0.0, 0.0);

         glColor3f(0.0, 1.0, 0.0);
         glVertex3f(0.0, 0.0, 0.0);
         glVertex3f(0.0, 1.0, 0.0);

         glColor3f(0.0, 0.0, 1.0);
         glVertex3f(0.0, 0.0, 0.0);
         glVertex3f(0.0, 0.0, 1.0);
      glEnd();
  glPopAttrib();
}

void SmokeSim::grabScreen()  // Code adapted from asst#1
{
	unsigned int image;
   ilGenImages(1, &image);
	ilBindImage(image);

	ILenum error = ilGetError();
	assert(error == IL_NO_ERROR);

	ilTexImage(640, 480, 1, 3, IL_RGB, IL_UNSIGNED_BYTE, NULL);

	error = ilGetError();
	assert(error == IL_NO_ERROR);

	unsigned char* data = ilGetData();

	error = ilGetError();
	assert(error == IL_NO_ERROR);

	for (int i=479; i>=0; i--) 
	{
		glReadPixels(0,i,640,1,GL_RGB, GL_UNSIGNED_BYTE, 
			data + (640 * 3 * i));
	}

	char anim_filename[2048];
	sprintf_s(anim_filename, 2048, "output/smoke_%04d.png", mFrameNum++); 

	ilSave(IL_PNG, anim_filename);

	error = ilGetError();
	assert(error == IL_NO_ERROR);

	ilDeleteImages(1, &image);

	error = ilGetError();
	assert(error == IL_NO_ERROR);
}