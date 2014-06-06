// ==========================================================================
// Copyright (C) 2009 Aline Normoyle
// ==========================================================================
#ifndef smokeSim_H_
#define smokeSim_H_

#include "MACGrid.h"

class Camera;
class SmokeSim
{
public:
   SmokeSim();
   virtual ~SmokeSim();

   virtual void reset();
   virtual void step();
   virtual void draw(const Camera& c);
   virtual void setGridDimensions(int x, int y, int z);
   virtual void setRecording(bool on);
   virtual bool isRecording();

protected:
   virtual void drawAxes();
   virtual void grabScreen();

protected:
   MACGrid mGrid;
   bool mRecordEnabled;
   int mFrameNum;
};

#endif