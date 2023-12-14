#include "myFace.h"
#include "myvector3d.h"
#include "myHalfedge.h"
#include "myVertex.h"
#include <GL/glew.h>

myFace::myFace(void)
{
	adjacent_halfedge = NULL;
	normal = new myVector3D(1.0, 1.0, 1.0);
}

myFace::~myFace(void)
{
	if (normal) delete normal;
}

void myFace::computeNormal() // do with M. Pluta in TP
{
	myHalfedge* edge = this->adjacent_halfedge;
	myHalfedge* next = edge->next;
	myHalfedge* prev = edge->prev;

	myVertex* v1 = edge->source;
	myVertex* vNext = next->source;
	myVertex* vPrev = prev->source;

	myPoint3D* p1 = v1->point;
	myPoint3D* pNext = vNext->point;
	myPoint3D* pPrev = vPrev->point;

	myVector3D vect1ToNext = *pNext - *p1;
	myVector3D vectPrevTo1 = *p1 - *pPrev;

	myVector3D v3 = vectPrevTo1.crossproduct(vect1ToNext);
	v3.normalize();

	normal->dX = v3.dX;
	normal->dY = v3.dY;
	normal->dZ = v3.dZ;
}