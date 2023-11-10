#include "myVertex.h"
#include "myvector3d.h"
#include "myHalfedge.h"
#include "myFace.h"

myVertex::myVertex(void)
{
	point = NULL;
	originof = NULL;
	normal = new myVector3D(0.0, 0.0, 0.0);
}

myVertex::~myVertex(void)
{
	if (normal) delete normal;
}

void myVertex::computeNormal()
{
	myHalfedge* startEdge = originof;
	myHalfedge* currentEdge = startEdge;
	float cpt = 0;

	do {
		cpt++;
		*normal += *(currentEdge->adjacent_face->normal);
		currentEdge = currentEdge->twin->next;
	} while (currentEdge != startEdge);

	*normal = *normal / (float)cpt;
}