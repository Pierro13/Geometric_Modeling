#pragma once
#include "myFace.h"
#include "myHalfedge.h"
#include "myVertex.h"
#include <vector>
#include <string>

class myMesh
{
public:
	std::vector<myVertex *> vertices;
	std::vector<myHalfedge *> halfedges;
	std::vector<myFace *> faces;
	std::string name;

	void checkMesh();
	void checkVertices();
	void checkFaces();
	void checkEdgeConnectivity();
	void checkNormalsConsistency();
	double longueur(myHalfedge* he);
	void simplification();
	bool readFile(std::string filename);
	void computeNormals();
	void normalize();

	void subdivisionCatmullClark();

	void splitFaceTRIS(myFace *, myPoint3D *);

	void splitEdge(myHalfedge *, myPoint3D *);
	void splitFaceQUADS(myFace *, myPoint3D *);

	myVertex* gravityCenter(myFace* f);
	void triangulate();
	bool triangulate(myFace *);

	void simplify();

	void clear();

	myMesh(void);
	~myMesh(void);
};

