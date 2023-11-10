#include "myMesh.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>
#include <GL/glew.h>
#include "myvector3d.h"

using namespace std;

myMesh::myMesh(void)
{
	/**** TODO ****/
}


myMesh::~myMesh(void)
{
	/**** TODO ****/
}

void myMesh::clear()
{
	for (unsigned int i = 0; i < vertices.size(); i++) if (vertices[i]) delete vertices[i];
	for (unsigned int i = 0; i < halfedges.size(); i++) if (halfedges[i]) delete halfedges[i];
	for (unsigned int i = 0; i < faces.size(); i++) if (faces[i]) delete faces[i];

	vector<myVertex*> empty_vertices;    vertices.swap(empty_vertices);
	vector<myHalfedge*> empty_halfedges; halfedges.swap(empty_halfedges);
	vector<myFace*> empty_faces;         faces.swap(empty_faces);
}

void myMesh::checkMesh()
{
	vector<myHalfedge*>::iterator it;
	for (it = halfedges.begin(); it != halfedges.end(); it++)
	{
		if ((*it)->twin == NULL)
			break;
	}
	if (it != halfedges.end())
		cout << "Error! Not all edges have their twins!\n";
	else cout << "Each edge has a twin!\n";
}


bool myMesh::readFile(std::string filename)
{
	string s, t, u;
	vector<int> faceids;
	myHalfedge** hedges;

	ifstream fin(filename);
	if (!fin.is_open()) {
		cout << "Unable to open file!\n";
		return false;
	}
	name = filename;

	map<pair<int, int>, myHalfedge*> twin_map;
	map<pair<int, int>, myHalfedge*>::iterator it;

	while (getline(fin, s))
	{
		stringstream myline(s);
		myline >> t;
		if (t == "g") {}
		else if (t == "v") //todo
		{
			float x, y, z;
			myline >> x >> y >> z;
			//cout << "v " << x << " " << y << " " << z << endl;

			myVertex* v = new myVertex();
			v->point = new myPoint3D(x, y, z);
			vertices.push_back(v);
		}
		else if (t == "mtllib") {}
		else if (t == "usemtl") {}
		else if (t == "s") {}
		else if (t == "f") //todo
		{
			//cout << endl;
			//cout << "f";
			while (myline >> u) {
				int vertexIndex = atoi((u.substr(0, u.find("/"))).c_str());
				//cout << " " << vertexIndex;
				faceids.push_back(vertexIndex); // list of the vertex index
			}

			if (faceids.size() < 3) continue;

			hedges = new myHalfedge * [faceids.size()];
			for (int i = 0; i < faceids.size(); i++) hedges[i] = new myHalfedge();

			myFace* face = new myFace();
			face->adjacent_halfedge = hedges[0];

			for (unsigned int i = 0; i < faceids.size(); i++) {
				int iplusone = (i + 1) % faceids.size();
				int iminusone = (i - 1 + faceids.size()) % faceids.size();

				myHalfedge* next = hedges[iplusone]; //new myHalfedge();
				myHalfedge* prev = hedges[iminusone];

				//hedges[i]->index = faceids[i];
				hedges[i]->source = vertices[faceids[i] - 1];
				hedges[i]->next = next;
				hedges[i]->prev = prev;
				hedges[i]->adjacent_face = face;
				hedges[i]->source->originof = hedges[i];

				it = twin_map.find(make_pair(faceids[iplusone], faceids[i]));
				if (it == twin_map.end()) {
					// This means there was no myHalfedge* present at location(a, b)
					twin_map[make_pair(faceids[i], faceids[iplusone])] = hedges[i];
				}
				else {
					// It was found.The variable it->second is of type myHalfedge* is the halfedge present at location(a, b)
					hedges[i]->twin = it->second;
					it->second->twin = hedges[i];
					twin_map.erase(it);
				}

				halfedges.push_back(hedges[i]);
			}
			faces.push_back(face);

			faceids.clear();
		}
	}

	checkMesh();
	normalize();

	return true;
}


void myMesh::computeNormals()
{
	/**** TODO ****/
	// compute the normal of each face
	for (myFace* f : faces) {
		f->computeNormal();
	}

	// compute the normal of each vertex
	for (myVertex* v : vertices) {
		v->computeNormal();
	}
}

void myMesh::normalize()
{
	if (vertices.size() < 1) return;

	int tmpxmin = 0, tmpymin = 0, tmpzmin = 0, tmpxmax = 0, tmpymax = 0, tmpzmax = 0;

	for (unsigned int i = 0; i < vertices.size(); i++) {
		if (vertices[i]->point->X < vertices[tmpxmin]->point->X) tmpxmin = i;
		if (vertices[i]->point->X > vertices[tmpxmax]->point->X) tmpxmax = i;

		if (vertices[i]->point->Y < vertices[tmpymin]->point->Y) tmpymin = i;
		if (vertices[i]->point->Y > vertices[tmpymax]->point->Y) tmpymax = i;

		if (vertices[i]->point->Z < vertices[tmpzmin]->point->Z) tmpzmin = i;
		if (vertices[i]->point->Z > vertices[tmpzmax]->point->Z) tmpzmax = i;
	}

	double xmin = vertices[tmpxmin]->point->X, xmax = vertices[tmpxmax]->point->X,
		ymin = vertices[tmpymin]->point->Y, ymax = vertices[tmpymax]->point->Y,
		zmin = vertices[tmpzmin]->point->Z, zmax = vertices[tmpzmax]->point->Z;

	double scale = (xmax - xmin) > (ymax - ymin) ? (xmax - xmin) : (ymax - ymin);
	scale = scale > (zmax - zmin) ? scale : (zmax - zmin);

	for (unsigned int i = 0; i < vertices.size(); i++) {
		vertices[i]->point->X -= (xmax + xmin) / 2;
		vertices[i]->point->Y -= (ymax + ymin) / 2;
		vertices[i]->point->Z -= (zmax + zmin) / 2;

		vertices[i]->point->X /= scale;
		vertices[i]->point->Y /= scale;
		vertices[i]->point->Z /= scale;
	}
}


void myMesh::splitFaceTRIS(myFace* f, myPoint3D* p)
{
	/**** TODO ****/
}

void myMesh::splitEdge(myHalfedge* e1, myPoint3D* p)
{

	/**** TODO ****/
}

void myMesh::splitFaceQUADS(myFace* f, myPoint3D* p)
{
	/**** TODO ****/
}


void myMesh::subdivisionCatmullClark()
{
	/**** TODO ****/
}

myVertex* myMesh::gravityCenter(myFace* f)
{
	myVertex* centroid = new myVertex();
	centroid->point = new myPoint3D();

	myHalfedge* startEdge = f->adjacent_halfedge;
	myHalfedge* currentEdge = startEdge;
	int cpt = 0;
	do {
		*(centroid->point) += *(currentEdge->source->point);
		cpt++;

		currentEdge = currentEdge->next;

	} while (currentEdge != startEdge);

	*(centroid->point) /= cpt;

	return centroid;
}

void myMesh::triangulate()
{
	// save all the faces in a vector
	std::vector<myFace*> savedFaces;
	for (myFace* f : faces) {
		savedFaces.push_back(f);
	}

	// interate over the saved vector and triangulate each face
	// the face will then be deleted in the original faces vector
	for (myFace* f : savedFaces) {
		triangulate(f);
	}

	cout << "triangulate done" << endl;
	checkMesh();
}

//return false if already triangle, true othewise.
bool myMesh::triangulate(myFace* f)
{
	// check if the face is already a triangle, if yes return false
	if (f->adjacent_halfedge->next->next->next == f->adjacent_halfedge) {
		return false;
	}

	myVertex* centroid = gravityCenter(f);
	vertices.push_back(centroid);

	myHalfedge* startEdge = f->adjacent_halfedge;
	myHalfedge* currentEdge = startEdge;
	myHalfedge* savedEdge = currentEdge;
	myHalfedge* twinEdge;
	bool running = true;
	do {
		savedEdge = savedEdge->next;
		//myHalfedge* savedNextEdge = currentEdge->next->next;

		// create a new face
		myFace* triangleFace = new myFace();
		triangleFace->adjacent_halfedge = currentEdge;

		//// create a new halfedge
		myHalfedge* centerToCurrent = new myHalfedge();
		myHalfedge* NextToCenter = new myHalfedge();

		// fill the new halfedges comming from the centroid to the current edge
		centerToCurrent->source = centroid;
		centerToCurrent->next = currentEdge;
		centerToCurrent->prev = NextToCenter;
		centerToCurrent->adjacent_face = triangleFace;
		centerToCurrent->source->originof = centerToCurrent;

		// fill the new halfedges comming from the current edge to the centroid
		NextToCenter->source = currentEdge->next->source;
		NextToCenter->next = centerToCurrent;
		NextToCenter->prev = currentEdge;
		NextToCenter->adjacent_face = triangleFace;
		NextToCenter->source->originof = NextToCenter;

		// update the twin of the current edge
		if (currentEdge != startEdge) {
			centerToCurrent->twin = twinEdge;
			twinEdge->twin = centerToCurrent;
		}
		twinEdge = NextToCenter;


		// update the current edge
		currentEdge->next = NextToCenter;
		currentEdge->prev = centerToCurrent;
		currentEdge->adjacent_face = triangleFace;

		halfedges.push_back(centerToCurrent);
		halfedges.push_back(NextToCenter);
		faces.push_back(triangleFace);

		running = false;
		currentEdge = savedEdge;
	} while (/*running*/currentEdge != startEdge);

	startEdge->prev->twin = twinEdge;
	twinEdge->twin = startEdge->prev;

	// delete the old face
	faces.erase(std::remove(faces.begin(), faces.end(), f), faces.end());
	delete f;

	return true;
}