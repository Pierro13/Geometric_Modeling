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

void myMesh::checkVertices() {
	for (auto vertex : vertices) {
		if (vertex->originof == NULL) {
			cout << "Error! Vertex does not have an originating halfedge!\n";
			return;
		}
	}
	cout << "All vertices have an originating halfedge!\n";
}

void myMesh::checkFaces() {
	for (auto face : faces) {
		if (face->adjacent_halfedge == NULL) {
			cout << "Error! Face does not have an adjacent halfedge!\n";
			return;
		}
	}
	cout << "All faces have an adjacent halfedge!\n";
}

void myMesh::checkEdgeConnectivity() {
	for (auto edge : halfedges) {
		if (edge->next == NULL || edge->prev == NULL) {
			cout << "Error in edge connectivity!\n";
			return;
		}
	}
	cout << "Edge connectivity is correct!\n";
}

void myMesh::checkNormalsConsistency() {

	for (auto face : faces) {
		myHalfedge* startEdge = face->adjacent_halfedge;
		myHalfedge* currentEdge = startEdge;
		myVector3D normal = myVector3D(0.0, 0.0, 0.0);
		int cpt = 0;

		do {
			cpt++;
			normal += *(currentEdge->source->normal);
			currentEdge = currentEdge->next;
		} while (currentEdge != startEdge);

		normal = normal / (float)cpt;

		if (normal != *(face->normal)) {
			cout << "Error in face normal!\n";
			return;
		}
		else {
			cout << "Face normal is correct!\n";
		}
	}
}

double myMesh::longueur(myHalfedge* he)
{
	myPoint3D* p1 = he->source->point;
	myPoint3D* p2 = he->next->source->point;
	return p1->dist(*p2);
}

void myMesh::simplification()
{
	if (faces.size() <= 4) { return; } // regarde si on a au moins 4 faces

	// check de l halfedge le plus petit
	myHalfedge* h1 = halfedges[0];
	for (myHalfedge* heCheck : halfedges) { if (longueur(heCheck) < longueur(h1)) { h1 = heCheck; } }

	// calcul du point du milieu de lhalfedges
	myVertex* v2 = h1->next->source;
	myPoint3D* ve = new myPoint3D();
	*ve = *(h1->source->point) + *(v2->point);
	*ve /= 2;

	h1->source->point = ve;

	// Changer le point d'accroche de chaque halfedges autour de v2->originOf
	// en utiisant twin->next
	myHalfedge* vSS = h1;
	myHalfedge* vPP = v2->originof;
	do
	{
		vPP->source = vSS->source;
		vPP = vPP->twin->next;
	} while (vPP != v2->originof);

	// changer les next et prev pour mes next et prev et les next et prev de mon twin
	h1->source->originof = h1->next;

	h1->next->prev = h1->prev;
	h1->prev->next = h1->next;

	h1->twin->next->prev = h1->twin->prev;
	h1->twin->prev->next = h1->twin->next;

	h1->adjacent_face->adjacent_halfedge = h1->next;
	h1->twin->adjacent_face->adjacent_halfedge = h1->twin->next;

	// supprimer le halfedge et son twin et son 2eme point
	halfedges.erase(std::remove(halfedges.begin(), halfedges.end(), h1->twin), halfedges.end());
	halfedges.erase(std::remove(halfedges.begin(), halfedges.end(), h1), halfedges.end());
	vertices.erase(std::remove(vertices.begin(), vertices.end(), v2), vertices.end());

	// Gestion du cas ou on a une figure inconcevable cest a dire une face de 2 halfedges
	// initialiser un vecteur pour stocker les halfedges qu'on doit supprimer
	vector<myHalfedge*> ttt;


	for (myHalfedge* h : halfedges)
	{
		// Si en faisant next -> next je retombe sur moi meme alors 2 halfedges
		if (h->next->next == h)
		{
			for (myVertex* v : vertices)
			{
				// pour les vertices qui ont comme originOf moi meme ou mon next,
				// vous prenez quelqu'un d'autre
				if (v->originof == h)
				{
					v->originof = h->twin->next;
				}
				else if (v->originof == h->next)
				{
					v->originof = h->next->twin->next;
				}
			}

			// Etape de transformation de deux halfedges en 1 (assimilation)
			h->twin->twin = h->next->twin;
			h->next->twin->twin = h->twin;
			ttt.push_back(h);
			ttt.push_back(h->next);

			// suppression de la face entre les deux halfedges
			faces.erase(std::remove(faces.begin(), faces.end(), h->adjacent_face), faces.end());
		}
	}

	// suppression des halfedges qui sont dans la liste des halfedges a supprimer
	halfedges.erase(std::remove_if(halfedges.begin(), halfedges.end(), [ttt](myHalfedge* hhh) {if (count(ttt.begin(), ttt.end(), hhh))
	{
		return true;
	}
	else
	{
		return false;
	}}), halfedges.end());

	checkMesh();
	checkVertices();
	checkFaces();
	checkEdgeConnectivity();

}


bool myMesh::readFile(string filename) {
	string s, t, u;
	vector<int> faceids;
	myHalfedge** hedges;

	ifstream fin(filename);
	if (!fin.is_open()) {
		cout << "Unable to open file!\n";
		return false;
	}
	name = filename;

	string obj_name;

	size_t lastSlash = filename.find_last_of("\\/");
	if (lastSlash != string::npos)
		obj_name = filename.substr(lastSlash + 1);
	else obj_name = filename;

	cout << "Reading file " << obj_name << "\n\n";

	map<pair<int, int>, myHalfedge*> twin_map;
	map<pair<int, int>, myHalfedge*>::iterator it;

	while (getline(fin, s)) {
		stringstream myline(s);
		myline >> t;
		if (t == "g") {}
		else if (t == "v") {
			float coords[3];
			for (int i = 0; i < 3; ++i) myline >> coords[i];
			myVertex* v = new myVertex();
			v->point = new myPoint3D(coords[0], coords[1], coords[2]);
			vertices.push_back(v);
		}
		else if (t == "mtllib") {}
		else if (t == "usemtl") {}
		else if (t == "s") {}
		else if (t == "f") {
			while (myline >> u) {
				int vertexIndex = atoi((u.substr(0, u.find("/"))).c_str());
				faceids.push_back(vertexIndex);
			}

			int taille = faceids.size();
			if (taille >= 3) { 
				hedges = new myHalfedge * [taille];
				int a = 0;
				while (a < taille) {
					hedges[a] = new myHalfedge();
					a++;
				}

				myFace* face = new myFace();
				face->adjacent_halfedge = hedges[0];

				for (unsigned int i = 0; i < taille; i++) {
					int i_plus = (i + 1) % taille;
					int i_moins = (i - 1 + taille) % taille;

					myHalfedge* next = hedges[i_plus];
					myHalfedge* prev = hedges[i_moins];

					hedges[i]->next = next;
					hedges[i]->prev = prev;
					hedges[i]->adjacent_face = face;
					hedges[i]->source = vertices[faceids[i] - 1];
					hedges[i]->source->originof = hedges[i];

					it = twin_map.find(make_pair(faceids[i_plus], faceids[i]));
					if (it != twin_map.end()) {
						hedges[i]->twin = it->second;
						it->second->twin = hedges[i];
						twin_map.erase(it);
					}
					else {
						twin_map[make_pair(faceids[i], faceids[i_plus])] = hedges[i];
					}

					halfedges.push_back(hedges[i]);
				}
				faces.push_back(face);
				faceids.clear();
			}
		}
	}

	checkMesh();
	checkVertices();
	checkFaces();
	checkEdgeConnectivity();
	normalize();

	cout << "\nReading file done!\n";
	cout << "------------------\n";

	return true;
}


void myMesh::computeNormals() // do with M. Pluta in TP
{
	for (myFace* f : faces) {
		f->computeNormal();
	}

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
	myVertex* centre_de_gravite = new myVertex();
	centre_de_gravite->point = new myPoint3D();

	myHalfedge* start_Edge = f->adjacent_halfedge;
	myHalfedge* current_Edge = start_Edge;
	int cpt = 0;
	do {
		*(centre_de_gravite->point) += *(current_Edge->source->point);
		cpt++;

		current_Edge = current_Edge->next;

	} while (current_Edge != start_Edge);

	*(centre_de_gravite->point) /= cpt;

	return centre_de_gravite;
}

void myMesh::triangulate()
{
	// save all the faces in a vector
	vector<myFace*> savedFaces;
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
	checkVertices();
	checkFaces();
	checkEdgeConnectivity();
}

//return false if already triangle, true othewise.
bool myMesh::triangulate(myFace* f)
{
	if (f->adjacent_halfedge->next->next->next != f->adjacent_halfedge) {
		myVertex* centre_de_gravite = gravityCenter(f);
		vertices.push_back(centre_de_gravite);

		myHalfedge* start_Edge = f->adjacent_halfedge;
		myHalfedge* current_Edge = start_Edge;
		myHalfedge* old_edge = current_Edge;
		myHalfedge* twin;

		bool running = true;
		do {
			old_edge = old_edge->next;

			myFace* triangleFace = new myFace();
			triangleFace->adjacent_halfedge = current_Edge;

			myHalfedge* center_actual = new myHalfedge();
			myHalfedge* next_center = new myHalfedge();

			center_actual->source = centre_de_gravite;
			center_actual->next = current_Edge;
			center_actual->prev = next_center;
			center_actual->adjacent_face = triangleFace;
			center_actual->source->originof = center_actual;

			next_center->source = current_Edge->next->source;
			next_center->next = center_actual;
			next_center->prev = current_Edge;

			next_center->adjacent_face = triangleFace;
			next_center->source->originof = next_center;

			if (!(current_Edge == start_Edge)) {
				center_actual->twin = twin;
				twin->twin = center_actual;
			}

			twin = next_center;

			current_Edge->next = next_center;
			current_Edge->prev = center_actual;
			current_Edge->adjacent_face = triangleFace;

			halfedges.push_back(center_actual);
			halfedges.push_back(next_center);
			faces.push_back(triangleFace);

			running = false;
			current_Edge = old_edge;

		} while (!(current_Edge == start_Edge));

		start_Edge->prev->twin = twin;
		twin->twin = start_Edge->prev;

		faces.erase(remove(faces.begin(), faces.end(), f), faces.end());
		delete f;

		return true;
	}

	return false;
}


void myMesh::simplify() {
	
	for (myVertex* v : vertices) {
		myVertex* v2 = v->originof->next->source;

		// millieux des deux points
		myPoint3D* ve = new myPoint3D();
		*ve = *(v->point) + *(v2->point);
		*ve /= 2;

		// on bouge les deux points vers le millieu
		v->point = ve;

		// change les orgiineof de v2 pour les mettre sur ve
		for (myHalfedge* h : halfedges) {
			if (h->source == v2) {
				h->source = v;
			}
		}



		break;
	}

}