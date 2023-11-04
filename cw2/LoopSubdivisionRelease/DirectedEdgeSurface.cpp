///////////////////////////////////////////////////
//
//  Hamish Carr
//  September, 2020
//
//  ------------------------
//  DirectedEdgeSurface.cpp
//  ------------------------
//  
//  Base code for rendering assignments.
//
//  Minimalist (non-optimised) code for reading and 
//  rendering an object file
//  
//  We will make some hard assumptions about input file
//  quality. We will not check for manifoldness or 
//  normal direction, &c.  And if it doesn't work on 
//  all object files, that's fine.
//
//  While I could set it up to use QImage for textures,
//  I want this code to be reusable without Qt, so I 
//  shall make a hard assumption that textures are in 
//  ASCII PPM and use my own code to read them
//  
///////////////////////////////////////////////////

// include the header file
#include "DirectedEdgeSurface.h"

// include the C++ standard libraries we want
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <math.h>
#include <unordered_map>
#include <algorithm>
#include <set>

// include the Cartesian 3- vector class
#include "Cartesian3.h"
#include "SphereVertices.h"

#define MAXIMUM_LINE_LENGTH 1024
#define PI 3.141593
using namespace std;

// constructor will initialise to safe values
DirectedEdgeSurface::DirectedEdgeSurface()
    : centreOfGravity(0.0,0.0,0.0)
    { // DirectedEdgeSurface()
    // force arrays to size 0
    vertices.resize(0);
    normals.resize(0);
	firstDirectedEdge.resize(0);
	faceVertices.resize(0);
	otherHalf.resize(0);
    } // DirectedEdgeSurface()

// read routine returns true on success, failure otherwise
bool DirectedEdgeSurface::ReadObjectStream(std::istream &geometryStream)
    { // ReadObjectStream()
    cout<< "start reading file..." << endl;
    // create a read buffer
    char readBuffer[MAXIMUM_LINE_LENGTH];
    
    // the rest of this is a loop reading lines & adding them in appropriate places
    while (true)
        { // not eof
		// token for identifying meaning of line
		std::string token;

        // character to read
        geometryStream >> token;
        
        // check for eof() in case we've run out
        if (geometryStream.eof())
            break;

        // otherwise, switch on the token we read
		if (token == "#")
			{ // comment 
			// read and discard the line
			geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
            } // comment
		else if (token == "Vertex")
			{ // vertex
			// variables for the read
			unsigned int vertexID;
			geometryStream >> vertexID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (vertexID != vertices.size())
				{ // bad vertex ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad vertex ID				
			
			// read in the new vertex position
			Cartesian3 newVertex;
			geometryStream >> newVertex;
			
			// and add it to the vertices
			vertices.push_back(newVertex);
			} // vertex
		else if (token == "Normal")
			{ // normal
			// variables for the read
			unsigned int normalID;
			geometryStream >> normalID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (normalID != normals.size())
				{ // bad ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad ID				
			
			// read in the new normal
			Cartesian3 newNormal;
			geometryStream >> newNormal;
			
			// and add it to the vertices
			normals.push_back(newNormal);
			} // normal
		else if (token == "FirstDirectedEdge")
			{ // first directed edge
			// variables for the read
			unsigned int FDEID;
			geometryStream >> FDEID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (FDEID != firstDirectedEdge.size())
				{ // bad ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad ID				
			
			// read in the new FDE
			unsigned int newFDE;
			geometryStream >> newFDE;
			
			// and add it to the vertices
			firstDirectedEdge.push_back(newFDE);
			} // first directed edge
		else if (token == "Face")
			{ // face
			// variables for the read
			unsigned int faceID;
			geometryStream >> faceID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (faceID != faceVertices.size()/3)
				{ // bad face ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad face ID				
			
			// read in the new face vertex (3 times)
			unsigned int newFaceVertex;
			geometryStream >> newFaceVertex;
			faceVertices.push_back(newFaceVertex);
			geometryStream >> newFaceVertex;
			faceVertices.push_back(newFaceVertex);
			geometryStream >> newFaceVertex;
			faceVertices.push_back(newFaceVertex);
			} // face
		else if (token == "OtherHalf")
			{ // other half
			// variables for the read
			unsigned int otherHalfID;
			geometryStream >> otherHalfID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (otherHalfID != otherHalf.size())
				{ // bad ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad ID				
			
			// read in the new face vertex (3 times)
			unsigned int newOtherHalf;
			geometryStream >> newOtherHalf;
			otherHalf.push_back(newOtherHalf);
			} // other half
        } // not eof

    // compute centre of gravity
    // note that very large files may have numerical problems with this
    centreOfGravity = Cartesian3(0.0, 0.0, 0.0);

    // if there are any vertices at all
    if (vertices.size() != 0)
        { // non-empty vertex set
        // sum up all of the vertex positions
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            centreOfGravity = centreOfGravity + vertices[vertex];
        
        // and divide through by the number to get the average position
        // also known as the barycentre
        centreOfGravity = centreOfGravity / vertices.size();

        // start with 0 radius
        objectSize = 0.0;

        // now compute the largest distance from the origin to a vertex
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            { // per vertex
            // compute the distance from the barycentre
            float distance = (vertices[vertex] - centreOfGravity).length();         
            
            // now test for maximality
            if (distance > objectSize)
                objectSize = distance;
            } // per vertex
        } // non-empty vertex set


	LoopSubdivision();
	

    // return a success code
	cout << "finished reading! "<<endl;
    return true;
    } // ReadObjectStream()

// write routine
void DirectedEdgeSurface::WriteObjectStream(std::ostream &geometryStream)
    { // WriteObjectStream()
	geometryStream << "#" << std::endl; 
	geometryStream << "# Created for Leeds COMP 5821M Autumn 2020" << std::endl; 
	geometryStream << "#" << std::endl; 
	geometryStream << "#" << std::endl; 
	geometryStream << "# Surface vertices=" << vertices.size() << " faces=" << faceVertices.size()/3 << std::endl; 
	geometryStream << "#" << std::endl; 

	// output the vertices
    for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
        geometryStream << "Vertex " << vertex << " " << std::fixed << vertices[vertex] << std::endl;

    // and the normal vectors
    for (unsigned int normal = 0; normal < normals.size(); normal++)
        geometryStream << "Normal " << normal << " " << std::fixed << normals[normal] << std::endl;

	// and the first directed edges
    for (unsigned int vertex = 0; vertex < firstDirectedEdge.size(); vertex++)
        geometryStream << "FirstDirectedEdge " << vertex<< " " << std::fixed << firstDirectedEdge[vertex] << std::endl;

    // and the faces - increment is taken care of internally
    for (unsigned int face = 0; face < faceVertices.size(); )
        { // per face
        geometryStream << "Face " << face/3 << " ";
        
        // read in three vertices
        geometryStream << faceVertices[face++] << " ";
        geometryStream << faceVertices[face++] << " ";
        geometryStream << faceVertices[face++];
            
        geometryStream << std::endl;
        } // per face

	// and the other halves
	for (unsigned int dirEdge = 0; dirEdge < otherHalf.size(); dirEdge++)
		geometryStream << "OtherHalf " << dirEdge << " " << otherHalf[dirEdge] << std::endl;
    } // WriteObjectStream()


// private function for calculating barycentic coordinate
Cartesian3 calculateBC(Cartesian3 &v0, Cartesian3 &v1, Cartesian3 &v2, Cartesian3 &vert)
{
	float area_total = ((v0-v1).cross(v1-v2)).length();
	float alpha, beta, gamma;
	alpha = ((v2-v1).cross(vert-v1)).length() / area_total;
    beta  = ((v0-v2).cross(vert-v2)).length() / area_total;
    gamma = ((v1-v0).cross(vert-v0)).length() / area_total;

	return Cartesian3(alpha, beta, gamma);
}



void DirectedEdgeSurface::LoopSubdivision()
{
	unsigned int originalVertsNum = vertices.size();
	unsigned int originalFacetNum = faceVertices.size();
	unordered_map<int,int> myHash;
	unordered_map<int,vector<int>> neighbours_dictionary;


	// for (auto vert : faceVertices){
	// 		cout << vert << endl;
			
	// 	}
    // cout <<"========================================================================================="<<endl;

	// loop every face
	for (unsigned int i = 0; i < originalFacetNum; i += 3)
	{
		// 3 verts in this face
		unsigned int v0 = faceVertices[i];
		unsigned int v1 = faceVertices[i + 1];
		unsigned int v2 = faceVertices[i + 2];
		
		// this will record the newly calculated vertices
		unsigned int v3 = -1;

		int newVert_id_1 = -1;
		int newVert_id_2 = -1;
		int newVert_id_3 = -1;

		// push the neighbours of the old vertex into there neighbour vector.
		neighbours_dictionary[faceVertices[i]].push_back(faceVertices[i+1]);		
		neighbours_dictionary[faceVertices[i]].push_back(faceVertices[i+2]);
		neighbours_dictionary[faceVertices[i+1]].push_back(faceVertices[i]);		
		neighbours_dictionary[faceVertices[i+1]].push_back(faceVertices[i+2]);
		neighbours_dictionary[faceVertices[i+2]].push_back(faceVertices[i]);		
		neighbours_dictionary[faceVertices[i+2]].push_back(faceVertices[i+1]);
		

		// -----------------handle the edge v2-v0:
		if (!myHash[i])
		{
			// otherHalf of edge v2-v0:
			int other_v2v0 = otherHalf[i];

			// triangle shares the edge v2-v0:
			int tri_1 = ceil(other_v2v0 / 3);

			// calculate new point between v2-v0:
			int index_mod = (otherHalf[i]) % 3;
			v3 = -1;
			switch(index_mod){
				case 0:
					v3 = faceVertices[tri_1 * 3 + 1];
					break;
				case 1:
					v3 = faceVertices[tri_1 * 3 + 2];
					break;
				case 2:
					v3 = faceVertices[tri_1 * 3];
					break;
			}
			// position of newVert
			Cartesian3 newVert_1 = 3.f/8.f * vertices[v2] + 3.f/8.f * vertices[v0] + 1.f/8.f * vertices[v1] + 1.f/8.f * vertices[v3];
			// normal of newvert
			Cartesian3 bc_1 = calculateBC(vertices[v0], vertices[v1], vertices[v2], newVert_1);
			Cartesian3 newNorm_1 = (normals[v0] * bc_1.x + normals[v1] * bc_1.y + normals[v2] * bc_1.z).unit();

			// mark this newVert and prevent repeating adding this vertex
			newVert_id_1 = vertices.size();
			myHash[i] = newVert_id_1;
			myHash[other_v2v0] = newVert_id_1;
			
			vertices.push_back(newVert_1);
			normals.push_back(newNorm_1);

		}
		else
		{
			// prevent repeating adding this vertex
			newVert_id_1 = myHash[i];
		}



		// -----------------handle the edge v0-v1:
		if (!myHash[i + 1])
		{
			// otherHalf of edge v0-v1:
			int other_v0v1 = otherHalf[i+1];

			// triangle shares the edge v0-v1:
			int tri_2 = ceil(other_v0v1 / 3);

			// calculate new point between v0-v1:
			int index_mod = (otherHalf[i + 1]) % 3;
			v3 = -1;
			switch(index_mod){
				case 0:
					v3 = faceVertices[tri_2 * 3 + 1];
					break;
				case 1:
					v3 = faceVertices[tri_2 * 3 + 2];
					break;
				case 2:
					v3 = faceVertices[tri_2 * 3];
					break;
			}
			// position of newVert
			Cartesian3 newVert_2 = 3.f/8.f * vertices[v0] + 3.f/8.f * vertices[v1] + 1.f/8.f * vertices[v2] + 1.f/8.f * vertices[v3];

			// normal of newVert
			Cartesian3 bc_2 = calculateBC(vertices[v0], vertices[v1], vertices[v2], newVert_2);
			Cartesian3 newNorm_2 = (normals[v0] * bc_2.x + normals[v1] * bc_2.y + normals[v2] * bc_2.z).unit();
			
			// mark this newVert and prevent repeating adding this vertex
			newVert_id_2 = vertices.size();
			myHash[i+1] = newVert_id_2;
			myHash[other_v0v1] = newVert_id_2;
			
			vertices.push_back(newVert_2);
			normals.push_back(newNorm_2);


		}
		else
		{
			// prevent repeating adding this vertex
			newVert_id_2 = myHash[i+1];
		}




		// -----------------handle the edge v1-v2:
		if (!myHash[i + 2])
		{
			// otherHalf of edge v1-v2:
			int other_v1v2 = otherHalf[i+2];

			// triangle shares the edge v1-v2:
			int tri_3 = ceil(other_v1v2 / 3);

			// calculate new point between v2-v0:
			int index_mod = (otherHalf[i + 2]) % 3;
			v3 = -1;
			switch(index_mod){
				case 0:
					v3 = faceVertices[tri_3 * 3 + 1];
					break;
				case 1:
					v3 = faceVertices[tri_3 * 3 + 2];
					break;
				case 2:
					v3 = faceVertices[tri_3 * 3];
					break;
			}
			// position of newVert
			Cartesian3 newVert_3 = 3.f/8.f * vertices[v1] + 3.f/8.f * vertices[v2] + 1.f/8.f * vertices[v0] + 1.f/8.f * vertices[v3];
			// normal of newVert
			Cartesian3 bc_3 = calculateBC(vertices[v0], vertices[v1], vertices[v2], newVert_3);
			Cartesian3 newNorm_3 = (normals[v0] * bc_3.x + normals[v1] * bc_3.y + normals[v2] * bc_3.z).unit();

			// mark this newVert and prevent repeating adding this vertex
			newVert_id_3 = vertices.size();
			myHash[i+2] = newVert_id_3;
			myHash[other_v1v2] = newVert_id_3;
			
			vertices.push_back(newVert_3);
			normals.push_back(newNorm_3);


		}
		else
		{
			// prevent repeating adding this vertex
			newVert_id_3 = myHash[i+2];
		}

		// push new "top" face
		faceVertices.push_back(newVert_id_2);
		faceVertices.push_back(newVert_id_1);
		faceVertices.push_back(faceVertices[i]);

		// push new "left" face
		faceVertices.push_back(faceVertices[i+1]);
		faceVertices.push_back(newVert_id_3);
		faceVertices.push_back(newVert_id_2);

		// push new "right" face
		faceVertices.push_back(newVert_id_3);
		faceVertices.push_back(faceVertices[i+2]);
		faceVertices.push_back(newVert_id_1);
		

		// change the old face
		faceVertices[i] = newVert_id_1;
		faceVertices[i+1] = newVert_id_2;
		faceVertices[i+2] = newVert_id_3;
		
		
		

		
	}

	// update the old vertices	
	
	vector<Cartesian3> oldVerts_update;

	// for each old vertex
	for (unsigned int i = 0; i < originalVertsNum; i++)
	{
		// get its neighbours
		vector<int> neighbours = neighbours_dictionary[i];

		// remove the repeating neighbours
    	sort(neighbours.begin(), neighbours.end());
		auto it = unique(neighbours.begin(), neighbours.end());
		neighbours.erase(it, neighbours.end());

		// count of neighbours
		float n = neighbours.size();
		// use the formula: 1-beta(n) * self.pos + sum( beta * neighbours.pos )
		float beta = 1.f/n *(5.f/8.f - pow(3.f/8.f + 1.f/4.f * cos(2.f * PI / n), 2));
		Cartesian3 beta_sumVert;
		for (auto neighbour : neighbours)
		{
			beta_sumVert  = beta_sumVert + beta * vertices[neighbour];
			//cout << neighbour<<endl;
		}

		//record the updated values
		oldVerts_update.push_back(vertices[i] * (1.f -  n * beta) + beta_sumVert);
		
	}

	
	for (unsigned int i = 0; i < oldVerts_update.size(); i++)
	{
		// update position of old vertices
		vertices[i] = oldVerts_update[i];
	}

	
	//=============================== for test ================================

	// for (auto id : firstDirectedEdge){
	// 	cout << id << endl;
	// }

	// for (auto vert : faceVertices){
	// 		cout << vert << endl;
	// 	}
	// 	cout << originalVertsNum << " | " << vertices.size()<< endl;
	// 	cout << originalFacetNum <<  " | "<< faceVertices.size() <<endl;
	// 	cout << normals.size();
	//=============================== for test ================================
	updateFirstDiredges();
	//updateOldNormals(originalVertsNum);
	updateOtherhalves();

	return;

	
}

void DirectedEdgeSurface::updateFirstDiredges()
{
	// clear the vector
	firstDirectedEdge.clear();
	unsigned int vertexCount = vertices.size();
	unordered_map<int, int> myHash;

	// for each vertex of the faces
	for (unsigned int faceVertex = 0; faceVertex < faceVertices.size() || myHash.size() < vertexCount; faceVertex++)
	{
		int vertexID = faceVertices[faceVertex];

		// if there are still some firstDirEdges to find...
		if (myHash.find(vertexID) == myHash.end())
		{	
			int mod = faceVertex % 3;
			switch (mod)
			{
				case 0:
					myHash[vertexID] = faceVertex + 1;
					break;
				case 1:
					myHash[vertexID] = faceVertex + 1;
					break;
				case 2:
					myHash[vertexID] = faceVertex - 2;
					break;
			}
		}
	}

	// push the results
	for (unsigned int vertex = 0; vertex < vertices.size(); vertex ++)
	{
		if (myHash.find(vertex) == myHash.end())
		{
			firstDirectedEdge.push_back(-1);
		}

		else 
		{
			firstDirectedEdge.push_back(myHash[vertex]);
		}
	}
}



/**===================================
 * 
 * this function attempt to update the normals of the old vertices. 
 * But it seems the normals of the old vertices should remain unchanged...
 * so no need to call the function below
 * 
===================================**/


// void DirectedEdgeSurface::updateOldNormals(int originalVertNum)
// {
// 	// for each old vertex
// 	for (unsigned int i = 0; i < originalVertNum; i++)
// 	{
// 		// find which triangle includes this vertex
// 		int faceID = int(firstDirectedEdge[i] / 3);
// 		// get 3 vertices of this triangle
// 		Cartesian3 v0 = vertices[faceVertices[faceID * 3]];
// 		Cartesian3 v1 = vertices[faceVertices[faceID * 3 + 1]];
// 		Cartesian3 v2 = vertices[faceVertices[faceID * 3 + 2]];
// 		// get the two edges of from & to this vertex
// 		int mod = firstDirectedEdge[i] % 3;
// 		Cartesian3 prev, next, curr;
// 		if (mod == 0)
// 		{
// 			prev = v2;
// 			curr = v0;
// 			next = v1;
// 		}
// 		else if (mod == 1)
// 		{
// 			prev = v0;
// 			curr = v1;
// 			next = v2;
// 		}
// 		else
// 		{
// 			prev = v1;
// 			curr = v2;
// 			next = v0;
// 		}
// 		// calculate normal
// 		normals[i] = (curr - prev).cross(next - curr).unit();
// 	}
// }

void DirectedEdgeSurface::updateOtherhalves()
{
	// clear otherHalf vector
	otherHalf.clear();
	unsigned int edgeCount = faceVertices.size();
	vector<int> myOtherHalf;
	myOtherHalf.resize(edgeCount, -1);

	// for each point of the faces
	for (unsigned int vertex = 0; vertex < faceVertices.size(); vertex ++)
	{
		// if this edge does not have otherHalf yet
		if (myOtherHalf[vertex] != -1) continue;
		unsigned int preID = vertex % 3 == 0 ? vertex + 2 : vertex - 1;

		// the edge is from preId to currID
		unsigned int v1 = faceVertices[preID];
		unsigned int v2 = faceVertices[vertex];

		// find its otherHalf
		for (unsigned int other = vertex + 1; other < faceVertices.size(); other++)
		{
			if (faceVertices[other] == v2)
			{
				unsigned int otherNext = other % 3 == 2 ? other - 2 : other + 1;
				if (faceVertices[otherNext] == v1)
				{
					// record this pair
					myOtherHalf[vertex] = otherNext;
					myOtherHalf[otherNext] = vertex;
				}
			}
		}
		
	}

	// push results
	for (auto& other : myOtherHalf)
	{
		otherHalf.push_back(other);
	}

}


// routine to render
void DirectedEdgeSurface::Render(RenderParameters *renderParameters)
    { // Render()
    // Ideally, we would apply a global transformation to the object, but sadly that breaks down
    // when we want to scale things, as unless we normalise the normal vectors, we end up affecting
    // the illumination.  Known solutions include:
    // 1.   Normalising the normal vectors
    // 2.   Explicitly dividing the normal vectors by the scale to balance
    // 3.   Scaling only the vertex position (slower, but safer)
    // 4.   Not allowing spatial zoom (note: sniper scopes are a modified projection matrix)
    //
    // Inside a game engine, zoom usually doesn't apply. Normalisation of normal vectors is expensive,
    // so we will choose option 2.  

    // Scale defaults to the zoom setting
    float scale = renderParameters->zoomScale;
    scale /= objectSize;
        
    //  now scale everything
    glScalef(scale, scale, scale);

    // apply the translation to the centre of the object if requested
    glTranslatef(-centreOfGravity.x, -centreOfGravity.y, -centreOfGravity.z);

    // start rendering
    glBegin(GL_TRIANGLES);

	// set colour for pick render - ignored for regular render
	glColor3f(1.0, 1.0, 1.0);

    // loop through the faces
	for (unsigned int face = 0; face < faceVertices.size(); face +=3)
		{ // per face
		// if we want flat normals, compute them here
		if (renderParameters->useFlatNormals)
			{ // flat normals
			// find two vectors along edges of the triangle
			Cartesian3 pq = vertices[faceVertices[face+1]] - vertices[faceVertices[face]];
			Cartesian3 pr = vertices[faceVertices[face+2]] - vertices[faceVertices[face]];

			// take their cross product and normalise
			Cartesian3 faceNormal = pq.cross(pr).unit();

			// and use it to set the glNormal
			glNormal3f(faceNormal.x * scale, faceNormal.y * scale, faceNormal.z * scale);
			} // flat normals

		// we have made a HARD assumption that we have enough normals
		for (unsigned int vertex = face; vertex < face+3; vertex++)
			{ // per vertex
		
			// if we are using smooth normals
			if (!renderParameters->useFlatNormals)
				// set the normal vector
				glNormal3f
					(
					normals[faceVertices[vertex]].x * scale,
					normals[faceVertices[vertex]].y * scale,
					normals[faceVertices[vertex]].z * scale
					);
			
			// and set the vertex position
			glVertex3f
				(
				vertices[faceVertices[vertex]].x,
				vertices[faceVertices[vertex]].y,
				vertices[faceVertices[vertex]].z
				);

			} // per vertex

		} // per face

    // close off the triangles
    glEnd();
    
    // now we add a second loop to render the vertices if desired
    if (!renderParameters->showVertices)
    	return;

	glDisable(GL_LIGHTING);

	// loop through the vertices
	for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
		{ // per vertex
		// use modelview matrix (not most efficient solution, but quickest to code)
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glTranslatef(vertices[vertex].x, vertices[vertex].y, vertices[vertex].z);
		glScalef(0.1 * renderParameters->vertexSize, 0.1 * renderParameters->vertexSize, 0.1 * renderParameters->vertexSize);
		renderTriangulatedSphere();
		glPopMatrix();
		} // per vertex 
    
    } // Render()

