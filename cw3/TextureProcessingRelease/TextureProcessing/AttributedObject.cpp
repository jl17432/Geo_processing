///////////////////////////////////////////////////
//
//  Hamish Carr
//  September, 2020
//
//  ------------------------
//  AttributedObject.cpp
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
//	Variant on TexturedObject that stores explicit RGB
//	values for each vertex
//
///////////////////////////////////////////////////

// include the header file
#include "AttributedObject.h"

// include the C++ standard libraries we want
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <algorithm>

// include the Cartesian 3- vector class
#include "Cartesian3.h"

#define MAXIMUM_LINE_LENGTH 1024
#define REMAP_TO_UNIT_INTERVAL(x) (0.5 + (0.5*(x)))
#define REMAP_FROM_UNIT_INTERVAL(x) (-1.0 + (2.0*(x)))

#define N_ITERATIONS 100000

// constructor will initialise to safe values
AttributedObject::AttributedObject()
    : centreOfGravity(0.0,0.0,0.0)
    { // AttributedObject()
    // force arrays to size 0
    vertices.resize(0);
    colours.resize(0);
    normals.resize(0);
    textureCoords.resize(0);
    firstDirectedEdge.resize(0);
    faceVertices.resize(0);
    otherHalf.resize(0);
    } // AttributedObject()

// read routine returns true on success, failure otherwise
bool AttributedObject::ReadObjectStream(std::istream &geometryStream)
    { // ReadObjectStream()

    // create a read buffer
    char readBuffer[MAXIMUM_LINE_LENGTH];

    // the rest of this is a loop reading lines & adding them in appropriate places
    while (true)
        { // not eof
        // character to read
        char firstChar = geometryStream.get();

//         std::cout << "Read: " << firstChar << std::endl;

        // check for eof() in case we've run out
        if (geometryStream.eof())
            break;

        // otherwise, switch on the character we read
        switch (firstChar)
            { // switch on first character
            case '#':       // comment line
                // read and discard the line
                geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
                break;

            case 'v':       // vertex data of some type
                { // some sort of vertex data
                // retrieve another character
                char secondChar = geometryStream.get();

                // bail if we ran out of file
                if (geometryStream.eof())
                    break;

                // now use the second character to choose branch
                switch (secondChar)
                    { // switch on second character
                    case ' ':       // space - indicates a vertex
                        { // vertex read
                        Cartesian3 vertex;
                        geometryStream >> vertex;
                        vertices.push_back(vertex);
//                         std::cout << "Vertex " << vertex << std::endl;
                        break;
                        } // vertex read
                    case 'c':       // c indicates colour
                        { // normal read
                        Cartesian3 colour;
                        geometryStream >> colour;
                        colours.push_back(colour);
//                         std::cout << "Colour " << colour << std::endl;
                        break;
                        } // normal read
                    case 'n':       // n indicates normal vector
                        { // normal read
                        Cartesian3 normal;
                        geometryStream >> normal;
                        normals.push_back(normal);
//                         std::cout << "Normal " << normal << std::endl;
                        break;
                        } // normal read
                    case 't':       // t indicates texture coords
                        { // tex coord
                        Cartesian3 texCoord;
                        geometryStream >> texCoord;
                        textureCoords.push_back(texCoord);
//                         std::cout << "Tex Coords " << texCoord << std::endl;
                        break;
                        } // tex coord
                    default:
                        break;
                    } // switch on second character
                break;
                } // some sort of vertex data

            case 'f':       // face data
                { // face
                // make a hard assumption that we have a single triangle per line
                unsigned int vertexID;

                // read in three vertices
                for (unsigned int vertex = 0; vertex < 3; vertex++)
                    { // per vertex
                    // read a vertex ID
                    geometryStream >> vertexID;

                    // subtract one and store them (OBJ uses 1-based numbering)
                    faceVertices.push_back(vertexID-1);
                    } // per vertex
                break;
                } // face

            // default processing: do nothing
            default:
                break;

            } // switch on first character

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

// 	std::cout << "Centre of Gravity: " << centreOfGravity << std::endl;
// 	std::cout << "Object Size:       " << objectSize << std::endl;

    EntranceFunc();

    // return a success code
    return true;
    } // ReadObjectStream()

// write routine
void AttributedObject::WriteObjectStream(std::ostream &geometryStream)
    { // WriteObjectStream()
    geometryStream << "# " << (faceVertices.size()/3) << " triangles" << std::endl;
    geometryStream << std::endl;

    // output the vertex coordinates
    geometryStream << "# " << vertices.size() << " vertices" << std::endl;
    for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
        geometryStream << "v  " << std::fixed << vertices[vertex] << std::endl;

    // output the vertex colours
    geometryStream << "# " << colours.size() << " vertex colours" << std::endl;
    for (unsigned int vertex = 0; vertex < colours.size(); vertex++)
        geometryStream << "vc " << std::fixed << colours[vertex] << std::endl;

    // output the vertex normals
    geometryStream << "# " << normals.size() << " vertex normals" << std::endl;
    for (unsigned int vertex = 0; vertex < normals.size(); vertex++)
        geometryStream << "vn " << std::fixed << normals[vertex] << std::endl;

    // output the vertex coords
    geometryStream << "# " << textureCoords.size() << " vertex tex coords" << std::endl;
    for (unsigned int vertex = 0; vertex < textureCoords.size(); vertex++)
        geometryStream << "vt " << std::fixed << textureCoords[vertex] << std::endl;

    // and the faces
    for (unsigned int face = 0; face < faceVertices.size(); face+=3)
        { // per face
        geometryStream << "f";

        // loop through # of vertices
        for (unsigned int vertex = 0; vertex < 3; vertex++)
            { // per vertex
            geometryStream << " ";
            geometryStream << faceVertices[face+vertex] + 1;
            } // per vertex
        // end the line
        geometryStream << std::endl;
        } // per face

    } // WriteObjectStream()

// routine to render
void AttributedObject::Render(RenderParameters *renderParameters)
    { // Render()
    // make sure that textures are disabled
    glDisable(GL_TEXTURE_2D);

    float scale = renderParameters->zoomScale;
    scale /= objectSize;
    // Scale defaults to the zoom setting
    glTranslatef(-centreOfGravity.x * scale, -centreOfGravity.y * scale, -centreOfGravity.z * scale);

    if (renderParameters->useWireframe)
        glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    else
        glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

    // start rendering
    glBegin(GL_TRIANGLES);

    // loop through the faces: note that they may not be triangles, which complicates life
    for (unsigned int face = 0; face < faceVertices.size(); face+=3)
        { // per face

        // now do a loop over three vertices
        for (unsigned int vertex = 0; vertex < 3; vertex++)
            { // per vertex
            // set colour using vertex ID

            // use texCoords as colours
            if(renderParameters->useTexCoords)
            {
                glColor3f
                    (
                    textureCoords[faceVertices[face+vertex]].x,
                    textureCoords[faceVertices[face+vertex]].y,
                    textureCoords[faceVertices[face+vertex]].z
                    );
            }

            // use normals as colours
            else if (renderParameters->useNormal || renderParameters->renderNormalMap)
            {
                glColor3f
                    (
                    normals[faceVertices[face+vertex]].x,
                    normals[faceVertices[face+vertex]].y,
                    normals[faceVertices[face+vertex]].z
                    );
            }
            else
            {
                // use original rgb colours
                glColor3f
                    (
                    colours[faceVertices[face+vertex]].x,
                    colours[faceVertices[face+vertex]].y,
                    colours[faceVertices[face+vertex]].z
                    );
            }


                // render original colours in uvw plane
                if (renderParameters->renderTexture)
                {
                    glVertex3f
                        (
                        scale * textureCoords[faceVertices[face+vertex]].x,
                        scale * textureCoords[faceVertices[face+vertex]].y,
                        scale * textureCoords[faceVertices[face+vertex]].z
                        );
                }

                // render normals on the uvw plane
                else if (renderParameters->renderNormalMap)
                {
                    glVertex3f
                        (
                        scale * textureCoords[faceVertices[face+vertex]].x,
                        scale * textureCoords[faceVertices[face+vertex]].y,
                        scale * textureCoords[faceVertices[face+vertex]].z
                        );
                }

                //use scaled xyz for vertex position
                else
                {
                    glVertex3f
                        (
                        scale * vertices[faceVertices[face+vertex]].x,
                        scale * vertices[faceVertices[face+vertex]].y,
                        scale * vertices[faceVertices[face+vertex]].z
                        );
                }
            } // per vertex
        } // per face

    // close off the triangles
    glEnd();

    // revert render mode
    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

    } // Render()



void AttributedObject::EntranceFunc()
{
    ConstructDirEdge();
    GetEdges();
    VerticesClassification();
    GenerateNeighbourVec();
    SetBoundaryToSquare();
    SetInteriorVertsToCentre();
    Floater();
}


void AttributedObject::ConstructDirEdge()
{
    // =============================looking for otherhalves...==============================
    cout << "looking for otherHalves..." << endl << endl;

    // clear and resize the otherhalf vector
    otherHalf.clear();
    otherHalf.resize(faceVertices.size(), -1);

    // loop over all the edges:
    for (unsigned int vertex = 0; vertex < faceVertices.size(); vertex++)
    {
        // if the edge does not have the otherhalf yet
        if(otherHalf[vertex] != -1) continue;
        unsigned int preID = vertex % 3 == 0 ? vertex + 2 : vertex - 1;

        // get the two vertices on the edge:
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
                    otherHalf[vertex] = otherNext;
                    otherHalf[otherNext] = vertex;
                }
            }

        }

    }

    cout << "FINISHED looking for otherhalves! " << endl << endl;
    cout << "===========================================================" << endl;
    cout << "===========================================================" << endl;

    // =============================end of looking for otherhalves...==============================


    // =============================looking for firstDirEdges...==============================
//    cout << "looking for first directed edges..."<< endl << endl;

//    // clear the vector
//        firstDirectedEdge.clear();
//        unsigned int vertexCount = vertices.size();
//        unordered_map<int, int> myHash;

//        // for each vertex of the faces
//        for (unsigned int faceVertex = 0; faceVertex < faceVertices.size() || myHash.size() < vertexCount; faceVertex++)
//        {
//            int vertexID = faceVertices[faceVertex];

//            // if there are still some firstDirEdges to find...
//            if (myHash.find(vertexID) == myHash.end())
//            {
//                int mod = faceVertex % 3;
//                switch (mod)
//                {
//                    case 0:
//                        myHash[vertexID] = faceVertex + 1;
//                        break;
//                    case 1:
//                        myHash[vertexID] = faceVertex + 1;
//                        break;
//                    case 2:
//                        myHash[vertexID] = faceVertex - 2;
//                        break;
//                }
//            }
//        }

//        // push the results
//        for (unsigned int vertex = 0; vertex < vertices.size(); vertex ++)
//        {
//            if (myHash.find(vertex) == myHash.end())
//            {
//                firstDirectedEdge.push_back(-1);
//            }

//            else
//            {
//                firstDirectedEdge.push_back(myHash[vertex]);
//            }
//        }

    cout << "FINISHED looking for firstDirEdges! " << endl << endl;
    cout << "===========================================================" << endl;
    cout << "===========================================================" << endl;
    // =============================end of looking for firstDirEdges...==============================

}




void AttributedObject::GetEdges()
{
    dirEdges.resize(0);

    unordered_map<unsigned int, int> has_normal;

    for (unsigned int i = 0; i < faceVertices.size(); i += 3)
    {
        // generate half-edges for every triangles
        DirEdge dirEdge0 = DirEdge(faceVertices[i + 2], faceVertices[i]);
        DirEdge dirEdge1 = DirEdge(faceVertices[i], faceVertices[i + 1]);
        DirEdge dirEdge2 = DirEdge(faceVertices[i + 1], faceVertices[i + 2]);

        dirEdges.push_back(dirEdge0);
        dirEdges.push_back(dirEdge1);
        dirEdges.push_back(dirEdge2);


        //===================

        // calculate normals

        normals.resize(vertices.size());
        Cartesian3 v0 = vertices[dirEdge0.from];
        Cartesian3 v1 = vertices[dirEdge1.from];
        Cartesian3 v2 = vertices[dirEdge2.from];

        Cartesian3 n0 = (v0-v2).cross(v1-v0);
        Cartesian3 n1 = (v1-v0).cross(v2-v1);
        Cartesian3 n2 = (v2-v1).cross(v0-v2);

        normals[dirEdge0.from] = normals[dirEdge0.from] + n0;
        has_normal[dirEdge0.from] ++;

        normals[dirEdge1.from] = normals[dirEdge1.from] + n1;
        has_normal[dirEdge1.from] ++;

        normals[dirEdge2.from] = normals[dirEdge2.from] + n2;
        has_normal[dirEdge2.from] ++;

//        normals[dirEdge0.from] = n0.unit();
//        normals[dirEdge1.from] = n1.unit();
//        normals[dirEdge2.from] = n2.unit();


    }

    // interpolating the normals
    for (unsigned int i = 0; i < normals.size(); i++)
    {
        normals[i] = (normals[i] / has_normal[i]).unit();
    }

}

// separate interior vertices and boundary vertices
void AttributedObject::VerticesClassification()
{
    unordered_map<unsigned int, unsigned int> dictionary;

    cout << "starting Vertices Classification..." << endl << endl;

    // build a dictionary, pairing the two vertices of all the boundary edges
    cout << "building dictionary..." << endl << endl;

    // this is the logic for finding the same starting vertex as the one in the result Hamish shows,
    // for making the head lying at the diagonal line.
    unsigned int begin = 0;
    float currDis = std::numeric_limits<double>::infinity();
    for (unsigned int edge = 0; edge < otherHalf.size(); edge++)
    {
        if (otherHalf[edge] == -1)
        {
            unsigned int preID = (edge % 3 == 0) ? edge + 2 : edge - 1;
            dictionary[faceVertices[preID]] = faceVertices[edge];
            float dis = vertices[faceVertices[edge]].x - centreOfGravity.x;
            if (dis < currDis)
            {
                currDis = dis;
                begin = faceVertices[edge];
            }
        }
    }
    cout << "FINISHED building dictionary!" << endl << endl;
    cout << "===========================================================" << endl;
    cout << "===========================================================" << endl;


//    auto pair = dictionary.begin();
//    unsigned int start = pair->first;
//    unsigned int next  = pair->second;
//    boundary.push_back(start);

    // ordering the boundary vertices
    unsigned int start = begin;
    unsigned int next = dictionary[start];
    boundary.push_back(start);

    while (next != start)
    {
        if (next == boundary[-1]) break;
        boundary.push_back(next);
        next = dictionary[next];
    }

    cout << "FINISH finding the boundaries!" << endl << endl;

    // finding all the interior vertices
    for (unsigned int i = 0; i < vertices.size(); i++)
    {
        if (dictionary.count(i) == 0)
        {
            interiorVertices.push_back(i);
        }
    }

    cout << "FINISH finding interiorVertices!" << endl << endl;
    cout << "===========================================================" << endl;
    cout << "===========================================================" << endl;


}

void AttributedObject::GenerateNeighbourVec()
{
    // finding the neighbours of every vertices
    neighbourVec.resize(vertices.size());
    for (auto& dirEdge : dirEdges)
    {
        // to is a neighbour of from
        neighbourVec[dirEdge.from].push_back(dirEdge.to);
        sort(neighbourVec[dirEdge.from].begin(), neighbourVec[dirEdge.from].end());
        neighbourVec[dirEdge.from].erase(unique(neighbourVec[dirEdge.from].begin(), neighbourVec[dirEdge.from].end()), neighbourVec[dirEdge.from].end());

        // from is a neighbour of to
        neighbourVec[dirEdge.to].push_back(dirEdge.from);
        sort(neighbourVec[dirEdge.to].begin(), neighbourVec[dirEdge.to].end());
        neighbourVec[dirEdge.to].erase(unique(neighbourVec[dirEdge.to].begin(), neighbourVec[dirEdge.to].end()), neighbourVec[dirEdge.to].end());
    }
}


// mapping the boundary vertices onto the edges of the square
void AttributedObject::SetBoundaryToSquare()
{
    textureCoords.resize(vertices.size());

    int mod = boundary.size() % 4;
    int vNumOnEachEdge = boundary.size() / 4;
    float gapLength = 1.0f / vNumOnEachEdge;


    // 1st edge:
    for (unsigned int i = 0; i < vNumOnEachEdge; i++)
    {
        textureCoords[boundary[i]] = Cartesian3(0 + i * gapLength, 0, 0);
    }
    
    // 2nd edge:
    for (unsigned int i = vNumOnEachEdge; i < 2 * vNumOnEachEdge; i++)
    {
        textureCoords[boundary[i]] = Cartesian3(1, 0 + (i - vNumOnEachEdge) * gapLength, 0);
    }

    // 3rd edge
    for (unsigned int i = 2 * vNumOnEachEdge; i < 3 * vNumOnEachEdge; i++)
    {
        textureCoords[boundary[i]] = Cartesian3(1 - (i - 2 * vNumOnEachEdge) * gapLength, 1, 0);
    }

    // 4th edge
    gapLength = 1.0f / (vNumOnEachEdge + mod);
    for (unsigned int i = 3 * vNumOnEachEdge; i < boundary.size(); i++)
    {
        textureCoords[boundary[i]] = Cartesian3(0, 1 - (i - 3 * vNumOnEachEdge)* gapLength, 0);
    }

}

// put every interior vertices on to the centre of the square
void AttributedObject::SetInteriorVertsToCentre()
{
    for (auto& interiorVert : interiorVertices)
    {
        textureCoords[interiorVert] = Cartesian3(0.5f, 0.5f, 0);
    }
}

// implement the floater-tutte algorithm
void AttributedObject::Floater()
{
    int n = 0;
    vector<Cartesian3> temp;
    temp.resize(interiorVertices.size());

    // this is the maxDisplacement amongst all interior vertices
    float maxDisplacement = 0;

    do{
        // calculate the updated position for all interior vertices
        maxDisplacement = 0;
        for (unsigned int i = 0; i < interiorVertices.size(); i++)
        {
            int numOfNeighbours = neighbourVec[interiorVertices[i]].size();
            Cartesian3 sum;
            for (auto& neighbour : neighbourVec[interiorVertices[i]])
            {
                sum = sum + textureCoords[neighbour];
            }
            sum = sum / numOfNeighbours;
            temp[i] = sum;

        }

        // update the new position symulteniously for all interior vertices
        for (unsigned int i = 0; i < interiorVertices.size(); i++)
        {
            maxDisplacement = max(maxDisplacement, (temp[i] - textureCoords[interiorVertices[i]]).length());
            float test = textureCoords[interiorVertices[i]].length();
            textureCoords[interiorVertices[i]] = temp[i];
        }
        cout << "iteration time: " << n << endl;
        n++;
    } while (maxDisplacement > 0.000001f); // compare the maxDisplacement to the tolerance. break if maxDisplacement is smaller
}
