#include "face2faceindex.h"
#include "faceindex2directedge.h"

/*
//=============================================================================================



This file format begins with a header block with identifying information, and basic information about the file. It is similar in nature to the standard .obj file, but much simpler, as it leaves out normals, texture coordinates, &c.

You should have a header and two blocks in the following order:
	1.	Vertex block
	2.	Face block

In the vertex block, each vertex is listed by number with its xyz coordinates on a single line.

In the face block, each face is listed by number with its face indices on a single line in CCW order.



//=============================================================================================
*/

vector<Cartesian3> myVertsVec;

bool face2faceindex(char *fileName)
{
    cout << endl <<"reading file path: " << fileName << endl;

    ifstream inFile(fileName);
	if (inFile.bad()) 
	{
		cout << "cannot read file: " << fileName << ", please double check the file name." << endl;
		return false;
	}
	// set the number of vertices and faces
	long nTriangles = 0, nVertices = 0;

	// read in the number of vertices
	inFile >> nTriangles;
	nVertices = nTriangles * 3;

    // now loop to read the vertices in, and hope nothing goes wrong
	for (int vertex = 0; vertex < nVertices; vertex++)
		{ // for each vertex
			float x,y,z;
			inFile >> x >> y >> z;
			myVertsVec.push_back(Cartesian3(x,y,z));

		//=======================for test=========================== to print out raw vertex vector
        // cout<<"\n";
        // cout<<"size: "<< nVertices << " now: "<< vertex <<"\n";
        // cout<< myVertsVec[vertex];
        // cout<<"\n";
		//=======================for test===========================


		} // for each vertex



	//now loop to set vertex IDs for each vertex
	vector<int> vertexIdVec(myVertsVec.size(), -1);
	//cout<< "vertexID.size(): " <<vertexIdVec.size() << "\n";

	//set initial vertex ID to be 0
	int next_vertex_ID = 0;

	//for each vertex in myVertsVec
	for (long unsigned int vertex = 0; vertex < vertexIdVec.size(); vertex++)
	{
		//loop through to see if the vertex is already exist
		for (long unsigned int other = 0; other < vertex; other++)
		{
			//already exist, grab the vertexID previously set
			if (myVertsVec[vertex] == myVertsVec[other])
			{
				vertexIdVec[vertex] = vertexIdVec[other];
			}
		}

		//is a new point, set a new ID
		if (vertexIdVec[vertex] == -1)
		{
			vertexIdVec[vertex] = next_vertex_ID++;
		}
	}

	//=======================for test=========================== to print out the vertexId vector
	// for (auto id : vertexIdVec){
	// 	cout<<" "<<id;
	// }
	//=======================for test===========================
    

	cout << "Finish reading!" << endl << endl;





    //=============================now write everything into the output file======================================




	cout << "Start generating .face file..." << endl;

	//construct the name of the file
	string output_filename = fileName;		//grab infile name
	output_filename = output_filename.substr(0, output_filename.length() - 4);		//remove ".tri"
	output_filename += ".face";		//add suffix ".face"
	string obj_name = output_filename.substr(18, output_filename.length() - 23);

	//output file
	ofstream output_file;
	//creat a new file with clearing all its content
	output_file.open(output_filename, ios::trunc);

	//something wrong when opening the file
	if (!output_file)
	{
		cout << "\n" << "can not open this file : " << output_filename << "\n";
		return false;
	}
	
	//nothing wrong when opening the file
	else
	{
		//title and intro
		output_file << "# University of Leeds 2022-2023" << endl;
		output_file << "# COMP 5812M Assignment 1" << endl;
		output_file << "# Name: Jingxuan Liu" << endl;
		output_file << "# Student Number: 201606569" << endl;
		output_file << "# " << endl;
		output_file << "# Object Name: " << obj_name << endl;
		output_file << "# Vertice=" << next_vertex_ID << " " << "Faces=" << myVertsVec.size()/3 << endl;
		output_file << "# " << endl;

		//writing vertex info
		int writingID = 0;

		//loop to write vertex:
		for (long unsigned int vertex = 0; vertex < myVertsVec.size(); vertex++)
		{
			//if first time found
			if (vertexIdVec[vertex] == writingID)
			{
				//write a vertex format
				cout.setf(ios::fixed);
				output_file << "Vertex " << writingID << " " << 
					setw(7) << setfill('0') << setiosflags(ios::fixed) << setprecision(6) <<
					myVertsVec[vertex].x << "  " <<
					myVertsVec[vertex].y << "  " <<
					myVertsVec[vertex].z << "  " <<
					endl;
				
				writingID++;
			}
		}

		//loop to write face:
		for (long unsigned int face = 0; face < myVertsVec.size() / 3; face++)
		{
			//write a face format
			output_file << "Face " << face << "  " <<
				vertexIdVec[3 * face]     << "  " <<
				vertexIdVec[3 * face + 2] << "  " <<
				vertexIdVec[3 * face + 1] << "  " <<
				endl;

		}

		//finish writing, close file
		output_file.close();
		cout << obj_name <<".face file has been generated!" << endl << endl;
	}


	// call the function in faceindex2directedge.cpp, to generate a .directedge file.

	//pass the filename
	char *passing_file_name = const_cast<char*>(output_filename.c_str());
	faceindex2directedge(passing_file_name, next_vertex_ID, myVertsVec.size()/3);


    return true;
}
