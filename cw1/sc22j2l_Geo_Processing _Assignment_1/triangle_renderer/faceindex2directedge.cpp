#include "faceindex2directedge.h"

/*
//=============================================================================================



This file format is an extended version of the .face file format. It begins with the same header block, and adds additional blocks to the file.

You should have a header and four blocks in the following order:
	1.	Vertex block (as for the .face file)
	2.	First Directed Edge block
	3.	Face block (as for the .face file)
	4.	Other Half block

The reason why we insert 2. before 3. rather than after it is that the Vertex & FirstDirectedEdge blocks are the same size.

In the vertex block, each vertex is listed by number with its xyz coordinates on a single line.

In the first directed edge block, the ID of the vertex is listed, with the ID of the first directed edge at each vertex. Note that this is not canonical, as any directed edge from the vertex may be chosen.

In the face block, each face is listed by number with its face indices on a single line in CCW order.

In the other half block, each directed edge is listed by number, with the other directed edge pointing in the opposite direction that pairs with it


//=============================================================================================
*/


bool faceindex2directedge(char *fileName, int num_verts, int num_faces)
{
	
	cout<< "reading file path: "<< fileName << endl;

	myVerts.resize(num_verts);
	myFirsts.resize(num_verts);
	vector<Cartesian3> myFaces;
	myFaces.resize(num_faces);





	




	//======================================start reading the file==================================================

	ifstream inFile(fileName);
	if (inFile.bad())
	{
		cout << "cannot read file: " << fileName << ", please double check the file name." << endl;
		return false;
	}
		

	//reading a line and storing the line into a string
	string str;
	while (!inFile.eof())
	{
		inFile >> str;
		// cout << str<< endl;

		//reading a vertex:
		long unsigned int vertex = 0;
		if (str == "Vertex")
		{
			inFile >> vertex;
			// exceeding the length
			if (vertex > myVerts.size()-1)
			{
				cout << "num of vertice exceeding the num that has declaired in the file"<< endl;
				return false;
 			}
			//copy vertex info into the vector
			inFile >> myVerts[vertex].x >> myVerts[vertex].y >> myVerts[vertex].z;
			vertex++;
		}
		
		//reading a face:
		long unsigned int face = 0;
		if (str == "Face")
		{
			inFile >> face;
			// exceeding the length
			if (face > myFaces.size()-1)
			{
				cout << "num of faces exceeding the num that has declaired in the file"<< endl;
				return false;
 			}
			//copy faces info into the vector
			inFile >> myFaces[face].x >> myFaces[face].y >> myFaces[face].z;
		}

		//construct directed edge vector, each item is a pair of vertexFrom and vertexTo


	}

	// now construct the edge vector

	//for each face
	for (auto face : myFaces){

		//grab the three vertices
		int v0 = int (face.x);
		int v1 = int (face.y);
		int v2 = int (face.z);

		//push edges into the vector
		myEdges.push_back(pair<int, int>(v0, v1));
		myEdges.push_back(pair<int, int>(v1, v2));
		myEdges.push_back(pair<int, int>(v2, v0));
	}

	inFile.close();

	cout << "Finish reading!" << endl << endl;

	//======================================finished reading the file==================================================

















	//======================================start writing the file==================================================



	cout << "Start generating .directedge file..." << endl;

	//construct the name of the output file
	string output_filename = fileName;		//grab infile name
	output_filename = output_filename.substr(0, output_filename.length() - 5);		//remove ".face"
	output_filename += ".diredge";		//add suffix ".directedge"
	string obj_name = output_filename.substr(18, output_filename.length() - 26);

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

	else{

		



	//title and intro
		output_file << "# University of Leeds 2022-2023" << endl;
		output_file << "# COMP 5812M Assignment 1" << endl;
		output_file << "# Name: Jingxuan Liu" << endl;
		output_file << "# Student Number: 201606569" << endl;
		output_file << "# " << endl;
		output_file << "# Object Name: " << obj_name << endl;
		output_file << "# Vertice=" << num_verts << " " << "Faces=" << num_faces << endl;
		output_file << "# " << endl;

//=====================================VERTEX BLOCK===============================================
	
	// the vertices block, write vertices in
	for (long unsigned int vertex = 0; vertex < myVerts.size(); vertex++)
		{
			// writing the vertex info
			cout.setf(ios::fixed);
				output_file << "Vertex " << vertex << " " << 
					setw(7) << setfill('0') << setiosflags(ios::fixed) << setprecision(6) <<
					myVerts[vertex].x << "  " <<
					myVerts[vertex].y << "  " <<
					myVerts[vertex].z << "  " <<
					endl;
		}

	}

	




//=====================================FIRST DIRECTED EDGE BLOCK===============================================
	
	// the first directed edge block, find the first dir-edge FROM the vertex, write it in


	//for each vertex
	for (int vertex = 0; vertex < num_verts; vertex++)
	{
		//this var will record the index of the first found edge, if not found, its value will be -1
		int firstDirEdge = -1;

		//for each edge
		for (int edge = 0; edge < num_faces * 3; edge++)
		{
			//found the first dir-edge
			if (myEdges[edge].first == vertex)
			{
				//update the value from -1 to the index
				firstDirEdge = edge;
				break;
			}
		}

		//record all first dir-edges, for testing manifolds
		myFirsts[vertex] = firstDirEdge;

		// has an edge from this vertex, write it into the output file
		output_file << "FirstDirectedEdge " << vertex << "  " << firstDirEdge << endl;
		
	}


	//================for test================

	// int testID = 0;
	// for (auto item : myEdges)
	// {
	// 	output_file << "test! edge: " << testID << " " << item.first << " " << item.second << endl;
	// 	testID++;
	// }

	//================for test================






//=====================================FACE BLOCK===============================================


	// face block, write in all faces info
	for (long unsigned int face = 0; face < myFaces.size(); face++)
	{
		output_file << "Face " << face << "  " << 
					int(myFaces[face].x) << "  " <<
					int(myFaces[face].y) << "  " <<
					int(myFaces[face].z) << "  " <<
					endl;
	}



//=====================================OTHERHALF BLOCK===============================================

	
	//otherhalf block, match all the directed edges, then write in
	
	
	myOtherhalves.resize(myEdges.size());

	// for every edge, looks for its otherhalf
	for (long unsigned int edge  = 0; edge < myEdges.size(); edge++)
	{
		//this is what the otherhalf should be like
		pair<int, int> otherHalf(myEdges[edge].second, myEdges[edge].first);
		//initialize the index to be -1
		int otherhalf_id = -1;
		//start looking up
		for (long unsigned int other = 0; other < myEdges.size(); other++){
			//found the otherhalf
			if (myEdges[other].first == otherHalf.first && myEdges[other].second == otherHalf.second)
			{
				//record the index
				otherhalf_id = other;
				//end the loop
				break;
			}
		}

		//record all the otherhalves,  for testing manifolds
		myOtherhalves[edge] = otherhalf_id;
		//write in the output file
		output_file << "otherHalf " << edge << "  " << otherhalf_id << endl;
		
		
	}

	cout << obj_name <<".directedge file has been generated!" << endl << endl;


//=====================================DETERMINE IF THE OBJECT IS MANIFOLDS BLOCK===============================================

	cout << "Applying manifold tests and calculating the amount of genus..." << endl << endl;

	int bad_vertex_id = -1;
	
	// logic: doing a BFS, for every vertex it visit, grab its 1st-ring then apply the test
	//
	// the following code could be simplifide as following:
	// 		do:
	//		{ 
	//			grab the first vertex in the queue as entry point, visit its neighbours;
	//				if (the neibour hasn't been visited){
	//					push the neibour into the queue;
	//				}
	//			then apply the test on this entry test, return the bad point ID or negative value
	// 		} while (!queue.empty);


	//first entry point to be vertex 0
	myQueue.push(0);

	do
	{
		int test_result = manifoldTest(myQueue.front());
		if (test_result > 0)
		{
			bad_vertex_id = test_result;
			break;
		}
		else
		{
			myQueue.pop();
		}
	} while (!myQueue.empty());

	// if the vertex[0] cannot reach all the point in the fill, 
	// suggest there may be more than 1 object in the file.
	if ( int (myHash.size()) < num_verts - 1 )
	{
		bad_vertex_id = -2;
	}

	//===================================backup=======================================



	// int num_genus = (num_verts - int(myEdges.size() / 2) + num_faces - 2) / -2;
	
	if (bad_vertex_id == -2)
	{
		// if we got here, means there may be multiple objects inside the file
		// print out the conclusion
		cout << endl << "This object is manifold =  NO!" << endl;
		cout << "The first vertex cannot reach all the other vertices in the file, there may be multiple objects inside the file. It is not a single object. The bad point is 0" << endl;
	}

	else if (bad_vertex_id == -1)
	{
		//if we got here, means there are no bad points, aka object is a manifold

		//calculate the number of genus of the object:


		//=================================for testing====================================
		// cout <<"verts "<<num_verts<<endl;
		// cout<<"edges "<<myEdges.size()<<endl;
		// cout<<"face "<<num_faces<<endl;
		// cout<<"test" << num_verts - (myEdges.size() / 2)<<endl;
		//=================================for testing====================================


		int num_genus = (num_verts - int(myEdges.size() / 2) + num_faces - 2) / -2;
		//print out result 
		cout << endl << "This object is manifold =  YES!" << endl;
		cout <<	"The object has " << num_genus << " genus" << endl << endl;
	}
	else 
	{
		//if we got here, means there is at least 1 bad point, aka the object is not a manifold

		//print out result
		cout << endl << "This object is manifold =  NO!" << endl;
		cout << "The bad vertex is: " << bad_vertex_id << endl;
	}

	output_file.close();

    return true;
}


//the manifold test, return the bad point ID, or negative value if all the vertex passed
int manifoldTest(int vertex)
{

	int num_cycle = 0;

			//visit the vertex's first-ring
			//count the amount of edges FROM this vertex
			int num_edge_from_vertex = 0;
			for (auto edge : myEdges)
			{
				if (edge.first == vertex)
				{
					num_edge_from_vertex ++;
					//record the vertex that has been visit into hash map
					//if first time visit this vertex, push it into the queue. 
					if (myHash[edge.second] == 0)
					{
						myHash[edge.second]++;
						myQueue.push(edge.second);
					}
				}
			}
			
			

			//first edge From the vertex
			int edge_index = myFirsts[vertex];

			//check if the vertex has first dir-edge
			if (edge_index == -1)
			{
				return vertex;
			}

			int edge_count = 0;

			int tempIndex;

			// using do{...}while{...};, so that the code inside the do{} part 
			// will be excuted at least once
			do 
			{
				//record its current index, incase it does not have otherhalf
				tempIndex = edge_index;
				//get otherhalf
				edge_index = myOtherhalves[edge_index];

				if (edge_index == -1)
				{
					return myEdges[tempIndex].second;
				}

				//get the next edge of the otherhalf
				edge_index = ( (edge_index + 1) % 3 ==0 ) ? edge_index - 2 : edge_index + 1;

				//increase the edge count
				edge_count++;

				//going back to the initial edge, counts as a cycle
				if (edge_index == myFirsts[vertex])
				{
					// first cycle found
					if (num_cycle < 1)
					{
						num_cycle++;
					}

					// more than 1 cycle for one vertex, not manifold
					else
					{
						return vertex;
						break;
					}
				}

			} while(edge_count < num_edge_from_vertex);

	return -1;
}


