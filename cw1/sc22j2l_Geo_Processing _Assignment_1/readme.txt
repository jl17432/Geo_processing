# University of Leeds 2022-23
# COMP 5812M Assignment 1 Triangle Soup Renderer

=================== Makefile Information ========================================================================================
I had provided a .makefile which should be able build the code. 

I had also tested on the University Linux machines and the code was built successfully

However, just in case if it does not build, the following code could help:

[userid@machine Assignment2MeshDSHandout]$ module add legacy-eng
[userid@machine Assignment2MeshDSHandout]$ module add qt/5.15.2
[userid@machine Assignment2MeshDSHandout]$ qmake -project QT+=opengl
[userid@machine Assignment2MeshDSHandout]$ qmake
[userid@machine Assignment2MeshDSHandout]$ make


Once the code is built, to execute the renderer, pass the file name on the command line:

---[userid@machine Assignment2MeshDSHandout]$ ./triangle_renderer ../handout_models/FILENAME.tri

please replace "FILENAME" with the file need to be tested, if the file is good, the program 

will automatically:

	first render the object, this is to visualize the geometry shape
	
	then generate its .face file and the .diredge file in the same folder with the original .tri file




================= Introduction ===========================================================================================

This is my submission to the first assignment for COMP 5812M, including:

    1) a very simple renderer for triangle soup. It makes no attempt to be efficient, 
    which keeps the code clearer. ----This part was credit to and provided by Hamish, 
    I did not contribute to that part and I did not make any change in that part

    2) face2faceindex.cpp and its header file. -----Wrote by me
    this file will read a .tri file, compute the face index format and 
    output a .face file.

    3) faceindex2directedge.cpp and its header file. -----Wrote by me
    this file will do:

        a. read a face file, compute the direct edges and 
           output a .diredge file  

        b. determine the file provided whether it is manifold or not,
           print out the vertex that fails if the file is not a manifold

        c. if the file is tested to be a manifold, using the 
           Euler formula to determine the number of genus


=================== Analysis of the Algorithm Complexity ===================================================================

-----note: n stands for num of vertex data in the file

face2faceindex.cpp:

    The main algorithm in this file are
        
        1. reading all vertex in -----this will cost O(n) of time 

        2. generating a vertexIndex vector, the method will loop all the vertex
           at the mean time it loops all the vertices before the current vertex.
           -----this will cost O(n^2) of time. 

        Time Complexity
        
        The space Complexity of the algorithm is ~O(2n)
        

faceindex2directedge.cpp:

    The main algorithm in this file are

        1. reading all the vertex in ----- O(n)

        2. construct direct edge vector ----- O(n/3) ~O(n)

        3. write vertex bolck -----O(n)

        4. look for the first directedge: for each vertex, 
           loop the direct edge vector ----- ~O(n^2)

        5. write face block -----O(n)

        6. write otherHalf block: for every edge, loop 
           the direct edge vector ----- ~O(n^2)

            -----finishes task b, total time complexity is O(2n^2 + 4n) ~ O(2n^2), 

		spatial complexity is:

                myVerts       -----O(n)
                myFaces       -----O(n)
                myFirst       -----O(n)
                myOtherhalves -----O(n)

                O(4n) max in total at runtime


        7. manifold test:

            BFS: for the vertex and all its neighbours, compute all 
            the edges From this vertex, loop them -----O(2n^2)

        8. calculate genus -----O(1)

            -----finished manifold testing and genus galculating
            total time complexity is O(2n^2) ~ O(n^2), 

            space complexity is:

            myEdge         -----O(n)
            myFirst        -----O(n)
            myOtherhalves  -----O(n)
            myQueue        -----O(n)
            myHash         -----O(n)

            O(5n) max in total at run time
