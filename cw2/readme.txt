To compile on feng-linux / feng-gps:

module add legacy-eng
module add qt/5.13.0
qmake -project QT+=opengl
qmake
make

To run on feng-linux / feng-gps:

./LoopSubdivisionRelease ../LoopSubdivisionRelease/diredgenormals/cube.diredgenormal

# Please change the "cube" to any other object name in the "diredgenormal" folder.


===========================================================================================

---This is a program that can read a .diredgenormal file, apply loop-subdivision,

and then render the result of the smoothed object. 

===========================================================================================

---The reading and writing .diredgenormal parts are credit to Hamish Carr, 

I am credite to the loop-subdivision part.

===========================================================================================

---The logic of my subdivision is detailedly illustrated inside the code by comment, 

The method can be concluded as  follows:

		1) for each face: calculate the 3 new vertice and their normalson the 3 edge, 
			so that it will generate 4 new faces
		2) change the original face to the new face in the middle
		3) add the rest 3 new faces to the end of face vector
		4) update the old vertices wrt. their own neighbours
		5) calculate first-dirEdge and otherhalf

===========================================================================================

---There is a diagram (order_logic.jpg) shows the logic of the order I push the new vertices

and new faces.

===========================================================================================

---The program will read a .diredgenormal file and generate a new .diredgenormal file with 

the name ends with "_new.diredgenormal" as a result. Then this newly generated file can be 

use as the next input file to reach the next loop-subdivision iteration.


