To compile:

module add legacy-eng
module add qt/5.13.0
qmake -project QT+=opengl
qmake
make

./TextureProcessing ../TextureProcessing/models/hamishbig.obj

(please replace 'hamishbig.obj' to any other obj files when testing)

===================================================

I had attempt to finished all 7 questions. 


For the first 6 questions, my result looks very similar to the 

result Hamish provided in the doc. The programme took 4804 

iterations while applying the floater and tutte algorithm on 

'hamishbigtrunc.obj'and setting the tolerance to be 0.000001.



However, there was a significant difference in the 7th question

while calculating the normals.

===================================================

I have built a class and its header file, namely 'DirEdge.cpp' 

and 'DirEdge.h'. Apart from these two new files, all the other 

code was added in 'AttributedObject.cpp'



Thanks, really enjoyed playing with the picture of the lecturer's head! XD
