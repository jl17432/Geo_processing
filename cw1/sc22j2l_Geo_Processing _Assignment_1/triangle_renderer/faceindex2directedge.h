//   University of Leeds 2022-2023
//   COMP 5812M Assignment 1
//   Jingxuan Liu
//   ID : 201606569
//=====================================================================
#include <iostream>
#include <fstream>
#include <cstring>
#include <math.h>
#include <iomanip>
#include <vector>
#include <queue>
#include <unordered_map>
#include "Cartesian3.h"


using namespace std;

static vector<Cartesian3> myVerts;

static vector<pair<int, int>> myEdges;

static vector<int> myFirsts;

static vector<int> myOtherhalves;

static queue<int> myQueue;

static unordered_map<int, int> myHash;

bool faceindex2directedge(char *fileName, int num_verts, int num_faces);

int manifoldTest(int vertex);

