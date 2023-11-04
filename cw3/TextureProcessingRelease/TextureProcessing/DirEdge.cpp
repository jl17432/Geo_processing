#include "DirEdge.h"

DirEdge::DirEdge()
    : from(-1), to(-1)
    {}

DirEdge::DirEdge(long v1, long v2)
    : from(v1), to(v2)
    {}
