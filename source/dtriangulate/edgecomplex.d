module dtriangulate.edgecomplex;

/* Stores a set of edges describing a planar subdivision

 */


//These names are from the Guibas and Stolfi paper more or less
struct HalfEdge {
  size_t org, dest; //vertex indices

  //index of an half edge whose "org" == this.dest, who share the same left face
  //IOW, the index of the next half edge on a counterclockwise traversal of the polygon
  //on the left side of this half edge
  size_t nextL;

  size_t sym; //index of the half edge going from dest to org

  //TODO use different types for vertex indices and edge indices
}



struct EdgeComplex{
  import std.array;
  private HalfEdge[] edges; //TODO use appender if this is too slow
  private size_t[] freeEdgeList; //indices of entries in edges that are unused


  ///Create a new edge from a to b which is 2 half edges.  Return index of the half-edge from a to b
  size_t newEdge(size_t a, size_t b){
    const e1Index = allocEdge();
    const e2Index = allocEdge();

    HalfEdge e1 = {a, b, e2Index, e2Index};
    HalfEdge e2 = {b, a, e1Index, e1Index};

    edges[e1Index] = e1;
    edges[e2Index] = e2;
    return e1Index;
  }

  size_t newDanglingEdge(size_t e1Index, size_t newVertexIndex){
    const ret = allocEdge();
    const retSym = allocEdge();

    const eStart = edges[e1Index];
    HalfEdge newEdge = {eStart.dest, newVertexIndex, retSym, retSym};
    HalfEdge newEdgeSym = {newVertexIndex, eStart.dest, eStart.nextL, ret};

    edges[ret] = newEdge;
    edges[retSym] = newEdgeSym;

    edges[e1Index].nextL = ret;

    return ret;

  }


  size_t newEdgeBetween(size_t e1Index, size_t e2Index){

    auto e1 = edges[e1Index];
    auto e2 = edges[e2Index];

    const newIndex1 = allocEdge();
    const newIndex2 = allocEdge();

    HalfEdge newE1 = {e1.dest, e2.org, e2Index, newIndex2};
    HalfEdge newE2 = {e2.org, e1.dest, e1.sym, newIndex1};

    edges[newIndex1] = newE1;
    edges[newIndex2] = newE2;

    //clean nLeft for affected edges
    edges[edges[e1.nextL].sym].nextL = newIndex1;
    edges[prevL(e2Index)].nextL = newIndex2;

    return newIndex1;
  }

  //CCW a -> b -> c
  size_t newTriangle(size_t a, size_t b, size_t c){
    auto e1 = newEdge(a, b);
    auto e2 = newDanglingEdge(e1, c);
    newEdgeBetween(e2, e1);
    return e1;
  }


  size_t prevL (size_t eIndex) const {
    size_t ret = edges[eIndex].sym;
    size_t limit = 0;
    while(edges[ret].nextL != eIndex  && ++limit < 1000){
      ret = edges[edges[ret].nextL].sym; //rotate around the vertex
    }
    assert(limit < 1000);
    return ret;
  }

  //read only public indexing
  HalfEdge opIndex(size_t index) const{
    return edges[index];
  }

  string toString() const {
    import std.conv : to;
    string ret;
    foreach(i, he; edges){
      ret ~= to!string(i) ~ ": " ~ to!string(he) ~ '\n';
    }
    return ret;
  }

  private size_t allocEdge(){
    import std.range: empty, popBack, back;
    if(!freeEdgeList.empty){
      auto ret = freeEdgeList.back;
      freeEdgeList.popBack();
      return ret;
    }
    auto ret = edges.length;
    edges ~= HalfEdge.init;
    return ret;
  }
}


unittest {

  import std.stdio;

  EdgeComplex ec;

  auto e1 = ec.newEdge(10, 11);
  auto e2 = ec.newDanglingEdge(e1, 12);
  auto e3 = ec.newEdgeBetween(e2, e1);

  ec.newTriangle(13, 14, 15);

  writeln("EDGECOMPLEX TEST!!!");
  writeln(ec);

  writeln();
}
