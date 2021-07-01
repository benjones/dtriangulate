module dtriangulate.edgecomplex;

import std.stdio;
import std.range: iota;
import std.algorithm.iteration : filter, map;
import std.array;
import std.container.rbtree;

/* Stores a set of edges describing a planar subdivision

 */


//These names are from the Guibas and Stolfi paper more or less
struct HalfEdge {
  size_t org, dest; //vertex indices

  //index of an half edge whose "org" == this.dest, who share the same left face
  //IOW, the index of the next half edge on a counterclockwise traversal of the polygon
  //on the left side of this half edge
  EdgeId nextL;

  EdgeId sym; //index of the half edge going from dest to org

  //TODO use different types for vertex indices and edge indices
}

///Opaque index wrapper to use type safety to make sure we don't use a vertex index
///when we mean an edge index
struct EdgeId{
  private size_t index;
  @safe: nothrow:
  public bool opEquals(EdgeId other) const {
    return index == other.index;
  }

  public int opCmp(EdgeId other) const {
    import std.typecons: tuple;
    return tuple(this.tupleof).opCmp(tuple(other.tupleof));
  }

  public size_t toHash() const @safe pure nothrow{
    return index;
  }
}


EdgeId getIndex(const ref EdgeComplex ec, HalfEdge edge){
  return ec[edge.sym].sym;
}

struct EdgeComplex{
  private HalfEdge[] edges; //TODO use appender if this is too slow
  private EdgeId[] freeEdgeList; //indices of entries in edges that are unused


  ///Create a new edge from a to b which is 2 half edges.  Return index of the half-edge from a to b
  EdgeId newEdge(size_t a, size_t b){
    const e1Index = allocEdge();
    const e2Index = allocEdge();

    HalfEdge e1 = {a, b, e2Index, e2Index};
    HalfEdge e2 = {b, a, e1Index, e1Index};

    getEdge(e1Index) = e1;
    getEdge(e2Index) = e2;
    return e1Index;
  }

  EdgeId newDanglingEdge(EdgeId e1Index, size_t newVertexIndex){
    const ret = allocEdge();
    const retSym = allocEdge();

    const eStart = this[e1Index];
    HalfEdge newEdge = {eStart.dest, newVertexIndex, retSym, retSym};
    HalfEdge newEdgeSym = {newVertexIndex, eStart.dest, eStart.nextL, ret};

    getEdge(ret) = newEdge;
    getEdge(retSym) = newEdgeSym;

    getEdge(e1Index).nextL = ret;

    return ret;

  }

  ///Create a new edge going from e1.dest to e2.org
  EdgeId newEdgeBetween(EdgeId e1Index, EdgeId e2Index){
    auto e1 = this[e1Index];
    auto e2 = this[e2Index];

    const originalE2Prev = prevL(e2Index);

    writeln("new edge between  e1: ", edgeToString(e1), " e2: ", edgeToString(e2));
    const newIndex1 = allocEdge();
    const newIndex2 = allocEdge();

    HalfEdge newE1 = {e1.dest, e2.org, e2Index, newIndex2};
    HalfEdge newE2 = {e2.org, e1.dest, e1.nextL, newIndex1};

    getEdge(newIndex1) = newE1;
    getEdge(newIndex2) = newE2;
    writeln("new edges: ", edgeToString(newE1), edgeToString(newE2));

    //clean nLeft for affected edges

    getEdge(e1Index).nextL = newIndex1;


    //MODIFYING THE WRONG EDGE HERE.  Need something like the previous
    writeln("prevl: ", edgeToString(this[prevL(e2Index)].nextL));
    getEdge(originalE2Prev).nextL = newIndex2;
    writeln("modified :", edgeToString(this[e1Index]), " and ", edgeToString(this[originalE2Prev]));

    writeln("ec after neb: \n", this.toString);
    return newIndex1;
  }


  void deleteEdge(EdgeId ei){

  }


  //CCW a -> b -> c
  EdgeId newTriangle(size_t a, size_t b, size_t c){
    auto e1 = newEdge(a, b);
    auto e2 = newDanglingEdge(e1, c);
    newEdgeBetween(e2, e1);
    return e1;
  }

  ///Makes a new triangle from ccw vertices a, b, c, d.  Returns helf-edge a->b
  EdgeId newQuad(size_t a, size_t b, size_t c, size_t d){
    auto e1 = newEdge(a, b);
    auto e2 = newDanglingEdge(e1, c);
    auto e3 = newDanglingEdge(e2, d);
    newEdgeBetween(e3, e1);
    return e1;
  }

  EdgeId prevL (EdgeId eIndex) const {
    auto ret = this[eIndex].sym;
    size_t limit = 0;
    while(this[ret].nextL != eIndex  && ++limit < 1000){
      ret = this[this[ret].nextL].sym; //rotate around the vertex
    }
    assert(limit < 1000);
    return ret;
  }

  ///If the half edge's nextL may have been modified, get the up to date version by sym-ing twice
  deprecated HalfEdge sync(HalfEdge he) {
    return this[this[he.sym].sym];
  }

  //read only public indexing
  HalfEdge opIndex(EdgeId index) const{
    return edges[index.index];
  }

  private ref HalfEdge getEdge(EdgeId index){
    return edges[index.index];
  }

  string toString() const {
    import std.conv : to;
    string ret;
    foreach(i, he; edges){
      ret ~= to!string(i) ~ ": " ~ to!string(he) ~ '\n';
    }
    return ret;
  }

  string edgeToString(EdgeId i){
    return edgeToString(this[i]);
  }

  string edgeToString(HalfEdge e){
    import std.format;
    return format!"HalfEdge(org: %s, dest: %s, nextL: (%s, %s))"(e.org, e.dest, this[e.nextL].org, this[e.nextL].dest);
  }

  ///gives a list of all half edges in use.  It has to filter unused edges and make a copy, so this is slow-ish
  const(HalfEdge)[] getEdgeList() const{


    auto freeEdgeSet = redBlackTree(freeEdgeList);
    return iota(edges.length).map!(i => EdgeId(i)).filter!(i => i !in freeEdgeSet).map!(i => this[i]).array;
  }

  size_t[][] polygons() const {
    size_t[][] ret;
    auto freeEdgeSet = freeEdgeList.redBlackTree();
    auto usedEdgeList = iota(edges.length).map!(i => EdgeId(i)).filter!(i => i !in freeEdgeSet).redBlackTree();

    while(!usedEdgeList.empty){
      size_t[] polygon;
      const startIndex = usedEdgeList.front;
      EdgeId edgeIndex = startIndex;
      for(; this[edgeIndex].nextL != startIndex; edgeIndex = this[edgeIndex].nextL){
        assert(edgeIndex in usedEdgeList);
        usedEdgeList.removeKey(edgeIndex);
        polygon ~= this[edgeIndex].org;
      }
      polygon ~= this[edgeIndex].org; //finish the polygon
      usedEdgeList.removeKey(edgeIndex);
      ret ~= polygon;
    }
    return ret;
  }

  ///return the length of the polygon this half edge is part of.
  ///This will assert that the the polygon is convex and nondegenerate (only 2-gons can be 0 volume)
  size_t polygonSize(Vec)(const Vec[] points, HalfEdge he){
    import dtriangulate.predicates;
    size_t ret = 1;
    const startIndex = getIndex(this, he);
    for(auto e = he; e.nextL != startIndex && ret < 1000; e = this[e.nextL]){
      ++ret;
      const nextVert = this[e.nextL].dest;
      if(nextVert == e.org){
        assert(ret == 2);
        assert(e.nextL == startIndex);
        return ret;
      }
      auto o2d = orient2D(points[e.org], points[e.dest], points[nextVert]);
      assert(o2d > 0);

    }
    return ret;
  }


  private EdgeId allocEdge(){
    import std.range: empty, popBack, back;
    if(!freeEdgeList.empty){
      auto ret = freeEdgeList.back;
      freeEdgeList.popBack();
      return ret;
    }
    auto ret = EdgeId(edges.length);
    edges ~= HalfEdge.init;
    return ret;
  }


}



void writeSVG(Vec)(string filename, const Vec[] points, const ref EdgeComplex ec){
  import dtriangulate.svg;


  writefln("dumping %s", filename);

  bool[EdgeId] freeEdges;
  foreach(fe; ec.freeEdgeList){
    freeEdges[fe] = true;
  }

  auto edges = ec.getEdgeList();

  /*  size_t[] activePoints;
  foreach(edge; ec.edges){
    if(edge.org < edge.dest){
      activePoints ~= edge.org;
      activePoints ~= edge.dest;
    }
    }*/

  size_t[] activePoints = iota(points.length).array;
  auto f = File(filename, "w");
  auto scaledPoints = prepareSVG(f, points, activePoints);
  foreach(edge; edges){
    if(edge.org < edge.dest){
      svgLine(f, scaledPoints[edge.org], scaledPoints[edge.dest], 0.02);
    }
  }

  svgLabeledPoints(f, scaledPoints, activePoints);
  f.writeln("</g>\n</svg>");
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
