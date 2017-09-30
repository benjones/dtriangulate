module dtriangulate.triangulate;

import dtriangulate.ssoVector;
import dtriangulate.predicates;

import gl3n.linalg : vec2;

import std.stdio;

struct Pair{ int first, second; }

bool leftOf(Vec)(Vec p, Vec e1, Vec e2){
  return orient2D(e1, e2, p) > 0;
}

bool rightOf(Vec)(Vec p, Vec e1, Vec e2){
  return orient2D(e1, e2, p) < 0;
}


struct Triangle{

  this(int a, int b, int c){
	v[] = [a, b, c];
  }

  int opIndex(int i){
	return v[i > 2 ? i -3 : i]; //faster than modulo?
  }
	  
  bool opEquals(Triangle rhs){
	return v[0] == rhs.v[0] && v[1] == rhs.v[1] && v[2] == rhs.v[2];
  }
  
  private int[3] v;
}


struct TriDB{

  enum GHOST = -1;
  this(ulong numPoints){
	triangles = new SSOVector!(Pair, 9)[numPoints];
	ghostEdges = new int[numPoints];
	ghostEdges[] = GHOST;
  }

  void addTriangle(Triangle t){
	foreach(i; 0..3){
	  if(t[i] != GHOST){
		triangles[t[i]] ~= Pair(t[i+1], t[i+2]);
	  } else {
		ghostEdges[t[i+1]] = t[i+2];
	  }
	}
  }

  void deleteTriangle(Triangle t){
	foreach(i; 0..3){
	  int u = t[i];
	  auto vw = Pair(t[i+1], t[i+2]);
	  if(u != GHOST){
		foreach(ref pr ; triangles[u]){
		  if(pr == vw){
			pr = triangles[u][$-1];
			triangles[u].popBack();
			break;
		  }
		}
	  } else {
		ghostEdges[vw.first] = GHOST;
	  }
	}
  }

  int adjacent(int u, int v) const{
	if(u != GHOST){
	  foreach(const ref pr; triangles[u]){
		if(pr.first == v){
		  return pr.second;
		}
	  }
	} else {
	  return ghostEdges[v];
	}
	assert(false);
	//return GHOST;
  }
  
  bool adjacentExists(int u, int v) const{
	if(u != GHOST){
	  foreach(const ref pr ; triangles[u]){
		if(pr.first == v){
		  return true;
		}
	  }
	} else {
	  return ghostEdges[v] != GHOST;
	}
	return false;
  }

  Pair adjacentTriangle(int u) const{
	if(u != GHOST){
	  return triangles[u].length == 0 ? Pair(GHOST, GHOST) : triangles[u][0];
	} else {
	  foreach(int i, ge; ghostEdges){
		if(ge != GHOST){ return Pair(i, ge); }
	  }
	}
	assert(false);
  }

  Triangle[] getTriangles() const{
	Triangle[] ret;
	foreach(int i, const ref svec; triangles){
	  foreach(const ref pr; svec){
		if((pr.first == GHOST || i < pr.first) && (pr.second == GHOST || i < pr.second)){
		  ret ~= Triangle(i, pr.first, pr.second);
		}
	  }
	}
	return ret;
  }


  void dump() const {
	import std.stdio;
	writeln("normal entries");
	foreach(set ; triangles){
	  writeln(set);
	}
	writeln("ghost edges");
	writeln(ghostEdges);
	
  }

private:

  SSOVector!(Pair, 9)[] triangles;
  int[] ghostEdges;
}


import std.typecons;

alias CCWEdge = Typedef!(Pair, Pair.init, "ccw");
alias CWEdge = Typedef!(Pair, Pair.init, "cw");

CWEdge reverse(CCWEdge e){
  return CWEdge(Pair(e.second, e.first));
}
CCWEdge reverse(CWEdge e){
  return CCWEdge(Pair(e.second, e.first));
}

void partitionPoints(Vec, bool ByX)(const Vec[] points, int[] indices){
  import std.algorithm;
  static if(ByX){
	indices.topN!(delegate bool(int a, int b){return points[a].x < points[b].x; })(indices.length/2);
  } else {
	indices.topN!(delegate bool(int a, int b){return points[a].y < points[b].y; })(indices.length/2);
  }
}

CWEdge hullAdvance(const ref TriDB triDB, CWEdge e){
  return CWEdge(Pair(e.second, triDB.adjacent(TriDB.GHOST, e.second)));
}

CCWEdge hullAdvance(const ref TriDB triDB, CCWEdge e){
  return CCWEdge(Pair(e.second, triDB.adjacent(e.second, TriDB.GHOST)));
}


CWEdge hullMaxXCW(Vec)(const ref TriDB triDB, const ref Vec[] points, CWEdge e){
  while(points[e.second].x >= points[e.first].x){
	e = hullAdvance(triDB, e);
  }
  int prev = triDB.adjacent(e.first, TriDB.GHOST);
  while(points[prev].x > points[e.first].x){
	e = CWEdge(Pair(prev, e.first));
	prev = triDB.adjacent(e.first, TriDB.GHOST);
  }
  return e;
}

CCWEdge hullMinXCCW(Vec)(const ref TriDB triDB, const ref Vec[] points, CCWEdge e){

  while(points[e.second].x <= points[e.first].x){
	e = hullAdvance(triDB, e);
  }
  int prev = triDB.adjacent(TriDB.GHOST, e.first);
  while(points[prev].x < points[e.first].x){
	e = CCWEdge(Pair(prev, e.first));
	prev = triDB.adjacent(TriDB.GHOST, e.first);
  }
  return e;
}

CCWEdge hullMaxYCCW(Vec)(const ref TriDB triDB, const ref Vec[] points, CCWEdge e){
  while(points[e.second].y >= points[e.first].y){
	e = hullAdvance(triDB, e);
  }
  int prev = triDB.adjacent(TriDB.GHOST, e.first);
  while(points[prev].y > points[e.first].y){
	e = CCWEdge(Pair(prev, e.first));
	prev = triDB.adjacent(TriDB.GHOST, e.first);
  }
  return e;
}

CWEdge hullMinYCW(Vec)(const ref TriDB triDB, const ref Vec[] points, CWEdge e){
  while(points[e.second].y <= points[e.first].y){
	e = hullAdvance(triDB, e);
  }
  int prev = triDB.adjacent(e.first, TriDB.GHOST);
  while(points[prev].y < points[e.first].y){
	e = CWEdge(Pair(prev, e.first));
	prev = triDB.adjacent(e.first, TriDB.GHOST);
  }
  return e;
}


//The CW edge will start at the bottom, and the CCW will start at the top
//or the CW edge will be the right-most, and the CCW edge will be the left-most
struct EdgePair{ CWEdge cwEdge; CCWEdge ccwEdge;}


EdgePair delaunayBaseCase(Vec, bool ByX)(ref TriDB triDB, const ref Vec[] points, int[] indices){
  import std.algorithm : minElement, maxElement;

  float yMap(int a){ return points[indices[a]].y; }
  float xMap(int a){ return points[indices[a]].x; }
  
  if(indices.length == 2){
	triDB.addTriangle(Triangle(indices[0], indices[1], TriDB.GHOST));
	triDB.addTriangle(Triangle(indices[1], indices[0], TriDB.GHOST));
	static if(ByX){
	  //returning to vertical, so bottommost CW, topmost CCW
	  if(points[indices[0]].y < points[indices[1]].y){
		return EdgePair(CWEdge(Pair(indices[0], indices[1])),
						CCWEdge(Pair(indices[1], indices[0])));
	  } else {
		return EdgePair(CWEdge(Pair(indices[1], indices[0])),
						CCWEdge(Pair(indices[0], indices[1])));
	  }
	} else {
	  //returning to horizontal, so return leftmost CCW, rightmost CW
	  if(points[indices[0]].x < points[indices[1]].x){
		return EdgePair(CWEdge(Pair(indices[1], indices[0])),
						CCWEdge(Pair(indices[0], indices[1])));
	  } else {
		return EdgePair(CWEdge(Pair(indices[0], indices[1])),
						CCWEdge(Pair(indices[1], indices[0])));
	  }
	  
	}
	
  } else {
	assert(indices.length == 3);
	auto o2d = orient2D(points[indices[0]], points[indices[1]], points[indices[2]]);

	if(o2d > 0){
	  //positively oriented, add the triangle as is
	  triDB.addTriangle(Triangle(indices[0], indices[1], indices[2]));
	  //and the appropriate ghosts
	  triDB.addTriangle(Triangle(indices[1], indices[0], TriDB.GHOST));
	  triDB.addTriangle(Triangle(indices[2], indices[1], TriDB.GHOST));
	  triDB.addTriangle(Triangle(indices[0], indices[2], TriDB.GHOST));

	  static if(ByX){
		//returning to vertical

		int minIndex = [0,1,2].minElement!yMap;
		int maxIndex = [0,1,2].maxElement!yMap;
		
		
		//int minIndex = minElement!"points[a].y"([0,1,2]);
		//		int maxIndex = maxElement!"points[a].y"([0,1,2]);
		//CW from bottom, CCW from top
		return EdgePair(CWEdge(Pair(indices[minIndex], indices[(minIndex + 2)%3])),
						CCWEdge(Pair(indices[maxIndex], indices[(maxIndex + 1)%3])));
		
	  } else {
		
		int minIndex = [0,1,2].minElement!xMap;
		int maxIndex = [0,1,2].maxElement!xMap;
		//CW from right, CCW from left
		return EdgePair(CWEdge(Pair(indices[maxIndex], indices[(maxIndex + 2)%3])),
						CCWEdge(Pair(indices[minIndex], indices[(minIndex + 1)%3])));

	  }
	  
	} else if(o2d < 0){
	  //negatively oriented, swap 1 and 2
	  triDB.addTriangle(Triangle(indices[0], indices[2], indices[1]));
	  triDB.addTriangle(Triangle(indices[2], indices[0], TriDB.GHOST));
	  triDB.addTriangle(Triangle(indices[1], indices[2], triDB.GHOST));
	  triDB.addTriangle(Triangle(indices[0], indices[1], triDB.GHOST));

	  static if(ByX){
		//returning to vertical
			int minIndex = [0,1,2].minElement!yMap;//minElement!"points[a].y"([0,1,2]);
		int maxIndex = [0,1,2].maxElement!yMap;//"points[a].y"([0,1,2]);
		//CW from Bottom, CCW from top
		return EdgePair(CWEdge(Pair(indices[minIndex], indices[(minIndex +1)%3])),
						CCWEdge(Pair(indices[maxIndex], indices[(maxIndex + 2)%3])));
	  } else {
		//return to horizontal
		int minIndex = [0,1,2].minElement!xMap;
		int maxIndex = [0,1,2].maxElement!xMap;
		//CW from right, CCW from left
		return EdgePair(CWEdge(Pair(indices[maxIndex], indices[(maxIndex + 1)%3])),
						CCWEdge(Pair(indices[minIndex], indices[(minIndex + 2)%3])));
		
	  }
	} else {
	  //colinear, but unsorted!
	  //sort by opposite of ByX so we know what to return
	  import std.algorithm.sorting : sort;
	  static if(ByX){
		indices.sort!(delegate bool(int a, int b){return points[a].y < points[b].y;});
	  } else {
		indices.sort!(delegate bool(int a, int b){return points[b].x < points[a].x;});
	  }
	  triDB.addTriangle(Triangle(indices[0], indices[1], TriDB.GHOST));
	  triDB.addTriangle(Triangle(indices[1], indices[0], TriDB.GHOST));
	  triDB.addTriangle(Triangle(indices[1], indices[2], TriDB.GHOST));
	  triDB.addTriangle(Triangle(indices[2], indices[1], TriDB.GHOST));

	  static if(ByX){
		//return to vertical, CW from bottom, CCW from top
		return EdgePair(CWEdge(Pair(indices[0], indices[1])),
						CCWEdge(Pair(indices[2], indices[1])));
	  } else {
		//return to horizontal, CW from right, CCW from left
		return EdgePair(CWEdge(Pair(indices[2], indices[1])),
						CCWEdge(Pair(indices[0], indices[1])));
	  }
	}
	
  }
}






//ldi, inner clockwise pointing edge from left half
//rdi, inner ccw pointing edge from right half
//connecting the two hulls will be CW, either right to left for horizontal,
//or bottom to top for vertical
CWEdge getLowerHullEdge(Vec)(const ref TriDB triDB, const ref Vec[] points,
							 CWEdge ldi, CCWEdge rdi){

  while(true){
	if(leftOf(points[rdi.first], points[ldi.first], points[ldi.second])){
	  ldi = hullAdvance(triDB, ldi);
	} else if(rightOf(points[ldi.first], points[rdi.first], points[rdi.second])){
	  rdi = hullAdvance(triDB, rdi);
	} else {
	  break;
	}
	
  }
  return CWEdge(Pair(rdi.first, ldi.first));
}

EdgePair zipHulls(Vec, bool ByX)(ref TriDB triDB, const ref Vec[] points, CWEdge rToL){
  SSOVector!(int, 4) leftHull;
  SSOVector!(int, 4) rightHull;

  leftHull ~= triDB.adjacent(rToL.second, TriDB.GHOST);
  rightHull ~= triDB.adjacent(TriDB.GHOST, rToL.first);

  //add an external edge for the lower hull
  triDB.deleteTriangle(Triangle(leftHull.back(), rToL.second, TriDB.GHOST));
  triDB.deleteTriangle(Triangle(rToL.first, rightHull.back(), TriDB.GHOST));
  triDB.addTriangle(Triangle(rToL.first, rToL.second, TriDB.GHOST));

  CWEdge start = rToL; //save this for later
  
  //are there ghost triangles that need to be removed?
  //not to start with, since we just removed them above
  bool leftGhost = false, rightGhost = false;
  while(true){

	int lCand;
	if(!leftHull.empty){
	  lCand = leftHull.back;
	  leftHull.popBack();
	} else {
	  lCand = triDB.adjacent(rToL.second, TriDB.GHOST);
	  leftGhost = true;
	}

	bool lValid = rightOf(points[lCand], points[rToL.first], points[rToL.second]);

	if(lValid){
	  int nextCand = triDB.adjacent(rToL.second, lCand);
	  while(nextCand != TriDB.GHOST &&
			isInCircle(points[rToL.second], points[rToL.first], points[lCand], points[nextCand])){
		//nextCand would be hosed by triangle rtoL.second, rToL.first, lCand, so delete edge
		leftHull ~= lCand;
		if(leftGhost){
		  triDB.deleteTriangle(Triangle(lCand, rToL.second, TriDB.GHOST));
		}
		triDB.deleteTriangle(Triangle(rToL.second, lCand, nextCand));

		leftGhost = false; //just deleted the ghost, if there was one
		lCand = nextCand;
		nextCand = triDB.adjacent(rToL.second, lCand);
	  }
	}

	int rCand;
	if(!rightHull.empty){
	  rCand = rightHull.back;
	  rightHull.popBack();
	} else {
	  rCand = triDB.adjacent(TriDB.GHOST, rToL.first);
	  rightGhost = true;
	}

	bool rValid = rightOf(points[rCand], points[rToL.first], points[rToL.second]);

	if(rValid){
	  int nextCand = triDB.adjacent(rCand, rToL.first);
	  while(nextCand != TriDB.GHOST &&
			isInCircle(points[rToL.second], points[rToL.first], points[rCand], points[nextCand])){
		rightHull ~= rCand;

		if(rightGhost){
		  triDB.deleteTriangle(Triangle(rToL.first, rCand, TriDB.GHOST));
		}
		triDB.deleteTriangle(Triangle(rCand, rToL.first, nextCand));

		rightGhost = false;
		rCand = nextCand;
		nextCand = triDB.adjacent(rCand, rToL.first);
	  }
	}

	if(!lValid && !rValid){ break; }
	if(!lValid ||
	   (rValid && isInCircle(points[lCand], points[rToL.second],
							 points[rToL.first], points[rCand]))){
	  //pick rCand
	  if(rightGhost){
		triDB.deleteTriangle(Triangle(rToL.first, rCand, TriDB.GHOST));
		rightGhost = false;
	  }
	  triDB.addTriangle(Triangle(rToL.second, rToL.first, rCand));
	  rToL.first = rCand;
	  leftHull ~= lCand;
	} else {
	  if(leftGhost){
		triDB.deleteTriangle(Triangle(lCand, rToL.second, TriDB.GHOST));
		leftGhost = false;
	  }
	  triDB.addTriangle(Triangle(rToL.second, rToL.first, lCand));
	  rToL.second = lCand;
	  rightHull ~= rCand;
	}
  }

  //add ghost for the top
  triDB.addTriangle(Triangle(rToL.second, rToL.first, TriDB.GHOST));

  CCWEdge stop = CCWEdge(Pair(rToL.first, rToL.second));
  static if(ByX){
	//return to vertical, CW from bottom, CCW from top
	//believe the hull calls will be no-ops, but maybe not
	return EdgePair(hullMinYCW!Vec(triDB, points, start),
					hullMaxYCCW!Vec(triDB, points, stop));
  } else {
	//return to horizontal, CW from right, CCW from left
	return EdgePair(hullMaxXCW!Vec(triDB, points, reverse(stop)),
					hullMinXCCW!Vec(triDB, points, reverse(start)));
	
  }
}


EdgePair delaunayRecurse(Vec, bool ByX)(ref TriDB triDB, const  Vec[] points, int[] indices){

  assert(indices.length > 1);

  if(indices.length < 4){
	return delaunayBaseCase!(Vec,ByX)(triDB, points, indices);
  } else {
	partitionPoints!(Vec,ByX)(points, indices);

	ulong middle = indices.length/2;

	auto ep1 = delaunayRecurse!(Vec, !ByX)(triDB, points, indices[0..middle]);
	auto ep2 = delaunayRecurse!(Vec, !ByX)(triDB, points, indices[middle..$]);

	//because ep2 will be below, if byY, and ep2 will be to the right if byX
	static if(ByX){
	  CWEdge ldi = ep1.cwEdge;
	  CCWEdge rdi = ep2.ccwEdge;
	} else {
	  CWEdge ldi = ep2.cwEdge;
	  CCWEdge rdi = ep1.ccwEdge;
	}
	
	CWEdge rToL = getLowerHullEdge(triDB, points, ldi, rdi);

	return zipHulls!(Vec,ByX)(triDB, points, rToL);
	
  }
}


TriDB delaunayTriangulate(Vec)(const Vec[] points){
  import std.array;
  import std.range: iota;
  
  TriDB triDB = TriDB(points.length);
  if(points.length < 2){ return triDB; }
  
  int[] indices = array(iota(0, cast(int)points.length));

  delaunayRecurse!(Vec,true)(triDB, points, indices);
  return triDB;
  
}


void writeSVG(Vec)(string filename, const Vec[] points, const Triangle[] tris){
  import std.stdio : File;
  File f = File(filename, "w");

  Vec minP = Vec(points.map("a.x").min(), points.map("a.y").min());
  Vec maxP = Vec(points.map("a.x").max(), points.map("a.y").max());

  int[] activePoints = triangles.map!(a => [a[0], a[1], a[2]])()
	.joiner()
	.filter!(a => a != TriDB.GHOST)()
	.sort()
	.uniq()
	.array();

  Vec size = maxP - minP;
  minP -= .03*size;
  maxP += .03*size;
  size = maxP - minP;

  
}

unittest{

  vec2[] points = [ vec2(0,0), vec2(1,1) ];
  auto triDB = delaunayTriangulate(points);
  auto tris = triDB.getTriangles();

  import std.stdio;

  triDB.dump();
  
  writefln("numtris: %d",  tris.length);
  foreach(Triangle  tri ; tris){
	writeln(tri);
  }
  assert(tris.length == 2);
}


unittest{

  vec2[] points  = [vec2(0,0), vec2(1,1), vec2(0.2, 0.9), vec2(0.9, 1)];
  auto triDB = delaunayTriangulate(points);
  auto tris = triDB.getTriangles();

  triDB.dump();
  foreach(const ref tri; tris){
	writeln(tri);
  }

  
  assert(tris.length == 6);

  
}

  
