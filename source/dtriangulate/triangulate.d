module dtriangulate.triangulate;
public import dtriangulate.tridb;

import std.typecons;
import std.stdio;
import std.conv;
import std.range;


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
  //write("partition, num points: ");
  //writeln(indices.length);
  static if(ByX){
	indices.topN!(delegate bool(int a, int b){
		//lexicographically
		return points[a].x == points[b].x ? points[a].y < points[b].y : points[a].x < points[b].x;
	  })(indices.length/2);
  } else {
	indices.topN!(delegate bool(int a, int b){
		return points[a].y == points[b].y ? points[a].x < points[b].x : points[a].y < points[b].y;
	  })(indices.length/2);
  }
}



CWEdge hullAdvance(const ref TriDB triDB, CWEdge e){
  return CWEdge(Pair(e.second,
					 TriDB.unGhost(triDB.adjacent(e.first, e.second))));
}

CCWEdge hullAdvance(const ref TriDB triDB, CCWEdge e){
  return CCWEdge(Pair(e.second,
					  triDB.adjacent(e.second, TriDB.makeGhost(e.first))));
}

//go backwards along the hull
CWEdge hullReverse(const ref TriDB triDB, CWEdge e){
  return reverse(hullAdvance(triDB, reverse(e)));
	//  return CWEdge(Pair(TriDB.adjacent(triDB.adjacentGhost(e.second, e.first)),
	//					 e.first));
}

 CCWEdge hullReverse(const ref TriDB triDB, CCWEdge e){
   return reverse(hullAdvance(triDB, reverse(e)));
   //   return CCWEdge(Pair( triDB.adjacentReal(e.first, TriDB.makeGhost(e.second)),
   //						e.first));
}


void hullCheck(Vec)(const ref TriDB triDB, const ref Vec[] points, CWEdge start){
  //writeln("hull check");
  CWEdge e = hullAdvance(triDB, start);
  int cMax = 100000;
  int count = 0;
  while(count < cMax && e != start){
	CWEdge next = hullAdvance(triDB, e);
	assert(!leftOf(points[e.first], points[e.second], points[next.second]));
	e = next;
	++count;
  }
  if(count == cMax){
	writeln("didn't make it back to start");
	assert(false);
  }

  /*
  foreach(const ref tri ; triDB.triangles){

	foreach (const  ref pr ; tri){
	  if(TriDB.isGhost(pr.first)){
		if(pr.second == TriDB.unGhost(pr.first)){
		  writeln("wtf: ghost " , TriDB.unGhost(pr.first), "...", pr.second);
		}
		assert(pr.second != TriDB.unGhost(pr.first));
	  }
	  if(TriDB.isGhost(pr.second)){
		assert(pr.first != TriDB.unGhost(pr.second));
	  }
	  assert(pr.first != pr.second);
	  assert(!TriDB.isGhost(pr.first) || !TriDB.isGhost(pr.second));
	}
	
	}*/

  
}


/*
Hulls are convex, so if the x or y coordinates of adjacent points are equal,
they must be extreme values (or else the second one wouldn't be on the hull)
so getLowerHullEdge will work for either choice.

HOWEVER,
For the pathological case of collinear points, we want these edges to be
correct for BOTH ByX AND ByY splits BC getHullLowerEdge will only chase
in one direction.

 */

//want the rightmost edge.  If many edges have the same x coord
//pick the one with the smallest y coord

CWEdge hullMaxXCW(Vec)(const ref TriDB triDB, const ref Vec[] points, CWEdge e){

  while(points[e.second].x > points[e.first].x ||
		(points[e.second].x == points[e.first].x && points[e.second].y < points[e.first].y)){
	
	e = hullAdvance(triDB, e);
  }
  
  int prev = hullReverse(triDB, e).first;
  while(points[prev].x > points[e.first].x ||
		(points[prev].x == points[e.first].x && points[prev].y < points[e.first].y)){
	e = CWEdge(Pair(prev, e.first));
	prev = hullReverse(triDB, e).first;
  }
  return e;
}

//want the smallest x edge.  If there's a choice, pick the one with the topmost Y
CCWEdge hullMinXCCW(Vec)(const ref TriDB triDB, const ref Vec[] points, CCWEdge e){

  while(points[e.second].x < points[e.first].x ||
		(points[e.second].x ==  points[e.first].x && points[e.second].y > points[e.first].y)){
	e = hullAdvance(triDB, e);
  }
  
  int prev = hullReverse(triDB, e).first;
  while(points[prev].x < points[e.first].x ||
		(points[prev].x == points[e.first].x && points[prev].y > points[e.first].y)){
	e = CCWEdge(Pair(prev, e.first));
	prev = hullReverse(triDB, e).first;
  }
  return e;
}


//want the largest Y value.  If we get to pick, chose the smallest X
CCWEdge hullMaxYCCW(Vec)(const ref TriDB triDB, const ref Vec[] points, CCWEdge e){
  
  while(points[e.second].y > points[e.first].y ||
		(points[e.second].y == points[e.first].y && points[e.second].x < points[e.first].x)){
	e = hullAdvance(triDB, e);
  }
  
  int prev = hullReverse(triDB, e).first;
  while(points[prev].y > points[e.first].y ||
		(points[prev].y == points[e.first].y && points[prev].x < points[e.first].x)){
	e = CCWEdge(Pair(prev, e.first));
	prev = hullReverse(triDB, e).first;
  }
  return e;
}
//want smallest Y.  If we get to pick, choose the biggest X
CWEdge hullMinYCW(Vec)(const ref TriDB triDB, const ref Vec[] points, CWEdge e){
  while(points[e.second].y < points[e.first].y ||
		(points[e.second].y == points[e.first].y && points[e.second].x > points[e.first].x)){
	e = hullAdvance(triDB, e);
  }
  
  int prev = hullReverse(triDB, e).first;
  while(points[prev].y < points[e.first].y ||
		(points[prev].y == points[e.first].y && points[prev].x > points[e.first].x)){
	e = CWEdge(Pair(prev, e.first));
	prev = hullReverse(triDB, e).first;
  }
  
  return e;
}


//The CW edge will start at the bottom, and the CCW will start at the top
//or the CW edge will be the right-most, and the CCW edge will be the left-most
struct EdgePair{ CWEdge cwEdge; CCWEdge ccwEdge;}


//replace CCW edges ab, and bc with triangle a, c, b, fixing all necessary ghost edges
void makeTriFromTwoEdges(Vec)(ref TriDB triDB, int a, int b, int c, const ref Vec[] points){
  auto prev = hullReverse(triDB, CCWEdge(Pair(a, b))).first;
  auto next = hullAdvance(triDB, CCWEdge(Pair(b, c))).second;
  //delete the 3 ghost triangles that are now bad
  triDB.deleteTriangle(Triangle(b, a, TriDB.makeGhost(prev)));
  triDB.deleteTriangle(Triangle(c, b, TriDB.makeGhost(a)));
  triDB.deleteTriangle(Triangle(next, c, TriDB.makeGhost(b)));
  //and add the appropriate new triangles:
  triDB.addTriangle(Triangle(a, c, b));
  triDB.addTriangle(Triangle(c, a, TriDB.makeGhost(prev)));
  triDB.addTriangle(Triangle(next, c, TriDB.makeGhost(a)));
  
}



void eraseEdge(Vec)(ref TriDB triDB, CCWEdge e, const ref Vec[] points){
  //delete the two edges:
  auto prev = hullReverse(triDB, e).first;
  auto next = hullAdvance(triDB, e).second;

  //the real value MUST exist if we're calling this function
  auto newPoint = triDB.adjacent(e.first, e.second);
  assert(!TriDB.isGhost(newPoint));
  
  triDB.deleteTriangle(Triangle(e.second, e.first, TriDB.makeGhost(prev)));
  triDB.deleteTriangle(Triangle(next, e.second, TriDB.makeGhost(e.first)));
  triDB.deleteTriangle(Triangle(e.first, e.second, newPoint));
  //and add the 3 new ones
  triDB.addTriangle(Triangle(newPoint, e.first, TriDB.makeGhost(prev)));
  triDB.addTriangle(Triangle(e.second, newPoint,  TriDB.makeGhost(e.first)));
  triDB.addTriangle(Triangle(next, e.second, TriDB.makeGhost(newPoint)));
}



int meshNumber = 0;

EdgePair delaunayBaseCase(Vec, bool ByX)(ref TriDB triDB, const ref Vec[] points, int[] indices){
  import std.algorithm : minElement, maxElement;

  //writefln("base case byX: %s with %d points", ByX, indices.length);

  //foreach(ind; indices){
  //	writefln("%d: %.8f, %.8f", ind, points[ind].x, points[ind].y);
  //}
  
  auto yMap(int a){ return Tuple!(float, float)(points[indices[a]].y, points[indices[a]].x); }
  auto xMap(int a){ return Tuple!(float, float)(points[indices[a]].x, points[indices[a]].y);}

  if(indices.length == 2){
	triDB.addTriangle(Triangle(indices[0], indices[1], TriDB.makeGhost(indices[0])));
	triDB.addTriangle(Triangle(indices[1], indices[0], TriDB.makeGhost(indices[1])));
	static if(ByX){
	  //returning to vertical, so bottommost CW, topmost CCW
	  if(points[indices[0]].y == points[indices[1]].y){
		//ys are equal, CW should start at rightmost
		if(points[indices[0]].x < points[indices[1]].x){
		  return EdgePair(CWEdge(Pair(indices[1], indices[0])),
						  CCWEdge(Pair(indices[0], indices[1])));
		} else {
		  return EdgePair(CWEdge(Pair(indices[0], indices[1])),
						  CCWEdge(Pair(indices[1], indices[0])));
		}
	  
	  } else if(points[indices[0]].y < points[indices[1]].y){
		//the y's are different, so pick based on that
		return EdgePair(CWEdge(Pair(indices[0], indices[1])),
						CCWEdge(Pair(indices[1], indices[0])));
	  } else {
		return EdgePair(CWEdge(Pair(indices[1], indices[0])),
						CCWEdge(Pair(indices[0], indices[1])));
	  }
	} else {
	  //returning to horizontal, so return leftmost CCW, rightmost CW
	  if(points[indices[0]].x == points[indices[1]].x){
		//the x's are equal, so CW should be bottom, CCW should be top
		if(points[indices[0]].y < points[indices[1]].y){
		  return EdgePair(CWEdge(Pair(indices[0], indices[1])),
						  CCWEdge(Pair(indices[1], indices[0])));
		} else {
		  return EdgePair(CWEdge(Pair(indices[1], indices[0])),
						  CCWEdge(Pair(indices[0], indices[1])));
		}
	  }
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
	  //	  writeln("oriented correctly");
	  //positively oriented, add the triangle as is
	  triDB.addTriangle(Triangle(indices[0], indices[1], indices[2]));
	  //and the appropriate ghosts
	  triDB.addTriangle(Triangle(indices[0], indices[2], TriDB.makeGhost(indices[1])));
	  triDB.addTriangle(Triangle(indices[2], indices[1], TriDB.makeGhost(indices[0])));
	  triDB.addTriangle(Triangle(indices[1], indices[0], TriDB.makeGhost(indices[2])));

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
	  //	  writeln("oriented incorrectly");
	  //negatively oriented, swap 1 and 2
	  triDB.addTriangle(Triangle(indices[0], indices[2], indices[1]));
	  //ghosts
	  triDB.addTriangle(Triangle(indices[0], indices[1], TriDB.makeGhost(indices[2])));
	  triDB.addTriangle(Triangle(indices[1], indices[2], triDB.makeGhost(indices[0])));
	  triDB.addTriangle(Triangle(indices[2], indices[0], triDB.makeGhost(indices[1])));

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
	  //colinear, but possibly unsorted!
	  //
	  //	  writeln("collinear!");
	  import std.algorithm.sorting : sort;
	  static if(ByX){

		//sort by Y descending, breaking ties with x ascending
		indices.sort!(
					  delegate bool(int a, int b){
						return points[a].y == points[b].y ?
						  points[a].x < points[b].x : points[a].y > points[b].y;
					  });
	  } else {
		//sort by X ascending, breaking ties by Y descending
		indices.sort!(
					  delegate bool(int a, int b){
						return points[a].x == points[b].x ?
						  points[a].y > points[b].y : points[a].x < points[b].x;
					  });
	  }
	  
	  
	  triDB.addTriangle(Triangle(indices[0], indices[1], TriDB.makeGhost(indices[2])));
	  triDB.addTriangle(Triangle(indices[1], indices[0], TriDB.makeGhost(indices[1])));
	  triDB.addTriangle(Triangle(indices[1], indices[2], TriDB.makeGhost(indices[1])));
	  triDB.addTriangle(Triangle(indices[2], indices[1], TriDB.makeGhost(indices[0])));

	  //same for both since we sorted intelligently
	  return EdgePair(CWEdge(Pair(indices[2], indices[1])),
					  CCWEdge(Pair(indices[0], indices[1])));
	  
	}
	
  }
}



//ldi, inner clockwise pointing edge from left half
//rdi, inner ccw pointing edge from right half
//connecting the two hulls will be CW, either right to left for horizontal,
//or bottom to top for vertical
//if ldi/rdi are on the same vertical line (or parts of them are)
//remember the other coordinate must be greater for RDI
void advanceToLowerHullEdge(Vec)(const ref TriDB triDB, const ref Vec[] points,
							 ref CWEdge ldi, ref CCWEdge rdi){

  struct Segment{Vec first, second;}
  while(true){
	auto lo2d = orient2D(points[ldi.first], points[ldi.second], points[rdi.first]);
	if(lo2d > 0){ //rdi.first is to the left of edge ldi, advance ldi
	  ldi = hullAdvance(triDB, ldi);
	  continue;
	} else if(lo2d == 0 &&
			  inFrontOf(Segment(points[ldi.first], points[ldi.second]), points[rdi.first])){

	  ldi = hullAdvance(triDB, ldi);
	  continue;
	}

	//since rdi must be "above" ldi, if ldi's point is collinear with rdi,
	//we want to keep advancing it.  
	auto ro2d = orient2D(points[rdi.first], points[rdi.second], points[ldi.first]);
	if(ro2d < 0){ //ldi.first is right of or collinear with edge rdi
	  rdi = hullAdvance(triDB, rdi);
	  continue;
	} else if(ro2d == 0 &&
			  inFrontOf(Segment(points[rdi.first], points[rdi.second]), points[ldi.first])){

	  rdi = hullAdvance(triDB, rdi);
	  continue;
	}
	
	//no more possible advancements
	break;
  }
}
 
EdgePair zipHulls(Vec)(ref TriDB triDB, const ref Vec[] points, CWEdge ldi, CCWEdge rdi, bool byX){
  
  //  writefln("zipping, ByX: %s, ldi: %s, rdi: %s", byX, ldi, rdi);
  //  writeln("ldi: ", points[ldi.first], points[ldi.second]);
  //  writeln("rdi: ", points[rdi.first], points[rdi.second]);
  advanceToLowerHullEdge(triDB, points, ldi, rdi);
  
  CCWEdge rToL = CCWEdge(Pair(rdi.first, ldi.first)); //CCW when viewed from the top
  //it's CW if we're on the bottom
  CWEdge start = CWEdge(Pair(rToL.first, rToL.second)); //save this for later
  
  //  writeln("updated ldi: ", ldi, " rdi: " , rdi);
  //  writefln("rtol to start: %d, %d", rToL.first, rToL.second);

  auto preRdi = hullReverse(triDB, rdi).first;
  //writeln("preRdi: ", preRdi);
  //update the ghost triangle above rdi
  triDB.deleteTriangle(Triangle(rdi.second, rdi.first, TriDB.makeGhost(preRdi)));
  triDB.addTriangle(Triangle(rToL.second, rToL.first, TriDB.makeGhost(preRdi)));

  //update the ghost triangle on the ldi side
  auto preLdi = hullReverse(triDB, ldi).first;
  //  writeln("preLdi: ", preLdi);
  triDB.deleteTriangle(Triangle(preLdi, ldi.first, TriDB.makeGhost(ldi.second)));
  triDB.addTriangle(Triangle(rToL.first, rToL.second, TriDB.makeGhost(ldi.second)));

  //go from l to r to rdi.second for the closing edge of the hull
  triDB.addTriangle(Triangle(rdi.second, rToL.first, TriDB.makeGhost(rToL.second)));
  //and from r to l on top
  triDB.addTriangle(Triangle(preLdi, rToL.second, TriDB.makeGhost(rToL.first)));

  
  //are there ghost triangles that need to be removed?
  //yes, to start with
  while(true){
	//	writeln("while true loop.  RToL: ", rToL);
	//	writeln("computing next lCand");
	auto lCand = hullAdvance(triDB, rToL).second;

	//	writeln("lcand: ", lCand);
	bool lValid = rightOf(points[lCand], points[rToL.first], points[rToL.second]);

	if(lValid){
	  int nextCand = triDB.adjacent(rToL.second, lCand);
	  //	  writeln("left nextCand, ", nextCand);
	  while(!TriDB.isGhost(nextCand) &&
			isInCircle(points[rToL.second], points[rToL.first], points[lCand], points[nextCand])){
		//		writefln("erasing edge left side: %d, %d", rToL.second, lCand);
		eraseEdge(triDB, CCWEdge(Pair(rToL.second, lCand)), points);
		lCand = nextCand;
		nextCand = triDB.adjacent(rToL.second, lCand);
		//		writeln("left nextCand, ", nextCand);
	  }
	}
	//	writeln("computing next rCand");
	auto rCand = hullAdvance(triDB, reverse(rToL)).second;
	//	writeln("rcand: ", rCand);

	bool rValid = rightOf(points[rCand], points[rToL.first], points[rToL.second]);

	if(rValid){
	  int nextCand = triDB.adjacent(rCand, rToL.first);
	  //	  writeln("right nextCand, ", nextCand);
	  while(!TriDB.isGhost(nextCand) &&
			isInCircle(points[rToL.second], points[rToL.first], points[rCand], points[nextCand])){

		//		writefln("erasing edge right side: %d, %d", rCand, rToL.first);
		eraseEdge(triDB, CCWEdge(Pair(rCand, rToL.first)), points);
		rCand = nextCand;
		nextCand = triDB.adjacent(rCand, rToL.first);
		//		writeln("right nextCand, ", nextCand);
	  }
	}
	//	writefln("candidates: lcand: %d, rcand: %d", lCand, rCand);
	if(!lValid && !rValid){ break; }
	
	if(!lValid ||
	   (rValid && isInCircle(points[lCand], points[rToL.second],
							 points[rToL.first], points[rCand]))){

	  makeTriFromTwoEdges(triDB, rCand, rToL.first, rToL.second, points);
	  rToL.first = rCand;

	} else {

	  auto nextLeft = hullAdvance(triDB, rToL).second;
	  makeTriFromTwoEdges(triDB, rToL.first, rToL.second, nextLeft, points);
						
	  rToL.second = lCand;
	}

  }

  //add ghost for the bottom and top


  
  //triDB.addTriangle(Triangle(rToL.first, rToL.second, TriDB.makeGhost(lCand)));
  //  writeln("computing return edges");
  
  CCWEdge stop = CCWEdge(Pair(rToL.first, rToL.second));
  //  writeHulls(to!string("hull"~to!string(meshNumber)~".svg"), points, triDB);
  //  writeSVG(to!string("mesh"~to!string(meshNumber++)~".svg"), points, triDB);
  //writeln("zipped. Checking hulls");
  hullCheck(triDB, points, start);
  if(byX){
	//return to vertical, CW from bottom, CCW from top
	//believe the hull calls will be no-ops, but maybe not
	
	return EdgePair(hullMinYCW!Vec(triDB, points, start),
					hullMaxYCCW!Vec(triDB, points, stop));
	//}
  } else {
	//return to horizontal, CW from right, CCW from left
	return EdgePair(hullMaxXCW!Vec(triDB, points, reverse(stop)),
					hullMinXCCW!Vec(triDB, points, reverse(start)));
  }
}

int svgCount = 0;
EdgePair delaunayRecurse(Vec, bool ByX)(ref TriDB triDB, const  Vec[] points, int[] indices){

  assert(indices.length > 1);
  if(indices.length < 4){
	return delaunayBaseCase!(Vec,ByX)(triDB, points, indices);
  } else {
	//	writeln("recurse.  Indices:");
	//	writeln(indices);

	partitionPoints!(Vec,ByX)(points, indices);

	ulong middle = indices.length/2;

	auto ep1 = delaunayRecurse!(Vec, !ByX)(triDB, points, indices[0..middle]);
	//	writeSVG(to!string("mesh"~to!string(meshNumber++)~".svg"), points, triDB.getTriangles());
	auto ep2 = delaunayRecurse!(Vec, !ByX)(triDB, points, indices[middle..$]);
	//	writeSVG(to!string("mesh"~to!string(meshNumber++)~".svg"), points, triDB.getTriangles());

	//because ep2 will be below, if byY, and ep2 will be to the right if byX
	bool zipByX = ByX;
	CWEdge ldi;
	CCWEdge rdi;

	//writeln("ep1: ", ep1);
	//writeln("ep2: ", ep2);

	
	static if(ByX){
	  //if the x coordinates are equal, ep2's point's y coordinate must be greater thatn ep1's
	  if(points[ep1.ccwEdge.first].x == points[ep2.cwEdge.first].x){
		zipByX = false;
		writefln("swapped zipByX from %s  to %s", ByX, zipByX);
		ldi = ep2.cwEdge;
		rdi = ep1.ccwEdge;
		//all the points must be collinear!  Merge by Y!
	  } else {
		ldi = ep1.cwEdge;
		rdi = ep2.ccwEdge;
	  }
	} else {

	  //??? tricky and scary!!!
	  if(points[ep1.cwEdge.second].y == points[ep2.ccwEdge.first].y){
		//collinear, merge by X
		zipByX = true;
		writefln("swapped zipByX from %s to %s", ByX, zipByX);
		ldi = ep1.cwEdge;
		rdi = ep2.ccwEdge;
	  } else {
		ldi = ep2.cwEdge;
		rdi = ep1.ccwEdge;
	  }
	}


	//	writefln("getting lower hull edge, byX: %s", zipByX);
	//	writefln("ep1: %s", ep1);
	//	writefln("ep2: %s", ep2);
	//	writefln("ldi: %s,   rdi:  %s", ldi, rdi);


	//	writefln("about to zip hulls with %d total points", indices.length);
	auto retEdgePair = zipHulls!(Vec)(triDB, points, ldi, rdi, zipByX);

	//	writeHulls(to!string("hulls" ~ to!string(svgCount) ~ ".svg"), points, triDB); 
	//	writeSVG(to!string("step" ~ to!string(svgCount++) ~ ".svg"), points, triDB);

	//ohullCheck(triDB, points, retEdgePair.cwEdge);
	
	//since we are either correct or collinear, we should be good.
	//zip should make sure that things work correctly for the collinear case
	//	if(ByX == zipByX){
	  return retEdgePair;
	  //	} else {
	  //if ByX , then we want to give back the 
	  
	  //	}
	
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


void makeConstrainedDelaunay(Vec)(const Vec[] points, ref TriDB triDB, bool[Pair] segmentSet){
  foreach(seg; segmentSet.byKey){
	if(!triDB.edgeExists(seg.first, seg.second)){
	  addSegment(points, triDB, seg);
	}
  }
}

void cutOffScraps(Vec)(const Vec[] points, ref TriDB triDB, bool[Pair] segmentSet){

  
  auto outsideEdge = triDB.findOutsideEdge();
  Triangle[] toDelete = [ outsideEdge];

  while(toDelete.length > 0){
	auto tri = toDelete[$-1];
	toDelete.popBack();

	if(!triDB.containsTriangle(tri)){
	  continue;
	}
	foreach(i; 0..3){
	  auto v = tri[i];
	  auto w = tri[i+1];

	  
	  if((Pair(v, w) in segmentSet ) || (Pair(w, v) in segmentSet) ){
		continue;
	  } //don't cross real edges.  if either is a ghost, we DO want to delete those triangles

	  //3 cases depending on which if any of the 3 vertices is a ghost (at most 1 can be)
	  if(TriDB.isGhost(v)){
		//delete the ghost triangle in the CCW direction of us, if it exists
		
		auto x = tri[i+2];
		if(triDB.adjacentExists(w, TriDB.makeGhost(x)) ){
		  auto prev = triDB.adjacent(w, TriDB.makeGhost(x));
		  toDelete ~= Triangle(prev, w, triDB.makeGhost(x));
		}



	  } else if(TriDB.isGhost(w)){
		//delete the ghost triangle in the CW direction of us, if it exists
		if(triDB.adjacentExists(v, TriDB.unGhost(w)) ){
			toDelete ~= Triangle(v, TriDB.unGhost(w), triDB.adjacent(v, TriDB.unGhost(w)));
		  }
	  } else { //all reals, easy case
		if(triDB.adjacentExists(w, v)){
		  toDelete ~= Triangle(w, v, triDB.adjacent(w, v));
		}
	  }
	}
	triDB.deleteTriangle(tri);
	
  }
  
}



void addSegment(Vec)(const Vec[] points, ref TriDB triDB, Pair s){
  auto holes = clearCavity(points, triDB, s);
  fillCavity(points, triDB, holes.left);
  fillCavity(points, triDB, holes.right);
}


struct ClearCavityList{ int[] left; int[] right;}

//erase all triangles that cross the segment s
//return the points on the polygons on other side of the segment s
//once we add s to the triangulation, we'll need to retriangulate those holes
ClearCavityList clearCavity(Vec)(const Vec[] points, ref TriDB triDB, Pair s){

  writeln("clearing cavity for segment ", s);
  
  int[] leftPoints = [s.second];
  int[] rightPoints = [s.first];

  auto s1 = points[s.first];
  auto s2 = points[s.second];

  auto u = s.first;
  auto vw = triDB.adjacentRealTriangle(u);
  struct Segment{ Vec first, second; }

  //find the triangle at s that contains the segment
  while(!segmentsCross(Segment(s1, s2),
					   Segment(points[vw.first], points[vw.second]))){
	int adj = triDB.adjacent(u, vw.second);
	vw = Pair(vw.second, adj);
  }
  //and delete it
  writeln("about to delete first edge: ", u, vw);
  triDB.deleteTriangle(Triangle(u, vw.first, vw.second));

  if(leftOf(points[vw.first], s1, s2)){
	leftPoints ~= vw.first;
	rightPoints ~= vw.second;
  } else {
	leftPoints ~= vw.second;
	rightPoints ~= vw.first;
  }

  //it had better be real
  auto adj = triDB.adjacent(vw.second, vw.first);

  while(adj != s.second){
	int newAdj;
	//	writeln("in while loop, about to delete: ", vw.second, ", ", vw.first, ", ", adj);
	triDB.deleteTriangle(Triangle(vw.second, vw.first, adj));

	if(leftOf(points[adj], s1, s2)){
	  leftPoints ~= adj;
	} else {
	  rightPoints ~= adj;
	}
	if(segmentsCross(Segment(s1, s2), Segment(points[vw.first], points[adj]))){
	  newAdj = triDB.adjacent(adj, vw.first);
	  vw = Pair(vw.first, adj);
	} else {
	  newAdj = triDB.adjacent(vw.second, adj);
	  vw = Pair(adj, vw.second);
	}
	adj = newAdj;
  }
  //  writeln("deleting ", Triangle(adj, vw.second, vw.first));
  triDB.deleteTriangle(Triangle(adj, vw.second, vw.first));

  rightPoints ~= s.second;
  leftPoints ~= s.first;
  import std.algorithm : reverse;
  reverse(leftPoints[1..$-1]);
  return ClearCavityList(leftPoints, rightPoints);
}

void cavityInsertVertex(Vec)(const Vec[] points, ref TriDB triDB, int[] poly, int u, int v, int w){

  if(triDB.adjacentExists(w,v)){

	auto x = triDB.adjacent(w, v);
	if(orient2D(points[poly[u]], points[poly[v]], points[poly[w]]) > 0 &&
	   !inCircle(points[poly[u]], points[poly[v]], points[poly[w]], points[poly[x]])){
	  //uvw is constrained delaunay because the point on the other side is far enough away
	  writeln("adding ", Triangle(u, v, w));
	  triDB.addTriangle(Triangle(u, v, w));
	} else {
	  triDB.deleteTriangle(Triangle(w, v, x));
	  cavityInsertVertex(points, triDB, poly, u, v, x);
	  cavityInsertVertex(points, triDB, poly, u, x, w);
	}
  } else {
	//uvw is constrained delaunay, because there's nothing on the other side of vw
	writeln("adding ", Triangle(u, v, w));
	triDB.addTriangle(Triangle(u, v, w));
  }
  
}



//fill a cavity cleared by clearCavity.
void fillCavity(Vec)(const Vec[] points, ref TriDB triDB, int[] poly){
  import std.range;
  import std.array;
  import std.random;
  import std.algorithm;

  
  auto first = poly[0];
  auto last = poly[$-1];
  int[] perm = iota(1, to!int(poly.length - 1)).array;
  randomShuffle(perm);

  int[] prev = new int[poly.length];
  int[] next = new int[poly.length];
  foreach(i; 0..poly.length){
	prev[i] = to!int((poly.length + i -1) % poly.length);
	next[i] = to!int((i + 1) % poly.length);
  }


  auto o2ds = poly.map!(delegate(int i){
	  return orient2D(points[first], points[last], points[i]);
	}).array;

  foreach(i; iota(perm.length -1, -1, -1)){
	while(o2ds[perm[i]] < o2ds[prev[perm[i]]]
		  && o2ds[perm[i]] < o2ds[next[perm[i]]]){
	  auto j = uniform(0, -1);
	  swap(perm[i], perm[j]);
	}
	next[prev[perm[i]]] = next[perm[i]];
	prev[next[perm[i]]] = prev[perm[i]];
  }

  TriDB cavityTriangles = TriDB(poly.length);

  //  auto cavityPoints = poly.map!(i => Vec(points[i])).array;
  
  cavityTriangles.addTriangle(Triangle(0, perm[0], to!int(poly.length -1)));
  
  foreach(i; 1..perm.length){
	cavityInsertVertex(points, cavityTriangles, poly, perm[i], next[perm[i]], prev[perm[i]]);
  }

  
  
  triDB.addTriangulatedPolygon(cavityTriangles, poly);
}

bool isEncroached(Vec)(const ref Vec[] points, const ref TriDB triDB, Pair segment){
  struct Segment{ Vec first, second; }
  //is the point on either side of the segment inside the diametrical circle of this segment?
  if(triDB.adjacentExists(segment.first, segment.second)){
	if(pointInDiametricCircle(Segment(points[segment.first], points[segment.second]),
							  points[triDB.adjacent(segment.first, segment.second)])){
	  return true;
	}
  }
  if(triDB.adjacentExists(segment.second, segment.first)){
	if(pointInDiametricCircle(Segment(points[segment.second], points[segment.first]),
							  points[triDB.adjacent(segment.second, segment.first)])){
	  return true;
	}
  }
  return false;

  
}


//returns true if things changed
bool refinementStep(Vec, FP)(ref Vec[] points, ref TriDB triDB,  ref bool[Pair] segmentSet,
						 FP minSegmentLength, FP minAngleDegrees_, FP maxArea){

  import std.algorithm;
  import std.array;
  import std.container.binaryheap;
  
  auto encroachedSegments = segmentSet.byKey.filter!((seg) => isEncroached(points, triDB, seg)).array;

  auto encroachedTriangles = triDB.getTriangles.filter!((tri) =>
														minAngleDegrees(tri, points) < minAngleDegrees_ ||
														area(tri, points) > maxArea).array;

  if(encroachedSegments.empty && encroachedTriangles.empty){
	return false; //early out
  }
  
  auto segmentHeap = heapify!((a, b) =>
							  (points[a.first] - points[a.second]).magnitude <
							  (points[b.first] - points[b.second]).magnitude
							  )(encroachedSegments);

  //max heap, compare greater by min angle so the max heap pulls the smallest angle triangle first
  auto triangleHeap = heapify!( (a, b) =>
								minAngleDegrees(a, points) > minAngleDegrees(b, points)
								)(encroachedTriangles);
  
  bool modifiedMesh = false;

  writeln("encroached segs: ", segmentHeap.length, " enrcoached tris: ", triangleHeap.length);
  
  while(!segmentHeap.empty || !triangleHeap.empty){
	//prefer to split based on segments first
	if(!segmentHeap.empty){
	  auto s = segmentHeap.front;
	  segmentHeap.popFront;

	  //don't shrink too small
	  if( (points[s.first] - points[s.second]).magnitude < 2*minSegmentLength){
		continue;
	  }

	  auto newIndex = boyerWatsonSplitEdge(points, triDB, segmentSet, s);
	  modifiedMesh = true;

	  if(isEncroached(points, triDB, Pair(s.first, newIndex))){
		segmentHeap.insert(Pair(s.first, newIndex));
	  }

	  if(isEncroached(points, triDB, Pair(newIndex, s.second))){
		segmentHeap.insert(Pair(newIndex, s.second));
	  }

	  
	} else { //split a big or bad triangle
	  auto tri = triangleHeap.front;
	  triangleHeap.popFront;

	  if(triDB.containsTriangle(tri)){
		auto newIndex = boyerWatsonSplitTriangle(points, triDB, segmentSet,
												 segmentHeap, tri,
												 minSegmentLength, minAngleDegrees_, maxArea);
		if(!TriDB.isGhost(newIndex)){
		  modifiedMesh = true;
		}
		
	  }
	}
	
	
  }

  return modifiedMesh;
}


int boyerWatsonSplitEdge(Vec)(ref Vec[] points, ref TriDB triDB, ref bool[Pair] segmentSet, Pair s){
  
  import std.stdio;
  writeln("spltting ", s);
  int newIndex = to!int(points.length);
  Vec midpoint = 0.5*(points[s.first] + points[s.second]);
  points ~= midpoint;

  triDB.addPoint();

  
  segmentSet.remove(s);
  segmentSet[Pair(s.first, newIndex)] = true;
  segmentSet[Pair(newIndex, s.second)] = true;

  Triangle[] triangleStack;

  if(triDB.adjacentExists(s.first, s.second)){
	triangleStack ~= Triangle(s.first, s.second, triDB.adjacent(s.first, s.second));
  }

  if(triDB.adjacentExists(s.second, s.first)){
	triangleStack ~= Triangle(s.second, s.first, triDB.adjacent(s.second, s.first));
  }

  bool[Pair] cavitySegments;

  while(!triangleStack.empty()){

	Triangle tri = triangleStack.back;
	triangleStack.popBack();

	if(!triDB.adjacentExists(tri[0], tri[1])){
	  continue;
	}

	//if this triangle has midpoint in its circumcircle, remove it and
	//fill the hole eventually
	if(isInCircle(points[tri[0]], points[tri[1]], points[tri[2]], midpoint)){
	  //triangle is not constrained delaunay
	  foreach(i; 0..3){
		//we need the edges of the cavity.  Add them the first time they're seen
		//remove them the second (triangles on both sides deleted, so edge isn't on the border)
		auto e = Pair(tri[i], tri[i+1]);
		auto eReverse = Pair(e.second, e.first);
		//don't step across segments
		if( (e in segmentSet) || (eReverse in segmentSet) ){
		  cavitySegments[e] = true; //triangulate this edge
		  continue;
		}

		if(eReverse in cavitySegments){
		  cavitySegments.remove(eReverse);
		} else {
		  cavitySegments[e] = true;
		}

		if(triDB.adjacentExists(e.second, e.first)){
		  triangleStack ~= Triangle(e.second, e.first, triDB.adjacent(e.second, e.first));
		}
	  }
	  triDB.deleteTriangle(tri);
	}

  }

  cavitySegments.remove(s);
  cavitySegments.remove(Pair(s.second, s.first));

  foreach(const ref cs; cavitySegments.byKey){
	triDB.addTriangle(Triangle(cs.first, cs.second, newIndex));
  }
  return newIndex;
}


int boyerWatsonSplitTriangle(Vec, Heap, FP)(ref Vec[] points, ref TriDB triDB, ref bool[Pair] segmentSet,
											ref Heap segmentHeap, Triangle triToSplit,
											FP minSegmentLength, FP minAngleDegrees, FP maxArea){

  struct Segment{Vec first, second;}

  import std.stdio;
  writeln("splitting ", triToSplit);

  
  Vec circumcenter = getCircumcenter(triToSplit, points);
  if(!circumcenter.isFinite()){
	return TriDB.makeGhost(1);
  }

  bool pointInTriangle = false;

  bool[Triangle] toDelete;
  bool[Pair] cavityEdges;
  Pair[] encroachedSegments;


  //clear out the cavity
  Triangle[] triangleStack = [triToSplit];
  while(!triangleStack.empty()){

	Triangle tri = canonical(triangleStack.back);
	triangleStack.popBack();

	if(tri in toDelete){ //already deleted
	  continue;
	}

	//point makes triangle non-delaunay
	if(isInCircle(points[tri[0]], points[tri[1]], points[tri[2]], circumcenter)){

	  //have we found a triangle that contains the cirumcenter yet?
	  //idaelly it's the triToSplit itself
	  if(inside(tri, points, circumcenter)){
		pointInTriangle = true;
	  }

	  foreach(i; 0..3){

		auto u = tri[i];
		auto v = tri[i+1];
		auto e = Pair(u,v);
		auto eReverse = Pair(v, u);

		//is this edge a segment
		if( (e in segmentSet) || (eReverse in segmentSet) ) {
		  if(pointInDiametricCircle(Segment(points[e.first], points[e.second]), circumcenter)){
			//add this to encroached segments
			encroachedSegments ~= (e in segmentSet) ? e : eReverse;
		  }
		  
		} else {
		  //cross the segment
		  auto other = triDB.adjacent(v, u);
		  auto neighbor = canonical(Triangle(v, u, other));
		  if(! (neighbor in toDelete)){
			triangleStack ~= neighbor;
		  }
		}

		//add or remove this from cavity edges so we know what to fill
		//TODO, probably only one of these can happen)
		if( (e in cavityEdges) || (eReverse in cavityEdges)){
		  cavityEdges.remove(e);
		  cavityEdges.remove(eReverse);
		} else {
		  cavityEdges[e] = true;
		}
		
		
	  }
	  toDelete[tri] = true;
	}
  }

  if(!encroachedSegments.empty){
	//don't split the triangle, split the segments
	foreach(seg; encroachedSegments){
	  segmentHeap.insert(seg);
	}
	return TriDB.makeGhost(1);
  } else if(!pointInTriangle){
	//	writeln("split triangle, no encroached segments, but circumcenter outside searched triangles");
	return TriDB.makeGhost(1);
  } else {
	foreach(t; toDelete.byKey()){
	  triDB.deleteTriangle(t);
	}

	int newIndex = to!int(points.length);
	points ~= circumcenter;
	triDB.addPoint();
	foreach(seg; cavityEdges.byKey()){
	  triDB.addTriangle(Triangle(seg.first, seg.second, newIndex));
	}
	return newIndex;
  }

}




//write the header and return points scaled appropriately
Vec[] prepareSVG(Vec)(ref File f, const Vec[] points, const int[] activePoints){
  import std.array;
  import std.stdio : File;
  import std.algorithm;
  
  //activePoints.sort().uniq();
  if(activePoints.empty){
	f.writeln("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"800\" height=\"800\" ></svg>");
	return [];
  }
  //using a delecate to store a copy of the points, not a reference
  Vec[] usedPoints = activePoints.map!(delegate Vec(int a){ return points[a]; }).array;
  
  
  Vec minP = Vec(usedPoints.map!("a.x").minElement(), usedPoints.map!("a.y").minElement());
  Vec maxP = Vec(usedPoints.map!("a.x").maxElement(), usedPoints.map!("a.y").maxElement());

  Vec center = 0.5*(minP + maxP);
  
  Vec size = maxP - minP;
  minP -= .03*size;
  maxP += .03*size;

  size = maxP - minP;
  auto maxDim = max(size.x, size.y);

  auto radius = max(size.x, size.y)/200.0f;

  auto aspectRatio = size.y/size.x;

  auto width = 1600;
  auto height = 1600;//*aspectRatio;


  Vec[] scaledPoints = points.map!( a => (a - center)*(10/maxDim))
	//	.map!(a => Vec(a.x*2.5, a.y))
	.array;
  
  f.write("<svg xmlns=\"http://www.w3.org/2000/svg\"  ");
  f.writefln("width=\"%d\" height=\"%d\" viewBox=\"-5 -5 10 10\" >", width, height);
  
  f.writeln("<g transform=\"scale(1, -1)\" >");

  return scaledPoints;
  
}




void writeSVG(Vec)(string filename, const Vec[] points, const ref TriDB triDB){
  import std.array;
  import std.stdio : File;
  import std.algorithm;

  writefln("dumping %s", filename);
  File f = File(filename, "w");
  auto activePoints = triDB.getActiveVertices();
  auto scaledPoints = prepareSVG(f, points, activePoints);
  Triangle[] tris = triDB.getTriangles();
  
  foreach(const ref tri ; tris){
	if(!TriDB.isGhost(tri[0])  && !TriDB.isGhost(tri[1]) && !TriDB.isGhost(tri[2])){
	  foreach(i; [0,1,2]){
		f.writefln("<line x1=\"%.8f\" y1=\"%.8f\" x2=\"%.8f\" y2=\"%.8f\" " ~
				   "stroke=\"black\" stroke-width=\"%.8f\" />",
				   scaledPoints[tri[i]].x, scaledPoints[tri[i]].y,
				   scaledPoints[tri[i+1]].x, scaledPoints[tri[i+1]].y, .0002);
	  }
	} else {
	  auto p1 = TriDB.isGhost(tri[0]) ? tri[1] : tri[0];
	  auto p2 = TriDB.isGhost(tri[2]) ? tri[1] : tri[2];
	  f.writefln("<line x1=\"%.8f\" y1=\"%.8f\" x2=\"%.8f\" y2=\"%.8f\" " ~
				 "stroke=\"black\" stroke-width=\"%.8f\" stroke-dasharray=\"2%%, 10%%\"/>",
				 scaledPoints[p1].x, scaledPoints[p1].y, scaledPoints[p2].x, scaledPoints[p2].y, .002);
	  
	}
  }

  foreach(i ; activePoints){
	const auto p = scaledPoints[i];
	f.writefln( "<circle cx=\"%.8f\" cy=\"%.8f\" r=\"%.8f\" fill=\"black\" />", p.x, p.y, .01);
	f.writefln( "<g transform=\"translate(%.8f, %.8f) scale(1, -1)\" >" ~
				"<text x=\"0\" y=\"0\" font-family=\"Verdana\" font-size=\"%.8f\" fill=\"red\" >%d</text></g>",
				p.x, p.y, .02, i);

  }
  f.writeln("</g>\n</svg>");
}

void writeHulls(Vec)(string filename, const Vec[] points, const ref TriDB triDB){
  enum colors = ["red", "green", "blue", "purple", "orange", "brown", "pink", "gold", "crimson"];

  writefln("dumping %s", filename);
  File f = File(filename, "w");
  auto activePoints = triDB.getActiveVertices();
  auto scaledPoints = prepareSVG(f, points, activePoints);

  bool[int] neededPoints;
  foreach(i ; activePoints){ neededPoints[i] = true;}

  int color = 0;
  foreach(i; activePoints){
	if(!neededPoints[i]){ continue; }
	int cwOfI = i;
	foreach(vw; triDB.getTriangles(i)){
	  if(!TriDB.isGhost(vw.first) && TriDB.isGhost(vw.second)){
		cwOfI = vw.first;
		break;
	  }
	}
	//not on the hull, skip this
	if(cwOfI == i){ continue; }
	CWEdge e = CWEdge(Pair(i, cwOfI));
	CWEdge start = e;
	do{
	  f.writefln("<line x1=\"%.8f\" y1=\"%.8f\" x2=\"%.8f\" y2=\"%.8f\" stroke=\"%s\" stroke-width=\"%.8f\" />",
			   scaledPoints[e.first].x, scaledPoints[e.first].y,
			   scaledPoints[e.second].x, scaledPoints[e.second].y,
			   colors[color % colors.length], .001);
	  neededPoints[e.first] = false;
	  e = hullAdvance(triDB, e);

	  auto p = scaledPoints[e.first];
	  f.writefln( "<circle cx=\"%.8f\" cy=\"%.8f\" r=\"%.8f\" fill=\"black\" />", p.x, p.y, .002);
	  f.writefln( "<g transform=\"translate(%.8f, %.8f) scale(1, -1)\" >" ~
				  "<text x=\"0\" y=\"0\" font-family=\"Verdana\" font-size=\"%.8f\" fill=\"%s\" >%d</text></g>",
				  p.x, p.y, .015, colors[color % colors.length], e.first);
	  
	  
	  
	}while(e != start);
	++color;
  }

  /*foreach(i ; activePoints){
	const auto p = scaledPoints[i];
	f.writefln( "<circle cx=\"%.8f\" cy=\"%.8f\" r=\"%.8f\" fill=\"black\" />", p.x, p.y, .002);
	f.writefln( "<g transform=\"translate(%.8f, %.8f) scale(1, -1)\" >" ~
				"<text x=\"0\" y=\"0\" font-family=\"Verdana\" font-size=\"%.8f\" fill=\"red\" >%d</text></g>",
				p.x, p.y, .015, i);
	
				}*/

  
  writefln("hulls this frame: %d" , color);
  f.writeln("</g>\n</svg>");
}






unittest{

  writeln("\n\n\n2 points test");
  
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

  writeln("\n\n\n\n4 points test");
  vec2[] points  = [vec2(0,0), vec2(1,1), vec2(0.2, 0.9), vec2(0.9, 1)];
  auto triDB = delaunayTriangulate(points);
  auto tris = triDB.getTriangles();

  triDB.dump();
  foreach(const ref tri; tris){
	writeln(tri);
  }


  writeSVG("test4.svg", points, triDB);
  
  
  assert(tris.length == 6);

  
}


/*unittest{
  //horizontal/vertical line test

  foreach(np ; 2..20){
	writefln("\n\n\nhorizontal line test: %d", np);
	vec2[] points = new vec2[np];
	foreach(i ; 0..points.length){
	  points[i] = vec2(i, 0);
	}

	auto triDB = delaunayTriangulate(points);
	writeln("HLTest finished.  TriDB Contents:");
	triDB.dump();
  }

  foreach(np ; 2..20){
	writefln("\n\n\nvertical line test: %d", np);
	vec2[] points = new vec2[np];
	foreach(i ; 0..points.length){
	  points[i] = vec2(0, i);
	}

	auto triDB = delaunayTriangulate(points);
	writeln("VLTest finished.  TriDB Contents:");
	triDB.dump();
  }


  //combined horizontal/vertical:

  foreach(np ; 2..10){
	writefln("\n\n\npathalogical lines test: %d", np);
	vec2[] points;
	foreach(i ; 1..np){
	  points  ~= vec2(0, i);
	  points  ~= vec2(i, 0);
	}
	writeln(points);
	auto triDB = delaunayTriangulate(points);
	writeln("HLTest finished.  TriDB Contents:");
	triDB.dump();
	writeSVG(to!string("pathalogical" ~ to!string(np) ~ ".svg"), points, triDB);
  }

  
}
*/
  

unittest {
  writeln("clearCavityTest");
  //vec2[] points = [ vec2(0, 0), vec2(1, 1), vec2(.4, .5), vec2(.4, .3), vec2(.6, .7), vec2(.6, .5) ];
  vec2[] points = [vec2(-114.052963256836, 37.592784881592), vec2(-114.052474975586, 37.604774475098), vec2(-114.052474975586, 37.604778289795), vec2(-114.051727294922, 37.745998382568), vec2(-114.051788330078, 37.746250152588), vec2(-114.051666259766, 37.746959686279)];
  foreach( ref p; points){
	writefln("%0.12f, %0.12f", p.x, p.y);
  }
  auto triDB = delaunayTriangulate(points);
  writeln("triangulated");
  triDB.dump();

  foreach(i; 0..points.length){
	if(!triDB.edgeExists(cast(int)i, cast(int)((i+1)%points.length))){
	  writeln("clearing cavity: ", i, ' ', (i+1) % points.length);
	  clearCavity(points, triDB, Pair(cast(int)i, cast(int)((i+1)%points.length)));
	  triDB.dump();
	}
  }

  
}
