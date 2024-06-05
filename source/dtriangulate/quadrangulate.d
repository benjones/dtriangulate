module dtriangulate.quadrangulate;

import dtriangulate.edgecomplex;
import dtriangulate.predicates;

import std.stdio;
import std.typecons : Typedef;
import std.algorithm;
import std.range;

alias CCWEdge = Typedef!(HalfEdge, HalfEdge.init, "ccw");
alias CWEdge = Typedef!(HalfEdge, HalfEdge.init, "cw");

//The CW edge will start at the bottom, and the CCW will start at the top
//or the CW edge will be the right-most, and the CCW edge will be the left-most
struct EdgePair{ CCWEdge ccwEdge; CWEdge cwEdge;}


void dumpSvg(Vec)(const Vec[] points, const ref EdgeComplex ec){
  import std.format;
  static int frameNumber = 0;
  writeln("Dumping frame ", frameNumber);
  writeSVG(format!"quads_%s.svg"(frameNumber++), points, ec);
}

EdgeComplex quadrangulate(Vec)(const Vec[] points){

  import std.array;
  import std.range : iota;
  import std.algorithm.sorting;
  import std.algorithm.comparison;

  EdgeComplex ec;

  if(points.length < 2){
    return ec;
  }

  size_t[] indices = iota(points.length).array;

  indices.sort!((a, b) => cmp(points[a].vector[], points[b].vector[]) < 0 ); //sort points lexicographically

  foreach(i; indices){
    writeln(points[i]);
  }

  quadrangulateRecurse(ec, points, indices);

  return ec;
}



EdgePair quadrangulateRecurse(Vec)(ref EdgeComplex ec, const Vec[] points, const size_t[] indices){
  writeln("quadrangulating these indices: ", indices);

  if(indices.length <= 3){
    return quadrangulateBaseCase(ec, points, indices);
  } else {
    const midpoint = indices.length/2;
    auto leftEdgePair = quadrangulateRecurse(ec, points, indices[0..midpoint]);
    auto rightEdgePair = quadrangulateRecurse(ec, points, indices[midpoint..$]);

    auto topEdgeCcwEdge = zipHulls(ec, points, leftEdgePair.cwEdge, rightEdgePair.ccwEdge);


    //trace the top edge that's returned to the left/right since ldo and rdo might have changed during zipping
    return EdgePair();
  }
}


///ldi: CW edge whose origin is the rightmost point in the left hull
///rdi: CCW edge whose org is the leftmost point in the right hull
CCWEdge zipHulls(Vec)(ref EdgeComplex ec, const Vec[] points, CWEdge ldi, CCWEdge rdi){

  writeln("\nzipping, ldi: ", ldi, " rdi: ", rdi, "\n\n");
  //start by finding the lowest edge in the combined hull
  advanceToLowerHullEdge(ec, points, ldi, rdi);

  writeln("bottom edge of zip, ldi: ", ldi, " rdi: ", rdi);

  const ldiIndex = ec.getIndex(cast(HalfEdge)ldi);
  const ldiPrevIndex = ec.prevL(ldiIndex);
  writeln("ldiIndex, ldiPrevIndex, ldiPrev", ldiIndex, ldiPrevIndex, ec[ldiPrevIndex]);

  //for rdi, we want the "outside" edge, which is the CW edge
  const rightIndex = ec[rdi.sym].nextL;
  writeln("used rigth edge: ", ec.edgeToString(rightIndex));

  //this is the TOP side of the edge.
  //following nextL will go up the being-zipped side of the right hull
  auto bottomLToRIndex = ec.newEdgeBetween(ldiPrevIndex, rightIndex, false);
  writeln("bottom edge index: ", bottomLToRIndex);

  auto ilegdIndex = ec.prevL(bottomLToRIndex);
  auto ireguIndex = ec[bottomLToRIndex].nextL;

  writeln("dumping before entering zip loop");
  dumpSvg(points, ec);
  while(true){

    auto bottomEdge = ec[bottomLToRIndex];
    writeln("bottom edge: ", ec.edgeToString(bottomEdge));

    auto innerLeftEdgeGoingDown = ec[ilegdIndex];
    writeln("ilegd: ", ec.edgeToString(innerLeftEdgeGoingDown));

    auto innerRightEdgeGoingUp = ec[ireguIndex];
    writeln("iregu ", ec.edgeToString(innerRightEdgeGoingUp));




    const leftEdgeOK = canNeighborQuad(ec, points, innerLeftEdgeGoingDown);
    const rightEdgeOK = canNeighborQuad(ec, points, innerRightEdgeGoingUp);
    writeln("left can border quad? ", leftEdgeOK);
    writeln("right can border quad? ", rightEdgeOK);

    if(leftEdgeOK && rightEdgeOK){
      //if the quad is convex, let's use it
      if(isQuadConvex(points, [innerLeftEdgeGoingDown.org,
                                            innerLeftEdgeGoingDown.dest,
                                            innerRightEdgeGoingUp.org,
                               innerRightEdgeGoingUp.dest])){
        writeln("rdi/ldi sides can be joined to make a convex quad, checked 4 corners");

        const newIlegdIndex = ec.prevL(ilegdIndex);
        const newIreguIndex = innerRightEdgeGoingUp.nextL;

        //returned edge goes R to L
        bottomLToRIndex = ec[ec.newEdgeBetween(ireguIndex, ilegdIndex, true)].sym;
        ilegdIndex = newIlegdIndex;
        ireguIndex = newIreguIndex;

        dumpSvg(points, ec);
        continue;
      }

      writeln("could not make a convex quad from the 2 side edges");

    }

    //check to see if this could be a triangle (meaning this is the "top" of the zip)
    //try using iregu.dest as the top of the triangle
    const rightTriOrientation = orient2D(points[innerLeftEdgeGoingDown.dest],
                                         points[innerRightEdgeGoingUp.org],
                                         points[innerRightEdgeGoingUp.dest]);
    if(rightTriOrientation > 0){
      writeln("using top of iregu would make a valid triangle, let's see if it would be on the hull");

      //TODO: double check
      //does this new edge make a convex hull segment with the edges around it?
      const nextRightIndex = ec[innerRightEdgeGoingUp.nextL].dest;

      const hullOrientation = orient2D(points[nextRightIndex],
                                       points[innerRightEdgeGoingUp.dest],
                                       points[innerLeftEdgeGoingDown.org]);


      if(hullOrientation > 0){
        writeln("OK to use this triangle to close with");

        const newEdge = ec.newEdgeBetween(ireguIndex, bottomLToRIndex, true);
        dumpSvg(points, ec);
        return CCWEdge(ec[newEdge]);

      }

    }


    const leftTriOrientation = orient2D(points[innerLeftEdgeGoingDown.org],
                                        points[innerLeftEdgeGoingDown.dest],
                                        points[innerRightEdgeGoingUp.org]);
    if(leftTriOrientation > 0){
      writeln("using top of ilegd would make a valid triangle, let's see if it would be on the hull (TODO)");

    }


    writeln("stealing from left/right side");
    //check left first

    if(isInTriangle(ec, points,ec[innerLeftEdgeGoingDown.sym])){
      writeln("left side is a triangle");
      //let's see if we could make a convex quad using its vertices plus iregu.org
      if(isQuadConvex(points, [innerLeftEdgeGoingDown.org,
                           ec[ec[innerLeftEdgeGoingDown.sym].nextL].dest,
                           innerLeftEdgeGoingDown.dest,
                           innerRightEdgeGoingUp.org])){
        writeln("can make a convex quad from the left triangle.  Yay!");

        const triEdgeIndex = ec[innerLeftEdgeGoingDown.sym].nextL;
        writeln("edge from triangle to use: ", ec[triEdgeIndex]);

        ec.deleteEdge(ilegdIndex);
        writeln("triangle edge deleted, new triangle coming up next");
        writeln(ec);
        dumpSvg(points, ec);
        const newEdge = ec.newEdgeBetween(bottomLToRIndex, triEdgeIndex, true);
        dumpSvg(points, ec);
      }


    }


    const size_t[] toCheck = [15, 19, 23, 16];
    writeln("convex check: ", isQuadConvex(points, toCheck));

    assert(0);


    //Can detect that some edges CAN't be in the quadrangulation
    // (on the interior and can't find a point to make a convex quad using the neighboring edges


    // Seems impossible to guarantee that any triangles must have at least one EDGE on the boundary
    // What if we just require at least one vertex on the boundary?
    // Need to mark that vertex so that we don't zip a triangle in during a future pass
    // we'll need to make sure we delete anything outside of it and make it a quad




    /*if(lPolygonSize == 3){
      writeln("lcand part of a triangle, deleting it");
      ec.deleteEdge(ilegdIndex);
      writeln(ec);
      dumpSvg(points, ec);
    }
    if(rPolygonSize == 3){
      writeln("rcand part of a triangle, deleting it");
      ec.deleteEdge(ireguIndex);
      writeln(ec);
      dumpSvg(points, ec);
      }*/


  }
  return CCWEdge();
}


/// Is the edge either part of a quad, or
/// or does it have at least one edge on the boundary?
/// If so, then it's fine to make a new quad using this edge
bool canNeighborQuad(T)(const ref EdgeComplex ec, const T[] points, HalfEdge edge){

  const sym = ec[edge.sym];
  //boundary edges can be part of quads, no problem
  if(sym.isOutside){
    writeln("sym is outside");
    return true;
  }
  //OK, it must be part of a polygon.  If it's a quad, then we're good
  const size = ec.polygonSize(points, sym);
  if(size == 4){
    return true;
  }

  //todo, check if sym is part of a triangle on the edge
  assert(size == 3);

  return ec.triangleOnBoundary(sym);
}


/// is this edge in a triangle?

bool isInTriangle(T)(const ref EdgeComplex ec,T[] points, HalfEdge edge){
  if(edge.isOutside){
    return false;
  }
  return ec.polygonSize(points, edge) == 3;
}


bool isQuadConvex(T)(const ref T[] points, scope const size_t[] vertices){
  writeln("is quad convex? ", vertices);
  foreach( i; 0..4){
    const o2d = orient2D(points[vertices[i]], points[vertices[(i + 1) %4]], points[vertices[(i +2) %4]]);
    if(o2d <= 0){
      writeln("this corner is bad: ", vertices[i], " ", vertices[(i + 1) %4], " ", vertices[(i + 2) %4],
              " o2d: ", o2d);
      return false;
    }
  }
  return true;
}



auto hullAdvance(T)(const ref EdgeComplex ec, T edge)
     if(is(T == CWEdge) || is(T == CCWEdge)){
       return T(ec[edge.nextL]);
     }

// ldi, inner clockwise pointing edge from left half
// rdi, inner ccw pointing edge from right half
// connecting the two hulls will be CW, right to left
// if ldi/rdi are on the same vertical line (or parts of them are)
// remember the other coordinate must be greater for RDI
// IOW, since the points were sorted lexicographically, ldi's vertices
// must be lexicographically before rdi's
void advanceToLowerHullEdge(Vec)(const ref EdgeComplex ec, const ref Vec[] points,
                                 ref CWEdge ldi, ref CCWEdge rdi){

  struct Segment{Vec first, second;}
  while(true){
	auto lo2d = orient2D(points[ldi.org], points[ldi.dest], points[rdi.org]);
	if(lo2d > 0){ //rdi.first is to the left of edge ldi, advance ldi
	  ldi = hullAdvance(ec, ldi);
	  continue;
	} else if(lo2d == 0 &&
			  inFrontOf(Segment(points[ldi.org], points[ldi.dest]), points[rdi.org])){

	  ldi = hullAdvance(ec, ldi);
	  continue;
	}

	//since rdi must be "above" ldi, if ldi's point is collinear with rdi,
	//we want to keep advancing it.
	auto ro2d = orient2D(points[rdi.org], points[rdi.dest], points[ldi.org]);
	if(ro2d < 0){ //ldi.first is right of or collinear with edge rdi
	  rdi = hullAdvance(ec, rdi);
	  continue;
	} else if(ro2d == 0 &&
			  inFrontOf(Segment(points[rdi.org], points[rdi.dest]), points[ldi.org])){

	  rdi = hullAdvance(ec, rdi);
	  continue;
	}

	//no more possible advancements
	break;
  }
}


EdgePair quadrangulateBaseCase(Vec)(ref EdgeComplex ec, const Vec[] points, const size_t[] indices){

  scope(exit){
    //writeln("ec after triangulating", indices);
    //writeln(ec);
  }
  writeln("base case for ", indices);

  if(indices.length == 2){
    assert(points[indices[0]].x <= points[indices[1]].x);
    auto newEdge = ec[ec.newEdge(indices[0], indices[1])];
    return EdgePair(CCWEdge(newEdge), CWEdge(ec[newEdge.sym]));
  } else {
    assert(indices.length == 3);

    const orientResult = orient2D(points[indices[0]], points[indices[1]], points[indices[2]]);
    if(orientResult < 0){
      //make triangle a, c, b
      auto edgeAC = ec[ec.newTriangle(indices[0], indices[2], indices[1])];

      //AC is the CCW edge (a is leftmost, and AC goes CCW
      //CW edge is either ca, because c is furthest right.  B must be above segment AC, between them
      return EdgePair(CCWEdge(edgeAC), CWEdge(ec[edgeAC.sym]));

    } else if(orientResult > 0){
      //triangle a, b, c
      auto edgeAB = ec[ec.newTriangle(indices[0], indices[1], indices[2])];
      //AB is the CCW edge, CB is the CW edge
      auto edgeCB = ec[ec[edgeAB.nextL].sym];
      return EdgePair(CCWEdge(edgeAB), CWEdge(edgeCB));
    } else {
      //two line segments since the points are collinear
      auto edgeAB = ec.newEdge(indices[0], indices[1]);
      auto edgeBC = ec.newDanglingEdge(edgeAB, indices[2]);
      //AB leftmost CCW, CB rightmost CW edge
      return EdgePair(CCWEdge(ec[edgeAB]), CWEdge(ec[ec[edgeBC].sym]));
    }

  }
  // handle this in zipping, no need to add it as a base case
  /*else {
    assert(indices.length == 4);



    //can we make a convex quad?
    //start by finding the leftmost CCW edge in the convex hull of these for points (h1, h2)
    //and the leftmost CW edge (h1, h4)
    const h1 = indices[0];
    size_t h2 = indices[1];
    size_t h4 = indices[1];
    foreach(i; [indices[2], indices[3]]){
      auto orientationDown = orient2D(points[h1], points[h2], points[i]);
      if(orientationDown == 0){
        assert(false, "Colinear points in quad base case");
        //h1, h2, i are collinear, no way we can get a convex quad
      } else if(orientationDown < 0){
        //h2 isn't on this hull edge, replace it with i
        h2 = i;
      }

      //if this is positive for the other points, then they are the better choice for h4
      auto orientationUp = orient2D(points[h1], points[h4], points[i]);
      if(orientationUp == 0){
        assert(false, "Colinear points in quad base case");
      } else if(orientationUp > 0){
        h4 = i;
      }
    }
    //h4, h1, h2 should be on the convex hull.  If h2, h3, h4 is a left turn, then we've got a convex quad.  Otherwise it's concave
    auto h3 = indices.filter!(x => x != h1 && x != h2 && x != h4).front;
    writeln("hull points: h1: ", h1, " h2: ", h2, "h3: ", h3, "h4: ", h4);

    auto closingOrient = orient2D(points[h2], points[h3], points[h4]);
    if(closingOrient == 0){
      assert(false, "Colinear points in quad base case");
    } else if(closingOrient > 0){
      //it's convex, we're all good
      auto ccwEdgeIndex = ec.newQuad(h1, h2, h3, h4);

      auto tempEdge = ec[ccwEdgeIndex];
      while(tempEdge.dest != indices[3]){ //right most, and upper most since we did lexicographic ordering
        tempEdge = ec[tempEdge.nextL]; //go around the lower part of the hull
      }
      return EdgePair(CCWEdge(ec[ccwEdgeIndex]), CWEdge(ec[tempEdge.sym]));
    }
    assert(0, "unreachable?");

    }*/
}

unittest {

  import gl3n.linalg : vec2;

  //test collinear sets

  auto points = [vec2(1.0f + float.epsilon, 1.0f + float.epsilon),
                 vec2(1.0f + float.epsilon, 2.0f + float.epsilon),
                 vec2(1.0f + float.epsilon, 3.0f + float.epsilon),
                 vec2(1.0f + float.epsilon, 4.0f + float.epsilon)];

  foreach(i; 2..points.length + 1){
    writeln("testing ", i , " collinear points");
    auto ec = quadrangulate(points[0..i]);

    writeln("ec for first ", i, " points: ", ec);
    auto polygons = ec.polygons;
    writeln("polygons: ", polygons);
    assert(polygons.length == 1);
    // 2 points: 2 segments
    // 3 points: 4 segments
    // 4 points: 6 segments
    assert(polygons[0].length == 2*(i -1));
  }

}
