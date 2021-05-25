module dtriangulate.tridb;

public import dtriangulate.ssoVector;
public import dtriangulate.predicates;

import std.stdio;
import std.conv;
import std.algorithm;

struct Pair{
  int first, second;
  bool opEquals(Pair rhs) const{
	return first == rhs.first && second == rhs.second;
  }

  size_t toHash() const @safe nothrow{
	size_t ret = cast(size_t)(first);
	ret <<= 32;
	size_t ss = cast(size_t)(second);
	ret |= ss;
	return ret;
  }
}

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


  ref int opIndex(int i) return{
	return v[i > 2 ? i -3 : i]; //faster than modulo?
  }

  int opIndex(int i) const{
	return v[i > 2 ? i -3 : i]; //faster than modulo?
  }

  bool opEquals(Triangle rhs) const{
	return v[0] == rhs.v[0] && v[1] == rhs.v[1] && v[2] == rhs.v[2];
  }

  //todo: do better?
  size_t toHash() const @safe nothrow{
	size_t ret = cast(size_t)(v[0]);
	ret <<= 32;
	size_t b = cast(size_t)(v[1]);
	ret |= b;
	size_t c = cast(size_t)(v[2]);
	c <<= 16;
	ret |= c;
	return ret;
  }

  string toString() const{
	string ret =  "Triangle( ";
	foreach(i ; v){
	  if(TriDB.isGhost(i)){
		ret ~= "GHOST " ~ to!string(TriDB.unGhost(i)) ~ ", ";
	  } else {
		ret ~= to!string(i) ~ ", ";
	  }
	}
	return to!string(ret ~ " )");
  }

  private int[3] v;
}

//rotate the indices so t[0] is the smallest
Triangle canonical(Triangle t){
  foreach(i; 0..2){ //rotate at most twice
	if(t[1] < t[0] || t[2] < t[0]){
	  t = Triangle(t[1], t[2], t[0]);
	}
  }
  return t;
}

unittest {
  auto t1 = Triangle(0, 1, 2);
  auto t2 = Triangle(2, 0, 1);
  auto t3 = Triangle(1, 2, 0);

  assert(canonical(t1) == t1);
  assert(canonical(t2) == t1);
  assert(canonical(t3) == t1);


  auto t4 = Triangle(0, 2, 1);
  auto t5 = Triangle(2, 1, 0);
  auto t6 = Triangle(1, 0, 2);

  assert(canonical(t4) == t4);
  assert(canonical(t5) == t4);
  assert(canonical(t6) == t4);

}

auto minAngleDegrees(Vec)(Triangle tri, const ref Vec[] points){
  import std.traits : Unqual;
  import std.math: acos, PI;
  import std.algorithm: clamp;

  alias FP = Unqual!(typeof(points[0].x));
  FP ret = 180;
  foreach(i; 0..3){
	auto a = points[tri[i]];
	auto b = points[tri[i+1]];
	auto c = points[tri[i+2]];

	auto ab = b - a;
	auto ac = c - a;

	//times is dot product for gl3n, hopefully other libs?
	auto angle = acos(clamp(ab*ac/(ab.magnitude*ac.magnitude) , 0, 1));
	ret = min(ret, angle);
  }
  return ret*180/PI;

}

auto area(Vec)(Triangle tri, const ref Vec[] points){
  import std.traits : Unqual;
  alias FP = Unqual!(typeof(points[0].x));
  auto a = points[tri[0]];
  auto b = points[tri[1]];
  auto c = points[tri[2]];
  return FP(0.5)*( (b.x-a.x)*(c.y-a.y)- (c.x-a.x)*(b.y -a.y));

}


//todo does this need to be adaptive precision?
Vec getCircumcenter(Vec)(Triangle tri, const ref Vec[] points){

  auto a = points[tri[0]];
  auto b = points[tri[1]];
  auto c = points[tri[2]];

  const auto d = 2*(a.x*(b.y - c.y) + b.x*(c.y - a.y) + c.x*(a.y - b.y));

  return Vec (
			  (a.magnitude_squared()*(b.y - c.y) +
			   b.magnitude_squared()*(c.y - a.y) +
			   c.magnitude_squared()*(a.y - b.y) )/d,
			  (a.magnitude_squared()*(c.x - b.x) +
			   b.magnitude_squared()*(a.x - c.x) +
			   c.magnitude_squared*(b.x - a.x))/d
			  );
}

bool inside(Vec)( Triangle t, const ref Vec[] points, Vec p){
  foreach(i ; 0..3){
	if(orient2D(points[t[i]], points[t[i+1]], p) < 0){
	  return false;
	}
  }
  return true;
}



struct TriDB{

  this(ulong numPoints){
	triangles = new SSOVector!(Pair, 9)[numPoints];
  }

  //ghosts will have their MSB set to 1

  public static bool isGhost(int i){
	return i < 0; //is MSB set?
  }

  public static int unGhost(int i){
	assert(isGhost(i));
	enum int mask = 1 << (8*int.sizeof -1);
	return i ^ mask;
  }

  public static int makeGhost(int i){
	assert(i >= 0);
	enum int mask = 1 << (8*int.sizeof -1);
	return i | mask;
  }



  void addTriangle(Triangle t){
	foreach(i; 0..3){
	  if(!isGhost(t[i])){ //don't add ghost edges
		triangles[t[i]] ~= Pair(t[i+1], t[i+2]);
	  }
	}
  }

  void deleteTriangle(Triangle t){
	int deletedCount = 0;
	int realCount = 0;
	foreach(i; 0..3){
	  int u = t[i];
	  if(!isGhost(u)){
		++realCount;
		auto vw = Pair(t[i+1], t[i+2]);
		foreach(j, pr ; triangles[u]){
		  if(pr == vw){
			triangles[u][j] = triangles[u][$-1];
			triangles[u].popBack();
			++deletedCount;
			break;
		  }
		}
	  }
	}
	assert(realCount == deletedCount);
  }

  int adjacent(int u, int v) const{
	assert(!isGhost(u));
	assert(u != v);
	assert( !isGhost(v) || u != unGhost(v));
	foreach(const ref pr; triangles[u]){
	  if(pr.first == v){
		return pr.second;
	  }
	}
	assert(false);
  }

  bool adjacentExists(int u, int v) const{
	assert(!isGhost(u));
	assert(u != v);
	assert( !isGhost(v) || u != unGhost(v));
	foreach(const ref pr; triangles[u]){
	  if(pr.first == v){
		return true;
	  }
	}
	return false;
  }

  bool edgeExists(int u, int v) const{
	assert(!isGhost(u));
	assert(!isGhost(v));
	foreach(const ref pr ; triangles[u]){
	  //necessary because of my backwards-ish ghost convention...
	  if(pr.first == v || pr.second == v){
		return true;
	  }
	}
	return false;
  }

  Pair adjacentRealTriangle(int u) const{
	assert(!isGhost(u));
	foreach(ref pr; triangles[u]){
	  if(!isGhost(pr.first) && !isGhost(pr.second)){
		return pr;
	  }
	}
	assert(false);

  }

  Pair[] getTriangles(int i) const{
	return triangles[i].dup();
  }


  //todo, return a range instead of allocating an array
  Triangle[] getTriangles() const{
	Triangle[] ret;
	foreach(size_t i, const ref svec; triangles){
	  foreach(const ref pr; svec){
		if((isGhost(pr.first) || i < pr.first) && (isGhost(pr.second) || i < pr.second)){
		  ret ~= Triangle(to!int(i), pr.first, pr.second);
		}
	  }
	}
	return ret;
  }

  int[] getActiveVertices() const{
	int[] ret;
	foreach(size_t i, const ref svec; triangles){
	  if(!svec.empty){
		ret ~= to!int(i);
	  }
	}
	return ret;
  }

  void dumpVertex(int i) const{

	auto set = triangles[i];
	foreach(pr ; set){
	  write("( ");
	  if(isGhost(pr.first)){
		write("Ghost ", unGhost(pr.first));
	  } else {
		write(pr.first);
	  }
	  write(", ");
	  if(isGhost(pr.second)){
		write("Ghost ", unGhost(pr.second));
	  } else {
		write(pr.second);
	  }

	  write(" ) ");
	}
  }

  void dump() const {
	import std.stdio;
	writeln("normal entries");
	foreach(i; 0..triangles.length){
	  if(!triangles[i].empty){
		write(i, ": ");
		dumpVertex(to!int(i));
	  }
	  writeln();
	}

  }


  void addTriangulatedPolygon(const ref TriDB triDB, int[] polyIndices){
	foreach(i; 0..polyIndices.length){
	  foreach(const ref pr; triDB.triangles[i]){
		triangles[polyIndices[i]] ~= Pair(polyIndices[pr.first], polyIndices[pr.second]);
	  }
	}
  }

  //TODO: do this better via caching or something.
  //ghost vertex is first
  Triangle findOutsideEdge() const {
	foreach(i; 0..triangles.length){
	  foreach(const ref pr; triangles[i]){
		if(isGhost(pr.first)){
		  return Triangle(pr.first, pr.second, to!int(i));
		}
		if(isGhost(pr.second)){
		  return Triangle(to!int(i), pr.first, pr.second);
		}
	  }
	}
	assert(false);
  }


  bool containsTriangle(const ref Triangle tri){
	if(!isGhost(tri[0])){
	  return triangles[tri[0]].canFind(Pair(tri[1], tri[2]));
	} else {
	  return triangles[tri[1]].canFind(Pair(tri[2], tri[0]));
	}
  }

  //add a steiner point connected to no triangles
  void addPoint(){
	triangles ~= typeof(triangles[0])();
  }

private:

  SSOVector!(Pair, 9)[] triangles;
  //each elements of triangles[i] means there is a triangle i, pr.first, pr.second
  //if one of the elements in the pair has its MSB set, it means it is a GHOST vertex
  //That vertex is "virtual," and means that it would be the next vertex in a CW traveral of
  //the hull.
  //for example if triangle[i] has j, GHOST(k), then i, j, and j, k are consecutive edges
  //in a CW ordering

  // Intuition:  a Ghost triangle is really 2 external edges!!


}



void writeSVG(Vec)(string filename, const Vec[] points, const ref TriDB triDB){
  import std.array;
  import std.stdio : File;
  import std.algorithm;

  import dtriangulate.svg;


  writefln("dumping %s", filename);
  File f = File(filename, "w");
  auto activePoints = triDB.getActiveVertices();
  auto scaledPoints = prepareSVG(f, points, activePoints);
  Triangle[] tris = triDB.getTriangles();

  foreach(const ref tri ; tris){
	if(!TriDB.isGhost(tri[0])  && !TriDB.isGhost(tri[1]) && !TriDB.isGhost(tri[2])){
	  foreach(i; [0,1,2]){
        svgLine(f, scaledPoints[tri[i]],scaledPoints[tri[i + 1]], .0002);
	  }
	} else {
	  auto p1 = TriDB.isGhost(tri[0]) ? tri[1] : tri[0];
	  auto p2 = TriDB.isGhost(tri[2]) ? tri[1] : tri[2];
      svgLine(f, scaledPoints[p1], scaledPoints[p2], .002, true);
	}
  }

  foreach(i ; activePoints){
	const auto p = scaledPoints[i];
	f.writefln( "<circle cx=\"%.8f\" cy=\"%.8f\" r=\"%.8f\" fill=\"black\" />", p.x, p.y, .01);
	f.writefln( "<g transform=\"translate(%.8f, %.8f) scale(1, -1)\" >" ~
				"<text x=\"0\" y=\"0\" font-family=\"Verdana\" font-size=\"%.8f\" fill=\"red\" >%d</text></g>",
				p.x, p.y, .025, i);

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
      svgLine(f, scaledPoints[e.first], scaledPoints[e.second], .001, false, colors[color % colors.length]);
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
