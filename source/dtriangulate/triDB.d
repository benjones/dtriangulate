module dtriangulate.tridb;

public import dtriangulate.ssoVector;
public import dtriangulate.predicates;

import gl3n.linalg : vec2;

import std.stdio;
import std.conv;

struct Pair{
  int first, second;
  bool opEquals(Pair rhs) const{
	return first == rhs.first && second == rhs.second;
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


  ref int opIndex(int i){
	return v[i > 2 ? i -3 : i]; //faster than modulo?
  }
  
  int opIndex(int i) const{
	return v[i > 2 ? i -3 : i]; //faster than modulo?
  }
	  
  bool opEquals(Triangle rhs) const{
	return v[0] == rhs.v[0] && v[1] == rhs.v[1] && v[2] == rhs.v[2];
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


  
  void addTriangle(Vec)(Triangle t, const ref Vec[] points){
	//writeln("adding: ", t);
	if(!isGhost(t[0]) && !isGhost(t[1]) && !isGhost(t[2])){
	  assert(orient2D(points[t[0]], points[t[1]], points[t[2]]) > 0);
	}
	foreach(i; 0..3){
	  if(!isGhost(t[i])){ //don't add ghost edges
		triangles[t[i]] ~= Pair(t[i+1], t[i+2]);
	  }
	}
  }
  
  void deleteTriangle(Triangle t){
	//writeln("deleting: ", t);
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


  //need to distinguish real neighbors from ghosts!
  /*  int adjacent(int u, int v) const{
	write("searching adjacent: ");
	if(isGhost(u)){ write("Ghost: ", unGhost(u), ", "); }
	else {write(u, ", ");}
	if(isGhost(v)){ writeln("Ghost: ", unGhost(v), ", "); }
	else {writeln(v, ", ");}
	  
	assert(!isGhost(u));
	assert(u != v);
	foreach(const ref pr; triangles[u]){
	  if(pr.first == v){
		return pr.second;
	  }
	}
	assert(false);
	//return -1;
	}*/


  int adjacentGhost(int u, int v) const{
	assert(!isGhost(u));
	assert(u != v);
	foreach(const ref pr; triangles[u]){
	  if(pr.first ==v && isGhost(pr.second)){
		return pr.second;
	  }
	}
	assert(false);
	
  }
  int adjacentReal(int u, int v) const{
	assert(!isGhost(u));
	assert(u != v);
	foreach(const ref pr; triangles[u]){
	  if(pr.first ==v && !isGhost(pr.second)){
		return pr.second;
	  }
	}
	assert(false);
	
  }

  bool adjacentRealExists(int u, int v) const{
	assert(!isGhost(u));
	assert(!isGhost(v));
	foreach(const ref pr; triangles[u]){
	  if(pr.first == v && !isGhost(pr.second)){
		return true;
	  }
	}
	return false;
  }
  
  //return the real adjacent vertex if it exists
  //otherwise returns the ghost version
  int adjacentRealIfExists(int u, int v) const{
	assert(!isGhost(u));
	assert(u != v);
	int ret;
	bool found;
	foreach(const ref pr; triangles[u]){
	  if(pr.first ==v){
		if(!isGhost(pr.second)){
		  return pr.second;
		} else {
		  ret = pr.second;
		  found = true;
		}
	  }
	}
	assert(found);
	return ret;
  }
  

  //
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
  
private:

  SSOVector!(Pair, 9)[] triangles;
  //each elements of triangles[i] means there is a triangle i, pr.first, pr.second
  //if one of the elements in the pair has its MSB set, it means it is a GHOST edge
  //which and the negative value is the next vertex in a CCW ordering
  //for example if triangle[i] has j, -k, then i, j, and j, k are consecutive edges
  //in a CCW ordering


}

