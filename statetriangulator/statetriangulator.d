
int main(string[] args){
  import std.stdio;
  import std.string;
  import std.conv;

  import shapefil;
  import dtriangulate.triangulate;

  immutable string usage = "stateTriangulator /path/to/shapefileFolder/fileNoPrefix [shape number [part number] ]";

  if(args.length < 2){
	writeln("usage: " ~ usage);
	return 1;
  }

  string prefix = args[1];
  
  auto dbfHandle = DBFOpen(toStringz(prefix ~ ".dbf"), toStringz("rb"));

  
  if(args.length < 3){
	auto nRecords = DBFGetRecordCount(dbfHandle);
	writefln("shape in file: ");
	foreach(i; 0..nRecords){
	  auto name =fromStringz( DBFReadStringAttribute(dbfHandle, i, 5));
	  writefln("shape %d: %s", i, name);
	}
	return 0;
  }

  int shapeNumber = to!int(args[2]);
  auto shpHandle = SHPOpen(toStringz(prefix ~ ".shp"), toStringz("rb"));
  auto shape = SHPReadObject(shpHandle, shapeNumber);
  auto nParts = shape.nParts;
  auto nVerts = shape.nVertices;
  
  if(args.length < 4){
	writeln("parts in shape: ");
	foreach(i; 0..shape.nParts){
	  writefln("part %d has %d verts", i,
			   (i < (nParts -1) ? shape.panPartStart[i + 1] : nVerts ) - shape.panPartStart[i]);
	}
	writefln("shape has %d parts and %d verts", nParts, nVerts);
	return 0;
  }

  int partNumber = to!int(args[3]);

  import gl3n.linalg : vec2;
  vec2[] points = new vec2[(partNumber < (nParts -1) ? shape.panPartStart[partNumber + 1] : nVerts ) -
						   shape.panPartStart[partNumber]];
  foreach(i; 0..points.length){
	points[i].x = shape.padfX[i + shape.panPartStart[partNumber]];
	points[i].y = shape.padfY[i + shape.panPartStart[partNumber]];
  }
  writeln(points);

  import std.algorithm.iteration : uniq;
  import std.array : array;
  points = points.uniq.array;

  import std.algorithm.searching : findAdjacent;
  vec2[] rng = points[];
  while((rng = rng.findAdjacent).length > 0){
	write("duplicate: ");
	writeln(rng[0]);
	rng = rng[1..$];
  }

  auto triDB = delaunayTriangulate(points);
  //  auto tris = triDB.getTriangles();
  string svgFile = to!string(fromStringz(DBFReadStringAttribute(dbfHandle, shapeNumber, 5)) ~ ".svg");
  writeSVG(svgFile, points, triDB);
  

  return 0;
}

