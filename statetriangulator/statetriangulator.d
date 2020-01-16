
int main(string[] args){
  import std.stdio;
  import std.string;
  import std.conv;
  import std.algorithm;
  import std.array;
  import std.range;
  
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
  //  writeln(points);

  import std.algorithm.iteration : uniq;
  import std.array : array;
  auto sizeBefore = points.length;
  points = points.uniq.array;
  auto sizeAfter = points.length;
  if(sizeBefore != sizeAfter){
	writeln("unique deleted ", sizeBefore - sizeAfter, " points");
  }
  if(points[0] == points[$-1]){
	writeln("deleted duplicate start/end point");
	points = points[0..$-1]; 
  }

  /*  auto pointsCopy = points.dup;
  sort(pointsCopy);
  if(!pointsCopy.findAdjacent.empty){
	writeln("duplicate points, oh shit!");
	assert(0);
	}*/
  

  //  points ~= 0.5*(points[259] + points[260]);
  
  /*
  import std.algorithm.searching : findAdjacent;
  vec2[] rng = points[];
  while((rng = rng.findAdjacent).length > 0){
	write("duplicate: ");
	writeln(rng[0]);
	rng = rng[1..$];
	}*/
  writeln("number of unique points: ", points.length);
  //points = points[0..282];
  //  foreach( i, ref p; points[1412 .. 1419]){
  //	writefln("%d: vec2(%0.12f, %0.12f) ", i + 1412, p.x, p.y);
  //    }

  
  auto triDB = delaunayTriangulate(points);

  //  auto tris = triDB.getTriangles();
  string svgFile = to!string(fromStringz(DBFReadStringAttribute(dbfHandle, shapeNumber, 5)) ~ ".svg");
  writeSVG(svgFile, points, triDB);


  auto segments = iota(points.length).map!(i => Pair(to!int(i), to!int((i + 1)%points.length))).array;
  //make constrained delaunay next
  
  bool[Pair] segmentSet;
  foreach(seg; segments){
	segmentSet[seg] = true;
  }

  //  segmentSet.remove(Pair(259,260));
  //  segmentSet[Pair(259, 282)] = true;
  //  segmentSet[Pair(282, 260)] = true;
  
  makeConstrainedDelaunay(points, triDB, segmentSet);

  string svgFileConstrained = to!string(svgFile[0..$-4] ~ "constrained.svg");
  writeSVG(svgFileConstrained, points, triDB);


  cutOffScraps(points, triDB, segmentSet);
  
  string svgFileTrimmed = to!string(svgFile[0..$-4] ~ "trimmed.svg");
  writeSVG(svgFileTrimmed, points, triDB);

  foreach(i; 0..10){
	if(!refinementStep(points, triDB, segmentSet, .1, 30, .5)){
	  break;
	}
	writeln("mesh modified in iteration ", i);
  }


  string svgFileRefined = to!string(svgFile[0..$-4] ~ "refined.svg");
  writeSVG(svgFileRefined, points, triDB);

  import dtriangulate.predicates;
  writeln("adaptive O2ds: ", extendedO2dCount);
  
  /*  
  foreach(i; 0..points.length){
	if(!triDB.edgeExists(cast(int)i, cast(int)((i+1)%points.length))){
	  writeln("clearing cavity: ", i, ' ', (i+1) % points.length);

	  writeln("vertices for ", i);
	  triDB.dumpVertex(to!int(i));
	  writeln();
	  writeln("vertices for ", (i +1) % points.length);
	  triDB.dumpVertex(to!int( (i + 1) % points.length));
	  writeln();
	  
	  clearCavity(points, triDB, Pair(cast(int)i, cast(int)((i+1)%points.length)));
	  string filename = to!string("cleared" ~ to!string(i) ~ ".svg");
	  writeSVG(filename, points, triDB);
	}
	}*/
  

  return 0;
}

