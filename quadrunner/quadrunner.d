
int main(string[] args){

  import dtriangulate.quadrangulate;
  import dtriangulate.edgecomplex;

  import std.stdio;
  import std.random;
  import std.range;
  import std.algorithm;
  import std.array;
  import gl3n.linalg : vec2;

  writeln("quadrangulating");

  auto rnd = Random(42); //predictable seed

  auto points = iota(25).map!( _ => vec2(uniform(0.0f,1.0f,rnd),uniform(0.0f,1.0f,rnd))).array;

  writeln(points);

  auto ec = quadrangulate(points);

  writeSVG("quads.svg", points, ec);
  //  writeln("EC: ", ec);

  foreach(pg; ec.polygons){
    writeln(pg);
  }

  return 0;
}
