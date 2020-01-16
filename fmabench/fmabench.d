module fmabench.fmabench;

import dtriangulate.predicates;
import dtriangulate.fma;

void main(){
  import std.stdio;
  import std.random : uniform;
  import std.math : sgn, nextUp;
  
  struct Vec{ float x, y; }

  static Vec rv(){
    return Vec(uniform(-1000.0f, 1000.0f), uniform(-1000.0f, 1000.0f));
  }

  struct Vecs{ Vec a, b, c; }

  static Vecs uniformVecs(size_t i){
    return Vecs(rv(), rv(), rv());
  }


  static void runBenchmark(){
    size_t naiveFailures = 0;
    size_t fmaFailures = 0;

    float min = 1.1f, max = 2.2f;
    Vecs vecs = Vecs(Vec(min, min), Vec(max, max), Vec(min, min));
    float curr = min;
    
    //foreach(iter; 0..1000000000){
    size_t iter = 0;
    while(curr <= max){
      //      const vecs = makeVecs(iter);
      float scale = iter % 2 == 0 ? -1.0f : 1.0f;
      vecs.c = Vec(curr + scale*float.epsilon, curr);
      curr = nextUp(curr);
      ++iter;
      with(vecs){
        auto exact = orient2D(a, b, c);
        auto naive = orient2DNaive(a, b, c);
        auto fma = orient2DFMA(a, b, c);
        
        if(sgn(exact) != sgn(naive)){
          ++naiveFailures;
          writeln("Naive Failure: ", a, b, c, " true: ", exact, " naive: ", naive);
        }
        
        if(sgn(exact) != sgn(fma)){
          ++fmaFailures;
          //writeln("FMA Failure: ", a, b, c, " true: ", exact, " fma: ", fma);
        }
      }
      
      //orient2DNaive, orient2DFMA
    }
    writeln("total failures naive: ", naiveFailures, " FMA: ", fmaFailures, " in ", iter, " iters");
  }

  //runBenchmark();
  triangulateBenchmark();

}

void triangulateBenchmark(){
  import dtriangulate.triangulate;
  import dtriangulate.predicates;
  import std.random : uniform;
  import std.stdio;
  
  struct Vec{ float x, y; }
  bool[Vec] seen;
  Vec[] points = new Vec[3000000];
  foreach(i; 0 .. points.length){
    Vec v = Vec(uniform(1.0f, 2.0f), uniform(1.0f, 2.0f));
    while(v in seen){
      v = Vec(uniform(1.0f, 2.0f), uniform(1.0f, 2.0f));
    }
    points[i] = v;
    seen[v] = true;
  }

  auto db = delaunayTriangulate(points);
  writeln("Made ", db.getTriangles.length, " triangles with " , extendedO2dCount, " adaptive calcs\n");
}
