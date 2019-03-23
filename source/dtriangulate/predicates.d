module dtriangulate.predicates;

public import dtriangulate.apFloat;

import std.traits : Unqual;

private auto square(T)(auto ref T t){ return t*t; }

auto orient2D(Vec)(auto ref Vec a, auto ref Vec b, auto ref Vec c){

  alias FP = Unqual!(typeof(a.x));
  //det of [ [ (ax - cx) (ay - cy) ] [bx - cx, by - cy]]
  FP acx = a.x - c.x;
  FP acy = a.y - c.y;
  FP bcx = b.x - c.x;
  FP bcy = b.y - c.y;

  //if any of those are 0, the subtraction is exact
  FP axby = acx*bcy;
  FP aybx = acy*bcx;
  
  FP detSum;
  
  //if the parts have opposite signs, we're cool
  if(axby > 0){
	if(aybx <= 0){ //the sign of axby can't be wrong, so 0 for aybx is fine
	  return axby - aybx;
	} else {
	  detSum = axby + aybx;
	}
  } else if(axby < 0){
	if(aybx >= 0){
	  return axby - aybx;
	} else {
	  detSum = -axby - aybx;
	}
  } else {
	//axby is 0, which is exact, so return -aybx
	return -aybx;
  }

  import std.math : abs;
  //from Shewchuck, of course
  enum FP errorBound = (3*FP.epsilon + 16*square(FP.epsilon));
  FP det = axby - aybx;

  if(abs(det) < errorBound*detSum){
	return orient2DAdaptive(a, b, c);
  } else {
	return det;
  }
}


auto orient2DAdaptive(Vec)(auto ref Vec a, auto ref Vec b, auto ref Vec c){

  alias FP = Unqual!(typeof(a.x));
  auto ax = AdaptiveFloat!FP(a.x);
  auto ay = AdaptiveFloat!FP(a.y);
  auto bx = AdaptiveFloat!FP(b.x);
  auto by = AdaptiveFloat!FP(b.y);
  auto cx = AdaptiveFloat!FP(c.x);
  auto cy = AdaptiveFloat!FP(c.y);

  auto acx = ax - cx;
  auto acy = ay - cy;
  auto bcx = bx - cx;
  auto bcy = by - cy;

  auto axby = acx*bcy;
  auto aybx = acy*bcx;

  return (axby - aybx).asReal();
}

bool isInCircle(Vec)(auto ref Vec a, auto ref Vec b, auto ref Vec c, auto ref Vec d){
  return inCircle(a, b, c, d) > 0;
}

auto inCircle(Vec)(auto ref Vec a, auto ref Vec b, auto ref Vec c, auto ref Vec d){
  import std.math;

  alias FP = Unqual!(typeof(a.x));
  FP adx = a.x - d.x;
  FP ady = a.y - d.y;
  FP bdx = b.x - d.x;
  FP bdy = b.y - d.y;
  FP cdx = c.x - d.x;
  FP cdy = c.y - d.y;

  auto aLift = square(adx) + square(ady);
  auto bLift = square(bdx) + square(bdy);
  auto cLift = square(cdx) + square(cdy);

  auto bxcy = bdx*cdy;
  auto bycx = bdy*cdx;
  auto aCofactor = aLift*(bxcy - bycx);

  auto axcy = adx*cdy;
  auto aycx = ady*cdx;
  auto bCofactor = bLift*(aycx - axcy);

  auto axby = adx*bdy;
  auto aybx = ady*bdx;
  auto cCofactor = cLift*(axby - aybx);

  auto permanent = aLift*(abs(bxcy) + abs(bycx)) +
	bLift*(abs(axcy) + abs(aycx)) +
	cLift*(abs(axby) + abs(aybx));

  auto det = aCofactor + bCofactor + cCofactor;
  
  enum FP errorBound =
	10*FP.epsilon +	96*square(FP.epsilon);
  
  if(abs(det) > errorBound*permanent){
	return det; //I don't think we ever care if it's on the circle...
  } else {
	return inCircleAdaptive(a, b, c, d);
  }
}

auto inCircleAdaptive( Vec)(auto ref Vec a, auto ref Vec b, auto ref Vec c, auto ref Vec d){

  alias FP = Unqual!(typeof(a.x));
  alias AF1 = AdaptiveFloat!(FP, 1);
  AF1 ax = AF1(a.x);
  auto ay = AdaptiveFloat!FP(a.y);
  auto bx = AdaptiveFloat!FP(b.x);
  auto by = AdaptiveFloat!FP(b.y);
  auto cx = AdaptiveFloat!FP(c.x);
  auto cy = AdaptiveFloat!FP(c.y);
  auto dx = AdaptiveFloat!FP(d.x);
  auto dy = AdaptiveFloat!FP(d.y);


   //2-component AFs
  auto adx = ax - dx;
  auto ady = ay - dy;
  auto bdx = bx - dx;
  auto bdy = by - dy;
  auto cdx = cx - dx;
  auto cdy = cy - dy;

  //8 component AFs
  auto aLift = square(adx) + square(ady);
  auto bLift = square(bdx) + square(bdy);
  auto cLift = square(cdx) + square(cdy);

  //128 compoent AFS
  auto bxcy = bdx*cdy;
  auto bycx = bdy*cdx;

  //256
  auto aCofactor = aLift*(bxcy - bycx);
  
  auto axcy = adx*cdy;
  auto aycx = ady*cdx;
  auto bCofactor = bLift*(aycx - axcy);

  auto axby = adx*bdy;
  auto aybx = ady*bdx;
  auto cCofactor = cLift*(axby - aybx);
	
  //too many components
  return (aCofactor + bCofactor + cCofactor).asReal();
  
}

unittest{
  import std.random;
  import std.math;
  import std.algorithm;
  import std.stdio;

  int failures = 0;
  foreach(i ; 0..100000){
	float x = uniform(0.0f, 1.0f);
	float y = uniform(0.0f, 1.0f);

	float[3] angles;
	foreach(j; 0..3){
	  angles[j] = uniform(0.0f, float(2*PI));
	}
	angles[].sort();

	struct Vec{ float x, y; }

	Vec[3] pts;
	foreach(j; 0..3){
	  pts[j] = Vec(cos(angles[j]), sin(angles[j]));
	}

	auto aX = AdaptiveFloat!float(x);
	auto aY = AdaptiveFloat!float(y);
	auto r = AdaptiveFloat!float(1);

	auto dist = aX*aX + aY*aY - r;

	bool insideRad = dist.asReal() < 0;

	bool inPred = inCircle(pts[0], pts[1], pts[2], Vec(x, y)) > 0;
	
	if(insideRad != inPred){
	  ++failures;
	  writefln("failure: (%f, %f), (%f, %f), (%f, %f), (%f, %f), dist: %d, pred: %d",
			   pts[0].x, pts[0].y, pts[1].x, pts[1].y, pts[2].x, pts[2].y,x, y, insideRad, inPred);
	}
	
  }
  assert(failures < 10);
}


/*
  closest point(p0, p1, t);  Which of p0 and p1 is closer to t?
  
  want to compute (p0 - t)^2 - p(1 -t)^2 and compare it to 0

  expanded out:
  p1x(p1x - 2tx) - p2x(p2x - 2tx) + p1y(p1y - 2ty) - p2y(p2y - 2ty) >? 0

  shewchuck style analysis:
  t (true)
  t0 = p1x - 2tx -- note that 2tx is exact since multiplying by powers of 2 doesn't add error
  t1 = p2x - 2tx
  t2 = p1y - 2ty
  t3 = p2y - 2ty
  t4 = p1x*t0
  t5 = p2x*t1
  t6 = p1y*t2
  t7 - p2y*t3
  t8 = t4 - t5
  t9 = t6 - t7
  t10 = t8 + t9

  //similary for Xs (the FP approximations)
  x0 = p1x - 2tx, etc
  t0 = x0 +- eps |x0| //t0 is within eps*mag of the FP approximation
  t1 = x1 +-eps |x1| ...

  t4 = p1x*to = p1x(x0 +- eps|x0|) 
     = (p1x*x0) +- (p1x *eps |x0|)
     = (x4 +- eps|x4|) +- (eps* (|x4| +- eps |x4|) )
     = x4 +- eps |x4| (2 + eps)

  similar for t5, 6, 7

  t8 = t4 - t5 = x4 +- eps |x4| (2 + eps) - x5 +- eps |x5| (2 + eps)
     = x8 +- eps |x8| +- eps(2+eps)(|x4| + |x5|)
     = x8 +- eps  (  |x8| +(2+eps)(|x4| + |x5|)
  similarly for t9

  t10 = t8 + t9 = x8 +- eps  (  |x8| +(2+eps)(|x4| + |x5|) + 
                  x9 +- eps  (  |x9| +(2+eps)(|x6| + |x7|)

      = x10 +- eps |x10| +- ...

 */

//which of p0 or p1 is closer to t?
int closer(Vec)(auto ref Vec p0, auto ref Vec p1, auto ref Vec t){

  alias FP = Unqual!(typeof(p0.x));
  FP tx2 = 2*t.x;
  FP ty2 = 2*t.y;
  FP x0 = p0.x - tx2;
  FP x1 = p1.x - tx2;
  FP x2 = p0.y - ty2;
  FP x3 = p1.y - ty2;

  FP x4 = p0.x*x0;
  FP x5 = p1.x*x1;
  FP x6 = p0.y*x2;
  FP x7 = p1.y*x3;

  FP x8 = x4 - x5;
  FP x9 = x6 - x7;
  FP x10 = x8 + x9;
  import std.math : abs;
  //check the error bounds:
  FP minError = FP.epsilon*(abs(x10) + abs(x8) + abs(x9) +
							(2 + FP.epsilon)*(abs(x4) + abs(x5) + abs(x6) + abs(x7)));
  
  if(abs(x10) > minError){
	//less than 0 if p0 is closer
	return x10 < 0 ? 0 : 1; 
  } else {

	return closerAdaptive(p0, p1, t);
  }
  
}

int closerAdaptive(Vec)(auto ref Vec p0, auto ref Vec p1, auto ref Vec t){

  alias FP = Unqual!(typeof(p0.x));
  auto tx2 = AdaptiveFloat!FP(2*t.x);
  auto ty2 = AdaptiveFloat!FP(2*t.y);

  auto p0x = AdaptiveFloat!FP(p0.x);
  auto p0y = AdaptiveFloat!FP(p0.y);
  auto p1x = AdaptiveFloat!FP(p1.x);
  auto p1y = AdaptiveFloat!FP(p1.y);
  
  auto x0 = p0x - tx2;
  auto x1 = p1x - tx2;
  auto x2 = p0y - ty2;
  auto x3 = p1y - ty2;

  auto x4 = p0x*x0;
  auto x5 = p1x*x1;
  auto x6 = p0y*x2;
  auto x7 = p1y*x3;

  auto x8 = x4 - x5;
  auto x9 = x6 - x7;
  auto x10 = x8 + x9;

  return x10.asReal() < 0 ? 0 : 1;
  
}

unittest{
  struct Vec{ float x, y; }


  auto p0 = Vec(0,0);
  auto p1 = Vec(float.epsilon, 0);
  auto t = Vec(1e10, 0);
  assert(closer(p0, p1, t) == 1);
  assert(closer(p1, p0, t) == 0);


  p1 = Vec(0, float.epsilon);
  t = Vec(0, 1e10);
  assert(closer(p0, p1, t) == 1);
  assert(closer(p1, p0, t) == 0);


  p1 = Vec(-float.epsilon, -float.epsilon);
  t = Vec(-1e10, -1e10);
  assert(closer(p0, p1, t) == 1);
  assert(closer(p1, p0, t) == 0);
  

  
  
}

bool segmentsCross(Segment)(auto ref Segment a, auto ref Segment b){
  bool bCrossesA = (orient2D(a.first, a.second, b.first) > 0) !=
	(orient2D(a.first, a.second, b.second) > 0);
  bool aCrossesB = (orient2D(b.first, b.second, a.first) > 0) !=
	(orient2D(b.first, b.second, a.second) > 0);
  
  //both must be true
  return bCrossesA && aCrossesB;

}
