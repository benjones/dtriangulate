module dtriangulate.predicates;
import dtriangulate.apFloat;

private auto square(T)(auto ref T t){ return t*t; }

auto orient2D(Vec)(auto ref Vec a, auto ref Vec b, auto ref Vec c){
  alias FP = typeof(a.x);
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

  import std.math;
  //from Shewchuck, of course
  enum FP = errorBound = (3*FP.epsilon + 16*square(FP.epsilon));
  Real det = axby - aybx;

  if(abs(det) < errorBound*detSum){
	return orient2DAdaptive(a, b, c);
  } else {
	return det;
  }
}


auto orient2DAdaptive(FP, Vec)(auto ref Vec a, auto ref Vec b, auto ref Vec c){
  alias FP = typeof(a.x);
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

auto inCircle(Vec)(auto ref Vec a, auto ref Vec b, auto ref Vec c, auto ref Vec d){
  import std.math;
  alias FP = typeof(a.x);
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
  alias FP = typeof(a.x);
  
  auto ax = AdaptiveFloat!FP(a.x);
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
  pragma(msg, typeof(aCofactor));
  
  auto axcy = adx*cdy;
  auto aycx = ady*cdx;
  auto bCofactor = bLift*(aycx - axcy);

  auto axby = adx*bdy;
  auto aybx = ady*bdx;
  auto cCofactor = cLift*(axby - aybx);
  pragma(msg, typeof(cCofactor));
	
  //too many components
  pragma(msg, typeof(aCofactor + bCofactor));
  pragma(msg, typeof(aCofactor + bCofactor + cCofactor));
  return (aCofactor + bCofactor + cCofactor).asReal();
  
}

unittest{
  import std.random;
  import std.math;
  import std.algorithm;
  import std.stdio;

  int failures = 0;
  foreach(i ; 0..1000000){
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
