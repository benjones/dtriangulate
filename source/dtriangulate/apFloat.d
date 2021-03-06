module dtriangulate.apFloat;

import dtriangulate.fma;

//adaptive precision floating point number
//based on the shewchuck paper
//N == 0 means "dynamically sized."  Positive int sizes use an internal static array
//dynamic AF's store data on the heap and keep a slice
//the notation is backwards from shewchuck.  We store things in decreasing order
struct AdaptiveFloat(FP, int N = 1){

private:
  static if(N > 0){
	FP[N] data;
  } else {
	FP[] data;
  }
  auto size() const{ return data.length; }

  enum StaticSizeLimit = 8;
  
public:
  static if(N == 1){
	this(FP _fp){
	  data[0] = _fp;
	}
  } 
  
  FP asReal(){
	return data[0]; 
  }
  
  AdaptiveFloat!(FP, N) opUnary(string op)() const if(op == "-"){
	AdaptiveFloat!(FP, N) ret;
	ret.data[] = data[];
	foreach(ref val; ret.data){
	  val = -val;
	}
	return ret;
  }


  
  auto opBinary(string op, int M)(auto ref AdaptiveFloat!(FP, M) rhs) if(op == "+"){

	static if(M == 1 && N == 1){
	  return twoAdd(data[0], rhs.data[0]);
	} else static if(M ==1){
	  auto ret = incrementExpansion();
	  ret.unsafePlusEq(rhs.data[0], 0, size());
	  return ret;
	} else static if(N ==1){
	  auto ret = rhs.incrementExpansion();
	  ret.unsafePlusEq(data[0], 0, rhs.size());
	  return ret;
	} else {
	  auto ret = expandToAdd(rhs);
	  foreach(i; 0..rhs.size()){
		ret.unsafePlusEq(rhs.data[rhs.size() - i - 1], i, size());
	  }
	  static if(M == 0 || N == 0 || (M+N > StaticSizeLimit)){
		if(ret.size() > StaticSizeLimit){
		  ret.compress();
		}
	  }
	  return ret;
	}
	
  }

  auto opBinary(string op, int M)(auto ref AdaptiveFloat!(FP, M) rhs) if( op == "-"){
	static if(M == 1 && N == 1){
	  return twoSub(data[0], rhs.data[0]);

	} else static if(M == 1){
	  auto ret = incrementExpansion();
	  ret.unsafeMinusEq(rhs.data[0], 0, size());
	  return ret;
	  
	} else static if(N == 1){
	  auto ret = rhs.incrementExpansion();
	  ret.unsafeMinusEq(data[0], 0, rhs.size());
	  return ret;
	  
	} else {
	  auto ret = expandToAdd(rhs);
	  foreach(i; 0..rhs.size()){
		ret.unsafeMinusEq(rhs.data[rhs.size() - i - 1], i, size());
	  }
	  static if(M == 0 || N == 0  || (M+N > StaticSizeLimit)){
		if(ret.size() > StaticSizeLimit){
		  ret.compress();
		}
	  }
	  return ret;
	}
  }


  auto opBinary(string op, int M)(auto ref AdaptiveFloat!(FP, M) rhs) if(op == "*"){
	import std.conv;
	import std.stdio;
	static if(M == 1 && N == 1){

	  return twoMul(this, rhs);

	} else static if(M == 1){ //RHS is one float

	  static if(N > 0){
		AdaptiveFloat!(FP, 2*N) ret; 
	  } else {
		AdaptiveFloat!(FP, 0) ret;
		ret.data = new FP[2*size()];
	  }
	  auto temp = AdaptiveFloat!FP(data[size() - 1])*rhs; //T.N == 2

	  ret.data[2*size() -1] = temp.data[1];

	  foreach(i ; 1..size()){
		auto temp2 = AdaptiveFloat!FP(data[size() - 1 - i])*rhs; //T2.N == 2

		//T3.N == 2
		auto temp3 = AdaptiveFloat!FP(temp.data[0]) + AdaptiveFloat!FP(temp2.data[1]);
		ret.data[2*size() - 2*i] = temp3.data[1];

		//sum should have N == 2
		temp = AdaptiveFloat!FP(temp3.data[0]) + AdaptiveFloat!FP(temp2.data[0]);
		ret.data[2*size() - 2*i - 1] = temp.data[1];
	  }
	  ret.data[0] = temp.data[0];

	  return ret;
	  
	} else static if(N ==1){ //LHS is one float, use the above
	  return rhs*this;
	} else {
	  //  M, N > 1... so just do it?
	  auto partial = rhs*AdaptiveFloat!FP(data[0]); //FP 2*M
	  auto ret = partial.expandToMultiply(this);
	  //		auto ret = (rhs*AdaptiveFloat!FP(data[0])).expandBy!(2*M*N - 2*M);
	  
	  foreach(i; 1..size()){
		auto partialProduct = rhs*AdaptiveFloat!FP(data[i]); //N = 2*M
		//writefln("partial prod, %s", partialProduct);
		//TODO, we probably don't need to start at j,
		//since it's the product of something bigger, there's probably
		//a proof that it won't kick in until higher up the chain
		foreach(j; 0..(2*rhs.size())){
		  ret.unsafePlusEq(partialProduct.data[j], j, 2*rhs.size()*i);
		}
	  }
	  static if(M == 0 || N == 0  || (M+N > StaticSizeLimit) ){
		if(ret.size() > StaticSizeLimit){
		  ret.compress();
		}
	  }
	  return ret;
	}
  }
  
  string toString(){
	import std.conv;
	import std.format;
	string ret = "[ ";
	foreach(d ; data){
	  ret ~= format("%.10f", d) ~ ", ";
	}
	ret ~= " ]";
	return ret;
  }

  void checkRep(){
	import std.math : abs;
	auto lastNonzero = data[0];
	foreach(i; 1..size()){
	  assert(abs(data[i]) < abs(lastNonzero));
	  if(abs(data[i]) > 0){
		lastNonzero = data[i];
	  }
	}
  }
  
private:

  //primitive operations for 2 normal floats
  AdaptiveFloat!(FP, 2) twoAdd(FP a, FP b){

	AdaptiveFloat!(FP, 2) ret;
	ret.data[0] = a + b;
	FP bVirt = ret.data[0] - a;
	FP aVirt = ret.data[0] - bVirt;
	FP bRoundoff = b - bVirt;
	FP aRoundoff = a - aVirt;
	ret.data[1] = aRoundoff + bRoundoff;
	return ret;

  }

  AdaptiveFloat!(FP, 2) twoSub(FP a, FP b){
	AdaptiveFloat!(FP, 2) ret;
	ret.data[0] = a - b;
	FP bVirt = a - ret.data[0];
	FP aVirt = ret.data[0] + bVirt;
	FP aRoundoff = a - aVirt;
	FP bRoundoff = bVirt - b;
	ret.data[1] = aRoundoff + bRoundoff;
	return ret;
  }

  AdaptiveFloat!(FP, 2) twoMul(AdaptiveFloat!(FP,1) a, AdaptiveFloat!(FP, 1) b){
	
	AdaptiveFloat!(FP, 2) ret;
	ret.data[0] = a.data[0]*b.data[0];
	auto aSplit = a.split();
	auto bSplit = b.split();
	FP err1 = ret.data[0] - (aSplit.data[0]*bSplit.data[0]);
	FP err2 = err1 - (aSplit.data[1]*bSplit.data[0]);
	FP err3 = err2 - (aSplit.data[0]*bSplit.data[1]);
	ret.data[1] = (aSplit.data[1]*bSplit.data[1]) - err3;
	return ret;
  }
  

  import std.stdio;
  //add 1 more float's worth of precision to this
  //make this a template so it's only compiled when needed
  //otherwise this stack overflows the compiler!
  auto incrementExpansion()() const{
	static if(N == 0){
	  AdaptiveFloat!(FP, 0) ret;
	  ret.data = new FP[size() +1];
	} else {
	  AdaptiveFloat!(FP, N+1) ret; 
	}

	ret.data[1.. $] = data[];
	ret.data[0] = 0; //expanded by 1, so set new entry to 0;
	return ret;
  }

  auto expandToAdd(int M)(const ref AdaptiveFloat!(FP, M) other) const{
	static if( N==0 || M==0 || ( (M+N) > StaticSizeLimit) ){
	  //	  writeln("expand add made a dynamic");

	  AdaptiveFloat!(FP, 0) ret;
	  ret.data = new FP[size() + other.size()];
	} else {
	  AdaptiveFloat!(FP, M+N) ret;
	}

	ret.data[other.size()  .. $] = data[];
	ret.data[0..other.size()] = 0;

	return ret;
  }


  auto expandToMultiply(int M)(const ref AdaptiveFloat!(FP, M) other) const{
	static if(N==0 || M==0 || ((M*N) > StaticSizeLimit) ){
	  //	  writeln("expand multiply made a dynamic");
	  AdaptiveFloat!(FP, 0) ret;
	  ret.data = new FP[size()*other.size()];
	} else {
	  AdaptiveFloat!(FP, M*N) ret;
	}

	ret.data[(ret.size() - size()) .. $] = data[];
	ret.data[0..(ret.size() - size())] = 0;
	return ret;
  }

  //eliminate zeros, and consolidate entries that don't use much precision
  static if(N == 0){
	void compress(){

	  import std.conv : to;
	  //loop down, then up
	  FP Q = data[0]; //biggest element
	  auto bottom = 0;
	  foreach(i; 1..size()){
		auto exact = AdaptiveFloat!FP(Q) + AdaptiveFloat!FP(data[i]);
		Q = exact.data[0];
		if(exact.data[1] != 0){
		  data[bottom] = Q;
		  ++bottom;
		  Q = exact.data[1];
		}
	  }
	  data[bottom] = Q; //this is the smallest nonzero component
	  //at this point, the beginning of the array should have big stuff

	  //I believe there is a typo in the shewchuck paper...
	  auto top = size() -1;
	  for(int i = to!int(bottom -1); i >= 0; --i){
		auto exact = AdaptiveFloat!FP(data[i]) + AdaptiveFloat!FP(Q);
		Q = exact.data[0];
		if(exact.data[1] != 0){
		  data[top] = exact.data[1]; //shewchuck writes "Q" here, but it should be "q"
		  top--;
		}
	  }
	  data[top] = Q;
	  import std.stdio;
	  data = data[top..$];
	  //this compresses everything, moving it towards the end
	}
  }
  
  
  void unsafePlusEq(FP f, ulong first, ulong len){
	auto q = AdaptiveFloat!FP(f);
	foreach(i; first..(first + len)){
	  auto t1 = AdaptiveFloat!FP(data[size() - i - 1]);
	  auto temp = q + t1;
	  q.data[0] = temp.data[0];
	  data[size() - i - 1] = temp.data[1];
	}
	data[size() - first - len -1] = q.data[0];
  }

  void unsafeMinusEq(FP f, ulong first, ulong len){
	auto q = AdaptiveFloat!FP(f);
	auto t1 = AdaptiveFloat!FP(data[size() - first - 1]);
	auto temp = t1 - q;
	data[size() - first - 1] = temp.data[1];
	unsafePlusEq(temp.data[0], first + 1, len -1);
  }

  static if(N == 1){
	AdaptiveFloat!(FP, 2) split() {
	  enum FP splitConst = cast(FP) ( 2 << ((FP.mant_dig + 1)/2) + 1);
	  FP c = splitConst*data[0];
	  FP aBig = c - data[0];
	  AdaptiveFloat!(FP, 2) ret;
	  ret.data[0] = c - aBig;
	  ret.data[1] = data[0] - ret.data[0];
	  return ret;
	}
  }
  
}

unittest{
  //testing + and -
  AdaptiveFloat!float ap = AdaptiveFloat!float(3.14);
  auto ap2 = -ap;
  assert(ap.asReal() == -(ap2.asReal()));
  assert(ap.asReal() == 3.14f);
  assert(ap2.asReal() == -3.14f);
  assert((ap + ap2).asReal == 0);

  AdaptiveFloat!float a = 1.0f;
  AdaptiveFloat!float b = float.epsilon/2;
  auto c = a + b;
  AdaptiveFloat!float d = -1.0f;
  assert((c + d).asReal() == float.epsilon/2);
  assert((c - a).asReal() == float.epsilon/2);
  assert(((c - a) - b).asReal() == 0);
  
}

unittest{
  import std.stdio;
  import std.conv;
  //multiplication
  AdaptiveFloat!(float,1) a = 1;
  auto eps = AdaptiveFloat!float( float.epsilon/2);
  auto b = a + eps;


  //writeln("b " ~ to!string(b));
  auto c = a * b; //1 + eps   (1)*(2)
  //writeln("c " ~ to!string(c));
  auto cPrime = b * a;
  //writeln("cPrime " ~ to!string(cPrime));
  auto d = c - eps;
  assert(d.asReal() == 1.0f);

  
  auto e = c*eps; //eps + eps*eps
  auto f = e - eps - eps*eps;
  assert(f.asReal() == 0.0f);
  
  //writeln("computing g");
  auto g = b*b; //1 + 2*eps + eps*eps
  //writeln("g " ~ to!string(g));
  auto h = eps*eps;
  auto i = g - h;
  //writeln("i " ~ to!string(i));
  auto j = AdaptiveFloat!float(2)*eps;
  //writeln("i = j " ~ to!string(i - j));
  assert( (i - j).asReal() == 1.0f);
  
  assert( (i - eps - eps).asReal() == 1.0f);
}


