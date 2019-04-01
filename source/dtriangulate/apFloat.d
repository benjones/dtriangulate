module dtriangulate.apFloat;

//adaptive precision floating point number
//based on the shewchuck paper
//N == 0 means "dynamically sized."  Positive int sizes use an internal static array
//dynamic AF's store data on the heap and keep a slice
struct AdaptiveFloat(FP, int N = 1){

private:
  static if(N > 0){
	FP[N] data;
  } else {
	FP[] data;
  }
  auto size() const{ return data.length; }
  
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


  
  auto  opBinary(string op, int M)(auto ref AdaptiveFloat!(FP, M) rhs) if(op == "+"){
	
	static if(M == 1 && N == 1){
	  return twoAdd(data[0], rhs.data[0]);
	} else static if(M ==1){
	  auto ret = incrementExpansion();
	  ret.unsafePlusEq(rhs.data[0], 0, N);
	  return ret;
	} else static if(N ==1){
	  auto ret = rhs.incrementExpansion();
	  ret.unsafePlusEq(data[0], 0, M);
	  return ret;
	} else {
	  auto ret = expandToAdd(rhs);
	  foreach(i; 0..M){
		ret.unsafePlusEq(rhs.data[M - i - 1], i, N);
	  }
	  return ret;
	}
	
  }

  AdaptiveFloat!(FP, M + N) opBinary(string op, int M)(auto ref AdaptiveFloat!(FP, M) rhs) if( op == "-"){
	static if(M == 1 && N == 1){
	  return twoSub(data[0], rhs.data[0]);

	} else static if(M == 1){
	  auto ret = incrementExpansion();
	  ret.unsafeMinusEq(rhs.data[0], 0, N);
	  return ret;
	  
	} else static if(N == 1){
	  auto ret = rhs.incrementExpansion();
	  ret.unsafeMinusEq(data[0], 0, M);
	  return ret;
	  
	} else {
	  auto ret = expandToAdd(rhs);
	  foreach(i; 0..M){
		ret.unsafeMinusEq(rhs.data[M - i - 1], i, N);
	  }
	  return ret;
	}
  }


  AdaptiveFloat!(FP, 2*M*N) opBinary(string op, int M)(auto ref AdaptiveFloat!(FP, M) rhs) if(op == "*"){
	import std.conv;
	import std.stdio;
	static if(M == 1 && N == 1){

	  return twoMul(this, rhs);

	} else static if(M == 1){ //RHS is one float
	  
	  AdaptiveFloat!(FP, 2*N) ret;
	  auto temp = AdaptiveFloat!FP(data[N - 1])*rhs; //T.N == 2

	  ret.data[2*N -1] = temp.data[1];

	  foreach(i ; 1..N){
		AdaptiveFloat!(FP, 2*M) temp2 = AdaptiveFloat!FP(data[N - 1 - i])*rhs; //T2.N == 2

		//T3.N == 2
		auto temp3 = AdaptiveFloat!FP(temp.data[0]) + AdaptiveFloat!FP(temp2.data[1]);
		ret.data[2*N - 2*i] = temp3.data[1];

		//sum should have N == 2
		temp = AdaptiveFloat!FP(temp3.data[0]) + AdaptiveFloat!FP(temp2.data[0]);
		ret.data[2*N - 2*i - 1] = temp.data[1];
	  }
	  ret.data[0] = temp.data[0];
	  return ret;
	  
	} else static if(N ==1){ //LHS is one float, use the above
	  return rhs*this;
	} else static if(M < N){
	  return rhs*this;
	} else {
	  //N <= M
	  AdaptiveFloat!(FP, 2*M) partial = rhs*AdaptiveFloat!FP(data[0]); //FP 2*M
	  auto ret = partial.expandToMultiply(this);
	  //		auto ret = (rhs*AdaptiveFloat!FP(data[0])).expandBy!(2*M*N - 2*M);
	  
	  foreach(i; 1..N){
		auto partialProduct = rhs*AdaptiveFloat!FP(data[i]); //N = 2*M
		//writefln("partial prod, %s", partialProduct);
		//TODO, we probably don't need to start at j,
		//since it's the product of something bigger, there's probably
		//a proof that it won't kick in until higher up the chain
		foreach(j; 0..(2*M)){
		  ret.unsafePlusEq(partialProduct.data[j], j, 2*M*i);
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
  

  //add 1 more float's worth of precision to this
  //make this a template so it's only compiled when needed
  //otherwise this stack overflows the compiler!
  auto incrementExpansion()() const{
	static if(N == 0){
	  AdaptiveFloat!(FP, 0) ret;
	  ret.data = new FP[size() +1];
	  assert(0);
	} else {
	  AdaptiveFloat!(FP, N+1) ret; 
	}

	ret.data[1.. $] = data[];
	ret.data[0] = 0; //expanded by 1, so set new entry to 0;
	return ret;
  }

  auto expandToAdd(int M)(const ref AdaptiveFloat!(FP, M) other) const{
	static if(N==0 || M==0){
	  AdaptiveFloat!(FP, 0) ret;
	  ret.data = new FP[size() + other.size()];
	  assert(0);
	} else {
	  AdaptiveFloat!(FP, M+N) ret;
	}

	ret.data[other.size()  .. $] = data[];
	ret.data[0..other.size()] = 0;

	return ret;
  }


  auto expandToMultiply(int M)(const ref AdaptiveFloat!(FP, M) other) const{
	static if(N==0 || M==0){
	  AdaptiveFloat!(FP, 0) ret;
	  ret.data = new FP[size()*other.size()];
	  assert(0);
	} else {
	  AdaptiveFloat!(FP, M*N) ret;
	}

	ret.data[(ret.size() - size()) .. $] = data[];
	ret.data[0..(ret.size() - size())] = 0;
	return ret;
  }
  
  void unsafePlusEq(FP f, int first, int len){
	auto q = AdaptiveFloat!FP(f);
	foreach(i; first..(first + len)){
	  auto t1 = AdaptiveFloat!FP(data[N - i - 1]);
	  auto temp = q + t1;
	  q.data[0] = temp.data[0];
	  data[N - i - 1] = temp.data[1];
	}
	data[N - first - len -1] = q.data[0];
  }

  void unsafeMinusEq(FP f, int first, int len){
	auto q = AdaptiveFloat!FP(f);
	auto t1 = AdaptiveFloat!FP(data[N - first - 1]);
	auto temp = t1 - q;
	data[N - first - 1] = temp.data[1];
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


