module dtriangulate.apFloat;

struct AdaptiveFloat(FP, int N = 1) if (N != 0){
  
private:
	FP[N] data;

  
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

  
  AdaptiveFloat!(FP, M + N) opBinary(string op, int M)(auto ref AdaptiveFloat!(FP, M) rhs) if(op == "+"){

	static if(M == 1 && N == 1){
	  //shewchuck 2-add
	  AdaptiveFloat!(FP, 2) ret;
	  ret.data[0] = data[0] + rhs.data[0];
	  FP bVirt = ret.data[0] - data[0];
	  FP aVirt = ret.data[0] - bVirt;
	  FP bRoundoff = rhs.data[0] - bVirt;
	  FP aRoundoff = data[0] - aVirt;
	  ret.data[1] = aRoundoff + bRoundoff;
	  return ret;
	} else static if(M ==1){
	  auto ret = expand!(N+1); //AdaptiveFloat!(FP, N + 1)(this);
	  ret.unsafePlusEq(rhs.data[0], 0, N);
	  return ret;
	} else static if(N ==1){
	  auto ret = AdaptiveFloat!(FP, M + 1)(rhs);
	  ret.unsafePlusEq(data[0], 0, M);
	  return ret;
	} else {
	  AdaptiveFloat!(FP, M + N) ret = expand!(M+N); //AdaptiveFloat!(FP, M + N)(this);
	  foreach(i; 0..M){
		ret.unsafePlusEq(rhs.data[M - i - 1], i, N);
	  }
	  return ret;
	}
	
  }

  AdaptiveFloat!(FP, M + N) opBinary(string op, int M)(auto ref AdaptiveFloat!(FP, M) rhs) if( op == "-"){
	static if(M == 1 && N == 1){
	  AdaptiveFloat!(FP, 2) ret;
	  ret.data[0] = data[0] - rhs.data[0];
	  FP bVirt = data[0] - ret.data[0];
	  FP aVirt = ret.data[0] + bVirt;
	  FP aRoundoff = data[0] - aVirt;
	  FP bRoundoff = bVirt - rhs.data[0];
	  ret.data[1] = aRoundoff + bRoundoff;
	  return ret;
	} else static if(M == 1){
	  auto ret = expand!(N+1); //AdaptiveFloat!(FP, N + 1)(this);
	  ret.unsafeMinusEq(rhs.data[0], 0, N);
	  return ret;
	} else static if(N == 1){
	  auto ret = AdaptiveFloat!(FP, M + 1)(rhs);
	  ret.unsafeMinusEq(data[0], 0, M);
	  return ret;
	} else {
	  auto ret = expand!(M+N);//AdaptiveFloat!(FP, M + N)(this);
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
	  
	  AdaptiveFloat!(FP, 2) ret;
	  ret.data[0] = data[0]*rhs.data[0];
	  auto aSplit = split();
	  auto bSplit = rhs.split();
	  FP err1 = ret.data[0] - (aSplit.data[0]*bSplit.data[0]);
	  FP err2 = err1 - (aSplit.data[1]*bSplit.data[0]);
	  FP err3 = err2 - (aSplit.data[0]*bSplit.data[1]);
	  ret.data[1] = (aSplit.data[1]*bSplit.data[1]) - err3;
	  return ret;
	  
	} else static if(M == 1){ //RHS is one float
	  
	  AdaptiveFloat!(FP, 2*N) ret;
	  auto temp = AdaptiveFloat!FP(data[N - 1])*rhs; //T.N == 2

	  ret.data[2*N -1] = temp.data[1];

	  foreach(i ; 1..N){
		auto temp2 = AdaptiveFloat!FP(data[N - 1 - i])*rhs; //T2.N == 2

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
	} else {
	  static if(M < N){
		return rhs*this;
	  } else {
		//N <= M
		AdaptiveFloat!(FP, 2*M*N) ret = (rhs*AdaptiveFloat!FP(data[0])).expand!(2*M*N);

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

  AdaptiveFloat!(FP, M) expand(int M)() if(M > 1 && M >= N){
	AdaptiveFloat!(FP, M) ret;
	foreach(i; 0..N){
	  ret.data[M - i - 1] = data[N - i - 1];
	}
	foreach(i; 0..(M - N)){
	  ret.data[i] = 0;
	}
	return ret;
  }
  /*
  this(int M)(auto ref AdaptiveFloat!(FP, M) rhs) if(M > 1 && M <= N){
	foreach(i; 0..M){
	  data[N - i - 1] = rhs.data[M - i - 1];
	}
	foreach(i; 0..(N - M)){
	  data[i] = 0;
	}
	}*/
  
  void unsafePlusEq(FP f, int first, int len){
	//import std.stdio;
	//writefln("unsafe peq, before, %f, %d, %d, %s", f, first, len, this);
	auto q = AdaptiveFloat!FP(f);
	foreach(i; first..(first + len)){
	  auto t1 = AdaptiveFloat!FP(data[N - i - 1]);
	  auto temp = q + t1;
	  q.data[0] = temp.data[0];
	  data[N - i - 1] = temp.data[1];
	}
	data[N - first - len -1] = q.data[0];
	//writefln("after: %s", this);
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
  /*  static if(N > 1){
	AdaptiveFloat!(FP, N -1) dropLast(){
	  import std.conv;
	  pragma(msg, "dropLast with N = " ~ to!string(N));
	  AdaptiveFloat!(FP, N -1) ret;
	  foreach(i ; 0..(N -1)){
		ret.data[i] = data[i];
	  }
	  return ret;
	}
	}*/
  
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
  //writeln("eps " ~ to!string(eps));
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
