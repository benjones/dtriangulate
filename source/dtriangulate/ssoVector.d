module dtriangulate.ssoVector;

struct SSOVector(T, ulong N) if(__traits(isPOD,T)){

private:
  byte[N*T.sizeof] localStorage;

public:
  T[] data;  //possible unsafe? Should probably disable some ops?
  alias data this;
  
  this(ulong n){
	if(n < N){
	  T* dataPtr = cast(T*)localStorage.ptr;
	  data = dataPtr[0..n];
	} else {
	  data = new T[n];
	}
  }
  
  void opOpAssign(string op)(auto ref T t) if(op == "~"){
	if(data.length < N){
	  T* dataPtr = cast(T*)localStorage.ptr;
	  data = dataPtr[0..data.length + 1];
	  data[$-1] = t;
	} else {

	  data ~= t;
	}
  }

  @property
  ref T back(){ return data[$-1]; }

  @property bool empty() const { return data.length == 0; }
  
  void popBack(){
	data = data[0..$-1];
  }
}


unittest {
  import std.stdio;
  
  SSOVector!(int, 4) i4;

  assert(i4.length == 0);
  assert(i4.empty());
  
  foreach(i; 0..10){
	i4 ~= i;
  }
  assert(!i4.empty);
  foreach(i, j; i4){
	assert(i == j);
  }

  foreach(ref i; i4){
	i *= 2;
  }

  foreach(i, j ; i4){
	assert(i*2 == j);
  }

  i4.popBack();
  assert(i4.length == 9);
  assert(i4.back == 16);
  
}
