module dtriangulate.ssoVector;

//we won't call destructors
struct SSOVector(T, ulong N) if(__traits(isPOD,T)){
    @safe:
private:
  T[N] localStorage;

public:
  T[] data;  //possible unsafe? Should probably disable some ops?
  alias data this;
  
  @trusted this(ulong n){
      if(n < N){
          T* dataPtr = &localStorage[0];
          data = dataPtr[0..n];
      } else {
          data = new T[n];
      }
  }
  
  @trusted void opOpAssign(string op)(auto ref T t) if(op == "~"){
	T* dataPtr = localStorage.ptr;
	if(data.length == 0 || (dataPtr == data.ptr && data.length < N)){
	  //use the built in array
	  data = dataPtr[0..(data.length + 1)];
	  data[$-1] = t;
	} else {
	  data ~= t;
	}
  }

  ref T back() return{ return data[$-1]; }

  @property bool empty() const { return data.length == 0; }
  
  void popBack(){
	data = data[0..$-1];
  }
}


@safe unittest {
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


  SSOVector!(int, 5) i5;
  foreach(i; 0..10){
	i5 ~= 1;
  }
  foreach(ref i ; i5){
	i = 0;
  }
  foreach(i; 0..8){
	i5.popBack();
  }
  i5 ~= 0;
  foreach(i ; i5){
	assert(i5.back == 0);
  }
}
