module dtriangulate.fma;

version(DigitalMars){

    ///returns round(a*b + c)  -- computed as if in infinite precision, rounded at the end
    double fma(double a, double b, double c) @safe pure @nogc nothrow{
        asm @safe pure @nogc nothrow {
            naked;
            vfmadd231sd XMM0, XMM1, XMM2;
            ret;
        }
    }

    ///ditto
    float fma(float a, float b, float c) @safe pure @nogc nothrow{
        asm @safe pure @nogc nothrow {
            naked;
            vfmadd231ss XMM0, XMM1, XMM2;
            ret;
        }
    }

    ///returns round(a*b + c)  -- computed as if in infinite precision, rounded at the end
    double fms(double a, double b, double c) @safe pure @nogc nothrow{
        asm @safe pure @nogc nothrow {
            naked;
            vfmsub231sd XMM0, XMM1, XMM2;
            ret;
        }
    }

    ///ditto
    float fms(float a, float b, float c) @safe pure @nogc nothrow{
        asm @safe pure @nogc nothrow {
            naked;
            vfmsub231ss XMM0, XMM1, XMM2;
            ret;
        }
    }
}

version(LDC){
    import ldc.intrinsics;
    import std.traits : isFloatingPoint;
    alias fma = llvm_fma;
    
    T fms(T)(T a, T b, T c) if(isFloatingPoint!T) {
        return fma(a, b, -c);
    }
}

version(GNU){
    static assert(0, "GDC FMA not supported yet");
}

@safe unittest {
    assert(fma(1.0, 2.0, 3.0) == 5);
    assert(fma(10.0, 5.0, 7.0) == 57);

    double a = 1 + double.epsilon;
    double b = 1 + 2*double.epsilon;

    assert(fma(a, a, -b) == double.epsilon*double.epsilon);


    assert(fms(2.0, 3.0, 4.0) == 2);
    assert(fms(10.0, 5.0, 7.0) == 43);

    assert(fms(a, a, b) == double.epsilon*double.epsilon);

    float c = 1 + float.epsilon;
    float d = 1 + 2*float.epsilon;

    assert(fma(1.0f, 2.0f, 3.0f) == 5.0f);
    assert(fms(2.0f, 3.0f, 4.0f) == 2.0f);

    assert(fma(c, c, -d) == float.epsilon*float.epsilon);
    assert(fms(c, c, d) == float.epsilon*float.epsilon);
}


