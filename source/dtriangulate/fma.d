module dtriangulate.fma;

///returns round(a*b + c)  -- computed as if in infinite precision, rounded at the end
double fma(double a, double b, double c) @safe pure @nogc nothrow{
    asm @safe pure @nogc nothrow {
        naked;
        vfmadd231sd XMM0, XMM1, XMM2;
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


@safe unittest {
    assert(fma(1, 2, 3) == 5);
    assert(fma(10, 5, 7) == 57);

    double a = 1 + double.epsilon;
    double b = 1 + 2*double.epsilon;

    assert(fma(a, a, -b) == double.epsilon*double.epsilon);


    assert(fms(2, 3, 4) == 2);
    assert(fms(10, 5, 7) == 43);

    assert(fms(a, a, b) == double.epsilon*double.epsilon);
}


