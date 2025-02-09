#module SLATEC

export DRC, DRD, DRF, DRJ

#***BEGIN PROLOGUE  DRF
#***PURPOSE  Compute the incomplete or complete elliptic integral of the
#            1st kind.  For X, Y, and Z non-negative and at most one of
#            them zero, RF(X,Y,Z) = Integral from zero to infinity of
#                                -1/2     -1/2     -1/2
#                      (1/2)(t+X)    (t+Y)    (t+Z)    dt.
#            If X, Y or Z is zero, the integral is complete.
#***LIBRARY   SLATEC
#***CATEGORY  C14
#***TYPE      DOUBLE PRECISION (RF-S, DRF-D)
#***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
#             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE FIRST KIND,
#             TAYLOR SERIES
#***AUTHOR  Carlson, B. C.
#             Ames Laboratory-DOE
#             Iowa State University
#             Ames, IA  50011
#           Notis, E. M.
#             Ames Laboratory-DOE
#             Iowa State University
#             Ames, IA  50011
#           Pexton, R. L.
#             Lawrence Livermore National Laboratory
#             Livermore, CA  94550

function DRF(X::A, Y::B, Z::C) where {A,B,C}
    T = promote_type(A, B, C)

    ERRTOL = (4 * eps(T) / 2)^T(1 / 6)
    LOLIM = 5floatmin(T)
    UPLIM = floatmax(T) / 5
    C1 = T(1 / 24)
    C2 = T(3 / 44)
    C3 = T(1 / 14)

    ans = zero(T)
    min(X, Y, Z) < zero(T) && return (ans, 1)
    max(X, Y, Z) > UPLIM && return (ans, 3)
    min(X + Y, X + Z, Y + Z) < LOLIM && return (ans, 2)

    XN = X
    YN = Y
    ZN = Z
    MU = 0
    XNDEV = 0
    YNDEV = 0
    ZNDEV = 0

    while true
        MU = (XN + YN + ZN) / 3
        XNDEV = 2 - (MU + XN) / MU
        YNDEV = 2 - (MU + YN) / MU
        ZNDEV = 2 - (MU + ZN) / MU
        EPSLON = max(abs(XNDEV), abs(YNDEV), abs(ZNDEV))
        (EPSLON < ERRTOL) && break
        XNROOT = sqrt(XN)
        YNROOT = sqrt(YN)
        ZNROOT = sqrt(ZN)
        LAMDA = XNROOT * (YNROOT + ZNROOT) + YNROOT * ZNROOT
        XN = (XN + LAMDA) / 4
        YN = (YN + LAMDA) / 4
        ZN = (ZN + LAMDA) / 4
    end

    E2 = XNDEV * YNDEV - ZNDEV * ZNDEV
    E3 = XNDEV * YNDEV * ZNDEV
    S = 1 + (C1 * E2 - T(1 / 10) - C2 * E3) * E2 + C3 * E3
    ans = S / sqrt(MU)

    return (ans, 0)
end

#***BEGIN PROLOGUE  DRD
#***PURPOSE  Compute the incomplete or complete elliptic integral of
#            the 2nd kind. For X and Y nonnegative, X+Y and Z positive,
#            DRD(X,Y,Z) = Integral from zero to infinity of
#                                -1/2     -1/2     -3/2
#                      (3/2)(t+X)    (t+Y)    (t+Z)    dt.
#            If X or Y is zero, the integral is complete.
#***LIBRARY   SLATEC
#***CATEGORY  C14
#***TYPE      DOUBLE PRECISION (RD-S, DRD-D)
#***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
#             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE SECOND KIND,
#             TAYLOR SERIES
#***AUTHOR  Carlson, B. C.
#             Ames Laboratory-DOE
#             Iowa State University
#             Ames, IA  50011
#           Notis, E. M.
#             Ames Laboratory-DOE
#             Iowa State University
#             Ames, IA  50011
#           Pexton, R. L.
#             Lawrence Livermore National Laboratory
#             Livermore, CA  94550

function DRD(X::A, Y::B, Z::C) where {A,B,C}
    T = promote_type(A, B, C)

    ERRTOL = (eps(T) / 6)^T(1 / 6)
    LOLIM = 2 / (floatmax(T))^T(2 / 3)
    TUPLIM = floatmin(T)^T(1 / 3)
    TUPLIM = (ERRTOL / 10)^T(1 / 3) / TUPLIM
    UPLIM = TUPLIM^2
    C1 = T(3 / 14)
    C2 = T(1 / 6)
    C3 = T(9 / 22)
    C4 = T(3 / 26)

    ans = zero(T)
    min(X, Y) < 0 && return (ans, 1)
    max(X, Y, Z) > UPLIM && return (ans, 3)
    min(X + Y, Z) < LOLIM && return (ans, 2)

    XN = X
    YN = Y
    ZN = Z
    SIGMA = zero(T)
    POWER4 = one(T)
    MU = zero(T)
    XNDEV = zero(T)
    YNDEV = zero(T)
    ZNDEV = zero(T)

    while true
        MU = 2 * (XN + YN + 3 * ZN) / 10
        XNDEV = (MU - XN) / MU
        YNDEV = (MU - YN) / MU
        ZNDEV = (MU - ZN) / MU
        EPSLON = max(abs(XNDEV), abs(YNDEV), abs(ZNDEV))
        (EPSLON < ERRTOL) && break
        XNROOT = sqrt(XN)
        YNROOT = sqrt(YN)
        ZNROOT = sqrt(ZN)
        LAMDA = XNROOT * (YNROOT + ZNROOT) + YNROOT * ZNROOT
        SIGMA = SIGMA + POWER4 / (ZNROOT * (ZN + LAMDA))
        POWER4 = POWER4 / 4
        XN = (XN + LAMDA) / 4
        YN = (YN + LAMDA) / 4
        ZN = (ZN + LAMDA) / 4
    end

    EA = XNDEV * YNDEV
    EB = ZNDEV * ZNDEV
    EC = EA - EB
    ED = EA - 6 * EB
    EF = ED + EC + EC
    S1 = ED * (-C1 + C3 * ED / 4 - 3C4 * ZNDEV * EF / 2)
    S2 = ZNDEV * (C2 * EF + ZNDEV * (-C3 * EC + ZNDEV * C4 * EA))
    ans = 3 * SIGMA + POWER4 * (1 + S1 + S2) / (MU * sqrt(MU))

    return (ans, 0)
end

#***BEGIN PROLOGUE  DRC
#***PURPOSE  Calculate a double precision approximation to
#             DRC(X,Y) = Integral from zero to infinity of
#                              -1/2     -1
#                    (1/2)(t+X)    (t+Y)  dt,
#            where X is nonnegative and Y is positive.
#***LIBRARY   SLATEC
#***CATEGORY  C14
#***TYPE      DOUBLE PRECISION (RC-S, DRC-D)
#***KEYWORDS  DUPLICATION THEOREM, ELEMENTARY FUNCTIONS,
#             ELLIPTIC INTEGRAL, TAYLOR SERIES
#***AUTHOR  Carlson, B. C.
#             Ames Laboratory-DOE
#             Iowa State University
#             Ames, IA  50011
#           Notis, E. M.
#             Ames Laboratory-DOE
#             Iowa State University
#             Ames, IA  50011
#           Pexton, R. L.
#             Lawrence Livermore National Laboratory
#             Livermore, CA  94550

function DRC(X::A, Y::B) where {A,B}
    T = promote_type(A, B)

    ERRTOL = (eps(T) / 32)^T(1 / 6)
    LOLIM = T(5) * floatmin(T)
    UPLIM = T(floatmax(T) / 5)
    C1 = T(1 / 7)
    C2 = T(9 / 22)

    ans = zero(T)

    X < zero(T) || Y <= zero(T) && return (ans, 1)
    max(X, Y) > UPLIM && return (ans, 3)
    X + Y < LOLIM && return (ans, 2)

    XN = X
    YN = Y
    MU = zero(T)
    SN = zero(T)

    while true
        MU = (XN + YN + YN) / 3
        SN = (YN + MU) / MU - 2
        abs(SN) < ERRTOL && break
        LAMDA = 2 * sqrt(XN) * sqrt(YN) + YN
        XN = (XN + LAMDA) / 4
        YN = (YN + LAMDA) / 4
    end

    S = SN * SN * (T(3 / 10) + SN * (C1 + SN * (T(0.3750) + SN * C2)))
    ans = (1 + S) / sqrt(MU)

    return (ans, 0)
end

#***BEGIN PROLOGUE  DRJ
#***PURPOSE  Compute the incomplete or complete (X or Y or Z is zero)
#            elliptic integral of the 3rd kind.  For X, Y, and Z non-
#            negative, at most one of them zero, and P positive,
#             RJ(X,Y,Z,P) = Integral from zero to infinity of
#                              -1/2     -1/2     -1/2     -1
#                    (3/2)(t+X)    (t+Y)    (t+Z)    (t+P)  dt.
#***LIBRARY   SLATEC
#***CATEGORY  C14
#***TYPE      DOUBLE PRECISION (RJ-S, DRJ-D)
#***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
#             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE THIRD KIND,
#             TAYLOR SERIES
#***AUTHOR  Carlson, B. C.
#             Ames Laboratory-DOE
#             Iowa State University
#             Ames, IA  50011
#           Notis, E. M.
#             Ames Laboratory-DOE
#             Iowa State University
#             Ames, IA  50011
#           Pexton, R. L.
#             Lawrence Livermore National Laboratory
#             Livermore, CA  94550

function DRJ(X::A, Y::B, Z::C, P::D) where {A,B,C,D}
    T = promote_type(A, B, C, D)

    ERRTOL = (eps(T) / 6)^T(1 / 6)
    LOLIM = (5floatmin(T))^T(1 / 3)
    UPLIM = T(3 / 10) * (floatmax(T) / 5)^T(1 / 3)

    C1 = T(3 / 14)
    C2 = T(1 / 3)
    C3 = T(3 / 22)
    C4 = T(3 / 26)

    ans = zero(T)

    min(X, Y, Z) < zero(T) && return (ans, 1)
    max(X, Y, Z, P) > UPLIM && return (ans, 3)
    min(X + Y, X + Z, Y + Z, P) < LOLIM && return (ans, 2)

    IER = 0
    XN = X
    YN = Y
    ZN = Z
    PN = P
    SIGMA = zero(T)
    POWER4 = one(T)
    MU = zero(T)
    XNDEV = zero(T)
    YNDEV = zero(T)
    ZNDEV = zero(T)
    PNDEV = zero(T)

    while true
        MU = 2(XN + YN + ZN + PN + PN) / 10
        XNDEV = (MU - XN) / MU
        YNDEV = (MU - YN) / MU
        ZNDEV = (MU - ZN) / MU
        PNDEV = (MU - PN) / MU
        EPSLON = max(abs(XNDEV), abs(YNDEV), abs(ZNDEV), abs(PNDEV))
        EPSLON < ERRTOL && break
        XNROOT = sqrt(XN)
        YNROOT = sqrt(YN)
        ZNROOT = sqrt(ZN)
        LAMDA = XNROOT * (YNROOT + ZNROOT) + YNROOT * ZNROOT
        ALFA = PN * (XNROOT + YNROOT + ZNROOT) + XNROOT * YNROOT * ZNROOT
        ALFA = ALFA * ALFA
        BETA = PN * (PN + LAMDA) * (PN + LAMDA)
        drc, IER = DRC(ALFA, BETA)
        SIGMA = SIGMA + POWER4 * drc
        POWER4 = POWER4 / 4
        XN = (XN + LAMDA) / 4
        YN = (YN + LAMDA) / 4
        ZN = (ZN + LAMDA) / 4
        PN = (PN + LAMDA) / 4
    end

    EA = XNDEV * (YNDEV + ZNDEV) + YNDEV * ZNDEV
    EB = XNDEV * YNDEV * ZNDEV
    EC = PNDEV * PNDEV
    E2 = EA - 3 * EC
    E3 = EB + 2 * PNDEV * (EA - EC)
    S1 = 1 + E2 * (-C1 + 3C3 / 4 * E2 - 3C4 / 2 * E3)
    S2 = EB * (C2 / 2 + PNDEV * (-C3 - C3 + PNDEV * C4))
    S3 = PNDEV * EA * (C2 - PNDEV * C3) - C2 * PNDEV * EC
    ans = 3SIGMA + POWER4 * (S1 + S2 + S3) / (MU * sqrt(MU))

    return (ans, IER)
end

#end # module
