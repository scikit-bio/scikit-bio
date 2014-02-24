#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The bipy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import division

from numpy import sqrt, log, exp, abs

ROUND_ERROR = 1e-14    # fp rounding error: causes some tests to fail
                        # will round to 0 if smaller in magnitude than this
MAXNUM = 1.7976931348623158E308  # 2**1024
SQRT2 = 1.41421356237309504880  # sqrt(2)
SQRTH = 7.07106781186547524401E-1  # sqrt(2)/2
MINLOG = -7.08396418532264106224E2  # log(2**-1022)
LOGPI = 1.14472988584940017414
MAXLGM = 2.556348e305
MINLOG = -7.08396418532264106224E2
MAXLOG = 7.09782712893383996843E2
MAXGAM = 171.624376956302725
LS2PI = 0.91893853320467274178
MACHEP = 1.11022302462515654042E-16
MAXSTIR = 143.01608
SQTPI = 2.50662827463100050242E0
SQRTH = 7.07106781186547524401E-1

PI = 3.14159265358979323846  # pi

exp_minus_2 = 0.13533528323661269189
big = 4.503599627370496e15
biginv = 2.22044604925031308085e-16


# Coefficients for Gamma follow:
GA = [
    8.11614167470508450300E-4,
    -5.95061904284301438324E-4,
    7.93650340457716943945E-4,
    -2.77777777730099687205E-3,
    8.33333333333331927722E-2,
]

GB = [
    -1.37825152569120859100E3,
    -3.88016315134637840924E4,
    -3.31612992738871184744E5,
    -1.16237097492762307383E6,
    -1.72173700820839662146E6,
    -8.53555664245765465627E5,
]

GC = [
    1.00000000000000000000E0,
    -3.51815701436523470549E2,
    -1.70642106651881159223E4,
    -2.20528590553854454839E5,
    -1.13933444367982507207E6,
    -2.53252307177582951285E6,
    -2.01889141433532773231E6,
]

GP = [
    1.60119522476751861407E-4,
    1.19135147006586384913E-3,
    1.04213797561761569935E-2,
    4.76367800457137231464E-2,
    2.07448227648435975150E-1,
    4.94214826801497100753E-1,
    9.99999999999999996796E-1,
]

GQ = [
    -2.31581873324120129819E-5,
    5.39605580493303397842E-4,
    -4.45641913851797240494E-3,
    1.18139785222060435552E-2,
    3.58236398605498653373E-2,
    -2.34591795718243348568E-1,
    7.14304917030273074085E-2,
    1.00000000000000000320E0,
]

STIR = [
    7.87311395793093628397E-4,
    -2.29549961613378126380E-4,
    -2.68132617805781232825E-3,
    3.47222221605458667310E-3,
    8.33333333333482257126E-2,
]

LP = [
    4.5270000862445199635215E-5,
    4.9854102823193375972212E-1,
    6.5787325942061044846969E0,
    2.9911919328553073277375E1,
    6.0949667980987787057556E1,
    5.7112963590585538103336E1,
    2.0039553499201281259648E1,
]
LQ = [
    1,
    1.5062909083469192043167E1,
    8.3047565967967209469434E1,
    2.2176239823732856465394E2,
    3.0909872225312059774938E2,
    2.1642788614495947685003E2,
    6.0118660497603843919306E1,
]
EP = [
    1.2617719307481059087798E-4,
    3.0299440770744196129956E-2,
    9.9999999999999999991025E-1,
]

EQ = [
    3.0019850513866445504159E-6,
    2.5244834034968410419224E-3,
    2.2726554820815502876593E-1,
    2.0000000000000000000897E0,
]

P0 = [
    -5.99633501014107895267E1,
    9.80010754185999661536E1,
    -5.66762857469070293439E1,
    1.39312609387279679503E1,
    -1.23916583867381258016E0,
]

Q0 = [
    1.00000000000000000000E0,
    1.95448858338141759834E0,
    4.67627912898881538453E0,
    8.63602421390890590575E1,
    -2.25462687854119370527E2,
    2.00260212380060660359E2,
    -8.20372256168333339912E1,
    1.59056225126211695515E1,
    -1.18331621121330003142E0,
]

s2pi = 2.50662827463100050242E0

P1 = [
    4.05544892305962419923E0,
    3.15251094599893866154E1,
    5.71628192246421288162E1,
    4.40805073893200834700E1,
    1.46849561928858024014E1,
    2.18663306850790267539E0,
    -1.40256079171354495875E-1,
    -3.50424626827848203418E-2,
    -8.57456785154685413611E-4,
]

Q1 = [
    1.00000000000000000000E0,
    1.57799883256466749731E1,
    4.53907635128879210584E1,
    4.13172038254672030440E1,
    1.50425385692907503408E1,
    2.50464946208309415979E0,
    -1.42182922854787788574E-1,
    -3.80806407691578277194E-2,
    -9.33259480895457427372E-4,
]

P2 = [
    3.23774891776946035970E0,
    6.91522889068984211695E0,
    3.93881025292474443415E0,
    1.33303460815807542389E0,
    2.01485389549179081538E-1,
    1.23716634817820021358E-2,
    3.01581553508235416007E-4,
    2.65806974686737550832E-6,
    6.23974539184983293730E-9,
]

Q2 = [
    1.00000000000000000000E0,
    6.02427039364742014255E0,
    3.67983563856160859403E0,
    1.37702099489081330271E0,
    2.16236993594496635890E-1,
    1.34204006088543189037E-2,
    3.28014464682127739104E-4,
    2.89247864745380683936E-6,
    6.79019408009981274425E-9,
]


# Translations of functions from Cephes Math Library, by Stephen L. Moshier
def polevl(x, coef):
    """evaluates a polynomial y = C_0 + C_1x + C_2x^2 + ... + C_Nx^N

    Coefficients are stored in reverse order, i.e. coef[0] = C_N
    """
    result = 0
    for c in coef:
        result = result * x + c
    return result


def expm1(x):
    """Something to do with exp? From Cephes."""
    if (x < -0.5) or (x > 0.5):
        return (exp(x) - 1.0)
    xx = x * x
    r = x * polevl(xx, EP)
    r /= polevl(xx, EQ) - r
    return r + r


def log1p(x):
    """Log for values close to 1: from Cephes math library"""
    z = 1 + x
    if (z < SQRTH) or (z > SQRT2):
        return log(z)
    z = x * x
    z = -0.5 * z + x * (z * polevl(x, LP) / polevl(x, LQ))
    return x + z


def igam(a, x):
    """Left tail of incomplete gamma function: see Cephes docs for details"""
    if x <= 0 or a <= 0:
        return 0
    if x > 1 and x > a:
        return 1 - igamc(a, x)

    # Compute x**a * exp(x) / Gamma(a)

    ax = a * log(x) - x - lgam(a)
    if ax < -MAXLOG:  # underflow
        return 0.0
    ax = exp(ax)

    # power series
    r = a
    c = 1
    ans = 1
    while True:
        r += 1
        c *= x / r
        ans += c
        if c / ans <= MACHEP:
            break

    return ans * ax / a


def lgam(x):
    """Natural log of the gamma fuction: see Cephes docs for details"""
    sgngam = 1
    if x < -34:
        q = -x
        w = lgam(q)
        p = floor(q)
        if p == q:
            raise OverflowError("lgam returned infinity.")
        i = p
        if i & 1 == 0:
            sgngam = -1
        else:
            sgngam = 1
        z = q - p
        if z > 0.5:
            p += 1
            z = p - q
        z = q * sin(PI * z)
        if z == 0:
            raise OverflowError("lgam returned infinity.")
        z = LOGPI - log(z) - w
        return z
    if x < 13:
        z = 1
        p = 0
        u = x
        while u >= 3:
            p -= 1
            u = x + p
            z *= u
        while u < 2:
            if u == 0:
                raise OverflowError("lgam returned infinity.")
            z /= u
            p += 1
            u = x + p
        if z < 0:
            sgngam = -1
            z = -z
        else:
            sgngam = 1
        if u == 2:
            return log(z)
        p -= 2
        x = x + p
        p = x * polevl(x, GB) / polevl(x, GC)
        return log(z) + p
    if x > MAXLGM:
        raise OverflowError("Too large a value of x in lgam.")
    q = (x - 0.5) * log(x) - x + LS2PI
    if x > 1.0e8:
        return q
    p = 1 / (x * x)
    if x >= 1000:
        q += ((7.9365079365079365079365e-4 * p
               - 2.7777777777777777777778e-3) * p
              + 0.0833333333333333333333) / x
    else:
        q += polevl(p, GA) / x
    return q


def igamc(a, x):
    """Complemented incomplete Gamma integral: see Cephes docs."""
    if x <= 0 or a <= 0:
        return 1
    if x < 1 or x < a:
        return 1 - igam(a, x)
    ax = a * log(x) - x - lgam(a)
    if ax < -MAXLOG:  # underflow
        return 0
    ax = exp(ax)
    # continued fraction
    y = 1 - a
    z = x + y + 1
    c = 0
    pkm2 = 1
    qkm2 = x
    pkm1 = x + 1
    qkm1 = z * x
    ans = pkm1 / qkm1

    while True:
        c += 1
        y += 1
        z += 2
        yc = y * c
        pk = pkm1 * z - pkm2 * yc
        qk = qkm1 * z - qkm2 * yc
        if qk != 0:
            r = pk / qk
            t = abs((ans - r) / r)
            ans = r
        else:
            t = 1
        pkm2 = pkm1
        pkm1 = pk
        qkm2 = qkm1
        qkm1 = qk
        if abs(pk) > big:
            pkm2 *= biginv
            pkm1 *= biginv
            qkm2 *= biginv
            qkm1 *= biginv
        if t <= MACHEP:
            break
    return ans * ax


def fix_rounding_error(x):
    """If x is almost in the range 0-1, fixes it.

    Specifically, if x is between -ROUND_ERROR and 0, returns 0.
    If x is between 1 and 1+ROUND_ERROR, returns 1.
    """
    if -ROUND_ERROR < x < 0:
        return 0
    elif 1 < x < 1 + ROUND_ERROR:
        return 1
    else:
        return x


def ndtri(y0):
    """Inverse normal distribution function.

    This is here and not in distributions because igami depends on it..."""
    y0 = fix_rounding_error(y0)
    # handle easy cases
    if y0 <= 0.0:
        return -MAXNUM
    elif y0 >= 1.0:
        return MAXNUM
    code = 1
    y = y0
    if y > (1.0 - exp_minus_2):
        y = 1.0 - y
        code = 0

    if y > exp_minus_2:
        y -= 0.5
        y2 = y * y
        x = y + y * (y2 * polevl(y2, P0) / polevl(y2, Q0))
        x = x * s2pi
        return x

    x = sqrt(-2.0 * log(y))
    x0 = x - log(x) / x

    z = 1.0 / x
    if x < 8.0:  # y > exp(-32) = 1.2664165549e-14
        x1 = z * polevl(z, P1) / polevl(z, Q1)
    else:
        x1 = z * polevl(z, P2) / polevl(z, Q2)
    x = x0 - x1
    if code != 0:
        x = -x
    return x


def betai(aa, bb, xx):
    """Returns integral of the incomplete beta density function, from 0 to x.

    See Cephes docs for details.
    """
    if aa <= 0 or bb <= 0:
        raise ValueError("betai: a and b must both be > 0.")
    if xx == 0:
        return 0
    if xx == 1:
        return 1
    if xx < 0 or xx > 1:
        raise ValueError("betai: x must be between 0 and 1.")
    flag = 0
    if (bb * xx <= 1) and (xx <= 0.95):
        t = pseries(aa, bb, xx)
        return betai_result(t, flag)
    w = 1 - xx
    # reverse a and b if x is greater than the mean
    if xx > (aa / (aa + bb)):
        flag = 1
        a = bb
        b = aa
        xc = xx
        x = w
    else:
        a = aa
        b = bb
        xc = w
        x = xx
    if (flag == 1) and ((b * x) <= 1) and (x <= 0.95):
        t = pseries(a, b, x)
        return betai_result(t, flag)
    # choose expansion for better convergence
    y = x * (a + b - 2) - (a - 1)
    if y < 0:
        w = incbcf(a, b, x)
    else:
        w = incbd(a, b, x) / xc
    y = a * log(x)
    t = b * log(xc)
    if ((a + b) < MAXGAM) and (abs(y) < MAXLOG) and (abs(t) < MAXLOG):
        t = pow(xc, b)
        t *= pow(x, a)
        t /= a
        t *= w
        t *= Gamma(a + b) / (Gamma(a) * Gamma(b))
        return betai_result(t, flag)
    # resort to logarithms
    y += t + lgam(a + b) - lgam(a) - lgam(b)
    y += log(w / a)
    if y < MINLOG:
        t = 0
    else:
        t = exp(y)
    return betai_result(t, flag)


def betai_result(t, flag):
    if flag == 1:
        if t <= MACHEP:
            t = 1 - MACHEP
        else:
            t = 1 - t
    return t


def pseries(a, b, x):
    """Power series for incomplete beta integral.

    Use when b * x is small and x not too close to 1.

    See Cephes docs for details.
    """
    ai = 1 / a
    u = (1 - b) * x
    v = u / (a + 1)
    t1 = v
    t = u
    n = 2
    s = 0
    z = MACHEP * ai
    while abs(v) > z:
        u = (n - b) * x / n
        t *= u
        v = t / (a + n)
        s += v
        n += 1
    s += t1
    s += ai

    u = a * log(x)
    if ((a + b) < MAXGAM) and (abs(u) < MAXLOG):
        t = Gamma(a + b) / (Gamma(a) * Gamma(b))
        s = s * t * pow(x, a)
    else:
        t = lgam(a + b) - lgam(a) - lgam(b) + u + log(s)
        if t < MINLOG:
            s = 0
        else:
            s = exp(t)
    return(s)


def incbd(a, b, x):
    """Incomplete beta integral, second continued fraction representation.

    See Cephes docs for details."""

    k1 = a
    k2 = b - 1.0
    k3 = a
    k4 = a + 1.0
    k5 = 1.0
    k6 = a + b
    k7 = a + 1.0
    k8 = a + 2.0

    pkm2 = 0.0
    qkm2 = 1.0
    pkm1 = 1.0
    qkm1 = 1.0
    z = x / (1.0 - x)
    ans = 1.0
    r = 1.0
    n = 0
    thresh = 3 * MACHEP

    while True:
        xk = - (z * k1 * k2) / (k3 * k4)
        pk = pkm1 + pkm2 * xk
        qk = qkm1 + qkm2 * xk
        pkm2 = pkm1
        pkm1 = pk
        qkm2 = qkm1
        qkm1 = qk

        xk = (z * k5 * k6) / (k7 * k8)
        pk = pkm1 + pkm2 * xk
        qk = qkm1 + qkm2 * xk
        pkm2 = pkm1
        pkm1 = pk
        qkm2 = qkm1
        qkm1 = qk

        if qk != 0:
            r = pk / qk
        if r != 0:
            t = abs((ans - r) / r)
            ans = r
        else:
            t = 1.0
        if t < thresh:
            return ans
        k1 += 1
        k2 -= 1
        k3 += 2
        k4 += 2
        k5 += 1
        k6 += 1
        k7 += 2
        k8 += 2

        if (abs(qk) + abs(pk)) > big:
            pkm2 *= biginv
            pkm1 *= biginv
            qkm2 *= biginv
            qkm1 *= biginv

        if (abs(qk) < biginv) or (abs(pk) < biginv):
            pkm2 *= big
            pkm1 *= big
            qkm2 *= big
            qkm1 *= big
        n += 1
        if n >= 300:
            return ans


def Gamma(x):
    """Returns the gamma function, a generalization of the factorial.

    See Cephes docs for details."""

    sgngam = 1
    q = abs(x)
    if q > 33:
        if x < 0:
            p = floor(q)
            if p == q:
                raise OverflowError("Bad value of x in Gamma function.")
            i = p
            if (i & 1) == 0:
                sgngam = -1
            z = q - p
            if z > 0.5:
                p += 1
                z = q - p
            z = q * sin(PI * z)
            if z == 0:
                raise OverflowError("Bad value of x in Gamma function.")
            z = abs(z)
            z = PI / (z * stirf(q))
        else:
            z = stirf(x)
        return sgngam * z
    z = 1
    while x >= 3:
        x -= 1
        z *= x
    while x < 0:
        if x > -1e9:
            return Gamma_small(x, z)
    while x < 2:
        if x < 1e-9:
            return Gamma_small(x, z)
        z /= x
        x += 1
    if x == 2:
        return float(z)
    x -= 2
    p = polevl(x, GP)
    q = polevl(x, GQ)
    return z * p / q


def Gamma_small(x, z):
    if x == 0:
        raise OverflowError("Bad value of x in Gamma function.")
    else:
        return z / ((1 + 0.5772156649015329 * x) * x)


def stirf(x):
    """Stirling's approximation for the Gamma function.

    Valid for 33 <= x <= 162.

    See Cephes docs for details.
    """
    w = 1.0 / x
    w = 1 + w * polevl(w, STIR)
    y = exp(x)
    if x > MAXSTIR:
        # avoid overflow in pow()
        v = pow(x, 0.5 * x - 0.25)
        y = v * (v / y)
    else:
        y = pow(x, x - 0.5) / y
    return SQTPI * y * w


def incbcf(a, b, x):
    """Incomplete beta integral, first continued fraction representation.

    See Cephes docs for details."""

    k1 = a
    k2 = a + b
    k3 = a
    k4 = a + 1
    k5 = 1
    k6 = b - 1
    k7 = k4
    k8 = a + 2

    pkm2 = 0
    qkm2 = 1
    pkm1 = 1
    qkm1 = 1
    ans = 1
    r = 1
    n = 0
    thresh = 3 * MACHEP

    while True:
        xk = -(x * k1 * k2) / (k3 * k4)
        pk = pkm1 + pkm2 * xk
        qk = qkm1 + qkm2 * xk
        pkm2 = pkm1
        pkm1 = pk
        qkm2 = qkm1
        qkm1 = qk

        xk = (x * k5 * k6) / (k7 * k8)
        pk = pkm1 + pkm2 * xk
        qk = qkm1 + qkm2 * xk
        pkm2 = pkm1
        pkm1 = pk
        qkm2 = qkm1
        qkm1 = qk

        if qk != 0:
            r = pk / qk
        if r != 0:
            t = abs((ans - r) / r)
            ans = r
        else:
            t = 1
        if t < thresh:
            return ans
        k1 += 1
        k2 += 1
        k3 += 2
        k4 += 2
        k5 += 1
        k6 -= 1
        k7 += 2
        k8 += 2

        if (abs(qk) + abs(pk)) > big:
            pkm2 *= biginv
            pkm1 *= biginv
            qkm2 *= biginv
            qkm1 *= biginv
        if (abs(qk) < biginv) or (abs(pk) < biginv):
            pkm2 *= big
            pkm1 *= big
            qkm2 *= big
            qkm1 *= big
        n += 1
        if n >= 300:
            return ans

# Coefficients for zdist follow:
ZP = [
    2.46196981473530512524E-10,
    5.64189564831068821977E-1,
    7.46321056442269912687E0,
    4.86371970985681366614E1,
    1.96520832956077098242E2,
    5.26445194995477358631E2,
    9.34528527171957607540E2,
    1.02755188689515710272E3,
    5.57535335369399327526E2,
]

ZQ = [
    1.0,
    1.32281951154744992508E1,
    8.67072140885989742329E1,
    3.54937778887819891062E2,
    9.75708501743205489753E2,
    1.82390916687909736289E3,
    2.24633760818710981792E3,
    1.65666309194161350182E3,
    5.57535340817727675546E2,
]

ZR = [
    5.64189583547755073984E-1,
    1.27536670759978104416E0,
    5.01905042251180477414E0,
    6.16021097993053585195E0,
    7.40974269950448939160E0,
    2.97886665372100240670E0,
]
ZS = [
    1.00000000000000000000E0,
    2.26052863220117276590E0,
    9.39603524938001434673E0,
    1.20489539808096656605E1,
    1.70814450747565897222E1,
    9.60896809063285878198E0,
    3.36907645100081516050E0,
]
ZT = [
    9.60497373987051638749E0,
    9.00260197203842689217E1,
    2.23200534594684319226E3,
    7.00332514112805075473E3,
    5.55923013010394962768E4,
]
ZU = [
    1.00000000000000000000E0,
    3.35617141647503099647E1,
    5.21357949780152679795E2,
    4.59432382970980127987E3,
    2.26290000613890934246E4,
    4.92673942608635921086E4,
]


def erf(a):
    """Returns the error function of a: see Cephes docs."""
    if abs(a) > 1:
        return 1 - erfc(a)
    z = a * a
    return a * polevl(z, ZT) / polevl(z, ZU)


def erfc(a):
    """Returns the complement of the error function of a: see Cephes docs."""
    if a < 0:
        x = -a
    else:
        x = a

    if x < 1:
        return 1 - erf(a)

    z = -a * a
    if z < -MAXLOG:  # underflow
        if a < 0:
            return 2
        else:
            return 0
    z = exp(z)

    if x < 8:
        p = polevl(x, ZP)
        q = polevl(x, ZQ)
    else:
        p = polevl(x, ZR)
        q = polevl(x, ZS)

    y = z * p / q

    if a < 0:
        y = 2 - y

    if y == 0:  # underflow
        if a < 0:
            return 2
        else:
            return 0
    else:
        return y
