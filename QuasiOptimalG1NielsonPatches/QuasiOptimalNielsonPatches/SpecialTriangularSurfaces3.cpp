#include "./SpecialTriangularSurfaces3.h"

using namespace cagd;
using namespace std;

GLdouble fact1(GLuint n)
{
    GLdouble res = 1.0;
    if (n > 0){
        for (GLuint i = 1; i <= n; i++)
            res = res * i;
    }
    return res;
}

/*CUBIC POLINOMIAL SURFACE*/
CubicPolinomialSurface::CubicPolinomialSurface():TriangularSurface3(){}

CubicPolinomialSurface::CubicPolinomialSurface(const CubicPolinomialSurface& cpc):TriangularSurface3(cpc){}

CubicPolinomialSurface& CubicPolinomialSurface::operator =(const CubicPolinomialSurface& cpc)
{
    if (this != &cpc) {
        TriangularSurface3::operator=(cpc);
    }
    return *this;
}

GLboolean CubicPolinomialSurface::CalculatePartialDerivatives(GLuint maximum_order_of_partial_derivatives,
                GLdouble u, GLdouble v, PartialDerivatives& pd) const
{
    if (u < 0.0 || u > _beta || v < 0.0 || v > _beta) // u + v < _beta?
    {
        pd.ResizeRows(0);
        return GL_FALSE;
    }

    GLdouble w = 1 - u - v, w2 = w * w, w3 = w2 * w, u2 = u * u, u3 = u2 * u, v2 = v * v, v3 = v2 * v;

    TriangularMatrix<GLdouble> T(4), duT(4), dvT(4);

    T(3,3) = w3;
    T(3,2) = 3 * v * w2;
    T(3,1) = 3 * v2 * w;
    T(3,0) = v3;

    T(2,2) = 3 * u * w2;
    T(2,1) = 6 * u * v * w;
    T(2,0) = 3 * u * v2;

    T(1,1) = 3 * u2 * w;
    T(1,0) = 3 * u2 * v;

    T(0,0) = u3;

    duT(3,3) = -3 * w2;
    duT(3,2) = -6 * v * w;
    duT(3,1) = -3 * v2;
    duT(3,0) = 0;

    duT(2,2) = 3 * (w2 - 2 * u * w);
    duT(2,1) = 6 * v * (w - u);
    duT(2,0) = 3 * v2;

    duT(1,1) = 3 * u * (2 * w - u);
    duT(1,0) = 6 * u * v;

    duT(0,0) = 3 * u2;

    dvT(3,3) = -3 * w2;
    dvT(3,2) = 3 * (w2 - 2 * v * w);
    dvT(3,1) = 3 * v * ( 2 * w - v);
    dvT(3,0) = 3 * v2;

    dvT(2,2) = -6 * u * w;
    dvT(2,1) = 6 * u * (w - v);
    dvT(2,0) = 6 * u * v;

    dvT(1,1) = -3 * u2;
    dvT(1,0) = 3 * u2;

    dvT(0,0) = 0;

    // eq. (25)
    pd.ResizeRows(maximum_order_of_partial_derivatives + 1);
    pd.LoadNullVectors();

    DCoordinate3 &point = pd(0, 0), &diff1u = pd(1, 0), &diff1v = pd(1, 1);

    for (GLuint k = 0; k < 4; ++k)
    {
        for (GLuint l = 0; l <= k; ++l)
        {
            point  += _data(k, l) *   T(k, l);
            diff1u += _data(k, l) * duT(k, l);
            diff1v += _data(k, l) * dvT(k, l);
        }
    }

    return GL_TRUE;
}

/*QUARTIC POLINOMIAL SURFACE*/
QuarticPolinomialSurface::QuarticPolinomialSurface():TriangularSurface3(){}

QuarticPolinomialSurface::QuarticPolinomialSurface(const QuarticPolinomialSurface& cpc):TriangularSurface3(cpc){}

QuarticPolinomialSurface& QuarticPolinomialSurface::operator =(const QuarticPolinomialSurface& cpc)
{
    if (this != &cpc) {
        TriangularSurface3::operator=(cpc);
    }
    return *this;
}

GLboolean QuarticPolinomialSurface::CalculatePartialDerivatives(GLuint maximum_order_of_partial_derivatives,
                                                                GLdouble u, GLdouble v, PartialDerivatives &pd) const
{
    if (u < 0.0 || u > _beta || v < 0.0 || v > _beta) // u + v > _beta?
    {
        pd.ResizeRows(0);
        return GL_FALSE;
    }

    pd.ResizeRows(maximum_order_of_partial_derivatives + 1);
    pd.LoadNullVectors();

    GLdouble w = 1 - u -v, w2 = w * w, w3 = w2 * w, w4 = w3 * w, u2 = u * u, u3 = u2 * u, u4 = u3 * u,
            v2 = v * v, v3 = v2 * v, v4 = v3 * v;

    TriangularMatrix<GLdouble> T(4), duT(4), dvT(4);

    T(3,3) = w4;
    T(3,2) = 4.0 * v * w3 + 3.0 * v2 * w2;
    T(3,1) = 3.0 * v2 * w2 + 4.0 * v3 * w;
    T(3,0) = v4;

    T(2,2) = 4.0 * u * w3 + 3.0 * u2 * w2;
    T(2,1) = 12.0 * (u * v * w2 + u * v2 * w + u2 * v * w);
    T(2,0) = 4.0 * u * v3 + 3.0 * u2 * v2;

    T(1,1) = 3.0 * u2 * w2 + 4.0 * u3 * w;
    T(1,0) = 3.0 * u2 * v2 + 4.0 * u3 * v;

    T(0,0) = u4;

    duT(3,3) = -4.0 * w3;
    duT(3,2) = -12.0 * v * w2 - 6.0 * v2 * w;
    duT(3,1) = -6.0 * v2 * w - 4.0 * v3;
    duT(3,0) = 0.0;

    duT(2,2) = 4.0 * w2 * (w - 3.0 * u) + 6.0 * u * w * (w - u);
    duT(2,1) = 12.0 * (v * w2 + v2 * w - u * v2 - u2 * v);
    duT(2,0) = 4.0 * v3 + 6.0 * u * v2;

    duT(1,1) = 6.0 * u * w *(w - u) + 4.0 * u2 * (3.0 * w - u);
    duT(1,0) = 6.0 * u * v2 + 12.0 * u2 * v;

    duT(0,0) = 4.0 * u3;

    dvT(3,3) = -4.0 * w3;
    dvT(3,2) = 4.0 * w2 * (w - 3.0 * v) + 6.0 * v * w * (w - v);
    dvT(3,1) = 6.0 * v * w * (w - v) + 4.0 * v2 * (3.0 * w - v);
    dvT(3,0) = 4.0 * v3;

    dvT(2,2) = -12.0 * u * w2 - 6.0 * u2 * w;
    dvT(2,1) = 12.0 * (u * w2 - u * v2 + u2 * w - u2 * v);
    dvT(2,0) = 12.0 * u * v2 + 6.0 * u2 * v;

    dvT(1,1) = -6.0 * u2 * w - 4.0 * u3;
    dvT(1,0) = 6.0 * u2 * v + 4.0 * u3;

    dvT(0,0) = 0.0;

    DCoordinate3 &point = pd(0, 0), &diff1u = pd(1, 0), &diff1v = pd(1, 1);


    for (GLuint k = 0; k < 4; ++k)
    {
        for (GLuint l = 0; l <= k; ++l)
        {
            point  += _data(k, l) *   T(k, l);
            diff1u += _data(k, l) * duT(k, l);
            diff1v += _data(k, l) * dvT(k, l);
        }
    }

    return GL_TRUE;
}


/*SECOND ORDER TRIGONOMETRIC SURFACE*/
SecondOrderTrigonometricSurface::SecondOrderTrigonometricSurface(GLdouble beta):TriangularSurface3(beta), _beta(beta){}

SecondOrderTrigonometricSurface::SecondOrderTrigonometricSurface(const SecondOrderTrigonometricSurface &cpc):TriangularSurface3(cpc){}

SecondOrderTrigonometricSurface& SecondOrderTrigonometricSurface::operator=(const SecondOrderTrigonometricSurface &cpc)
{
    if (this != &cpc) {
        TriangularSurface3::operator=(cpc);
    }
    return *this;
}

GLboolean SecondOrderTrigonometricSurface::CalculatePartialDerivatives(GLuint maximum_order_of_partial_derivatives,
                                                                       GLdouble u, GLdouble v, PartialDerivatives &pd) const
{
    if (u < 0.0 || u > _beta || v < 0.0 || v > _beta) // u + v > _beta?
    {
        pd.ResizeRows(0);
        return GL_FALSE;
    }

    pd.ResizeRows(maximum_order_of_partial_derivatives + 1);
    pd.LoadNullVectors();

    GLdouble w = _beta - u - v;

    GLdouble sb = sin(_beta / 2.0), sb4 = pow(sb, (GLint)4), sb5 = sb4 * sb,
             cb = cos(_beta / 2.0), cb2 = cb * cb, cb3 = cb2 * cb;

    GLdouble su = sin(u / 2.0), su2 = su * su, su3 = su2 * su, su4 = su2 * su2,
             cu = cos(u / 2.0), cu2 = cu * cu;

    GLdouble sv = sin(v / 2.0), sv2 = sv * sv, sv3 = sv2 * sv, sv4 = sv2 * sv2,
             cv = cos(v / 2.0), cv2 = cv * cv;

    GLdouble sw = sin(w / 2.0), sw2 = sw * sw, sw3 = sw2 * sw, sw4 = sw2 * sw2,
             cw = cos(w / 2.0), cw2 = cw * cw;

    GLdouble r440 = 1.0 / sb4,
             r430 = 4.0 * cb / sb4,
             r420 = (2.0 + 4.0 * cb2) / sb4,
             r431 = (4.0 + 8.0 * cb2) / sb5,
             r421 = (16.0 * cb + 8.0 * cb3) / sb5,
             r422 = (10.0 + 20.0 * cb2) / sb4;

    TriangularMatrix<GLdouble> T(4), duT(4), dvT(4);
    //

    T(0, 0) = r440 * su4;

    T(1, 1) = r430 * su3 * sw * cv + 0.5 * r420 * su2 * sw2 * cv2;
    T(1, 0) = r430 * su3 * sv * cw + 0.5 * r420 * su2 * sv2 * cw2;

    T(2, 2) = 0.5 * r420 * su2 * sw2 * cv2 + r430 * su * sw3 * cv;
    T(2, 1) = r431 * (su3 * sw * sv + sv3 * su * sw + sw3 * sv * su) +
              r421 * (su2 * sw2 * cv * sv + su2 * sv2 * cw * sw + sw2 * sv2 * cu * su) +
              r422 * su2 * sv2 * sw2;
    T(2, 0) = 0.5 * r420 * su2 * sv2 * cw2 + r430 * su * sv3 * cw;

    T(3, 3) = r440 * sw4;
    T(3, 2) = r430 * sw3 * sv * cu + 0.5 * r420 * sw2 * sv2 * cu2;
    T(3, 1) = r430 * sv3 * sw * cu + 0.5 * r420 * sw2 * sv2 * cu2;
    T(3, 0) = r440 * sv4;

    //

    duT(0, 0) = 2.0 * r440 * su3 * cu;

    duT(1, 1) = r430 * (1.5 * su2 * cu * sw * cv - 0.5 * su3 * cw * cv) + 0.5 * r420 * (su * cu * sw2 * cv2 - su2 * sw * cw * cv2);
    duT(1, 0) = r430 * (1.5 * su2 * cu * sv * cw + 0.5 * su3 * sv * sw) + 0.5 * r420 * (su * cu * sv2 * cw2 + su2 * sv2 * cw * sw);

    duT(2, 2) = 0.5 * r420 * (su * cu * sw2 * cv2 - su2 * sw * cw * cv2) + r430 * (0.5 * cu * sw3 * cv - 1.5 * su * sw2 * cw * cv);

    duT(2, 1) = r431 * ((1.5 * su2 * cu * sw * sv - 0.5 * su3 * cw * sv) + (0.5 * sv3 * cu * sw - 0.5 * sv3 * su * cw) + (-1.5 * sw2 * cw * sv * su + 0.5 * sw3 * sv * cu)) +
                r421 * ((su * cu * sw2 * cv * sv - su2 * sw * cw * cv * sv) + (su * cu * sv2 * cw * sw + 0.5 * su2 * sv2 * sw2 - 0.5 * su2 * sv2 * cw2) + (-sw * cw * sv2 * cu * su - 0.5 * sw2 * sv2 * su2 + 0.5 * sw2 * sv2 * cu2)) +
                r422 * (su * cu * sv2 * sw2 - su2 * sv2 * sw * cw);

    duT(2, 0) = 0.5 * r420 * (su * cu * sv2 * cw2 + su2 * sv2 * cw * sw) + r430 * (0.5 * cu * sv3 * cw + 0.5 * su * sv3 * sw);

    duT(3, 3) = -2.0 * r440 * sw3 * cw;
    duT(3, 2) = r430 * (-1.5 * sw2 * cw * sv * cu - 0.5 * sw3 * sv * su) + 0.5 * r420 * (-sw * cw * sv2 * cu2 - sw2 * sv2 * cu * su);
    duT(3, 1) = r430 * (-0.5 * sv3 * cw * cu - 0.5 * sv3 * sw * su) + 0.5 * r420 * (-sw * cw * sv2 * cu2 - sw2 * sv2 * cu * su);
    duT(3, 0) = 0.0;

    //

    dvT(0, 0) = 0.0;

    dvT(1, 1) = r430 * (-0.5 * su3 * cw * cv - 0.5 * su3 * sw * sv) + 0.5 * r420 * (- su2 * sw * cw * cv2 - su2 * sw2 * cv * sv);
    dvT(1, 0) = r430 * (0.5 * su3 * cv * cw + 0.5 * su3 * sv * sw) + 0.5 * r420 * (su2 * sv * cv * cw2 + su2 * sv2 * cw * sw);

    dvT(2, 2) = 0.5 * r420 * (-su2 * sw * cw * cv2 - su2 * sw2 * cv * sv) + r430 * (-1.5 * su * sw2 * cw * cv - 0.5 * su * sw3 * sv);
    dvT(2, 1) = r431 * ((-0.5 * su3 * cw * sv + 0.5 * su3 * sw * cv) + (1.5 * sv2 * cv * su * sw - 0.5 * sv3 * su * cw) + (-1.5 * sw2 * cw * sv * su + 0.5 * sw3 * cv * su)) +
                r421 * ((-su2 * sw * cw * cv * sv - 0.5 * su2 * sw2 * sv2 + 0.5 * su2 * sw2 * cv2) + (su2 * sv * cv * cw * sw + 0.5 * su2 * sv2 * sw2 - 0.5 * su2 * sv2 * cw2) + (-sw * cw * sv2 * cu * su + sw2 * sv * cv * cu * su)) +
                r422 * (su2 * sv * cv * sw2 - su2 * sv2 * sw * cw);
    dvT(2, 0) = 0.5 * r420 * (su2 * sv * cv * cw2 + su2 * sv2 * cw * sw) + r430 * (1.5 * su * sv2 * cv * cw + 0.5 * su * sv3 * sw);

    dvT(3, 3) = - 2.0 * r440 * sw3 * cw;
    dvT(3, 2) = r430 * (-1.5 * sw2 * cw * sv * cu + 0.5 * sw3 * cv * cu) + 0.5 * r420 * (sw2 * sv * cv * cu2 - sw * cw * sv2 * cu2);
    dvT(3, 1) = r430 * (1.5 * sv2 * cv * sw * cu - 0.5 * sv3 * cw * cu) + 0.5 * r420 * (-sw * cw * sv2 * cu2 + sw2 * sv * cv * cu2);
    dvT(3, 0) = 2.0 * r440 * sv3 * cv;

    pd.LoadNullVectors();

    DCoordinate3 &point = pd(0, 0), &diff1u = pd(1, 0), &diff1v = pd(1, 1);

    for (GLuint k = 0; k < 4; ++k)
    {
        for (GLuint l = 0; l <= k; ++l)
        {
            point  += _data(k, l) *   T(k, l);
            diff1u += _data(k, l) * duT(k, l);
            diff1v += _data(k, l) * dvT(k, l);
        }
    }

    return GL_TRUE;
}


/*SECOND ORDER HYPERBOLIC SURFACE*/
SecondOrderHyperbolicSurface::SecondOrderHyperbolicSurface(GLdouble beta):TriangularSurface3(beta), _beta(beta){}

SecondOrderHyperbolicSurface::SecondOrderHyperbolicSurface(const SecondOrderHyperbolicSurface &cpc):TriangularSurface3(cpc){}

SecondOrderHyperbolicSurface& SecondOrderHyperbolicSurface::operator=(const SecondOrderHyperbolicSurface &cpc)
{
    if (this != &cpc) {
        TriangularSurface3::operator=(cpc);
    }
    return *this;
}

GLboolean SecondOrderHyperbolicSurface::CalculatePartialDerivatives(GLuint maximum_order_of_partial_derivatives,
                                                                    GLdouble u, GLdouble v, PartialDerivatives &pd) const
{
    if (u < 0.0 || u > _beta || v < 0.0 || v > _beta) // u + v > _beta?
    {
        pd.ResizeRows(0);
        return GL_FALSE;
    }

    pd.ResizeRows(maximum_order_of_partial_derivatives + 1);
    pd.LoadNullVectors();

    GLdouble w = _beta - u - v;

    GLdouble sb = sinh(_beta / 2.0), sb4 = pow(sb, (GLint)4), sb5 = sb4 * sb,
             cb = cosh(_beta / 2.0), cb2 = cb * cb, cb3 = cb2 * cb;

    GLdouble su = sinh(u / 2.0), su2 = su * su, su3 = su2 * su, su4 = su2 * su2,
             cu = cosh(u / 2.0), cu2 = cu * cu;

    GLdouble sv = sinh(v / 2.0), sv2 = sv * sv, sv3 = sv2 * sv, sv4 = sv2 * sv2,
             cv = cosh(v / 2.0), cv2 = cv * cv;

    GLdouble sw = sinh(w / 2.0), sw2 = sw * sw, sw3 = sw2 * sw, sw4 = sw2 * sw2,
             cw = cosh(w / 2.0), cw2 = cw * cw;

    GLdouble r440 = 1.0 / sb4,
             r430 = 4.0 * cb / sb4,
             r420 = (2.0 + 4.0 * cb2) / sb4,
             r431 = (4.0 + 8.0 * cb2) / sb5,
             r421 = (16.0 * cb + 8.0 * cb3) / sb5,
             r422 = (10.0 + 20.0 * cb2) / sb4;

    TriangularMatrix<GLdouble> T(4), duT(4), dvT(4);
    //

    T(0, 0) = r440 * su4;

    T(1, 1) = r430 * su3 * sw * cv + 0.5 * r420 * su2 * sw2 * cv2;
    T(1, 0) = r430 * su3 * sv * cw + 0.5 * r420 * su2 * sv2 * cw2;

    T(2, 2) = 0.5 * r420 * su2 * sw2 * cv2 + r430 * su * sw3 * cv;
    T(2, 1) = r431 * (su3 * sw * sv + sv3 * su * sw + sw3 * sv * su) +
              r421 * (su2 * sw2 * cv * sv + su2 * sv2 * cw * sw + sw2 * sv2 * cu * su) +
              r422 * su2 * sv2 * sw2;
    T(2, 0) = 0.5 * r420 * su2 * sv2 * cw2 + r430 * su * sv3 * cw;

    T(3, 3) = r440 * sw4;
    T(3, 2) = r430 * sw3 * sv * cu + 0.5 * r420 * sw2 * sv2 * cu2;
    T(3, 1) = r430 * sv3 * sw * cu + 0.5 * r420 * sw2 * sv2 * cu2;
    T(3, 0) = r440 * sv4;

    //

    duT(0, 0) = 2.0 * r440 * su3 * cu;

    duT(1, 1) = r430 * (1.5 * su2 * cu * sw * cv - 0.5 * su3 * cw * cv) + 0.5 * r420 * (su * cu * sw2 * cv2 - su2 * sw * cw * cv2);
    duT(1, 0) = r430 * (1.5 * su2 * cu * sv * cw + 0.5 * su3 * sv * sw) + 0.5 * r420 * (su * cu * sv2 * cw2 + su2 * sv2 * cw * sw);

    duT(2, 2) = 0.5 * r420 * (su * cu * sw2 * cv2 - su2 * sw * cw * cv2) + r430 * (0.5 * cu * sw3 * cv - 1.5 * su * sw2 * cw * cv);

    duT(2, 1) = r431 * ((1.5 * su2 * cu * sw * sv - 0.5 * su3 * cw * sv) + (0.5 * sv3 * cu * sw - 0.5 * sv3 * su * cw) + (-1.5 * sw2 * cw * sv * su + 0.5 * sw3 * sv * cu)) +
                r421 * ((su * cu * sw2 * cv * sv - su2 * sw * cw * cv * sv) + (su * cu * sv2 * cw * sw + 0.5 * su2 * sv2 * sw2 - 0.5 * su2 * sv2 * cw2) + (-sw * cw * sv2 * cu * su - 0.5 * sw2 * sv2 * su2 + 0.5 * sw2 * sv2 * cu2)) +
                r422 * (su * cu * sv2 * sw2 - su2 * sv2 * sw * cw);

    duT(2, 0) = 0.5 * r420 * (su * cu * sv2 * cw2 + su2 * sv2 * cw * sw) + r430 * (0.5 * cu * sv3 * cw + 0.5 * su * sv3 * sw);

    duT(3, 3) = -2.0 * r440 * sw3 * cw;
    duT(3, 2) = r430 * (-1.5 * sw2 * cw * sv * cu - 0.5 * sw3 * sv * su) + 0.5 * r420 * (-sw * cw * sv2 * cu2 - sw2 * sv2 * cu * su);
    duT(3, 1) = r430 * (-0.5 * sv3 * cw * cu - 0.5 * sv3 * sw * su) + 0.5 * r420 * (-sw * cw * sv2 * cu2 - sw2 * sv2 * cu * su);
    duT(3, 0) = 0.0;

    //

    dvT(0, 0) = 0.0;

    dvT(1, 1) = r430 * (-0.5 * su3 * cw * cv - 0.5 * su3 * sw * sv) + 0.5 * r420 * (- su2 * sw * cw * cv2 - su2 * sw2 * cv * sv);
    dvT(1, 0) = r430 * (0.5 * su3 * cv * cw + 0.5 * su3 * sv * sw) + 0.5 * r420 * (su2 * sv * cv * cw2 + su2 * sv2 * cw * sw);

    dvT(2, 2) = 0.5 * r420 * (-su2 * sw * cw * cv2 - su2 * sw2 * cv * sv) + r430 * (-1.5 * su * sw2 * cw * cv - 0.5 * su * sw3 * sv);
    dvT(2, 1) = r431 * ((-0.5 * su3 * cw * sv + 0.5 * su3 * sw * cv) + (1.5 * sv2 * cv * su * sw - 0.5 * sv3 * su * cw) + (-1.5 * sw2 * cw * sv * su + 0.5 * sw3 * cv * su)) +
                r421 * ((-su2 * sw * cw * cv * sv - 0.5 * su2 * sw2 * sv2 + 0.5 * su2 * sw2 * cv2) + (su2 * sv * cv * cw * sw + 0.5 * su2 * sv2 * sw2 - 0.5 * su2 * sv2 * cw2) + (-sw * cw * sv2 * cu * su + sw2 * sv * cv * cu * su)) +
                r422 * (su2 * sv * cv * sw2 - su2 * sv2 * sw * cw);
    dvT(2, 0) = 0.5 * r420 * (su2 * sv * cv * cw2 + su2 * sv2 * cw * sw) + r430 * (1.5 * su * sv2 * cv * cw + 0.5 * su * sv3 * sw);

    dvT(3, 3) = - 2.0 * r440 * sw3 * cw;
    dvT(3, 2) = r430 * (-1.5 * sw2 * cw * sv * cu + 0.5 * sw3 * cv * cu) + 0.5 * r420 * (sw2 * sv * cv * cu2 - sw * cw * sv2 * cu2);
    dvT(3, 1) = r430 * (1.5 * sv2 * cv * sw * cu - 0.5 * sv3 * cw * cu) + 0.5 * r420 * (-sw * cw * sv2 * cu2 + sw2 * sv * cv * cu2);
    dvT(3, 0) = 2.0 * r440 * sv3 * cv;

    DCoordinate3 &point = pd(0, 0), &diff1u = pd(1, 0), &diff1v = pd(1, 1);


    for (GLuint k = 0; k < 4; ++k)
    {
        for (GLuint l = 0; l <= k; ++l)
        {
            point  += _data(k, l) *   T(k, l);
            diff1u += _data(k, l) * duT(k, l);
            diff1v += _data(k, l) * dvT(k, l);
        }
    }

    return GL_TRUE;
}

/*FIRST ORDER ALGEBRAIC TRIGONOMETRIC SURFACE*/
FirstOrderAlgebraicTrigonometricSurface::FirstOrderAlgebraicTrigonometricSurface(GLdouble beta):TriangularSurface3(beta), _beta(beta){}

FirstOrderAlgebraicTrigonometricSurface::FirstOrderAlgebraicTrigonometricSurface(const FirstOrderAlgebraicTrigonometricSurface &cpc):TriangularSurface3(cpc){}

FirstOrderAlgebraicTrigonometricSurface& FirstOrderAlgebraicTrigonometricSurface::operator=(const FirstOrderAlgebraicTrigonometricSurface &cpc)
{
    if (this != &cpc) {
        TriangularSurface3::operator=(cpc);
    }
    return *this;
}

GLboolean FirstOrderAlgebraicTrigonometricSurface::CalculatePartialDerivatives(GLuint maximum_order_of_partial_derivatives,
                                                                               GLdouble u, GLdouble v, PartialDerivatives &pd) const
{
    if (u < 0.0 || u > _beta || v < 0.0 || v > _beta) // u + v > _beta?
    {
        pd.ResizeRows(0);
        return GL_FALSE;
    }

    pd.ResizeRows(maximum_order_of_partial_derivatives + 1);
    pd.LoadNullVectors();

    /*GLdouble w = _beta - v - u;
    GLdouble sb = sin(_beta), cb = cos(_beta), cb_2 = cos(_beta / 2.0),
            su = sin(u), sv = sin(v), sw = sin(w), su_2 = sin(u / 2.0), sv_2 = sin(v / 2.0), sw_2 = sin(w / 2.0),
            cu = cos(u), cv = cos(v), cw = cos(w), cu_2 = cos(u / 2.0), cv_2 = cos(v / 2.0), cw_2 = cos(w / 2.0),
            sbu = sin(_beta - u), sbv = sin(_beta - v), sbw = sin(_beta - w),
            cbu = cos(_beta - u), cbv = cos(_beta - v), cbw = cos(_beta - w),
            swv = sin((w - v) / 2.0), swu = sin((w - u) / 2.0),
            sub = sin(u - (_beta / 2.0)), sb_2_w = sin((_beta / 2.0) - w), svb = sin(v - (_beta / 2.0)),
            cub = cos(u - (_beta / 2.0)), cvb = cos(v - (_beta / 2.0)), cb_2_w = cos((_beta / 2.0) - w),
            d1 = _beta - sb, d2 = d1 * (2.0 * sb - _beta - _beta * cb);*/

    TriangularMatrix<GLdouble> T(4), duT(4), dvT(4);

    /*T(0,0) = (u - su) / d1;

    T(1,0) = (w + sw + su - sbv + u * cbv - (_beta - v) * cu) * sb / d2;
    T(1,1) = (v + sv + su - sbw + u * cbw - (_beta - w) * cu) * sb / d2;

    T(2,2) = (u + su + sw - sbv + w * cbv - (_beta - v) * cw) * sb / d2;
    T(2,1) = (4.0 * (3.0 * _beta + 4.0 * sb - _beta * cb) * cb_2) / d2 * su_2 * sv_2 * sw_2 -
            (4.0 * sb * cb_2) / d2 * (u * cu_2 * sv_2 * sw_2 + v * su_2 * cv_2 * sw_2 + w * su_2 * sv_2 * cw_2);
    T(2,0) = (u + su + sv - sbw + v * cbw - (_beta - w) * cv) * sb / d2;

    T(3,3) = (w - sw) / d1;
    T(3,2) = (v + sv + sw - sbu + w * cbu - (_beta - u) * cw) * sb / d2;
    T(3,1) = (w + sw + sv - sbu + v * cbu - (_beta - u) * cv) * sb / d2;
    T(3,0) = (v - sv) / d1;

    duT(0,0) = (1.0 - cu) / d1;

    duT(1,1) = (su * (_beta - v) - cw + cu + cbv - 1.0) * sb / d2;
    duT(1,0) = ((_beta - w) * su - u * sbw) * sb / d2;

    duT(2,2) = (1.0 - (_beta - v) * sw - cw + cu - cbv) * sb / d2;
    duT(2,1) = (4.0 * (3.0 * _beta + 4.0 * sb - _beta * cb) * cb_2) / d2 * 1.0 / 2.0 * sv_2 * (cu_2 * sw_2 - su_2 * cw_2) -
            (4.0 * sb * cb_2) / d2 * (cu_2 * sv_2 * sw_2 - 0.5 * u * su_2 * sv_2 * sw_2 - 0.5 * u * cu_2 * sv_2 * cw_2 +
                                      0.5 * v * cu_2 * cv_2 * sw_2 - 0.5 * v * su_2 * cv_2 * cw_2 - su_2 * sv_2 * cw_2 +
                                      0.5 * w * cu_2 * sv_2 * cw_2 + 0.5 * w * su_2 * sv_2 * sw_2);
    duT(2,0) = (-v * sbw - cbw + cu - cv + 1.0) * sb / d2;

    duT(3,3) = (cw - 1.0) / d1;
    duT(3,2) = ((u - _beta) * sw + w * sbu) * sb / d2;
    duT(3,1) = (cbu + v * sbu + cw + cv + 1.0) * sb / d2;
    duT(3,0) = 0;

    dvT(0,0) = 0;

    dvT(1,1) = (u * sbv - cw + cu + cbv - 1.0) * sb / d2;
    dvT(1,0) = (-u * sbw - cbw - cu + cv + 1.0) * sb / d2;

    dvT(2,2) = ((v - _beta) * sw + w * sbv) * sb / d2;
    dvT(2,1) = (4.0 * (3.0 * _beta + 4.0 * sb - _beta * cb) * cb_2) / d2 * 1.0 / 2.0 * su_2 * (cv_2 * sw_2 - sv_2 * cw_2) -
            (4.0 * sb * cb_2) / d2 * 1.0 / 4.0 * (0.5 * u * cu_2 * cv_2 * sw_2 - 0.5 * u * cu_2 * sv_2 * cw_2 +
                                                  su_2 * cv_2 * sw_2 - 0.5 * v * su_2 * sv_2 * sw_2 - 0.5 * v * su_2 * cv_2 * cw_2 -
                                                  su_2 * sv_2 * cw_2 + 0.5 * w * su_2 * cv_2 * cw_2 + 0.5 * w * su_2 * sv_2 * sw_2);
    dvT(2,0) = ((_beta - w) * sv - v * sw) * sb / d2;

    dvT(3,3) = (cw - 1.0) / d1;
    dvT(3,2) = (-cbu + (u - _beta) * sw - cw + cv + 1.0) * sb / d2;
    dvT(3,1) = (cbu + (_beta - u) * sv - cw + cv - 1.0) * sb / d2;
    dvT(3,0) = (1.0 - cv) / d1;*/

    GLdouble b = _beta, x = u, y = v, z = b - x - y,
                 bx = b - x, by = b - y, bz = x + y;

    GLdouble sb = sin(b), cb = cos(b),
             sbp2 = sin(b / 2.0), cbp2 = cos(b / 2.0),
             sx = sin(x), sy = sin(y), sz = sin(z),
             cx = cos(x), cy = cos(y), cz = cos(z),
             sxp2 = sin(x / 2.0), syp2 = sin(y / 2.0), szp2 = sin(z / 2.0),
             cxp2 = cos(x / 2.0), cyp2 = cos(y / 2.0), czp2 = cos(z / 2.0),
             sbx = sin(bx), cbx = cos(bx),
             sby = sin(by), cby = cos(by),
             sbz = sin(bz), cbz = cos(bz);

    GLdouble c[4] =
    {
      1.0 / (b - sb),
      sb / (2.0 * sb - b - b * cb) / (b - sb),
      4.0 * (3.0 * b + 4.0 * sb - b * cb) * cbp2 / (2.0 * sb - b - b * cb) / (b - sb),
      4.0 * sb * cbp2 / (2.0 * sb - b - b * cb) / (b - sb)
    };

    TriangularMatrix<GLdouble> B(4), dxB(4), dyB(4);

    //

    T(3, 3) = c[0] * (z - sz);
    T(3, 2) = c[1] * (y + sy + sz - sbx + z * cbx - bx * cz);
    T(3, 1) = c[1] * (z + sz + sy - sbx + y * cbx - bx * cy);
    T(3, 0) = c[0] * (y - sy);

    T(2, 2) = c[1] * (x + sx + sz - sby + z * cby - by * cz);
    T(2, 1) = c[2] * sxp2 * syp2 * szp2 - c[3] * (x * cxp2 * syp2 * szp2 + y * sxp2 * cyp2 * szp2 + z * sxp2 * syp2 * czp2);
    T(2, 0) = c[1] * (x + sx + sy - sbz + y * cbz - bz * cy);

    T(1, 1) = c[1] * (z + sz + sx - sby + x * cby - by * cx);
    T(1, 0) = c[1] * (y + sy + sx - sbz + x * cbz - bz * cx);

    T(0, 0) = c[0] * (x - sx);

    // du

    duT(3, 3) = c[0] * (-1.0 + cz);
    duT(3, 2) = c[1] * (z * sbx - bx * sz);
    duT(3, 1) = c[1] * (-1.0 - cz + cbx + y * sbx + cy);
    duT(3, 0) = 0.0;

    duT(2, 2) = c[1] * (1.0 + cx - cz - cby - by * sz);
    duT(2, 1) = 0.5 * c[2] * cxp2 * syp2 * szp2 - 0.5 * c[2] * sxp2 * syp2 * czp2 - c[3] * (cxp2 * syp2 * szp2 - 0.5 * x * sxp2 * syp2 * szp2 - 0.5 * x * cxp2 * syp2 * czp2 + 0.5 * y * cxp2 * cyp2 * szp2 - 0.5 * y * sxp2 * cyp2 * czp2 - sxp2 * syp2 * czp2 + 0.5 * z * cxp2 * syp2 * czp2 + 0.5 * z * sxp2 * syp2 * szp2);
    duT(2, 0) = c[1] * (1.0 + cx - cbz - y * sbz - cy);

    duT(1, 1) = c[1] * (-1.0 - cz + cx + cby + by * sx);
    duT(1, 0) = c[1] * (-x * sbz + bz * sx);

    duT(0, 0) = c[0] * (1.0 - cx);

    // dv

    dvT(3, 3) = c[0] * (-1.0 + cz);
    dvT(3, 2) = c[1] * (1.0 + cy - cz - cbx - bx * sz);
    dvT(3, 1) = c[1] * (-1.0 - cz + cy + cbx + bx * sy);
    dvT(3, 0) = c[0] * (1.0 - cy);

    dvT(2, 2) = c[1] * (z * sby - by * sz);
    dvT(2, 1) = 0.5 * c[2] * sxp2 * cyp2 * szp2 - 0.5 * c[2] * sxp2 * syp2 * czp2 - c[3] * (0.5 * x * cxp2 * cyp2 * szp2 - 0.5 * x * cxp2 * syp2 * czp2 + sxp2 * cyp2 * szp2 - 0.5 * y * sxp2 * syp2 * szp2 - 0.5 * y * sxp2 * cyp2 * czp2 - sxp2 * syp2 * czp2 + 0.5 * z * sxp2 * cyp2 * czp2 + 0.5 * z * sxp2 * syp2 * szp2);
    dvT(2, 0) = c[1] * (-y * sbz + bz * sy);

    dvT(1, 1) = c[1] * (-1.0 - cz + cby + x * sby + cx);
    dvT(1, 0) = c[1] * (1.0 + cy - cbz - x * sbz - cx);

    dvT(0, 0) = 0.0;

    DCoordinate3 &point = pd(0, 0), &diff1u = pd(1, 0), &diff1v = pd(1, 1);


    for (GLuint k = 0; k < 4; ++k)
    {
        for (GLuint l = 0; l <= k; ++l)
        {
            point  += _data(k, l) *   T(k, l);
            diff1u += _data(k, l) * duT(k, l);
            diff1v += _data(k, l) * dvT(k, l);
        }
    }

    return GL_TRUE;
}
