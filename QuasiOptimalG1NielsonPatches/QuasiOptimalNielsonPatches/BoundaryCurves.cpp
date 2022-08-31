#include "./BoundaryCurves.h"

using namespace cagd;
using namespace std;

GLuint fact(GLuint n)
{
    int res = 1;
    for (GLuint i = 2; i <= n; i++)
        res = res * i;
    return res;
}

CubicPolinomialCurve::CubicPolinomialCurve():LinearCombination3(0.0, 1.0, 4) {}

CubicPolinomialCurve::CubicPolinomialCurve(const CubicPolinomialCurve& cpc): LinearCombination3(cpc) {}

CubicPolinomialCurve& CubicPolinomialCurve::operator=(const CubicPolinomialCurve &cpc) {
    if (this != &cpc) {
        LinearCombination3::operator=(cpc);
    }
    return *this;
}

GLboolean CubicPolinomialCurve::BlendingFunctionValues(GLdouble u, RowMatrix<GLdouble> &values) const {
    if (u < _u_min || u > _u_max)
    {
        values.ResizeColumns(0);
        return GL_FALSE;
    }

    values.ResizeColumns(4);

    GLdouble u2 = u * u, u3 = u2 * u, w = 1.0 - u, w2 = w * w, w3 = w2 * w;

    values[0] = w3;
    values[1] = 3.0 * u * w2;
    values[2] = 3.0 * u2 * w;
    values[3] = u3;

    return GL_TRUE;
}

GLboolean CubicPolinomialCurve::CalculateDerivatives(GLuint max_order_of_derivatives, GLdouble u, Derivatives &d) const {

    if (u < _u_min || u > _u_max || max_order_of_derivatives >= 4)
    {
        d.ResizeRows(0);
        return GL_FALSE;
    }

    GLdouble u2 = u * u, u3 = u2 * u, w = 1.0 - u, w2 = w * w, w3 = w2 * w;

    Matrix<GLdouble> dF(max_order_of_derivatives + 1, 4);

    if (max_order_of_derivatives >= 0)
    {
        dF(0, 0) = w3;
        dF(0, 1) = 3.0 * u * w2;
        dF(0, 2) = 3.0 * u2 * w;
        dF(0, 3) = u3;
    }

    if (max_order_of_derivatives >= 1)
    {
        dF(1, 0) = -3.0 * w2;
        dF(1, 1) = 3.0 * (w2 - 2.0 * u * w);
        dF(1, 2) = 3.0 * (2.0 * u * w - u2);
        dF(1, 3) = 3.0 * u2;
    }

    if (max_order_of_derivatives >= 2)
    {
        dF(2, 0) = 6.0 * w;
        dF(2, 1) = -6.0 * (2.0 * w - u);
        dF(2, 2) = 6.0 * (w - 2.0 * u);
        dF(2, 3) = 6.0 * u;
    }

    if (max_order_of_derivatives >= 3)
    {
        dF(3, 0) = -6.0;
        dF(3, 1) = 18.0;
        dF(3, 2) = -18.0;
        dF(3, 3) = 6.0;
    }



    d.ResizeRows(max_order_of_derivatives + 1);
    d.LoadNullVectors();
    // c(u) = _data[0] * F_0(u) + _data[1] * F_1(u) + _data[2] * F_2(u) + _data[3] * F_3(u) = d[0]
    // c^{(r)}(u) = _data[0] * F_0^{(r)}(u) + _data[1] * F_1^{(r)}(u) + _data[2] * F_2^{(r)}(u) + _data[3] * F_3^{(r)}(u) = d[r], r = 0,1,2,3,...

    for (GLuint r = 0; r <= max_order_of_derivatives; r++){
       for (GLuint i = 0; i <= 3; i++)
       {
           d[r] += _data[i] * dF(r, i);
       }
    }

    return GL_TRUE;
}


QuarticPolinomialCurve::QuarticPolinomialCurve():LinearCombination3(0.0, 1.0, 4) {}

QuarticPolinomialCurve::QuarticPolinomialCurve(const QuarticPolinomialCurve& cpc): LinearCombination3(cpc) {}

QuarticPolinomialCurve& QuarticPolinomialCurve::operator=(const QuarticPolinomialCurve &cpc) {
    if (this != &cpc) {
        LinearCombination3::operator=(cpc);
    }
    return *this;
}

GLboolean QuarticPolinomialCurve::BlendingFunctionValues(GLdouble u, RowMatrix<GLdouble> &values) const {
    if (u < _u_min || u > _u_max)
    {
        values.ResizeColumns(0);
        return GL_FALSE;
    }
    values.ResizeColumns(4);

    GLdouble w = 1 - u, w2 = w * w, w3 = w2 * w, w4 = w2 * w2;
    GLdouble u2 = u * u, u3 = u2 * u, u4 = u2 * u2;

    values[0] = w4;
    values[1] = 4.0 * u * w3 + 3.0 * u2 * w2;
    values[2] = 4.0 * u3 * w + 3.0 * u2 * w2;
    values[3] = u4;

    return GL_TRUE;
}

GLboolean QuarticPolinomialCurve::CalculateDerivatives(GLuint max_order_of_derivatives, GLdouble u, Derivatives &d) const {
    if (u < _u_min || u > _u_max || max_order_of_derivatives >= 4)
    {
        d.ResizeRows(0);
        return GL_FALSE;
    }

    Matrix<GLdouble> dF(max_order_of_derivatives + 1, 4);
    d.ResizeRows(max_order_of_derivatives + 1);
    d.LoadNullVectors();

    GLdouble w = 1 - u, w2 = w * w, w3 = w2 * w, w4 = w2 * w2;
    GLdouble u2 = u * u, u3 = u2 * u, u4 = u2 * u2;

    if (max_order_of_derivatives >= 0){
        dF(0,0) = w4;
        dF(0,1) = 4.0 * u * w3 + 3.0 * u2 * w2;
        dF(0,2) = 4.0 * u3 * w + 3.0 * u2 * w2;
        dF(0,3) = u4;
    }

    if (max_order_of_derivatives >= 1){
        dF(1,0) = -4.0 * w3;
        dF(1,1) = 4.0 * (w3 - 3.0 * u * w2) + 6.0 * (u * w2 - u2 * w);
        dF(1,2) = 4.0 * (3.0 * u2 * w - u3) + 6.0 * (u * w2 - u2 * w);
        dF(1,3) = 4.0 * u3;
    }

    if (max_order_of_derivatives >= 2){
        dF(2,0) = 12.0 * w2;
        dF(2,1) = 24.0 * (u * w - w2) + 6.0 * (w2 - 4.0 * u * w + u2);
        dF(2,2) = 24.0 * (u * w - u2) + 6.0 * (w2 - 4.0 * u * w + u2);
        dF(2,3) = 12.0 * u2;
    }

    if (max_order_of_derivatives >= 3){
        dF(3,0) = -24.0 * w;
        dF(3,1) = 24.0 * (3.0 * w - u) + 36.0 * (u - w);
        dF(3,2) = 24.0 * (w - 3.0 * u) + 36.0 * (u - w);
        dF(3,3) = 24.0 * u;
    }


    for (GLuint r = 0; r <= max_order_of_derivatives; r++){
       for (GLuint i = 0; i <= 3; i++)
       {
           d[r] += _data[i] * dF(r, i);
       }
    }
    return GL_TRUE;
}

SecondOrderTrigonometricCurve::SecondOrderTrigonometricCurve(GLdouble beta):LinearCombination3(0.0, beta, 4), _beta(beta) {}

SecondOrderTrigonometricCurve::SecondOrderTrigonometricCurve(const SecondOrderTrigonometricCurve& cpc): LinearCombination3(cpc) {
    _beta = cpc._beta;
}

SecondOrderTrigonometricCurve& SecondOrderTrigonometricCurve::operator=(const SecondOrderTrigonometricCurve &cpc) {
    if (this != &cpc) {
        LinearCombination3::operator=(cpc);
        _beta = cpc._beta;
    }
    return *this;
}

GLboolean SecondOrderTrigonometricCurve::BlendingFunctionValues(GLdouble u, RowMatrix<GLdouble> &values) const
{
    if (u < _u_min || u > _u_max)
    {
        values.ResizeColumns(0);
        return GL_FALSE;
    }
    values.ResizeColumns(4);

    RowMatrix<GLdouble> c;
    c.ResizeColumns(5);

    GLdouble sin_beta = sin(_beta / 2.0), sin_beta2 = sin_beta * sin_beta, sin_beta4 = sin_beta2 * sin_beta2,
            cos_beta = cos(_beta / 2.0), cos_beta2 = cos_beta * cos_beta;

    c[0] = 1.0 / sin_beta4;
    c[1] = 4.0 * cos_beta / sin_beta4;
    c[2] = (4.0 * cos_beta2 + 2.0) / sin_beta4;
    c[3] = c[1];
    c[4] = c[0];

    GLdouble sin_beta_u = sin((_beta - u) / 2.0), sin_beta_u2 = sin_beta_u * sin_beta_u, sin_beta_u3 = sin_beta_u * sin_beta_u2,
            sin_beta_u4 = sin_beta_u2 * sin_beta_u2;
    GLdouble sin_u = sin(u / 2.0), sin_u2 = sin_u * sin_u, sin_u3 = sin_u * sin_u2, sin_u4 = sin_u2 * sin_u2;

    values[0] = c[0] * sin_beta_u4;
    values[1] = c[1] * sin_beta_u3 * sin_u + 1.0 / 2.0 * c[2] * sin_beta_u2 * sin_u2;
    values[2] = c[3] * sin_beta_u * sin_u3 + 1.0 / 2.0 * c[2] * sin_beta_u2 * sin_u2;
    values[3] = c[4] * sin_u4;

    /*GLdouble _b = _beta;

    values[0] = pow(sin((_beta - u) / 2.0), (GLint)4) / pow(sin(_beta / 2.0), (GLint)4);
    values[1] = pow(sin((_beta - u) / 2.0), (GLint)3) * sin(u / 2.0) * 4.0 * cos(_b / 2.0) / pow(sin(_beta / 2.0), (GLint)4) +
                0.5 * pow(sin((_beta - u) / 2.0), (GLint)2) * pow(sin(u / 2.0), (GLint)2) * (2.0 + 4.0 * pow(cos(_beta / 2.0), (GLint)2)) / pow(sin(_beta / 2.0), (GLint)4);
    values[2] = 0.5 * pow(sin((_beta - u) / 2.0), (GLint)2) * pow(sin(u / 2.0), (GLint)2) * (2.0 + 4.0 * pow(cos(_beta / 2.0), (GLint)2)) / pow(sin(_beta / 2.0), (GLint)4) +
                sin((_beta - u) / 2.0) * pow(sin(u / 2.0), (GLint)3) * 4.0 * cos(_beta / 2.0) / pow(sin(_beta / 2.0), (GLint)4);
    values[3] = pow(sin(u / 2.0), (GLint)4) / pow(sin(_beta / 2.0), (GLint)4);*/

    return GL_TRUE;
}

GLboolean SecondOrderTrigonometricCurve::CalculateDerivatives(GLuint max_order_of_derivatives, GLdouble u, Derivatives &d) const {
    if (u < _u_min || u > _u_max || max_order_of_derivatives >= 4)
    {
        d.ResizeRows(0);
        return GL_FALSE;
    }

    Matrix<GLdouble> dF(max_order_of_derivatives + 1, 4);

    d.ResizeRows(max_order_of_derivatives + 1);
    d.LoadNullVectors();

    RowMatrix<GLdouble> c;
    c.ResizeColumns(5);

    GLdouble sin_beta = sin(_beta / 2.0), sin_beta2 = sin_beta * sin_beta, sin_beta4 = sin_beta2 * sin_beta2,
            cos_beta = cos(_beta / 2.0), cos_beta2 = cos_beta * cos_beta;

    c[0] = 1.0 / sin_beta4;
    c[1] = 4.0 * cos_beta / sin_beta4;
    c[2] = (4.0 * cos_beta2 + 2.0) / sin_beta4;
    c[3] = c[1];
    c[4] = c[0];

    GLdouble sin_beta_u = sin((_beta - u) / 2.0), sin_beta_u2 = sin_beta_u * sin_beta_u, sin_beta_u3 = sin_beta_u * sin_beta_u2,
            sin_beta_u4 = sin_beta_u2 * sin_beta_u2,
            cos_beta_u = cos((_beta - u) / 2.0), cos_beta_u2 = cos_beta_u * cos_beta_u, cos_beta_u3 = cos_beta_u * cos_beta_u2;
    GLdouble sin_u = sin(u / 2.0), sin_u2 = sin_u * sin_u, sin_u3 = sin_u * sin_u2, sin_u4 = sin_u2 * sin_u2,
            cos_u = cos(u / 2.0), cos_u2 = cos_u * cos_u, cos_u3 = cos_u * cos_u2;

    if (max_order_of_derivatives >= 0){
        dF(0,0) = c[0] * sin_beta_u4;
        dF(0,1) = c[1] * sin_beta_u3 * sin_u + 1.0 / 2.0 * c[2] * sin_beta_u2 * sin_u2;
        dF(0,2) = c[3] * sin_beta_u * sin_u3 + 1.0 / 2.0 * c[2] * sin_beta_u2 * sin_u2;
        dF(0,3) = c[4] * sin_u4;
    }

    if (max_order_of_derivatives >= 1){
        dF(1,0) = -2.0 * c[0] * sin_beta_u3 * cos_beta_u;
        dF(1,1) = c[1] * 1.0/2.0 * (sin_beta_u3 * cos_u - 3.0 * sin_beta_u2 * cos_beta_u * sin_u) +
                1.0 / 2.0 * c[2] * (sin_beta_u2 * sin_u * cos_u - sin_beta_u * cos_beta_u * sin_u2);
        dF(1,2) = c[3] * 1.0 / 2.0 * (3.0 * sin_beta_u * sin_u2 * cos_u - cos_beta_u * sin_u3) +
                1.0 / 2.0 * c[2] * (sin_beta_u2 * sin_u * cos_u - sin_beta_u * cos_beta_u * sin_u2);
        dF(1,3) = 2.0 * c[4] * sin_u3 * cos_u;
    }

    if (max_order_of_derivatives >= 2){
        dF(2,0) = c[0] * (3.0 * sin_beta_u2 * cos_beta_u2 - sin_beta_u4);
        dF(2,1) = c[1] * 1.0 / 2.0 * (-3.0 * sin_beta_u2 * cos_beta_u * cos_u - 2.0 * sin_beta_u3 * sin_u + 3.0 * sin_beta_u * cos_beta_u2 * sin_u) +
                1.0 / 4.0 * c[2] * (- 4.0 * sin_beta_u * cos_beta_u * sin_u * cos_u - 2.0 * sin_beta_u2 * sin_u2 + sin_beta_u2 * cos_u2 + cos_beta_u2 * sin_u2);
        dF(2,2) = c[3] * 1.0 / 2.0 * (-3.0 * cos_beta_u * sin_u2 * cos_u + 3.0 * sin_beta_u * sin_u * cos_u2 - 2.0 * sin_beta_u * sin_u3) +
                1.0 / 4.0 * c[2] * (- 4.0 * sin_beta_u * cos_beta_u * sin_u * cos_u - 2.0 * sin_beta_u2 * sin_u2 + sin_beta_u2 * cos_u2 + cos_beta_u2 * sin_u2);
        dF(2,3) = c[4] * (3.0 * sin_u2 * cos_u2 - sin_u4);
    }

    if (max_order_of_derivatives >= 3){
        dF(3,0) = c[0] * (5.0 * sin_beta_u3 * cos_beta_u - 3.0 * sin_beta_u * cos_beta_u3);
        dF(3,1) = c[1] * 1.0 / 4.0 * (9.0 * sin_beta_u * cos_beta_u2 * cos_u - 5.0 * sin_beta_u3 * cos_u + 15.0 * sin_beta_u2 * cos_beta_u * sin_u - 3.0 * cos_beta_u3 * sin_u) +
                1.0 / 4.0 * c[2] * (3.0 * cos_beta_u2 * sin_u * cos_u - 5.0 * sin_beta_u2 * sin_u * cos_u - 3.0 * sin_beta_u * cos_beta_u * cos_u2 + 5.0 * sin_beta_u * cos_beta_u * sin_u2);
        dF(3,2) = c[3] * 1.0 / 4.0 * (-15.0 * sin_beta_u * sin_u2 * cos_u + 5.0 * cos_beta_u * sin_u3 - 9.0 * cos_beta_u * sin_u * cos_u2 + 3.0 * sin_beta_u * cos_u3) +
               1.0 / 4.0 * c[2] * (3.0 * cos_beta_u2 * sin_u * cos_u - 5.0 * sin_beta_u2 * sin_u * cos_u - 3.0 * sin_beta_u * cos_beta_u * cos_u2 + 5.0 * sin_beta_u * cos_beta_u * sin_u2);
        dF(3,3) = c[4] * (3.0 * sin_u * cos_u3 - 5.0 * sin_u3 * cos_u);
    }

    for (GLuint r = 0; r <= max_order_of_derivatives; r++){
       for (GLuint i = 0; i <= 3; i++)
       {
           d[r] += _data[i] * dF(r, i);
       }
    }

    return GL_TRUE;
}

SecondOrderHyperbolicCurve::SecondOrderHyperbolicCurve(GLdouble beta):LinearCombination3(0.0, beta, 4), _beta(beta) {}

SecondOrderHyperbolicCurve::SecondOrderHyperbolicCurve(const SecondOrderHyperbolicCurve& cpc): LinearCombination3(cpc) {
    _beta = cpc._beta;
}

SecondOrderHyperbolicCurve& SecondOrderHyperbolicCurve::operator=(const SecondOrderHyperbolicCurve &cpc) {
    if (this != &cpc) {
        LinearCombination3::operator=(cpc);
        _beta = cpc._beta;
    }
    return *this;
}

GLboolean SecondOrderHyperbolicCurve::BlendingFunctionValues(GLdouble u, RowMatrix<GLdouble> &values) const {
    if (u < _u_min || u > _u_max)
    {
        values.ResizeColumns(0);
        return GL_FALSE;
    }
    values.ResizeColumns(4);

    RowMatrix<GLdouble> c;
    c.ResizeColumns(5);

    GLdouble sin_beta = sinh(_beta / 2.0), sin_beta2 = sin_beta * sin_beta, sin_beta4 = sin_beta2 * sin_beta2,
            cos_beta = cosh(_beta / 2.0), cos_beta2 = cos_beta * cos_beta;

    c[0] = 1.0 / sin_beta4;
    c[1] = 1.0 / sin_beta4 * 4.0 * cos_beta;
    c[2] = 1.0 / sin_beta4 * (4.0 * cos_beta2 + 2.0);
    c[3] = c[1];
    c[4] = c[0];

    GLdouble sin_beta_u = sinh((_beta - u) / 2.0), sin_beta_u2 = sin_beta_u * sin_beta_u, sin_beta_u3 = sin_beta_u * sin_beta_u2,
            sin_beta_u4 = sin_beta_u2 * sin_beta_u2;
    GLdouble sin_u = sinh(u / 2.0), sin_u2 = sin_u * sin_u, sin_u3 = sin_u * sin_u2, sin_u4 = sin_u2 * sin_u2;

    values[0] = c[0] * sin_beta_u4;
    values[1] = c[1] * sin_beta_u3 * sin_u + 1.0 / 2.0 * c[2] * sin_beta_u2 * sin_u2;
    values[2] = c[3] * sin_beta_u * sin_u3 + 1.0 / 2.0 * c[2] * sin_beta_u2 * sin_u2;
    values[3] = c[4] * sin_u4;

    return GL_TRUE;
}

GLboolean SecondOrderHyperbolicCurve::CalculateDerivatives(GLuint max_order_of_derivatives, GLdouble u, Derivatives &d) const {
    if (u < _u_min || u > _u_max || max_order_of_derivatives >= 4)
    {
        d.ResizeRows(0);
        return GL_FALSE;
    }

    Matrix<GLdouble> dF(max_order_of_derivatives + 1, 4);

    d.ResizeRows(max_order_of_derivatives + 1);
    d.LoadNullVectors();

    RowMatrix<GLdouble> c;
    c.ResizeColumns(5);

    GLdouble sin_beta = sinh(_beta / 2.0), sin_beta2 = sin_beta * sin_beta, sin_beta4 = sin_beta2 * sin_beta2,
            cos_beta = cosh(_beta / 2.0), cos_beta2 = cos_beta * cos_beta;

    c[0] = 1.0 / sin_beta4;
    c[1] = 1.0 / sin_beta4 * 4.0 * cos_beta;
    c[2] = 1.0 / sin_beta4 * (4.0 * cos_beta2 + 2.0);
    c[3] = c[1];
    c[4] = c[0];

    GLdouble sin_beta_u = sinh((_beta - u) / 2.0), sin_beta_u2 = sin_beta_u * sin_beta_u, sin_beta_u3 = sin_beta_u * sin_beta_u2,
            sin_beta_u4 = sin_beta_u2 * sin_beta_u2,
            cos_beta_u = cosh((_beta - u) / 2.0), cos_beta_u2 = cos_beta_u * cos_beta_u, cos_beta_u3 = cos_beta_u * cos_beta_u2;
    GLdouble sin_u = sinh(u / 2.0), sin_u2 = sin_u * sin_u, sin_u3 = sin_u * sin_u2, sin_u4 = sin_u2 * sin_u2,
            cos_u = cosh(u / 2.0), cos_u2 = cos_u * cos_u, cos_u3 = cos_u * cos_u2;

    if (max_order_of_derivatives >= 0){
        dF(0,0) = c[0] * sin_beta_u4;
        dF(0,1) = c[1] * sin_beta_u3 * sin_u + 1.0 / 2.0 * c[2] * sin_beta_u2 * sin_u2;
        dF(0,2) = c[3] * sin_beta_u * sin_u3 + 1.0 / 2.0 * c[2] * sin_beta_u2 * sin_u2;
        dF(0,3) = c[4] * sin_u4;
    }

    if (max_order_of_derivatives >= 1){
        dF(1,0) = -2.0 * c[0] * sin_beta_u3 * cos_beta_u;
        dF(1,1) = c[1] * 1.0/2.0 * (sin_beta_u3 * cos_u - 3.0 * sin_beta_u2 * cos_beta_u * sin_u) +
                1.0 / 2.0 * c[2] * (sin_beta_u2 * sin_u * cos_u - sin_beta_u * cos_beta_u * sin_u2);
        dF(1,2) = c[3] * 1.0 / 2.0 * (3.0 * sin_beta_u * sin_u2 * cos_u - cos_beta_u * sin_u3) +
                1.0 / 2.0 * c[2] * (sin_beta_u2 * sin_u * cos_u - sin_beta_u * cos_beta_u * sin_u2);
        dF(1,3) = 2.0 * c[4] * sin_u3 * cos_u;
    }

    if (max_order_of_derivatives >= 2){
        dF(2,0) = c[0] * (3.0 * sin_beta_u2 * cos_beta_u2 - sin_beta_u4);
        dF(2,1) = c[1] * 1.0 / 2.0 * (-3.0 * sin_beta_u2 * cos_beta_u * cos_u - 2.0 * sin_beta_u3 * sin_u + 3.0 * sin_beta_u * cos_beta_u2 * sin_u) +
                1.0 / 4.0 * c[2] * (- 4.0 * sin_beta_u * cos_beta_u * sin_u * cos_u - 2.0 * sin_beta_u2 * sin_u2 + sin_beta_u2 * cos_u2 + cos_beta_u2 * sin_u2);
        dF(2,2) = c[3] * 1.0 / 2.0 * (-3.0 * cos_beta_u * sin_u2 * cos_u + 3.0 * sin_beta_u * sin_u * cos_u2 - 2.0 * sin_beta_u * sin_u3) +
                1.0 / 4.0 * c[2] * (- 4.0 * sin_beta_u * cos_beta_u * sin_u * cos_u - 2.0 * sin_beta_u2 * sin_u2 + sin_beta_u2 * cos_u2 + cos_beta_u2 * sin_u2);
        dF(2,3) = c[4] * (3.0 * sin_u2 * cos_u2 - sin_u4);
    }

    if (max_order_of_derivatives >= 3){
        dF(3,0) = c[0] * (5.0 * sin_beta_u3 * cos_beta_u - 3.0 * sin_beta_u * cos_beta_u3);
        dF(3,1) = c[1] * 1.0 / 4.0 * (9.0 * sin_beta_u * cos_beta_u2 * cos_u - 5.0 * sin_beta_u3 * cos_u + 15.0 *sin_beta_u2 * cos_beta_u * sin_u - 3.0 * cos_beta_u3 * sin_u) +
                1.0 / 4.0 * c[2] * (3.0 * cos_beta_u2 * sin_u * cos_u - 5.0 * sin_beta_u2 * sin_u * cos_u - 3.0 * sin_beta_u * cos_beta_u * cos_u2 + 5.0 * sin_beta_u * cos_beta_u * sin_u2);
        dF(3,2) = c[3] * 1.0 / 4.0 * (-15.0 * sin_beta_u * sin_u2 * cos_u + 5.0 * cos_beta_u * sin_u3 - 9.0 * cos_beta_u * sin_u * cos_u2 + 3.0 * sin_beta_u * cos_u3) +
               1.0 / 4.0 * c[2] * (3.0 * cos_beta_u2 * sin_u * cos_u - 5.0 * sin_beta_u2 * sin_u * cos_u - 3.0 * sin_beta_u * cos_beta_u * cos_u2 + 5.0 * sin_beta_u * cos_beta_u * sin_u2);
        dF(3,3) = c[4] * (3.0 * sin_u * cos_u3 - 5.0 * sin_u3 * cos_u);
    }

    for (GLuint r = 0; r <= max_order_of_derivatives; r++){
       for (GLuint i = 0; i <= 3; i++)
       {
           d[r] += _data[i] * dF(r, i);
       }
    }


    return GL_TRUE;
}

FirstOrderAlgebraicTrigonometricCurve::FirstOrderAlgebraicTrigonometricCurve(GLdouble beta):LinearCombination3(0.0, beta, 4), _beta(beta) {}

FirstOrderAlgebraicTrigonometricCurve::FirstOrderAlgebraicTrigonometricCurve(const FirstOrderAlgebraicTrigonometricCurve& cpc): LinearCombination3(cpc) {
    _beta = cpc._beta;
}

FirstOrderAlgebraicTrigonometricCurve& FirstOrderAlgebraicTrigonometricCurve::operator=(const FirstOrderAlgebraicTrigonometricCurve &cpc) {
    if (this != &cpc) {
        LinearCombination3::operator=(cpc);
        _beta = cpc._beta;
    }
    return *this;
}

GLboolean FirstOrderAlgebraicTrigonometricCurve::BlendingFunctionValues(GLdouble u, RowMatrix<GLdouble> &values) const {
    if (u < _u_min || u > _u_max)
    {
        values.ResizeColumns(0);
        return GL_FALSE;
    }
    values.ResizeColumns(4);

    GLdouble sin_u = sin(u), sin_beta = sin(_beta), sin_beta_u = sin(_beta - u),
            cos_u = cos(u), cos_beta = cos(_beta), cos_beta_u = cos(_beta - u),
            denom1 = _beta - sin_beta, denom2 = denom1 * (2.0 * sin_beta - _beta - _beta * cos_beta);

    values[0] = ((_beta - u) - sin_beta_u) / denom1;
    values[1] = ((u + sin_u + sin_beta_u - sin_beta + (_beta - u) * cos_beta - _beta * cos_beta_u) * sin_beta) / denom2;
    values[2] = ((_beta - u + sin_beta_u + sin_u - sin_beta + u * cos_beta - _beta * cos_u) * sin_beta) / denom2;
    values[3] = ((u) - sin_u) / denom1;

    return GL_TRUE;
}

GLboolean FirstOrderAlgebraicTrigonometricCurve::CalculateDerivatives(GLuint max_order_of_derivatives, GLdouble u, Derivatives &d) const {
    if (u < _u_min || u > _u_max || max_order_of_derivatives >= 4)
    {
        d.ResizeRows(0);
        return GL_FALSE;
    }

    Matrix<GLdouble> dF(max_order_of_derivatives + 1, 4);

    d.ResizeRows(max_order_of_derivatives + 1);
    d.LoadNullVectors();

    GLdouble sin_u = sin(u), sin_beta = sin(_beta), sin_beta_u = sin(_beta - u),
            cos_u = cos(u), cos_beta = cos(_beta), cos_beta_u = cos(_beta - u),
            denom1 = _beta - sin_beta, denom2 = denom1 * (2.0 * sin_beta - _beta - _beta * cos_beta);

    if (max_order_of_derivatives >= 0){
        dF(0,0) = (_beta - u - sin_beta_u) / denom1;
        dF(0,1) = (u + sin_u + sin_beta_u - sin_beta + (_beta - u) * cos_beta - _beta * cos_beta_u) * sin_beta / denom2;
        dF(0,2) = (_beta - u + sin_beta_u + sin_u - sin_beta + u * cos_beta - _beta * cos_u) * sin_beta / denom2;
        dF(0,3) = (u - sin_u) / denom1;
    }

    if (max_order_of_derivatives >= 1){
        dF(1,0) = (-1.0 + cos_beta_u) / denom1;
        dF(1,1) = (1.0 + cos_u - cos_beta_u - cos_beta - _beta * sin_beta_u) * sin(_beta) / denom2;
        dF(1,2) = (-1.0 - cos_beta_u + cos_u + cos_beta + _beta * sin_u) * sin(_beta) / denom2;
        dF(1,3) = (1.0 - cos_u) / denom1;
    }

    if (max_order_of_derivatives >= 2){
        dF(2,0) = sin_beta_u / denom1;
        dF(2,1) = (-sin_u - sin_beta_u + _beta * cos_beta_u) * sin(_beta) / denom2;
        dF(2,2) = (-sin_beta_u - sin_u + _beta * cos_u) * sin(_beta) / denom2;
        dF(2,3) = sin_u / denom1;
    }

    if (max_order_of_derivatives >= 3){
        dF(3,0) = - cos_beta_u / denom1;
        dF(3,1) = (-cos_u + cos_beta_u + _beta * sin_beta_u) * sin(_beta) / denom2;
        dF(3,2) = (cos_beta_u - cos_u - _beta * sin_u) * sin(_beta) / denom2;
        dF(3,3) = cos_u / denom1;
    }

    for (GLuint r = 0; r <= max_order_of_derivatives; r++){
       for (GLuint i = 0; i <= 3; i++)
       {
           d[r] += _data[i] * dF(r, i);
       }
    }

    return GL_TRUE;
}


