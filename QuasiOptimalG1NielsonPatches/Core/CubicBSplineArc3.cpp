#include "CubicBSplineArc3.h"
#include <algorithm>

using namespace cagd;
using namespace std;

CubicBSplineArc3::CubicBSplineArc3(): LinearCombination3(0.0, 1.0, 4)
{
}

GLboolean CubicBSplineArc3::BlendingFunctionValues(GLdouble u, RowMatrix<GLdouble> &values) const
{
    if (!Derivative0(u, values))
        return GL_FALSE;

    return GL_TRUE;
}

GLboolean CubicBSplineArc3::CalculateDerivatives(GLuint max_order_of_derivatives, GLdouble u, Derivatives &d) const
{
    if (u < 0.0 || u > 1.0)
        return GL_FALSE;

    d.ResizeRows(max_order_of_derivatives + 1);
    d.LoadNullVectors();

    for (GLuint i = 0; i < max_order_of_derivatives; i++) {
        RowMatrix<GLdouble> derivative = RowMatrix<GLdouble>(_data.GetRowCount());
        if (!CalculateDerivative(i, u, derivative))
            return GL_FALSE;

        for (GLuint j = 0; j < _data.GetRowCount(); j++)
            d(i) += _data[j] * derivative(j);
    }

    return GL_TRUE;
}

GLboolean CubicBSplineArc3::CalculateDerivative(GLuint order_of_derivative, GLdouble u, RowMatrix<GLdouble> &derivative) const
{
    switch (order_of_derivative) {
        case 0:
        {
            if (!Derivative0(u, derivative))
                return GL_FALSE;

            break;
        }
        case 1:
        {
            if (!Derivative1(u, derivative))
                return GL_FALSE;

            break;
        }
        case 2:
        {
            if (!Derivative2(u, derivative))
                return GL_FALSE;

            break;
        }
        default: return GL_FALSE;
    }

    return GL_TRUE;
}

GLboolean CubicBSplineArc3::Derivative0(GLdouble u, RowMatrix<GLdouble> &derivative) const
{
    if (u < 0.0 || u > 1.0)
        return GL_FALSE;

    derivative.ResizeColumns(_data.GetRowCount());

    derivative(0) = pow(1.0 - u, 3) / 6.0;
    derivative(1) = (3.0 * pow(u, 3) - 6.0 * pow(u, 2) + 4.0) / 6.0;
    derivative(2) = - (3.0 * pow(u, 3) - 3 * pow(u, 2) - 3 * u - 1) / 6.0;
    derivative(3) = pow(u, 3) / 6.0;

    return GL_TRUE;
}

GLboolean CubicBSplineArc3::Derivative1(GLdouble u, RowMatrix<GLdouble> &derivative) const
{
    if (u < 0.0 || u > 1.0)
        return GL_FALSE;

    derivative.ResizeColumns(_data.GetRowCount());

    derivative(0) = - pow(1.0 - u, 2) / 2.0;
    derivative(1) = (u * (3.0 * u - 4)) / 2.0;
    derivative(2) = (- 3.0 * pow(u, 2) + 2 * u + 1) / 2.0;
    derivative(3) = pow(u, 2) / 2.0;

    return GL_TRUE;
}

GLboolean CubicBSplineArc3::Derivative2(GLdouble u, RowMatrix<GLdouble> &derivative) const
{
    if (u < 0.0 || u > 1.0)
        return GL_FALSE;

    derivative.ResizeColumns(_data.GetRowCount());

    derivative(0) = 1 - u;
    derivative(1) = 3.0 * u - 2.0;
    derivative(2) = 1 - 3 * u;
    derivative(3) = u;

    return GL_TRUE;
}
