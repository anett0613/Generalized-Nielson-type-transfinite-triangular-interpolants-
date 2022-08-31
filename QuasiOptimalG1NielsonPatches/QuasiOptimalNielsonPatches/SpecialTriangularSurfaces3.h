#ifndef SPECIALTRIANGULARSURFACES3_H
#define SPECIALTRIANGULARSURFACES3_H

#include "TriangularSurfaces3.h"

namespace cagd
{
    class CubicPolinomialSurface: public TriangularSurface3 {
    public:
        CubicPolinomialSurface();
        CubicPolinomialSurface(const CubicPolinomialSurface& cpc);
        CubicPolinomialSurface& operator =(const CubicPolinomialSurface& cpc);
        GLboolean CalculatePartialDerivatives(
                        GLuint maximum_order_of_partial_derivatives,
                        GLdouble u, GLdouble v, PartialDerivatives& pd) const;

    };

    class QuarticPolinomialSurface: public TriangularSurface3 {
    public:
        QuarticPolinomialSurface();
        QuarticPolinomialSurface(const QuarticPolinomialSurface& cpc);
        QuarticPolinomialSurface& operator =(const QuarticPolinomialSurface& cpc);
        GLboolean CalculatePartialDerivatives(
                        GLuint maximum_order_of_partial_derivatives,
                        GLdouble u, GLdouble v, PartialDerivatives& pd) const;
    };

    class SecondOrderTrigonometricSurface: public TriangularSurface3 {
    protected:
        GLdouble _beta;
    public:
        SecondOrderTrigonometricSurface(GLdouble beta);
        SecondOrderTrigonometricSurface(const SecondOrderTrigonometricSurface& cpc);
        SecondOrderTrigonometricSurface& operator =(const SecondOrderTrigonometricSurface& cpc);
        GLboolean CalculatePartialDerivatives(
                        GLuint maximum_order_of_partial_derivatives,
                        GLdouble u, GLdouble v, PartialDerivatives& pd) const;
    };

    class SecondOrderHyperbolicSurface: public TriangularSurface3 {
    protected:
        GLdouble _beta;
    public:
        SecondOrderHyperbolicSurface(GLdouble beta);
        SecondOrderHyperbolicSurface(const SecondOrderHyperbolicSurface& cpc);
        SecondOrderHyperbolicSurface& operator =(const SecondOrderHyperbolicSurface& cpc);
        GLboolean CalculatePartialDerivatives(
                        GLuint maximum_order_of_partial_derivatives,
                        GLdouble u, GLdouble v, PartialDerivatives& pd) const;
    };

    class FirstOrderAlgebraicTrigonometricSurface: public TriangularSurface3 {
    protected:
        GLdouble _beta;
    public:
        FirstOrderAlgebraicTrigonometricSurface(GLdouble beta);
        FirstOrderAlgebraicTrigonometricSurface(const FirstOrderAlgebraicTrigonometricSurface& cpc);
        FirstOrderAlgebraicTrigonometricSurface& operator =(const FirstOrderAlgebraicTrigonometricSurface& cpc);
        GLboolean CalculatePartialDerivatives(
                        GLuint maximum_order_of_partial_derivatives,
                        GLdouble u, GLdouble v, PartialDerivatives& pd) const;
    };
}

#endif // SPECIALTRIANGULARSURFACES3_H
