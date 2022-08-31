#ifndef BOUNDARYCURVES_H
#define BOUNDARYCURVES_H

#include "../Core/LinearCombination3.h"
#include "../Core/Constants.h"

namespace cagd {
    class CubicPolinomialCurve: public LinearCombination3 {
    public:
        CubicPolinomialCurve();
        CubicPolinomialCurve(const CubicPolinomialCurve& cpc);
        CubicPolinomialCurve& operator =(const CubicPolinomialCurve& cpc);
        GLboolean BlendingFunctionValues(GLdouble u, RowMatrix<GLdouble> &values) const;
        GLboolean CalculateDerivatives(GLuint max_order_of_derivatives, GLdouble u, Derivatives &d) const;

    };

    class QuarticPolinomialCurve: public LinearCombination3 {
    public:
        QuarticPolinomialCurve();
        QuarticPolinomialCurve(const QuarticPolinomialCurve& cpc);
        QuarticPolinomialCurve& operator =(const QuarticPolinomialCurve& cpc);
        GLboolean BlendingFunctionValues(GLdouble u, RowMatrix<GLdouble> &values) const;
        GLboolean CalculateDerivatives(GLuint max_order_of_derivatives, GLdouble u, Derivatives &d) const;

    };

    class SecondOrderTrigonometricCurve: public LinearCombination3 {
    protected:
        GLdouble _beta;
        //nu = 4
        //mu = 2
    public:
        SecondOrderTrigonometricCurve(GLdouble beta);
        SecondOrderTrigonometricCurve(const SecondOrderTrigonometricCurve& cpc);
        SecondOrderTrigonometricCurve& operator =(const SecondOrderTrigonometricCurve& cpc);
        GLboolean BlendingFunctionValues(GLdouble u, RowMatrix<GLdouble> &values) const;
        GLboolean CalculateDerivatives(GLuint max_order_of_derivatives, GLdouble u, Derivatives &d) const;

    };

    class SecondOrderHyperbolicCurve: public LinearCombination3 {
    protected:
        GLdouble _beta;
    public:
        SecondOrderHyperbolicCurve(GLdouble beta);
        SecondOrderHyperbolicCurve(const SecondOrderHyperbolicCurve& cpc);
        SecondOrderHyperbolicCurve& operator =(const SecondOrderHyperbolicCurve& cpc);
        GLboolean BlendingFunctionValues(GLdouble u, RowMatrix<GLdouble> &values) const;
        GLboolean CalculateDerivatives(GLuint max_order_of_derivatives, GLdouble u, Derivatives &d) const;

    };

    class FirstOrderAlgebraicTrigonometricCurve: public LinearCombination3 {
    protected:
        GLdouble _beta;
    public:
        FirstOrderAlgebraicTrigonometricCurve(GLdouble beta);
        FirstOrderAlgebraicTrigonometricCurve(const FirstOrderAlgebraicTrigonometricCurve& cpc);
        FirstOrderAlgebraicTrigonometricCurve& operator =(const FirstOrderAlgebraicTrigonometricCurve& cpc);
        GLboolean BlendingFunctionValues(GLdouble u, RowMatrix<GLdouble> &values) const;
        GLboolean CalculateDerivatives(GLuint max_order_of_derivatives, GLdouble u, Derivatives &d) const;

    };
}

#endif // BOUNDARYCURVES_H
