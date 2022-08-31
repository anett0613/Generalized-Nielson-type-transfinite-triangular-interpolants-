#include <Core/Constants.h> //Ha van alappaparameter
#include <Core/LinearCombination3.h>

namespace cagd
{
    class CubicBSplineArc3: public LinearCombination3
    {
    protected:
        GLdouble _alpha; //Ha van parametere a kepletnek

    public:
        CubicBSplineArc3(/*const GLdouble &alpha = PI / 2.0 */);

        //values[i] = F_i(u); i = 0,1,2,3
        GLboolean BlendingFunctionValues(GLdouble u, RowMatrix<GLdouble>& values) const;

        //at least for max_order_of_derivatives = 2
        //d[0] = curve point
        //d[1] = tangent vector
        //d[2] = acceleration vector
        GLboolean CalculateDerivatives(GLuint max_order_of_derivatives, GLdouble u, Derivatives& d) const;

        GLboolean CalculateDerivative(GLuint order_of_derivative, GLdouble u, RowMatrix<GLdouble> &derivative) const;
        GLboolean Derivative0(GLdouble u, RowMatrix<GLdouble>& derivative) const;
        GLboolean Derivative1(GLdouble u, RowMatrix<GLdouble>& derivative) const;
        GLboolean Derivative2(GLdouble u, RowMatrix<GLdouble>& derivative) const;

        //GLboolean SetAlpha(const GLdouble &alpha;

        //other setters/getters...


    };

}
