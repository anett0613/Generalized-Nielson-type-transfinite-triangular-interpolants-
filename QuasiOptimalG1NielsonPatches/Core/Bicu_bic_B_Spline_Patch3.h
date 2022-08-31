#include <Core/TensorProductSurfaces3.h>

namespace cagd
{
    class BicubicBSplinePatch3: public TensorProductSurface3
    {
    protected:
        //GLdouble _u_alpha, _v_alpha;

    public:
        BicubicBSplinePatch3();

        GLboolean UBlendingFunctionValues(GLdouble u_knot,RowMatrix<GLdouble>& blending_values) const override;
        GLboolean VBlendingFunctionValues(GLdouble v_knot,RowMatrix<GLdouble>& blending_values) const override;
        GLboolean CalculatePartialDerivatives(GLuint maximum_order_of_partial_derivatives, GLdouble u, GLdouble v,
                                              PartialDerivatives& pd) const override;
    };
}
