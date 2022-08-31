#include "Bicu_bic_B_Spline_Patch3.h"
namespace cagd {

BicubicBSplinePatch3::BicubicBSplinePatch3()
  : TensorProductSurface3(0, 1, 0, 1)
{}

GLboolean BicubicBSplinePatch3::UBlendingFunctionValues( GLdouble u_knot,RowMatrix<GLdouble>& blending_values) const
{

  if (u_knot < 0.0 && u_knot > 1.0)
    return GL_FALSE;

  blending_values.ResizeColumns(4);
  GLdouble u = u_knot, u2 = u * u, u3 = u2 * u;
  GLdouble w1 = 1 - u, w2 = w1 * w1, w3 = w1 * w2;

  blending_values(0) = w3 / 6.0;
  blending_values(1) = (3 * u * w2 + 3 * w1 + 1) / 6.0;
  blending_values(2) = (3 * u2 * w1 + 3 * u + 1) / 6.0;
  blending_values(3) = u3 / 6.0;

  return GL_TRUE;
}

GLboolean BicubicBSplinePatch3::VBlendingFunctionValues(GLdouble v_knot, RowMatrix<GLdouble>& blending_values) const
{

  if (v_knot < 0.0 && v_knot > 1.0)
    return GL_FALSE;

  blending_values.ResizeColumns(4);
  GLdouble v = v_knot, v2 = v * v, v3 = v2 * v;
  GLdouble w1 = 1 - v, w2 = w1 * w1, w3 = w1 * w2;

  blending_values(0) = w3 / 6.0;
  blending_values(1) = (3 * v * w2 + 3 * w1 + 1) / 6.0;
  blending_values(2) = (3 * v2 * w1 + 3 * v + 1) / 6.0;
  blending_values(3) = v3 / 6.0;

  return GL_TRUE;
}

GLboolean BicubicBSplinePatch3::CalculatePartialDerivatives(GLuint maximum_order_of_partial_derivatives,
                                                  GLdouble u, GLdouble v, PartialDerivatives& pd) const
{

  if (u < 0.0 || u > 1.0 || v < 0.0 || v > 1.0 ||
      maximum_order_of_partial_derivatives > 1)
    return GL_FALSE;

  RowMatrix<GLdouble> u_blending_values(4), d1_u_blending_values(4);

  GLdouble u2 = u * u, u3 = u2 * u;
  GLdouble w1 = 1 - u, w2 = w1 * w1, w3 = w1 * w2;

  u_blending_values(0) = w3 / 6.0;
  u_blending_values(1) = (3 * u * w2 + 3 * w1 + 1) / 6.0;
  u_blending_values(2) = (3 * u2 * w1 + 3 * u + 1) / 6.0;
  u_blending_values(3) = u3 / 6.0;

  d1_u_blending_values(0) = -w2 / 2.0; //((1-u)^3/6)'
  d1_u_blending_values(1) = (w2 - 2.0 * u * w1 - 1.0) / 2.0; //(u * (3 * u - 4)) / 2.0;
  d1_u_blending_values(2) = -(u2 - 2.0 * w1 * u - 1.0) / 2.0; //(3.0 * u2 - 2.0 * u - 1) / -2.0;
  d1_u_blending_values(3) = u2 / 2.0;//(u^3/6)'

  RowMatrix<GLdouble> v_blending_values(4), d1_v_blending_values(4);

  GLdouble v2 = v * v, v3 = v2 * v;
  w1 = 1 - v;
  w2 = w1 * w1;
  w3 = w1 * w2;

  v_blending_values(0) = w3 / 6.0;
  v_blending_values(1) = (3 * v * w2 + 3 * w1 + 1) / 6.0;
  v_blending_values(2) = (3 * v2 * w1 + 3 * v + 1) / 6.0;
  v_blending_values(3) = v3 / 6.0;

  d1_v_blending_values(0) = -w2 / 2.0;
  d1_v_blending_values(1) = (w2 - 2.0 * v * w1 - 1.0) / 2.0; //(v * (3 * v - 4)) / 2.0;
  d1_v_blending_values(2) = -(v2 - 2.0 * w1 * v - 1.0) / 2.0; //(3.0 * v2 - 2.0 * v - 1) / -2.0;
  d1_v_blending_values(3) = v2 / 2.0;

  pd.ResizeRows(2);
  pd.LoadNullVectors();

  for (GLuint row = 0; row < 4; row++) {

    DCoordinate3 aux_d0_v, aux_d1_v;
    for (GLuint column = 0; column < 4; column++) {
      aux_d0_v += _data(row, column) * v_blending_values(column);
      aux_d1_v += _data(row, column) * d1_v_blending_values(column);
    }
    pd(0, 0) += aux_d0_v * u_blending_values(row);
    pd(1, 0) += aux_d0_v * d1_u_blending_values(row);
    pd(1, 1) += aux_d1_v * u_blending_values(row);
  }
  return GL_TRUE;
}
}
