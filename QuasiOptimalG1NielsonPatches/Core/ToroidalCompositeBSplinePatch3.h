#include "Bicu_bic_B_Spline_CompositePatch3.h"

namespace cagd{

class ToroidalCompositeBSplineSurface3: public BicubicBSplineCompositeSurface3
{
protected:
   GLdouble _r, _R; // a tórusz paraméteres képletében megjelenő kis és nagy sugár

public:
    ToroidalCompositeBSplineSurface3(GLdouble r, GLdouble R, GLuint n, GLuint m); // lásd a fent említett 37. oldalt, n+1 = row_count, m+1 = column_count
    GLvoid setControlPoints();
};
}
