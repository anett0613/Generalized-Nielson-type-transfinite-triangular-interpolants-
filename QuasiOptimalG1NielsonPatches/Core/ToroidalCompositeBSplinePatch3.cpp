#include "ToroidalCompositeBSplinePatch3.h"
#include "Constants.h"

namespace cagd{
ToroidalCompositeBSplineSurface3::ToroidalCompositeBSplineSurface3(GLdouble r, GLdouble R, GLuint n, GLuint m):
    BicubicBSplineCompositeSurface3(n,m)
{
    _r = r;
    _R = R;
    setControlPoints();
    generateImage();
    updateIsoLines();

}
GLvoid ToroidalCompositeBSplineSurface3::setControlPoints(){
    DCoordinate3 point;
    GLuint n = _control_net.GetRowCount(), m = _control_net.GetColumnCount();
    GLdouble uStep = TWO_PI / (n), vStep = TWO_PI / (m);
    GLdouble u = 0, v = 0;
    //GLdouble a = 0.5*(_R - _r), c = 0.5*(_r + _R);

    for (GLuint i = 0; i < n; i++) {
      for (GLuint j = 0; j < m; j++) {
          point[0] = -(_R + _r * sin(u)) * sin(v);
          point[1] =  (_R + _r * sin(u)) * cos(v);
          point[2] = _r * cos(u);
        _control_net(i, j) = point;
        v = j*vStep;
      }
      v = 0;
      u = i*uStep;
    }
}
}
