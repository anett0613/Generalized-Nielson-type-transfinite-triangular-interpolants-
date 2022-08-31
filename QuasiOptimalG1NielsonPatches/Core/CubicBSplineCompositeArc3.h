#include "CubicBSplineArc3.h"
#include <Core/Colors4.h>

namespace cagd
{
    class SomeCompositeCurve3
    {
    public:
        enum Direction{LEFT, RIGHT};

        class ArcAttributes
        {
        public:
            SomeArch3        *arc;
            GenericCurve3   *image;
            Color4          *color;
            ArcAttributes   *previous, *next; // they are initializes to nullptr by default

            // TO DO: constructor, copy_constructor(deep copy), operator = (deep copy), destructor

        };

    protected:
        RowMatrix<DCoordinate3>     _big_control_polygon; // 9eshez

        // GLdouble                    _alpha;
        std::vector<ArcAttributes>  _attributes;

        // _attributes[i].arc = new SomeArc3(alpha);
        // (*_attributes[i].arc)[j] // p_j, j = 0, 1, 2, 3 <==> LinearCombination3::cagd
        // _attributes[i].image = _attributes[i].arc->GenerateImage(..., ...);
        // !!! _attrbiutes.resize(_attributes.size() + 1);

        // _attribute.image = _attributes[i].arc->GenerateImage(..., ...);
        //  glColor4fv(_attrbiutes[i].color)[0];
        // _attribute.image->RenderDerivatives(..., ...);

    public:
        SomeCompositeCurve3(const GLdouble &alpha = PI / 2.0, const size_t &minimum_arc_count = 1000); // _attributes.reserve(minimum_arc_count) -> es nem resize !

        /* nem kell a 9es projekthez
        GLboolean insertNewIsolatedArc();
        GLboolean continueExistingArc(const GLuint &arc_index, Direction direction);
        GLboolean joinExistingArcs(const GLuint &arc_index_1, Direction direction_1,
                                   const GLuint &arc_index_2, Direction direction_2);
        GLboolean mergeExistingArcs(const GLuint &arc_index_1, Direction direction_1,
                                    const GLuint &arc_index_2, Direction direction_2);
        */

        //setters, getters
        //renderAllArcs
        //renderSelectedArc
        //update methods
    };

}
