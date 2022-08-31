#ifndef TRIANGULATEDHALFEDGEDATASTRUCTURES_H
#define TRIANGULATEDHALFEDGEDATASTRUCTURES_H

#include "./Core/Constants.h"
#include "../Core/DCoordinates3.h"
#include "../Core/TriangulatedMeshes3.h"
#include "../Core/LinearCombination3.h"
#include "BoundaryCurves.h"
#include "../Core/RealSquareMatrices.h"
#include "TriangularSurfaces3.h"
#include "SpecialTriangularSurfaces3.h"
#include "../Core/Matrices.h"
#include "../Core/ShaderPrograms.h"
#include "math.h"
#include <map>
#include <omp.h>



// other include files

namespace cagd
{
    class TriangulatedHalfEdgeDataStructure3
    {
        friend class GLWidget;
    public:


        // represents a face of the graph
        class Face
        {
            // output to stream
            friend std::ostream& operator <<(std::ostream& lhs, const Face &rhs);

            // input from stream
            friend std::istream& operator >>(std::istream& lhs, Face &rhs);

        public:
            GLint              node[3];
            TriangularSurface3 *C0;     // points to the C^{0} optimal triangular surface patch (25)
            TriangulatedMesh3  *img_C0; // references the image of the C^{0} optimal triangular surface patch:  energy functional optimization -> control net -> (25)->GenerateImage
            TriangulatedMesh3  *G1;     // points to the image of the final G^{1} quasi-optimal surface: (37) <- convex blending of (34), (35) and (36)

            // constructor
            Face();

            // get node id by value
            GLint operator [](GLint i) const;

            // get node id by reference
            GLint& operator [](GLint i);

            // destructor: deletes the memory areas referenced by the pointers C0, img_C0, G1
            ~Face();
        };

        // describes an oriented half-edge of the graph
        class HalfEdge
        {
            friend class TriangulatedHalfEdgeDataStructure3;
        public:
            DCoordinate3        *vertex, *normal; // simple copy
            HalfEdge            *opposite, *next; // simple copy
            DCoordinate3        tangent;          // simple copy
            LinearCombination3  *boundary = nullptr;        // deep copy
            GenericCurve3       *img_boundary = nullptr;    // deep copy
            bool                is_img_boundary_renderable = true;
            Face                *face;            // simple copy

            // constructor (sets every pointer to nullptr)
            HalfEdge();

            // you should also write copy constructor and assignment operator
            HalfEdge(const HalfEdge& half_edge);
            HalfEdge& operator =(const HalfEdge& rhs);

            // destructor: deletes the memory areas referenced by the pointers boundary and img_boundary
            ~HalfEdge();
        };
        //ezekre kulon osztalyok
        enum BoundaryType{
            CUBIC_POLYNOMIAL,
            QUARTIC_POLYNOMIAL,
            SECOND_ORDER_TRIGONOMETRIC,
            SECOND_ORDER_HYPERBOLIC,
            FIRST_ORDER_ALGEBRAIC_TRIGONOMETRIC};

    protected:
        std::vector<DCoordinate3>                       _vertex;  // input data loaded from OFF
        std::vector<Face>                               _face;    // input data loaded from OFF
        std::vector<DCoordinate3>                       _normal;  // avareged unit normals calculated
                                                                  // based on the input data

        // left- and rightmost corners of the bounding box (calculated based on input data)
        DCoordinate3                                    _leftmost_vertex;
        DCoordinate3                                    _rightmost_vertex;

        std::map< std::pair<GLint, GLint>, HalfEdge* >  _edges;  // half-edge data structure

    private:
        GLdouble fact(GLuint n);
        GLdouble combination(GLuint n, GLuint k);

    public:
        // its behavior is similar to that of the method TriangulatedMesh3::LoadFromOFF, but you should also
        // build the map _edges
        GLboolean LoadFromOFF(const std::string& file_name, GLboolean translate_and_scale_to_unit_cube = GL_FALSE);

        bool GenerateBoundaryCurves(
                BoundaryType type, GLdouble beta,
                GLdouble theta_1, GLdouble theta_2, GLdouble theta_3,
                GLuint div_point_count, GLenum usage_flag = GL_STATIC_DRAW);

        void GenerateC0TriangularSurfaces(
                BoundaryType type, GLdouble beta,
                GLdouble epsilon_1, GLdouble epsilon_2, GLdouble epsilon_3,
                GLuint div_point_count, GLenum usage_flag = GL_STATIC_DRAW);

        void GenerateG1TriangularSurfaces(
                BoundaryType type, GLdouble beta,
                GLdouble theta_1, GLdouble theta_2, GLdouble theta_3,
                GLuint div_point_count, GLenum usage_flag = GL_STATIC_DRAW);

        // rendering methods, you may define other ones as well
        GLboolean RenderCurves(GLuint order, GLenum render_mode = GL_LINE_STRIP);
        GLboolean RenderControlPolygons(GLenum render_mode = GL_LINE_STRIP);

        GLboolean RenderHalfEdges();
        GLboolean RenderFaces();

        GLboolean RenderC0Meshes(GLenum render_mode = GL_TRIANGLES);
        GLboolean RenderC0ControlNets();
        GLboolean RenderG1Meshes(GLenum render_mode = GL_TRIANGLES);
        GLboolean RenderVertexNormals(GLdouble scale = 1.0);
        GLboolean RenderC0BoundaryNormals(GLdouble scale = 1.0);
        GLboolean RenderG1BoundaryNormals(GLdouble scale = 1.0);
        GLboolean RenderVertices();
        GLboolean RenderApproximatedUnitTangents(GLdouble scale = 1.0);

        // other setters/getters
        // ...

        // clean-up methods
        GLvoid DeleteBoundaryCurves();
        GLvoid DeleteC0Surfaces();
        GLvoid DeleteG1Surfaces();
        GLvoid CleanAll();

        // destructor
        ~TriangulatedHalfEdgeDataStructure3();

    };
}

#endif // TRIANGULATEDHALFEDGEDATASTRUCTURES_H
