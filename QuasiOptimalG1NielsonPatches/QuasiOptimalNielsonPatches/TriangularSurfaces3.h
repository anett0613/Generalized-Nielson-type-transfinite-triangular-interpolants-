#ifndef TRIANGULARSURFACES3_H
#define TRIANGULARSURFACES3_H

#include "../Core/DCoordinates3.h"
#include <GL/glew.h>
#include <iostream>
#include "../Core/Matrices.h"
#include "../Core/TriangulatedMeshes3.h"
#include <vector>

namespace cagd
{
    class TriangularSurface3
    {
    public:
        // a nested class the stores the zeroth and higher order partial derivatives associated with a
        // surface point
        class PartialDerivatives: public TriangularMatrix<DCoordinate3>
        {
        public:
            PartialDerivatives(GLuint maximum_order_of_partial_derivatives = 1);

            // homework: initializes all partial derivatives to the origin
            GLvoid LoadNullVectors();
        };


    protected:
       // GLuint                         _vbo_data;            // vertex buffer object of the control net
        GLuint                         _vbo_indices;
        GLuint                         _vbo_vertices;
        GLdouble                       _beta;                // definition domain in direction u
        TriangularMatrix<DCoordinate3> _data;                // the control net (usually stores position vectors)

    public:
        // homework: special constructor
        TriangularSurface3(GLdouble beta = 1.0, GLuint row_count = 4);

        // homework: copy constructor
        TriangularSurface3(const TriangularSurface3& surface);

        // homework: assignment operator
        TriangularSurface3& operator =(const TriangularSurface3& surface);

        // homework: set/get the definition domain of the surface
        GLvoid SetBeta(GLdouble beta);

        GLdouble GetBeta(GLdouble& u_min, GLdouble& u_max) const;

        // homework: set coordinates of a selected data point
        GLboolean SetData(GLuint row, GLuint column, GLdouble x, GLdouble y, GLdouble z);
        GLboolean SetData(GLuint row, GLuint column, const DCoordinate3& point);

        // homework: get coordinates of a selected data point
        GLboolean GetData(GLuint row, GLuint column, GLdouble& x, GLdouble& y, GLdouble& z) const;
        GLboolean GetData(GLuint row, GLuint column, DCoordinate3& point) const;


        // homework: get data by value
        DCoordinate3 operator ()(GLuint row, GLuint column) const;

        // homework: get data by reference
        DCoordinate3& operator ()(GLuint row, GLuint column);

        // eq. (25)
        virtual GLboolean CalculatePartialDerivatives(
                GLuint maximum_order_of_partial_derivatives,
                GLdouble u, GLdouble v, PartialDerivatives& pd) const = 0;

        // generates a triangulated mesh that approximates the shape of the surface above
        virtual TriangulatedMesh3* GenerateImage(
                GLuint div_point_count,
                GLenum usage_flag = GL_STATIC_DRAW) const;

        // homework: VBO handling methods
        virtual GLvoid    DeleteVertexBufferObjectsOfData();
        virtual GLboolean RenderData(GLenum render_mode = GL_LINE_STRIP) const;
        virtual GLboolean UpdateVertexBufferObjectsOfData(GLenum usage_flag = GL_STATIC_DRAW);


        // homework: destructor
        virtual ~TriangularSurface3();
    };
}


#endif // TRIANGULARSURFACES3_H
