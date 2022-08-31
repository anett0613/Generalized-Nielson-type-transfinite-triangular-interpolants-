#include "TriangularSurfaces3.h"
#include <algorithm>

using namespace cagd;
using namespace std;

TriangularSurface3::PartialDerivatives::PartialDerivatives(GLuint maximum_order_of_partial_derivatives):
    TriangularMatrix<DCoordinate3>(maximum_order_of_partial_derivatives + 1)
{

}

GLvoid TriangularSurface3::PartialDerivatives::LoadNullVectors()
{
    for (GLuint i = 0; i < _row_count; i++) {
        for (GLuint j = 0; j <= i; j++) {
            for (GLuint k = 0; k < 3; k++) {
                _data[i][j][k] = 0;
            }
        }
    }
}

//special contructor
TriangularSurface3::TriangularSurface3(GLdouble beta, GLuint row_count): _beta(beta), _data(row_count)
{
    //_data = TriangularMatrix<DCoordinate3>(row_count);
}

//copy constructor
TriangularSurface3::TriangularSurface3(const TriangularSurface3& surface):
    _beta(surface._beta), _data(TriangularMatrix<DCoordinate3>(surface._data))
{
    _vbo_vertices = 0;
    _vbo_indices = 0;

    if (surface._vbo_vertices && surface._vbo_indices) {
        UpdateVertexBufferObjectsOfData();
    }
}

TriangularSurface3& TriangularSurface3::operator =(const TriangularSurface3& surface)
{
    if (this != &surface) {

        DeleteVertexBufferObjectsOfData();

        _beta = surface._beta;

        _data = surface._data;

        if (surface._vbo_vertices && surface._vbo_indices)
            UpdateVertexBufferObjectsOfData();
    }

    return *this;
}

GLvoid TriangularSurface3::SetBeta(GLdouble beta)
{
    if (_beta != beta)
    {
        _beta = beta;
    }
}

GLdouble TriangularSurface3::GetBeta(GLdouble& u_min, GLdouble& u_max) const
{
    if (_beta >= u_min && _beta <= u_max)
    {
        return _beta;
    }
    else
    {
      cout << _beta << " is not in the [ " << u_min << ", " << u_max << " ] interval" <<endl;
      return 0;
    }
}

GLboolean TriangularSurface3::SetData(GLuint row, GLuint column, GLdouble x, GLdouble y, GLdouble z)
{
    if (row < column || row >= _data.GetRowCount() || row < 0 || column < 0) {
        return GL_FALSE;
    }

    _data(row, column).x() = x;
    _data(row, column).y() = y;
    _data(row, column).z() = z;

    return GL_TRUE;
}

GLboolean TriangularSurface3::SetData(GLuint row, GLuint column, const DCoordinate3& point)
{
    if (row < column || row >= _data.GetRowCount() || row < 0 || column < 0) {
        return GL_FALSE;
    }

    _data(row, column) = point;

    return GL_TRUE;
}

GLboolean TriangularSurface3::GetData(GLuint row, GLuint column, GLdouble& x, GLdouble& y, GLdouble& z) const
{
    if (row < column || row >= _data.GetRowCount() || row < 0 || column < 0) {
        return GL_FALSE;
    }

    x = _data(row, column).x();
    y = _data(row, column).y();
    z = _data(row, column).z();

    return GL_TRUE;
}

GLboolean TriangularSurface3::GetData(GLuint row, GLuint column, DCoordinate3& point) const
{
    if (row < column || row >= _data.GetRowCount() || row < 0 || column < 0) {
        return GL_FALSE;
    }

    point = _data(row, column);

    return GL_TRUE;
}

DCoordinate3 TriangularSurface3::operator ()(GLuint row, GLuint column) const
{
    return _data(row,column);
}

DCoordinate3& TriangularSurface3::operator ()(GLuint row, GLuint column)
{
    return _data(row,column);
}

GLvoid TriangularSurface3::DeleteVertexBufferObjectsOfData()
{
    /*if (_vbo_data)
    {
        glDeleteBuffers(1, &_vbo_data);
        _vbo_data = 0;
    }*/

    if (_vbo_indices)
    {
        glDeleteBuffers(1, &_vbo_indices);
        _vbo_indices = 0;
    }

    if (_vbo_vertices)
    {
        glDeleteBuffers(1, &_vbo_vertices);
        _vbo_vertices = 0;
    }
}

GLboolean TriangularSurface3::RenderData(GLenum render_mode) const
{
    if (!_vbo_vertices || !_vbo_indices)
    {
        cout<<"render data"<<endl;
        return GL_FALSE;
    }
    if (render_mode != GL_POINTS && render_mode != GL_LINE_STRIP && render_mode != GL_LINE_LOOP && render_mode != GL_TRIANGLES)
    {
        return GL_FALSE;
    }

   /* for (GLuint r = 0; r <= 3; r++)
    {
        glBegin(GL_LINE_STRIP);
            for (GLuint c = 0; c <= r; c++)
            {
                glVertex3dv(&_data(r, c)[0]);
            }
        glEnd();
    }

    glBegin(GL_LINE_STRIP);
        for (GLuint r = 0; r <= 3; r++)
        {
            glVertex3dv(&_data(r, 0)[0]);
        }
    glEnd();

    glBegin(GL_LINE_STRIP);
        for (GLuint r = 0; r <= 3; r++)
        {
            glVertex3dv(&_data(r, r)[0]);
        }
    glEnd();

    return GL_TRUE;*/

    GLuint r = _data.GetRowCount();
    GLuint face_count = (r - 1) * (r - 1) + 2;

    glEnableClientState(GL_VERTEX_ARRAY);

        // activate the VBO of vertices
        glBindBuffer(GL_ARRAY_BUFFER, _vbo_vertices);
        // specify the location and data format of vertices
        glVertexPointer(3, GL_FLOAT, 0, (const GLvoid *)0);

        // activate the element array buffer for indexed vertices of triangular faces
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _vbo_indices);

        // render primitives
        glDrawElements(render_mode, 3 * (GLsizei)face_count, GL_UNSIGNED_INT, (const GLvoid *)0);

   /* GLuint offset = 0;
    for (GLuint i = 0; i < _data.GetRowCount(); i++) {
        glDrawElements(render_mode, offset,);
        //glDrawArrays(render_mode, offset, i + 1);
        offset += i + 1;
        cout<<offset<<endl;
    }*/

    glDisableClientState(GL_VERTEX_ARRAY);

    // unbind any buffer object previously bound and restore client memory usage
    // for these buffer object targets
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    return GL_TRUE;
}

GLboolean TriangularSurface3::UpdateVertexBufferObjectsOfData(GLenum usage_flag)
{
    if (usage_flag != GL_STREAM_DRAW && usage_flag != GL_STREAM_READ && usage_flag != GL_STREAM_COPY &&
        usage_flag != GL_DYNAMIC_DRAW && usage_flag != GL_DYNAMIC_READ && usage_flag != GL_DYNAMIC_COPY &&
        usage_flag != GL_STATIC_DRAW && usage_flag != GL_STATIC_READ && usage_flag != GL_STATIC_COPY) {

        return GL_FALSE;
    }

    DeleteVertexBufferObjectsOfData();

   /* glGenBuffers(1, &_vbo_data);

    if (!_vbo_data) {
        return GL_FALSE;
    }*/

    glGenBuffers(1, &_vbo_vertices);

    if (!_vbo_vertices)
        return GL_FALSE;

    glGenBuffers(1, &_vbo_indices);
    if (!_vbo_indices) {
        glDeleteBuffers(1, &_vbo_vertices);
        _vbo_vertices = 0;

        return GL_FALSE;
    }

    GLuint r = _data.GetRowCount();
    GLuint vertex_byte_size = 3 * r * (r + 1) / 2 * sizeof(GLfloat);

    glBindBuffer(GL_ARRAY_BUFFER, _vbo_vertices);
    glBufferData(GL_ARRAY_BUFFER, vertex_byte_size, 0, usage_flag);

    GLfloat* coordinate = (GLfloat*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    /*if (!coordinate) {
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        DeleteVertexBufferObjectsOfData();
        return GL_FALSE;
    }*/

    for (GLuint i = 0; i < r; i++) {
        for (GLuint j = 0; j <= i; j++) {
            for (GLuint component = 0; component < 3; component++) {
                *coordinate = (GLfloat)_data(i, j)[component];
                ++coordinate;
            }
        }
    }

    GLuint index_byte_size = 3 * ((r - 1) * (r - 1) + 2)* sizeof(GLuint);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _vbo_indices);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, index_byte_size, 0, usage_flag);
    GLuint* element = (GLuint*)glMapBuffer(GL_ELEMENT_ARRAY_BUFFER, GL_WRITE_ONLY);

    for (GLuint i = 0; i < r - 1; i++) {
        for (GLuint j = 0; j <= i; j++) {
            GLuint index[4];

            index[0] = i * (i + 1) / 2 + j;
            index[1] = index[0] + i + 1;
            index[2] = index[1] + 1;
            index[3] = index[0] + 1;

            *element = index[1];
            ++element;

            *element = index[2];
            ++element;

            *element = index[0];
            ++element;

            if (j < i) {
                *element = index[0];
                ++element;

                *element = index[2];
                ++element;

                *element = index[3];
                ++element;
            }
            if(i==j){
                *element = index[1];
                ++element;

                *element = index[1] - 1;
                ++element;
            }

        }
        /*if(i == 1) {
            *element = i * (i + 1) / 2;
            ++element;

            *element = (i + 1) * (i + 2) / 2;
            ++element;
        }*/
    }

    glBindBuffer(GL_ARRAY_BUFFER, _vbo_vertices);
    if (!glUnmapBuffer(GL_ARRAY_BUFFER))
        return GL_FALSE;

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _vbo_indices);
    if (!glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER))
        return GL_FALSE;

    // unbind any buffer object previously bound and restore client memory usage
    // for these buffer object targets
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    return GL_TRUE;
}

TriangulatedMesh3* TriangularSurface3::GenerateImage(GLuint div_point_count, GLenum usage_flag) const
{
    if (div_point_count <= 1)
        return GL_FALSE;

    // calculating number of vertices, unit normal vectors and texture coordinates
    GLuint vertex_count = div_point_count * (div_point_count + 1) / 2;

    // calculating number of triangular faces
    GLuint face_count = (div_point_count - 1) * (div_point_count - 1);

    TriangulatedMesh3 *result = nullptr;

    result = new TriangulatedMesh3(vertex_count, face_count, usage_flag);

    if (!result)
    {
        return nullptr;
    }

    // for face indexing
    GLuint current_face = 0;

    GLdouble ad = 1.0 / (div_point_count - 1);
    GLfloat sd = 1.0f / (div_point_count - 1);

    DCoordinate3 pi(_beta, 0, 0);
    DCoordinate3 pj(0, _beta, 0);
    DCoordinate3 pk(0, 0, _beta);

    GLboolean is_parallel = GL_TRUE;
    GLboolean aborted = GL_FALSE;

    if(is_parallel){

        #pragma omp parallel for
        for (GLint i_j = 0; i_j < (GLint)vertex_count; i_j++)
        {

            #pragma omp flush(aborted)
            if (!aborted)
            {
                //GLint i = i_j * 2 / (div_point_count + 1) + 1;
                GLint i = (GLint)sqrt(i_j * 2);
                if (i * (i + 1) / 2 > i_j)
                {
                    i = i - 1;
                }

                GLint j = i_j - i * (i + 1) / 2;

                PartialDerivatives pd;

                GLdouble a = min(i * ad, 1.0);
                GLfloat s = min(i * sd, 1.0f);
                GLdouble db = (i == 0) ? 0.0 : 1.0 / i;

                GLfloat t = min(j * sd, 1.0f);
                GLdouble b = min(db * j, 1.0);

                DCoordinate3 z((1.0 - b) * ((1.0 - a) * pi + a * pj) + b * ((1.0 - a) * pi + a * pk));

                if(!CalculatePartialDerivatives(1, z[0], z[1], pd))
                {
                    aborted = GL_TRUE;
                    #pragma omp flush(aborted)
                }


                /*
                 * 0-3
                 * |\|
                 * 1-2
                */
                GLuint index[4];

                index[0] = i * (i + 1) / 2 + j;
                index[1] = index[0] + i + 1;
                index[2] = index[1] + 1;
                index[3] = index[0] + 1;

                // surface point
                (*result)._vertex[index[0]] = pd(0, 0);


                // unit surface normal
                (*result)._normal[index[0]] = pd(1, 0);
                (*result)._normal[index[0]] ^= pd(1, 1);
                (*result)._normal[index[0]].normalize();

                // texture coordinates
                (*result)._tex[index[0]].s() = s;
                (*result)._tex[index[0]].t() = t;

                #pragma omp critical
                // faces
                if (i < (GLint)div_point_count - 1 && j < (GLint)div_point_count - 1)
                {
                    //#pragma omp critical
                    (*result)._face[current_face][0] = index[0];
                    (*result)._face[current_face][1] = index[1];
                    (*result)._face[current_face][2] = index[2];


                    current_face++;

                    if (j < i)
                    {
                        (*result)._face[current_face][0] = index[0];
                        (*result)._face[current_face][1] = index[2];
                        (*result)._face[current_face][2] = index[3];

                        //#pragma omp atomic
                        current_face++;
                    }
                }
            }
        }

        if (aborted)
        {
            delete result;
            result = nullptr;
        }

    } else
    {
        // partial derivatives of order 0, 1, 2, and 3
        PartialDerivatives pd;
        for (GLuint i = 0; i < div_point_count; i++)
        {
            GLdouble a = min(i * ad, 1.0);
            GLfloat s = min(i * sd, 1.0f);
            GLdouble db = (i == 0) ? 0.0 : 1.0 / i;
            for (GLuint j = 0; j < i + 1; j++)
            {
                GLfloat t = min(j * sd, 1.0f);
                GLdouble b = min(db * j, 1.0);

                DCoordinate3 z((1.0 - b) * ((1.0 - a) * pi + a * pj) + b * ((1.0 - a) * pi + a * pk));

                if(!CalculatePartialDerivatives(1, z[0], z[1], pd))
                {
                    delete result;
                    result = nullptr;
                    return result;
                }


                /*
                 * 0-3
                 * |\|
                 * 1-2
                */
                GLuint index[4];

                index[0] = i * (i + 1) / 2 + j;
                index[1] = index[0] + i + 1;
                index[2] = index[1] + 1;
                index[3] = index[0] + 1;

                // surface point
                (*result)._vertex[index[0]] = pd(0, 0);


                // unit surface normal
                (*result)._normal[index[0]] = pd(1, 0);
                (*result)._normal[index[0]] ^= pd(1, 1);
                (*result)._normal[index[0]].normalize();

                // texture coordinates
                (*result)._tex[index[0]].s() = s;
                (*result)._tex[index[0]].t() = t;

                // faces
                if (i < div_point_count - 1 && j < div_point_count - 1)
                {
                    (*result)._face[current_face][0] = index[0];
                    (*result)._face[current_face][1] = index[1];
                    (*result)._face[current_face][2] = index[2];
                    current_face++;

                    if (j < i)
                    {
                        (*result)._face[current_face][0] = index[0];
                        (*result)._face[current_face][1] = index[2];
                        (*result)._face[current_face][2] = index[3];
                        current_face++;
                    }
                }
            }
        }
    }


    return result;
}

TriangularSurface3::~TriangularSurface3()
{
    DeleteVertexBufferObjectsOfData();
}
