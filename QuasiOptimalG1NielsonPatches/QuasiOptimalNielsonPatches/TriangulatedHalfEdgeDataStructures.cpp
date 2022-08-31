#include "./TriangulatedHalfEdgeDataStructures.h"

#include <algorithm>
#include <fstream>

using namespace cagd;
using namespace std;

TriangulatedHalfEdgeDataStructure3::Face::Face()
{
    node[0] = node[1] = node[2] = -1; // non-existing vertex ids
    C0 = nullptr;
    img_C0 = nullptr;
    G1 = nullptr;
}

// + copy ctor, op=

GLint TriangulatedHalfEdgeDataStructure3::Face::operator[](GLint i) const
{
    return node[i];
}

GLint& TriangulatedHalfEdgeDataStructure3::Face::operator[](GLint i)
{
    return node[i];
}


// output to stream
std::ostream& cagd::operator <<(std::ostream& lhs, const TriangulatedHalfEdgeDataStructure3::Face& rhs)
{
    lhs << 3;
    for (GLuint i = 0; i < 3; ++i)
        lhs  << " " << rhs[i];
    return lhs;
}

// homework
std::istream& cagd::operator >>(std::istream& lhs, TriangulatedHalfEdgeDataStructure3::Face& rhs) {
    GLuint three;
    return lhs >> three >> rhs[0] >> rhs[1] >> rhs[2];
}


TriangulatedHalfEdgeDataStructure3::Face::~Face()
{
    if (C0)
    {
        delete C0; C0 = nullptr;
    }

    if (img_C0)
    {
        delete img_C0; img_C0 = nullptr;
    }

    if (G1)
    {
        delete G1; G1 = nullptr;
    }
}

TriangulatedHalfEdgeDataStructure3::HalfEdge::HalfEdge():
    vertex(nullptr), normal(nullptr), opposite(nullptr), next(nullptr), tangent(DCoordinate3()), boundary(nullptr), img_boundary(nullptr), face(nullptr)
{
}

TriangulatedHalfEdgeDataStructure3::HalfEdge::~HalfEdge()
{
    if (boundary)
    {
        delete boundary; boundary = nullptr;
    }

    if (img_boundary)
    {
        delete img_boundary; img_boundary = nullptr;
    }
}

TriangulatedHalfEdgeDataStructure3::HalfEdge::HalfEdge(const HalfEdge& half_edge)
{
    this->vertex = half_edge.vertex;
    this->normal = half_edge.normal;
    this->opposite = half_edge.opposite;
    this->next = half_edge.next;
    this->tangent = half_edge.tangent;
    this->boundary = half_edge.boundary; // !
    this->img_boundary = half_edge.img_boundary; // !
    this->face = half_edge.face;
}
TriangulatedHalfEdgeDataStructure3::HalfEdge& TriangulatedHalfEdgeDataStructure3::HalfEdge::operator =(const HalfEdge& rhs)
{
    if (this != &rhs)
    {
        this->vertex = rhs.vertex;
        this->normal = rhs.normal;
        this->opposite = rhs.opposite;
        this->next = rhs.next;
        this->tangent = rhs.tangent;
        this->boundary = rhs.boundary; // !
        this->img_boundary = rhs.img_boundary; // !
        this->face = rhs.face;
    }

    return *this;
}

GLboolean TriangulatedHalfEdgeDataStructure3::LoadFromOFF(
        const string &file_name, GLboolean translate_and_scale_to_unit_cube)
{
    fstream f(file_name.c_str(), ios_base::in);

    if (!f || !f.good())
        return GL_FALSE;

    // loading the header
    string header;

    f >> header;

    if (header != "OFF")
        return GL_FALSE;

    // loading number of vertices, faces, and edges
    GLuint vertex_count, face_count, edge_count;

    f >> vertex_count >> face_count >> edge_count;

    // allocating memory for vertices, unit normal vectors and faces
    _vertex.resize(vertex_count);
    _normal.resize(vertex_count);
    _face.resize(face_count);

    // initializing the leftmost and rightmost corners of the bounding box
    _leftmost_vertex.x() = _leftmost_vertex.y() = _leftmost_vertex.z() = numeric_limits<GLdouble>::max();
    _rightmost_vertex.x() = _rightmost_vertex.y() = _rightmost_vertex.z() = -numeric_limits<GLdouble>::max();

    // loading vertices and correcting the leftmost and rightmost corners of the bounding box
    for (vector<DCoordinate3>::iterator vit = _vertex.begin(); vit != _vertex.end(); ++vit)
    {
        f >> *vit;

        if (vit->x() < _leftmost_vertex.x())
            _leftmost_vertex.x() = vit->x();
        if (vit->y() < _leftmost_vertex.y())
            _leftmost_vertex.y() = vit->y();
        if (vit->z() < _leftmost_vertex.z())
            _leftmost_vertex.z() = vit->z();

        if (vit->x() > _rightmost_vertex.x())
            _rightmost_vertex.x() = vit->x();
        if (vit->y() > _rightmost_vertex.y())
            _rightmost_vertex.y() = vit->y();
        if (vit->z() > _rightmost_vertex.z())
            _rightmost_vertex.z() = vit->z();
    }

    // if we do not want to preserve the original positions and coordinates of vertices
    if (translate_and_scale_to_unit_cube)
    {
        GLdouble scale = 1.0 / max(_rightmost_vertex.x() - _leftmost_vertex.x(),
                                   max(_rightmost_vertex.y() - _leftmost_vertex.y(),
                                       _rightmost_vertex.z() - _leftmost_vertex.z()));

        DCoordinate3 middle(_leftmost_vertex);
        middle += _rightmost_vertex;
        middle *= 0.5;
        for (vector<DCoordinate3>::iterator vit = _vertex.begin(); vit != _vertex.end(); ++vit)
        {
            *vit -= middle;
            *vit *= scale;
        }
    }

    // loading faces
    for (vector<Face>::iterator fit = _face.begin(); fit != _face.end(); ++fit)
        f >> *fit;

    // calculating average unit normal vectors associated with vertices
    for (vector<Face>::const_iterator fit = _face.begin(); fit != _face.end(); ++fit)
    {
        DCoordinate3 n = _vertex[(*fit)[1]];
        n -= _vertex[(*fit)[0]];

        DCoordinate3 p = _vertex[(*fit)[2]];
        p -= _vertex[(*fit)[0]];

        n ^= p;

        for (GLint node = 0; node < 3; ++node)
            _normal[(*fit)[node]] += n;
    }

    for (vector<DCoordinate3>::iterator nit = _normal.begin(); nit != _normal.end(); ++nit)
        nit->normalize();

    // creating half-edge data structure
    for (vector<Face>::iterator fit = _face.begin(); fit != _face.end(); ++fit)
    {
        for (GLuint i = 0; i < 3; i++) {

            GLuint start = i;
            GLuint end = (i + 1) % 3;

            pair<GLuint, GLuint> index_pair((*fit)[start], (*fit)[end]);

            _edges[index_pair] = new HalfEdge();

            _edges[index_pair]->vertex = &_vertex[(*fit)[start]];
            _edges[index_pair]->normal = &_normal[(*fit)[start]];
            _edges[index_pair]->face = &(*fit);

        }

        for (GLuint i = 0; i < 3; i++) {

            GLuint start = i;
            GLuint end = (i + 1) % 3;
            GLuint next_start = end;
            GLuint next_end = (end + 1) % 3;

            pair<GLuint, GLuint> index_pair((*fit)[start], (*fit)[end]);
            pair<GLuint, GLuint> next_pair((*fit)[next_start], (*fit)[next_end]);
            pair<GLuint, GLuint> opposite_pair((*fit)[end], (*fit)[start]);

            _edges[index_pair]->next = _edges[next_pair];

            if (_edges.find(opposite_pair) != _edges.end()) {

                _edges[index_pair]->opposite = _edges[opposite_pair];
                _edges[opposite_pair]->opposite = _edges[index_pair];
                _edges[opposite_pair]->is_img_boundary_renderable = false;

            } else {

                _edges[index_pair]->opposite = nullptr;
            }

            DCoordinate3 p_i = _vertex[(*fit)[start]];
            DCoordinate3 p_j = _vertex[(*fit)[end]];
            DCoordinate3 n_i = _normal[(*fit)[start]];

            DCoordinate3 f_i = -n_i;

            DCoordinate3 b_i_j = ((p_j - p_i) ^ f_i); // (((p_j - p_i) ^ (-n_i)).length());
            b_i_j.normalize();

            _edges[index_pair]->tangent = (f_i ^ b_i_j);
            _edges[index_pair]->tangent.normalize();


        }
    }

    f.close();

    return GL_TRUE;
}


bool TriangulatedHalfEdgeDataStructure3::GenerateBoundaryCurves(BoundaryType type, GLdouble beta,
                                                                GLdouble theta_1, GLdouble theta_2, GLdouble theta_3,
                                                                GLuint div_point_count, GLenum usage_flag)
{
    // phi lookup tables
    GLuint rho = 3;

    RowMatrix<GLdouble> theta(rho + 1);

    theta[0] = 0;
    theta[1] = theta_1;
    theta[2] = theta_2;
    theta[3] = theta_3;

    RowMatrix< TriangularMatrix<GLdouble> > phi(rho + 1); // phi[0](k, l), phi[1](k, l), phi[2](k, l), ..., phi[rho](k, l)
    TriangularMatrix<GLdouble> PHI(4);

    for (GLuint r = 0; r <= rho; r++) {
        phi[r].ResizeRows(4);
    }

    RowMatrix<GLdouble> c, ch;
    c.ResizeColumns(5);
    ch.ResizeColumns(5);

    GLdouble sin_beta = sin(beta), sin_beta2 = sin_beta * sin_beta, sin_beta3 = sin_beta * sin_beta2,
            //sin_beta4 = sin_beta2 * sin_beta2,
            sin_2beta = sin(2.0 * beta), sin_3beta = sin(3.0 * beta), sin_4beta = sin(4.0 * beta),
            cos_beta = cos(beta), cos_beta2 = cos_beta * cos_beta, cos_beta3 = cos_beta * cos_beta2,
            cos_beta4 = cos_beta2 * cos_beta2,
            cos_2beta = cos(2.0 * beta);

    GLdouble sin_beta_2 = sin(beta / 2.0), sin_beta_2_2 = sin_beta_2 * sin_beta_2, sin_beta_2_3 = sin_beta_2 * sin_beta_2_2,
            sin_beta_2_4 = sin_beta_2_2 * sin_beta_2_2, sin_beta_2_8 = sin_beta_2_4 * sin_beta_2_4,
            sin_3beta_2 = sin(3.0 * beta / 2.0), sin_5beta_2 = sin(5.0 * beta / 2.0), sin_7beta_2 = sin(7.0 * beta / 2.0);
    GLdouble cos_beta_2 = cos(beta / 2.0), cos_beta_2_2 = cos_beta_2 * cos_beta_2,
            cos_3beta_2 = cos(3.0 * beta / 2.0);
    GLdouble denom = 96.0 * sin_beta_2_8;
    GLdouble denom2 = 6.0 * (1.0 - 4.0 * cos_beta + 6.0 * cos_beta2 - 4.0 * cos_beta3 + cos_beta4);
    GLdouble denom3 = denom2 / 2.0;
    GLdouble d1 = beta - sin_beta, d2 = d1 * (2.0 * sin_beta - beta - beta * cos_beta);

    c[0] = 1.0 / sin_beta_2_4;
    c[1] = 1.0 / sin_beta_2_4 * 4.0 * cos_beta_2;
    c[2] = 1.0 / sin_beta_2_4 * (4.0 * cos_beta_2_2 + 2.0);
    c[3] = c[1];
    c[4] = c[0];


    GLdouble cosh1 = cosh(beta), cosh2 = cosh1 * cosh1, cosh3 = cosh2 * cosh1,
            cosh4 = cosh2 * cosh2;
    GLdouble sinh1 = sinh(beta), sinh2 = sinh(2.0 * beta), sinh3 = sinh(3.0 * beta),
            sinh4 = sinh(4.0 * beta);
    GLdouble sinh_beta = sinh(beta/2), sinh_beta2 = sinh_beta * sinh_beta, sinh_beta4 = sinh_beta2 * sinh_beta2,
            sinh_beta8 = sinh_beta4 * sinh_beta4;
    GLdouble denomh = 1536.0 * sinh_beta8;
    GLdouble denomh2 = 1.0 - 4.0 * cosh1 + 6.0 * cosh2 - 4.0 * cosh3 + cosh4;

    GLdouble a = 2.0 * sin_beta - beta - beta * cos_beta,
            b = beta - sin_beta, b2 = b * b;
    GLdouble denom_a = a * b;

    switch (type)
    {
    case CUBIC_POLYNOMIAL:
        //phi[1](0,0) = phi[1](3,0) = phi[1](3,3) = 0.0;
        phi[1](1,0) = phi[1](3,2) = -9.0 / 10.0;
        phi[1](2,0) = phi[1](3,1) = -3.0 / 5.0;
        phi[1](1,1) = phi[1](2,2) = 6.0 / 5.0;
        phi[1](2,1) = 3.0 / 10.0;

        phi[2](1,0) = phi[2](3,2) = phi[2](2,1) = -18.0;
        phi[2](2,0) = phi[2](3,1) = 0.0;
        phi[2](1,1) = phi[2](2,2) = 36.0;
        //phi[2](3,3) = phi[2](3,0) = phi[2](0,0)

        phi[3](0,0) = phi[3](3,3) =  36.0;
        phi[3](1,0) = phi[3](3,2) = -108.0;
        phi[3](1,1) = phi[3](2,2) = 324.0;
        phi[3](2,0) = phi[3](3,1) = 108.0;
        phi[3](2,1) = -324.0;
        phi[3](3,0) = -36.0;
        break;

    case QUARTIC_POLYNOMIAL:
        //phi[1](0,0) = phi[1](3,0) = phi[1](3,3) = 0.0;
        phi[1](1,0) = phi[1](3,2) = -52.0 / 35.0;
        phi[1](2,0) = phi[1](3,1) = -24.0 / 35.0;
        phi[1](1,1) = phi[1](2,2) = 66.0 / 35.0;
        phi[1](2,1) = 2.0 / 7.0;

        phi[2](1,0) = phi[2](3,2) = -204.0 / 5.0;
        phi[2](2,0) = phi[2](3,1) = 36.0 / 5.0;
        phi[2](1,1) = phi[2](2,2) = 324.0 / 5.0;
        phi[2](2,1) = -156.0 / 5.0;
        //phi[2](3,3) = phi[2](3,0) = phi[2](0,0)

        phi[3](0,0) = phi[3](3,3) = 192.0;
        phi[3](1,0) = phi[3](3,2) = -336.0;
        phi[3](1,1) = phi[3](2,2) = 624.0;
        phi[3](2,0) = phi[3](3,1) = 240.0;
        phi[3](2,1) = - 528.0;
        phi[3](3,0) = - 96.0;
        break;

    case SECOND_ORDER_TRIGONOMETRIC:
        //phi[1](0,0) = phi[1](3,0) = phi[1](3,3) = 0.0;
        phi[1](1,0) = phi[1](3,2) = ((24.0 + 4.0 * cos_beta - 18.0 * cos_beta2 + 5.0 * cos_beta3) * sin_beta -
                3.0 * beta * (6.0 + 2.0 * cos_beta - 3.0 * cos_beta2)) / denom;
        phi[1](2,0) = phi[1](3,1) = ((-12.0 + 24.0 * cos_beta + 2.0 * cos_beta2 + cos_beta3) * sin_beta +
                3.0 * beta * (2.0 - 2.0 * cos_beta - 5.0 * cos_beta2)) / denom;
        phi[1](1,1) = phi[1](2,2) = ((-32.0 - 12.0 * cos_beta + 34.0 * cos_beta2 - 5.0 * cos_beta3) * sin_beta +
                                     3.0 * beta * (8.0 + 4.0 * cos_beta - 5.0 * cos_beta2 - 2.0 * cos_beta3)) / denom2;
        phi[1](2,1) = ((20.0 - 16.0 * cos_beta - 18 * cos_beta2 - cos_beta3) * sin_beta -
                       3.0 * beta * (4.0  - 7.0 * cos_beta2 - 2.0 * cos_beta3)) / denom2;

        phi[2](1,0) = phi[2](3,2) = ((12.0 + 5.0 * cos_beta + 9.0 * cos_beta2 - 14.0 * cos_beta3) * sin_beta -
                                     3.0 * beta * (6.0 + cos_beta - 3.0 * cos_beta2)) / denom3;
        phi[2](2,0) = phi[2](3,1) = ((-18.0 + 21.0 * cos_beta + 7.0 * cos_beta2 + 2.0 * cos_beta3) * sin_beta +
                                     3.0 * beta * (4.0 - cos_beta - 7.0 * cos_beta2)) / denom3;
        phi[2](1,1) = phi[2](2,2) = ((-16.0 -24.0 * cos_beta + 11.0 * cos_beta2 + 17.0 * cos_beta3) * sin_beta +
                                     3.0 * beta * (10.0 - 7.0 * cos_beta2 + 2.0 * cos_beta - cos_beta3)) / denom3;
        phi[2](2,1) = ((22.0 - 2.0 * cos_beta - 27.0 * cos_beta2 - 5.0 * cos_beta3) * sin_beta -
                       3.0 * beta * (8.0 - 11.0 * cos_beta2 - cos_beta3)) / denom3;
        //phi[2](3,3) = phi[2](3,0) = phi[2](0,0)

        /*phi[3](0,0) = phi[3](3,3) = c[0] * c[0] * 1.0 / 24.0 * (15 * beta - 16.0 * sin_beta3 -
                                                                3.0 * sin_4beta - 3.0 * sin_beta * cos_beta);
        phi[3](3,0) = c[0] * c[4] * 1.0 / 24.0 * (3.0 * beta * cos_beta + 12.0 * beta * cos_2beta + sin_beta * (13.0 - 28.0 * cos_beta));
        phi[3](1,0) = phi[3](3,2) = 1.0 / 4.0 * c[4] * (c[1] * 5.0 / 24.0  * (18.0 * sin_beta - 9.0 * sin_2beta + 4.0 * (sin_3beta - 3.0 * beta) * cos_beta_2)+
                c[2] * 1.0 / 48.0 * (6.0 * beta - 30.0 * sin_beta - 3.0 * sin_2beta - 8.0 * sin_3beta + 54.0 * beta * cos_beta));
        phi[3](3,1) = 1.0 / 4.0 * c[4] * (c[1] * 1.0 / 4.0 * (-7.0 * sin_beta_2 +
                                     17.0 / 3.0 * sin_3beta_2 + 2.0 * sin_5beta_2 - beta * cos_beta_2 * (17.0 * cos_beta - 7.0)) +
                c[2] * 1.0 / 48.0 * (6.0 * beta - 30.0 * sin_beta - 3.0 * sin_2beta - 8.0 * sin_3beta + 54.0 * beta * cos_beta));
        phi[3](2,0) = 1.0 / 4.0 * c[0] * (c[3] * 1.0 / 24.0 * (-42.0 * sin_beta_2 + 34.0 * sin_3beta_2 + 12.0 * sin_5beta_2 -
                                                               9.0 * beta * cos_beta_2 - 51.0 * beta * cos_3beta_2) +
                c[2] * 1.0 / 48.0 * (6.0 * beta - 30.0 * sin_beta - 3.0 * sin_2beta - 8.0 * sin_3beta + 54.0 * beta * cos_beta));
        phi[3](1,1) = phi[3](2,2) = c[3] * c[3] * 1.0 / 16.0 * 1.0 / 48.0 * (444.0 * beta - 435 * sin_beta + 78 * sin_2beta - 67 * sin_3beta + 36 * beta * cos_beta) +
                c[2] * c[3] * 1.0 / 8.0 * (-3.0/4.0 * sin_beta_2_3 * (9.0 * cos_beta + 11.0) + 5.0 / 96.0 * (82.0 * sin_beta_2 + 3.0 * sin_3beta_2 + sin_5beta_2
                                                                                                             - 48.0 * beta * cos_beta) +
                                           1.0 / 4.0 * sin_beta_2_2 * (2.0 * (6.0 * sin_beta_2 + sin_3beta_2) - 3* beta * cos_beta_2) -
                                           5.0 / 48.0 * (2.0 * (sin_beta_2 - 1.0 * sin_3beta_2 + sin_5beta_2) + 21.0 * beta * cos_beta_2 +
                                                         3.0 * beta * cos_3beta_2)) +
                c[2] * c[2] * 1.0 / 16.0 * 1.0 / 12.0 * (27.0 * beta - 19.0 * sin_beta + 3.0 * beta * cos_beta - 11.0 * sin_beta * cos_beta);

        phi[3](2,1) = c[2] * c[2] * 1.0 / 16.0 * 1.0 / 12.0 * (27.0 * beta - 19.0 * sin_beta + 3.0 * beta * cos_beta - 11.0 * sin_beta * cos_beta) +
                c[2] * c[3] * 1.0 / 4.0 * (-3.0/4.0 * sin_beta_2_3 * (9.0 * cos_beta + 11.0) + 5.0 / 96.0 * (82.0 * sin_beta_2 + 3.0 * sin_3beta_2 + sin_5beta_2
                                                                                                             - 48.0 * beta * cos_beta) +
                                           1.0 / 4.0 * sin_beta_2_2 * (2.0 * (6.0 * sin_beta_2 + sin_3beta_2) - 3* beta * cos_beta_2) -
                                           5.0 / 48.0 * (2.0 * (sin_beta_2 - 1.0 * sin_3beta_2 + sin_5beta_2) + 21.0 * beta * cos_beta_2 +
                                                         3.0 * beta * cos_3beta_2)) +
                c[1] * c[3] *1.0 / 16.0 * (5.0 / 32.0 * (9.0 * beta - 16.0 * sin_beta - 16.0 * sin_2beta + 36.0 * beta * cos_beta + 3.0 * beta * cos_2beta) +
                                           5.0 / 96.0 * (3.0 * beta * sin_beta2 - 3.0 * beta * (cos_beta - 9.0) * (cos_beta + 1) + 4 * sin_beta * (cos_beta-13.0)) -
                                           9.0 / 4.0 * sin_beta_2_3 * (beta * sin_beta_2 + 4.0 * cos_beta_2) -
                                           3.0 / 8.0 * sin_beta_2_2 * (5.0 * beta - 26.0 * sin_beta + beta * cos_beta));*/
        break;

    case SECOND_ORDER_HYPERBOLIC:
        //phi[1](0,0) = phi[1](3,0) = phi[1](3,3) = 0.0;
        phi[1](1,0) = phi[1](3,2) = (48.0 * beta * (6.0 + 2.0 * cosh1 - 3.0 * cosh2)  - 312.0 * sinh1 -
                                     52.0 * sinh2 + 72.0 * sinh3 - 10.0 * sinh4) / denomh;
        phi[1](2,0) = phi[1](3,1) = (-48.0 * beta * (2.0 - 2.0 * cosh1 - 5.0 * cosh2) + 184.0 * sinh1 -
                                     196.0 * sinh2 - 8.0 * sinh3 - 2.0 * sinh4) / denomh;
        phi[1](1,1) = phi[1](2,2) = (-48.0 * beta * (8.0 + 4.0 * cosh1 - 5.0 * cosh2 - 2.0 * cosh3) +
                                     376.0 * sinh1 + 116.0 * sinh2 - 136.0 * sinh3 + 10.0 * sinh4) / ( 96.0 * denomh2);
        phi[1](2,1) = (48.0 * beta * (4.0 - 7.0 * cosh2 - 2.0 * cosh3) - 248.0 * sinh1 + 132.0 * sinh2 +
                       72.0 * sinh3 + 2.0 * sinh4) / ( 96.0 * denomh2);


        phi[2](1,0) = phi[2](3,2) = (-24.0 * beta * (6.0 + cosh1 - 3.0 * cosh2) + 114.0 * sinh1 -
                                     8.0 * sinh2 + 18.0 * sinh3 - 14.0 * sinh4) / ( 24.0 * denomh2);
        phi[2](2,0) = phi[2](3,1) = (24.0 * beta * (4.0 - cosh1 - 7.0 * cosh2) -130.0 * sinh1 + 88.0 * sinh2 +
                                     14.0 * sinh3 + 2.0 * sinh4) / ( 24.0 * denomh2);
        phi[2](1,1) = phi[2](2,2) = (48.0 * beta * (10.0 + 2.0 * cosh1 - 7.0 * cosh2 - cosh3) -212.0 * sinh1 -
                                     124.0 * sinh2 + 44.0 * sinh3 + 34.0 * sinh4) / ( 48.0 * denomh2);
        phi[2](2,1) = (-48.0 * beta * (8.0 - 11.0 * cosh2 - cosh3) + 244.0 * sinh1 - 36.0 * sinh2 - 108.0 * sinh3 -
                       10.0 * sinh4) / ( 48.0 * denomh2);
        //phi[2](3,3) = phi[2](3,0) = phi[2](0,0)

        break;

    case FIRST_ORDER_ALGEBRAIC_TRIGONOMETRIC:
        //phi[1](0,0) = phi[1](3,0) = phi[1](3,3) = 0.0;
        GLdouble beta2 = beta * beta, beta3 = beta * beta2;
        phi[1](1,0) = phi[1](3,2) = (sin_beta * (-3.0 * beta + 6.0 * sin_beta + 2.0 * beta * cos_beta -
                                             3.0 * sin_2beta + beta * cos_2beta)) / (4.0 * a * b2);
        phi[1](2,0) = phi[1](3,1) = (sin_beta * (- beta - sin_beta + beta2 * sin_beta +
                                             beta * cos_beta + cos_beta * sin_beta)) / (2.0 * a * b2);
        phi[1](1,1) = phi[1](2,2) = (sin_beta2 * (2.0 * beta3 - 4.0 * sin_beta - 4.0 * beta2 * sin_beta +
                                             4.0 * beta * cos_beta + 2.0 * sin_2beta - beta2 * sin_2beta -
                                             4.0 * beta * cos_2beta)) / (4.0 * denom_a * denom_a);
        phi[1](2,1) = (sin_beta2 * (6.0 * beta - 2.0 * sin_beta - 3.0 * beta2 * sin_beta -
                               6.0 * beta * cos_beta + beta3 * cos_beta + sin_2beta)) /(2.0 * denom_a * denom_a);

        phi[2](1,0) = phi[2](3,2) = (sin_beta * (- beta - 2.0 * sin_beta + 2.0 * beta * cos_beta +
                                             sin_2beta - beta * cos_2beta)) / (4.0 * a * b2);
        phi[2](2,0) = phi[2](3,1) = (sin_beta * (- beta - sin_beta + beta2 * sin_beta + beta * cos_beta +
                                             cos_beta * sin_beta)) / (2.0 * a * b2);
        phi[2](1,1) = phi[2](2,2) = (sin_beta2 * (2.0 * beta + 2.0 * beta3 + 4.0 * sin_beta -
                                             4.0 * beta2 * sin_beta - 4.0 * beta * cos_beta - 2.0 * sin_2beta +
                                             beta2 * sin_2beta + 2.0 * beta * cos_2beta)) / (4.0 * denom_a * denom_a);
        phi[2](2,1) = (sin_beta2 * (beta + 2.0 * sin_beta - beta2 * sin_beta - 2.0 * beta * cos_beta +
                               beta3 * cos_beta - sin_2beta + beta * cos_2beta)) / (2.0 * denom_a * denom_a);
        //phi[2](3,3) = phi[2](3,0) = phi[2](0,0)


        phi[3](0,0) = phi[3](3,3) = 1.0 / (2.0 * d1 * d1) * (beta + sin_beta * cos_beta);
        phi[3](1,0) = phi[3](3,2) = sin_beta / (d1 * d2) * sin_beta_2_2 * (sin_beta - beta * (cos_beta + 2.0));
        phi[3](1,1) = phi[3](2,2) = sin_beta2 / (4.0 * d2 * d2) * (2.0 * (beta3 + 3.0 * beta - 2.0 * sin_beta + sin_2beta) -
                                                     beta * (beta * (4.0 * sin_beta + sin_2beta) + 4.0 * cos_beta + 2.0 * cos_2beta));
        phi[3](2,0) = phi[3](3,1) = sin_beta / (2.0 * d1 * d2) * ((beta2 + 1) * sin_beta - beta + (beta - sin_beta) * cos_beta);
        phi[3](2,1) = sin_beta2 / (d2 * d2) * (-1.0 / 2.0 * (3.0 * beta2 + 2.0) * sin_beta +
                                               (1.0 / 2.0 * beta * (beta2 - 2) + sin_beta) * cos_beta + beta + beta * sin_beta2);
        phi[3](3,0) = -1.0 / (2.0 * d1 * d1) * (sin_beta + beta * cos_beta);

        break;
    }

    for (GLuint r = 1 ; r <= rho; r++) {
        for (GLuint k = 0; k <= 3; k++) {
            for (GLuint l = 0; l <= k; l++)
                PHI(k, l) += theta[r] * phi[r](k, l);
        }
    }

    //double time = omp_get_wtime();

    GLboolean is_parallel = GL_TRUE;

    if (is_parallel){
        GLint num_of_faces = _face.size();
        //omp_set_num_threads(1);
        #pragma omp parallel for
        for (GLint i_j = 0; i_j< num_of_faces * 3; i_j++)
        {
            GLint i = i_j / 3;
            GLint j = i_j % 3;

            pair<GLint, GLint> index_pair(_face[i][j % 3], _face[i][(j+1)%3]);

            HalfEdge* it = _edges[index_pair];

            if (it->vertex && it->normal && it->next->vertex && it->next->normal)
            {
                DCoordinate3 p_i = *it->vertex;
                DCoordinate3 p_j = *it->next->vertex;

                DCoordinate3 n_i = *it->normal;

                DCoordinate3 f_i = -n_i;

                DCoordinate3 b_i_j = ((p_j - p_i) ^ f_i); // (((p_j - p_i) ^ (-n_i)).length());
                b_i_j.normalize();

                it->tangent = (f_i ^ b_i_j);
                it->tangent.normalize();

                DCoordinate3 n_j = *it->next->normal;
                DCoordinate3 f_j = -n_j;

                DCoordinate3 b_j_i = ((p_i - p_j) ^ f_j);
                b_j_i.normalize();

                it->next->tangent = (f_j ^ b_j_i);
                it->next->tangent.normalize();


                DCoordinate3 t_i_j = it->tangent;
                DCoordinate3 t_j_i = it->next->tangent;

                switch (type)
                {
                case CUBIC_POLYNOMIAL:
                    it->boundary = new (nothrow) CubicPolinomialCurve();
                    break;

                case QUARTIC_POLYNOMIAL:
                    it->boundary = new (nothrow) QuarticPolinomialCurve();
                    break;

                case SECOND_ORDER_TRIGONOMETRIC:
                    it->boundary = new (nothrow) SecondOrderTrigonometricCurve(beta);
                    break;

                case SECOND_ORDER_HYPERBOLIC:
                    it->boundary = new (nothrow) SecondOrderHyperbolicCurve(beta);
                    break;

                case FIRST_ORDER_ALGEBRAIC_TRIGONOMETRIC:
                    it->boundary = new (nothrow) FirstOrderAlgebraicTrigonometricCurve(beta);
                    break;
                }


                GLdouble dot_t_i_j_t_j_i = t_i_j * t_j_i;

                GLdouble delta = PHI(1, 1) * PHI(2, 2) - dot_t_i_j_t_j_i * dot_t_i_j_t_j_i * PHI(2, 1) * PHI(2, 1);

                RealSquareMatrix inverse(2);

                inverse(0, 0) = PHI(2, 2);
                inverse(0, 1) = -dot_t_i_j_t_j_i * PHI(2, 1);
                inverse(1, 0) = inverse(0, 1);
                inverse(1, 1) = PHI(1, 1);

                for (GLuint m = 0; m <= 1; m++)
                {
                    for (GLuint n = 0; n <= 1; n++)
                    {
                        inverse(m, n) /= -delta;
                    }
                }

                ColumnMatrix<GLdouble> B(2);

                B[0] = (p_i * (PHI(1, 0) + PHI(1, 1)) + p_j * (PHI(2, 1) + PHI(3, 1))) * t_i_j;
                B[1] = (p_i * (PHI(2, 0) + PHI(2, 1)) + p_j * (PHI(2, 2) + PHI(3, 2))) * t_j_i;

                LinearCombination3 *ptr_boundary = it->boundary;

                if (!ptr_boundary)
                {
                    cout<<"no boundary"<<endl;
                }

                (*ptr_boundary)[0] = p_i;
                (*ptr_boundary)[3] = p_j;

                // (24) --> lambda[0], lambda[1]
                ColumnMatrix<GLdouble> lambda(2);

                for (GLuint n = 0; n <= 1; n++)
                {
                    for (GLuint m = 0; m <= 1; m++)
                    {
                        lambda[n] += inverse(n, m) * B[m];
                    }
                }

                // solving the optimization problem
                (*ptr_boundary)[1] = p_i + lambda[0] * t_i_j;
                (*ptr_boundary)[2] = p_j + lambda[1] * t_j_i;

                /*if(!ptr_boundary->UpdateVertexBufferObjectsOfData())
                {
                    cout<< "Could not update the VBOs of the boundary"<<endl;
                }

                it->img_boundary = ptr_boundary->GenerateImage(2, div_point_count, usage_flag);

                if (!it->img_boundary)
                {
                    cout << "could not create boundary curve!" << endl;
                }

                if (!it->img_boundary->UpdateVertexBufferObjects(1, usage_flag))
                {
                    cout << "could not update the VBOs of the boundary's image!" << endl;
                }*/

            }
        }

        for(map<pair<GLint, GLint>, HalfEdge* >::iterator it = _edges.begin(); it != _edges.end(); it++)
        {
            if (it->second->vertex && it->second->normal && it->second->next->vertex && it->second->next->normal)
            {
                LinearCombination3* ptr_boundary = it->second->boundary;
                if(!ptr_boundary->UpdateVertexBufferObjectsOfData())
                {
                    cout<< "Could not update the VBOs of the boundary"<<endl;
                }

                it->second->img_boundary = ptr_boundary->GenerateImage(2, div_point_count, usage_flag);

                if (!it->second->img_boundary)
                {
                    cout << "could not create boundary curve!" << endl;
                }

                if (!it->second->img_boundary->UpdateVertexBufferObjects(1, usage_flag))
                {
                    cout << "could not update the VBOs of the boundary's image!" << endl;
                }

            }
        }

    }
    else
    {
        for ( map< pair<GLint, GLint>, HalfEdge* >::iterator it = _edges.begin(); it != _edges.end(); it++)
        {
            if (it->second->vertex && it->second->normal && it->second->next->vertex && it->second->next->normal)
            {
                DCoordinate3 p_i = *it->second->vertex;
                DCoordinate3 p_j = *it->second->next->vertex;

                DCoordinate3 n_i = *it->second->normal;

                DCoordinate3 f_i = -n_i;

                DCoordinate3 b_i_j = ((p_j - p_i) ^ f_i); // (((p_j - p_i) ^ (-n_i)).length());
                b_i_j.normalize();

                it->second->tangent = (f_i ^ b_i_j);
                it->second->tangent.normalize();

                DCoordinate3 n_j = *it->second->next->normal;
                DCoordinate3 f_j = -n_j;

                DCoordinate3 b_j_i = ((p_i - p_j) ^ f_j);
                b_j_i.normalize();

                it->second->next->tangent = (f_j ^ b_j_i);
                it->second->next->tangent.normalize();


                DCoordinate3 t_i_j = it->second->tangent;
                DCoordinate3 t_j_i = it->second->next->tangent;

                switch (type)
                {
                case CUBIC_POLYNOMIAL:
                    it->second->boundary = new (nothrow) CubicPolinomialCurve();
                    break;

                case QUARTIC_POLYNOMIAL:
                    it->second->boundary = new (nothrow) QuarticPolinomialCurve();
                    break;

                case SECOND_ORDER_TRIGONOMETRIC:
                    it->second->boundary = new (nothrow) SecondOrderTrigonometricCurve(beta);
                    break;

                case SECOND_ORDER_HYPERBOLIC:
                    it->second->boundary = new (nothrow) SecondOrderHyperbolicCurve(beta);
                    break;

                case FIRST_ORDER_ALGEBRAIC_TRIGONOMETRIC:
                    it->second->boundary = new (nothrow) FirstOrderAlgebraicTrigonometricCurve(beta);
                    break;
                }


                GLdouble dot_t_i_j_t_j_i = t_i_j * t_j_i;

                GLdouble delta = PHI(1, 1) * PHI(2, 2) - dot_t_i_j_t_j_i * dot_t_i_j_t_j_i * PHI(2, 1) * PHI(2, 1);

                RealSquareMatrix inverse(2);

                inverse(0, 0) = PHI(2, 2);
                inverse(0, 1) = -dot_t_i_j_t_j_i * PHI(2, 1);
                inverse(1, 0) = inverse(0, 1);
                inverse(1, 1) = PHI(1, 1);

                for (GLuint m = 0; m <= 1; m++)
                {
                    for (GLuint n = 0; n <= 1; n++)
                    {
                        inverse(m, n) /= -delta;
                    }
                }

                ColumnMatrix<GLdouble> B(2);

                B[0] = (p_i * (PHI(1, 0) + PHI(1, 1)) + p_j * (PHI(2, 1) + PHI(3, 1))) * t_i_j;
                B[1] = (p_i * (PHI(2, 0) + PHI(2, 1)) + p_j * (PHI(2, 2) + PHI(3, 2))) * t_j_i;

                LinearCombination3 *ptr_boundary = it->second->boundary;

                if (!ptr_boundary)
                {
                    cout<<"no boundary"<<endl;
                }

                (*ptr_boundary)[0] = p_i;
                (*ptr_boundary)[3] = p_j;

                // (24) --> lambda[0], lambda[1]
                ColumnMatrix<GLdouble> lambda(2);

                for (GLuint n = 0; n <= 1; n++)
                {
                    for (GLuint m = 0; m <= 1; m++)
                    {
                        lambda[n] += inverse(n, m) * B[m];
                    }
                }

                // solving the optimization problem

                (*ptr_boundary)[1] = p_i + lambda[0] * t_i_j;
                (*ptr_boundary)[2] = p_j + lambda[1] * t_j_i;


                ptr_boundary->UpdateVertexBufferObjectsOfData();

                it->second->img_boundary = ptr_boundary->GenerateImage(2, div_point_count, usage_flag);

                if (!it->second->img_boundary)
                {
    //                cout << "could not create boundary curve!" << endl;
                }

                if (!it->second->img_boundary->UpdateVertexBufferObjects(1, usage_flag))
                {
    //                cout << "could not update the VBOs of the boundary's image!" << endl;
                }

            }
        }
    }

    //cout<< omp_get_wtime() - time <<endl;
    return true;
}

GLdouble TriangulatedHalfEdgeDataStructure3::fact(GLuint n)
{
    GLdouble res = 1.0;
    if (n > 0){
        for (GLuint i = 1; i <= n; i++)
            res = res * i;
    }
    return res;
}

GLdouble TriangulatedHalfEdgeDataStructure3::combination(GLuint n, GLuint k)
{
    return(fact(n) / (fact(k) * fact(n - k)));
}

void TriangulatedHalfEdgeDataStructure3::GenerateC0TriangularSurfaces(BoundaryType type, GLdouble beta,
                                                                      GLdouble epsilon_1, GLdouble epsilon_2, GLdouble epsilon_3,
                                                                      GLuint div_point_count, GLenum usage_flag)
{
    RowMatrix<GLdouble> epsilon(4);
    epsilon[0] = 0;
    epsilon[1] = epsilon_1;
    epsilon[2] = epsilon_2;
    epsilon[3] = epsilon_3;

    Matrix<TriangularMatrix<GLdouble>> tau(4, 4);

    for (GLuint i = 0; i <= 3; i++)
    {
        for (GLuint j = 0; j <= 3; j++)
        {
            //tau(i,j).ResizeColumns(4);
            tau(i,j).ResizeRows(4);
        }
    }

    GLdouble beta2 = beta * beta, beta_2 = beta / 2.0, beta3 = beta * beta2, beta4 = beta2 * beta2,
            sin_beta = sin(beta), cos_beta = cos(beta), sin_beta_2 = sin(beta_2),
            sin2_beta_2 = sin_beta_2 * sin_beta_2, sin4_beta_2 = sin2_beta_2 * sin2_beta_2, sin8_beta_2 = sin4_beta_2 * sin4_beta_2,
            cos_beta_2 = cos(beta_2), cos2_beta_2 = cos_beta_2 * cos_beta_2, cos3_beta_2 = cos2_beta_2 * cos_beta_2,
            cos4_beta_2 = cos2_beta_2 * cos2_beta_2, cos6_beta_2 = cos4_beta_2 * cos2_beta_2, cos8_beta_2 = cos4_beta_2 * cos4_beta_2,
            b_sin_cos_beta_2 = beta * sin_beta_2 * cos_beta_2,
            d1 = beta - sin_beta, d2 = d1 * (2.0 * sin_beta - beta - beta * cos_beta),
            c1 = 1.0 / d1, c2 = sin_beta / d2, c3 = 4.0 * (3.0 * beta + 4.0 * sin_beta - beta * cos_beta) * cos_beta_2 / d2,
            c3_2 = c3 * c3, c4 = 4 * sin_beta * cos_beta_2 / d2, c4_2 = c4 * c4;

    switch(type)
    {
    case CUBIC_POLYNOMIAL:
        //Variations of products of first order partial derivatives
        tau(0,1)(3,0) = tau(1,0)(3,0) = tau(1,0)(3,1) = tau(0,1)(3,3) = tau(0,1)(2,0) =
                tau(0,1)(2,2) = tau(1,0)(1,1) = tau(1,0)(0,0) = - 0.1;
        tau(0,1)(3,1) = tau(0,1)(3,2) = tau(1,0)(2,0) = tau(1,0)(1,0) = 0.1;
        tau(0,1)(2,1) = tau(1,0)(2,1) = 1.0 / 5.0;

        //Variations of products of second order partial derivatives
        tau(0,2)(3,0) = tau(2,0)(3,0) = tau(1,1)(3,1) = tau(0,2)(3,3) = tau(1,1)(2,0) =
                tau(1,1)(2,2) = tau(1,1)(1,1) = tau(2,0)(0,0) = -3.0;
        tau(1,1)(3,0) = tau(2,0)(3,2) = tau(1,1)(3,3) = tau(2,0)(3,3) = tau(2,0)(2,2) =
                tau(0,2)(1,0) = tau(0,2)(1,1) = tau(0,2)(0,0) =  tau(1,1)(0,0) = 0.0;
        tau(0,2)(3,1) = tau(0,2)(3,2) = tau(2,0)(2,0) = tau(1,1)(3,2) = tau(1,1)(1,0) = tau(2,0)(1,0) = 3.0;
        tau(2,0)(3,1) = tau(0,2)(2,0) = tau(0,2)(2,2) = tau(2,0)(1,1) = -6.0;
        tau(0,2)(2,1) = tau(2,0)(2,1) = 12.0;
        tau(1,1)(2,1) = 6;

        break;
    case QUARTIC_POLYNOMIAL:
        tau(0,1)(3,0) = tau(1,0)(3,0) = tau(0,1)(3,3)= tau(1,0)(0,0) = - 6.0 / 35.0;
        tau(0,1)(3,1) = tau(0,1)(3,2) = tau(1,0)(2,0) = tau(1,0)(1,0) = 6.0 /35.0;
        tau(1,0)(3,1) = tau(0,1)(2,0) = tau(0,1)(2,2) = tau(1,0)(1,1) = - 11.0 / 35.0;
        tau(1,0)(3,2) = tau(1,0)(2,2) = tau(0,1)(1,0) = tau(0,1)(1,1) = - 3.0 / 35.0;
        tau(1,0)(3,3) = tau(0,1)(0,0) = 0.0;
        tau(0,1)(2,1) = tau(1,0)(2,1) = 4.0 / 5.0;

        tau(2,0)(3,0) = tau(0,2)(3,0) = tau(0,2)(3,3) = tau(2,0)(0,0) = - 24.0 / 5.0;
        tau(1,1)(3,0) = 12.0 / 5.0;
        tau(0,2)(3,1) = tau(0,2)(3,2) = tau(1,1)(3,2) = tau(1,1)(1,0) = tau(2,0)(2,0) = tau(2,0)(1,0) = 24.0 / 5.0;
        tau(1,1)(3,1) = tau(1,1)(2,0) = tau(2,0)(3,2) = tau(2,0)(2,2) = tau(0,2)(1,0) = tau(0,2)(1,1) = - 36.0 / 5.0;
        tau(2,0)(3,1) = tau(0,2)(2,0) = tau(0,2)(2,2) = tau(2,0)(1,1) = - 84.0 / 5.0;
        tau(1,1)(3,3) = tau(1,1)(0,0) = tau(2,0)(3,3) = tau(0,2)(0,0) = 0.0;
        tau(0,2)(2,1) = tau(1,1)(2,1) = tau(2,0)(2,1) = 48.0;
        tau(1,1)(2,2) = tau(1,1)(1,1) = - 54.0 / 5.0;
        break;
    case SECOND_ORDER_TRIGONOMETRIC:
        tau(0,1)(3,0) = tau(1,0)(3,0) = tau(0,1)(3,3) = tau(1,0)(0,0) =
                (1.0 / (1152 * sin8_beta_2)) * (384 - 81 * beta2 - 2330 * cos2_beta_2 + 288 * beta2 * cos2_beta_2 + 1770 * cos4_beta_2 -
                                                180 * beta2 * cos4_beta_2 + 312 * cos6_beta_2 - 136 * cos8_beta_2 -
                                                3 * b_sin_cos_beta_2 * (37 - 286 * cos2_beta_2));
        tau(1,0)(3,1) = tau(0,1)(2,0) = tau(0,1)(2,2) = tau(1,0)(1,1) =
                (1.0 / (2304 * sin8_beta_2)) * (-27 * beta2 + 2 * (890 + 18 * beta2) * cos2_beta_2 - 180 * (1 + beta2) * cos4_beta_2 +
                                                48 * (3 * beta2 - 31) * cos6_beta_2 - 112 * cos8_beta_2 - 12 * b_sin_cos_beta_2 *
                                                (37 - 6 * cos2_beta_2 + 130 * cos4_beta_2 - 20 * cos6_beta_2));
        tau(0,1)(3,1) = tau(0,1)(3,2) = tau(1,0)(2,0) = tau(1,0)(1,0) =
                (1.0 / (2304 * sin8_beta_2)) * (640 - 45 * beta2 - 12 * (143 + 3 * beta2) * cos2_beta_2 + 36 * (93 - 7 * beta2) *
                                                cos4_beta_2 - 16 * (169 - 9 * beta2) * cos6_beta_2 + 432 * cos8_beta_2 + 12 *
                                                b_sin_cos_beta_2 *(1 + 14 * cos2_beta_2 - 38 * cos4_beta_2 - 4 * cos6_beta_2));
        tau(1,0)(3,2) = tau(1,0)(2,2) = tau(0,1)(1,0) = tau(0,1)(1,1) =
                (1.0 / (2304 * sin8_beta_2)) * (-27 * beta2 + 4 * (113 - 27 * beta2) * cos2_beta_2 + 12 * (13 - 9 * beta2) *
                                                cos4_beta_2 - 1200 * cos6_beta_2 + 592 * cos8_beta_2 + 12 * b_sin_cos_beta_2 *
                                                (9 + 6 * cos2_beta_2 - 14 * cos4_beta_2 + 20 * cos6_beta_2));
        tau(0,1)(2,1) = tau(1,0)(2,1) =
                (1.0 / (1152 * sin8_beta_2)) * (-1408 + 261 * beta2 + 4 * (1036 - 117 * beta2) * cos2_beta_2 - 12 * (572 - 75 * beta2) *
                                                cos4_beta_2 + 32 * (149 - 9 * beta2) * cos6_beta_2 - 640 * cos8_beta_2 + 6 *
                                                b_sin_cos_beta_2 * (91 - 338 * cos2_beta_2 + 364 * cos4_beta_2 - 72 * cos6_beta_2));
        tau(0,1)(0,0) = tau(1,0)(3,3) = 0;

        tau(0,2)(3,0) = tau(2,0)(3,0) = tau(0,2)(3,3) = tau(2,0)(0,0) =
                (1.0 / (576 * sin8_beta_2)) * (-336 - 27 * beta2 - 4 * (319 - 63 * beta2) * cos2_beta_2 + 12 * (181 - 21 * beta2) *
                                               cos4_beta_2 - 864 * cos6_beta_2 + 304 * cos8_beta_2 + 12 * b_sin_cos_beta_2 *
                                               (1 + 74 * cos2_beta_2));
        tau(1,1)(3,0) = (1.0 / (576 * sin8_beta_2)) *
                (-144 - 27 * beta2 - 4 * (451 - 63 * beta2) * cos2_beta_2 + 84 * (23 - 3 * beta2) * cos4_beta_2 + 288 * cos6_beta_2 -
                 272 * cos8_beta_2 + 12 * b_sin_cos_beta_2 * (19 + 62 * cos2_beta_2));
        tau(0,2)(3,1) = tau(0,2)(3,2) = tau(2,0)(2,0) = tau(2,0)(1,0) =
                (1.0 / (1152 * sin8_beta_2)) * (224 - 9 * beta2 + 6 * (137 - 39 * beta2) * cos2_beta_2 + 6 * (43 + 6 * beta2) *
                                                cos4_beta_2 + 8 * (80 + 9 * beta2) * cos6_beta_2 - 1944 * cos8_beta_2 + 3 *
                                                b_sin_cos_beta_2 * (131 - 242 * cos2_beta_2 - 472 * cos4_beta_2 - 80 * cos6_beta_2));
        tau(1,1)(3,1) = tau(1,1)(2,0) =
                (1.0 / (2304 * sin8_beta_2)) * (-224 - 45 * beta2 + 8 * (541 - 45 * beta2) * cos2_beta_2 - 36 * (50 - 11 * beta2) *
                                                cos4_beta_2 - 16 * (202 - 9 * beta2) * cos6_beta_2 + 928 * cos8_beta_2 - 6 *
                                                b_sin_cos_beta_2 * (193 + 14 * cos2_beta_2 + 236 * cos4_beta_2 + 40 * cos6_beta_2));
        tau(2,0)(3,1) = tau(0,2)(2,0) = tau(0,2)(2,2) = tau(2,0)(1,1) =
                (1.0 / (1152 * sin8_beta_2)) * (-384 - 54 * beta2 + 2 * (1357 - 99 * beta2) * cos2_beta_2 - 6 * (115 - 48 * beta2) *
                                                cos4_beta_2 - 24 * (88 - 3 * beta2) * cos6_beta_2 + 472 * cos8_beta_2 - 3 *
                                                b_sin_cos_beta_2 * (341 - 150 * cos2_beta_2 + 176 * cos4_beta_2 + 224 * cos6_beta_2));
        tau(1,1)(3,2) = tau(1,1)(1,0) =
                (1.0 / (2304 * sin8_beta_2)) * (-224 - 45 * beta2 + 12 * (131 - 15 * beta2) * cos2_beta_2 - 12 * (79 + 15 * beta2) *
                                                cos4_beta_2 + 2432 * cos6_beta_2 - 2832 * cos8_beta_2 + 24 * b_sin_cos_beta_2 *
                                                (49 - 53 * cos2_beta_2 - 43 * cos4_beta_2 - 10 * cos6_beta_2));
        tau(2,0)(3,2) = tau(2,0)(2,2) = tau(0,2)(1,0) = tau(0,2)(1,1) =
                (1.0 / (576 * sin8_beta_2)) * (-192 -27 * beta2 + 4 * (143 - 27 * beta2) * cos2_beta_2 + 12 * (115 - 9 * beta2) *
                                               cos4_beta_2 - 1296 * cos6_beta_2 - 464 * cos8_beta_2 + 12 * b_sin_cos_beta_2 *
                                               (6 - 38 * cos4_beta_2 - 28 * cos6_beta_2));
        tau(0,2)(2,1) = tau(1,1)(2,1) = tau(2,0)(2,1) =
                (1.0 / (576 * sin8_beta_2)) * (1216 + 171 * beta2 - 16 * (133 - 9 * beta2) * cos2_beta_2 - 12 * (556 - 33 * beta2) *
                                               cos4_beta_2 + 16 * (362 - 9 * beta2) * cos6_beta_2 + 1792 * cos8_beta_2 + 6 *
                                               b_sin_cos_beta_2 * (77 - 250 * cos2_beta_2 + 476 * cos4_beta_2 + 264 * cos6_beta_2));
        tau(1,1)(2,2) = tau(1,1)(1,1) =
                (1.0 / (2304 * sin8_beta_2)) * (-480 - 27 * beta2 - 4 * (41 + 27 * beta2) * cos2_beta_2 + 12 * (463 - 9 * beta2) *
                                                cos4_beta_2 - 5568 * cos6_beta_2 + 656 * cos8_beta_2 - 24 * b_sin_cos_beta_2 *
                                                (39 - 57 * cos2_beta_2 + 17 * cos4_beta_2 + 46 * cos6_beta_2));
        tau(1,1)(3,3) = tau(2,0)(3,3) = tau(0,2)(0,0) = tau(1,1)(0,0) = 0;

        break;
    case SECOND_ORDER_HYPERBOLIC:
        tau(0,1)(3,0) = tau(1,0)(3,0) = tau(0,1)(3,3) = tau(1,0)(0,0) =
                (1.0 / (1152 * sin8_beta_2)) * (384 - 81 * beta2 - 2330 * cos2_beta_2 + 288 * beta2 * cos2_beta_2 + 1770 * cos4_beta_2 -
                                                180 * beta2 * cos4_beta_2 + 312 * cos6_beta_2 - 136 * cos8_beta_2 -
                                                3 * b_sin_cos_beta_2 * (37 - 286 * cos2_beta_2));
        tau(1,0)(3,1) = tau(0,1)(2,0) = tau(0,1)(2,2) = tau(1,0)(1,1) =
                (1.0 / (2304 * sin8_beta_2)) * (-27 * beta2 + 2 * (890 + 18 * beta2) * cos2_beta_2 - 180 * (1 + beta2) * cos4_beta_2 +
                                                48 * (3 * beta2 - 31) * cos6_beta_2 - 112 * cos8_beta_2 - 12 * b_sin_cos_beta_2 *
                                                (37 - 6 * cos2_beta_2 + 130 * cos4_beta_2 - 20 * cos6_beta_2));
        tau(0,1)(3,1) = tau(0,1)(3,2) = tau(1,0)(2,0) = tau(1,0)(1,0) =
                (1.0 / (2304 * sin8_beta_2)) * (640 - 45 * beta2 - 12 * (143 + 3 * beta2) * cos2_beta_2 + 36 * (93 - 7 * beta2) *
                                                cos4_beta_2 - 16 * (169 - 9 * beta2) * cos6_beta_2 + 432 * cos8_beta_2 + 12 *
                                                b_sin_cos_beta_2 *(1 + 14 * cos2_beta_2 - 38 * cos4_beta_2 - 4 * cos6_beta_2));
        tau(1,0)(3,2) = tau(1,0)(2,2) = tau(0,1)(1,0) = tau(0,1)(1,1) =
                (1.0 / (2304 * sin8_beta_2)) * (-27 * beta2 + 4 * (113 - 27 * beta2) * cos2_beta_2 + 12 * (13 - 9 * beta2) *
                                                cos4_beta_2 - 1200 * cos6_beta_2 + 592 * cos8_beta_2 + 12 * b_sin_cos_beta_2 *
                                                (9 + 6 * cos2_beta_2 - 14 * cos4_beta_2 + 20 * cos6_beta_2));
        tau(0,1)(2,1) = tau(1,0)(2,1) =
                (1.0 / (1152 * sin8_beta_2)) * (-1408 + 261 * beta2 + 4 * (1036 - 117 * beta2) * cos2_beta_2 - 12 * (572 - 75 * beta2) *
                                                cos4_beta_2 + 32 * (149 - 9 * beta2) * cos6_beta_2 - 640 * cos8_beta_2 + 6 *
                                                b_sin_cos_beta_2 * (91 - 338 * cos2_beta_2 + 364 * cos4_beta_2 - 72 * cos6_beta_2));
        tau(0,1)(0,0) = tau(1,0)(3,3) = 0;

        tau(0,2)(3,0) = tau(2,0)(3,0) = tau(0,2)(3,3) = tau(2,0)(0,0) =
                (1.0 / (576 * sin8_beta_2)) * (-336 - 27 * beta2 - 4 * (319 - 63 * beta2) * cos2_beta_2 + 12 * (181 - 21 * beta2) *
                                               cos4_beta_2 - 864 * cos6_beta_2 + 304 * cos8_beta_2 + 12 * b_sin_cos_beta_2 *
                                               (1 + 74 * cos2_beta_2));
        tau(1,1)(3,0) = (1.0 / (576 * sin8_beta_2)) *
                (-144 - 27 * beta2 - 4 * (451 - 63 * beta2) * cos2_beta_2 + 84 * (23 - 3 * beta2) * cos4_beta_2 + 288 * cos6_beta_2 -
                 272 * cos8_beta_2 + 12 * b_sin_cos_beta_2 * (19 + 62 * cos2_beta_2));
        tau(0,2)(3,1) = tau(0,2)(3,2) = tau(2,0)(2,0) = tau(2,0)(1,0) =
                (1.0 / (1152 * sin8_beta_2)) * (224 - 9 * beta2 + 6 * (137 - 39 * beta2) * cos2_beta_2 + 6 * (43 + 6 * beta2) *
                                                cos4_beta_2 + 8 * (80 + 9 * beta2) * cos6_beta_2 - 1944 * cos8_beta_2 + 3 *
                                                b_sin_cos_beta_2 * (131 - 242 * cos2_beta_2 - 472 * cos4_beta_2 - 80 * cos6_beta_2));
        tau(1,1)(3,1) = tau(1,1)(2,0) =
                (1.0 / (2304 * sin8_beta_2)) * (-224 - 45 * beta2 + 8 * (541 - 45 * beta2) * cos2_beta_2 - 36 * (50 - 11 * beta2) *
                                                cos4_beta_2 - 16 * (202 - 9 * beta2) * cos6_beta_2 + 928 * cos8_beta_2 - 6 *
                                                b_sin_cos_beta_2 * (193 + 14 * cos2_beta_2 + 236 * cos4_beta_2 + 40 * cos6_beta_2));
        tau(2,0)(3,1) = tau(0,2)(2,0) = tau(0,2)(2,2) = tau(2,0)(1,1) =
                (1.0 / (1152 * sin8_beta_2)) * (-384 - 54 * beta2 + 2 * (1357 - 99 * beta2) * cos2_beta_2 - 6 * (115 - 48 * beta2) *
                                                cos4_beta_2 - 24 * (88 - 3 * beta2) * cos6_beta_2 + 472 * cos8_beta_2 - 3 *
                                                b_sin_cos_beta_2 * (341 - 150 * cos2_beta_2 + 176 * cos4_beta_2 + 224 * cos6_beta_2));
        tau(1,1)(3,2) = tau(1,1)(1,0) =
                (1.0 / (2304 * sin8_beta_2)) * (-224 - 45 * beta2 + 12 * (131 - 15 * beta2) * cos2_beta_2 - 12 * (79 + 15 * beta2) *
                                                cos4_beta_2 + 2432 * cos6_beta_2 - 2832 * cos8_beta_2 + 24 * b_sin_cos_beta_2 *
                                                (49 - 53 * cos2_beta_2 - 43 * cos4_beta_2 - 10 * cos6_beta_2));
        tau(2,0)(3,2) = tau(2,0)(2,2) = tau(0,2)(1,0) = tau(0,2)(1,1) =
                (1.0 / (576 * sin8_beta_2)) * (-192 -27 * beta2 + 4 * (143 - 27 * beta2) * cos2_beta_2 + 12 * (115 - 9 * beta2) *
                                               cos4_beta_2 - 1296 * cos6_beta_2 - 464 * cos8_beta_2 + 12 * b_sin_cos_beta_2 *
                                               (6 - 38 * cos4_beta_2 - 28 * cos6_beta_2));
        tau(0,2)(2,1) = tau(1,1)(2,1) = tau(2,0)(2,1) =
                (1.0 / (576 * sin8_beta_2)) * (1216 + 171 * beta2 - 16 * (133 - 9 * beta2) * cos2_beta_2 - 12 * (556 - 33 * beta2) *
                                               cos4_beta_2 + 16 * (362 - 9 * beta2) * cos6_beta_2 + 1792 * cos8_beta_2 + 6 *
                                               b_sin_cos_beta_2 * (77 - 250 * cos2_beta_2 + 476 * cos4_beta_2 + 264 * cos6_beta_2));
        tau(1,1)(2,2) = tau(1,1)(1,1) =
                (1.0 / (2304 * sin8_beta_2)) * (-480 - 27 * beta2 - 4 * (41 + 27 * beta2) * cos2_beta_2 + 12 * (463 - 9 * beta2) *
                                                cos4_beta_2 - 5568 * cos6_beta_2 + 656 * cos8_beta_2 - 24 * b_sin_cos_beta_2 *
                                                (39 - 57 * cos2_beta_2 + 17 * cos4_beta_2 + 46 * cos6_beta_2));
        tau(1,1)(3,3) = tau(2,0)(3,3) = tau(0,2)(0,0) = tau(1,1)(0,0) = 0;

        break;
    case FIRST_ORDER_ALGEBRAIC_TRIGONOMETRIC:
        tau(0,1)(3,0) = tau(1,0)(3,0) = tau(0,1)(3,3) = tau(1,0)(0,0) =
                (c1 / 48.0) * (-3.0 * c3 * (beta * sin_beta_2 + 2.0 * cos_beta_2 - beta2 * cos_beta_2 - 2.0 * cos3_beta_2) +
                             c4 * (6.0 * beta * sin_beta_2 + beta3 * sin_beta_2 + 12.0 * cos_beta_2 - 9.0 * beta2 * cos_beta_2 -
                                   12.0 * cos3_beta_2 + 6.0 * b_sin_cos_beta_2 * cos_beta_2));
        tau(0,1)(3,1) = tau(0,1)(3,2) = tau(1,0)(2,0) = tau(1,0)(1,0) =
                (c2 / 96.0) * (c3 * (beta * sin_beta_2 * (45.0- 4.0 * beta2 - 12.0 * cos2_beta_2) + 3.0 * cos_beta_2 *
                                   (14.0 - 9.0 * beta2 - 14.0 * cos2_beta_2)) - c4 * (beta * sin_beta_2 * (39.0 + 18.0 * cos2_beta_2 - beta2) +
                                    cos_beta_2 * (78.0 - 36.0 * beta2 - beta4 - 6.0 * (2.0  * beta2 + 13.0 ) * cos2_beta_2)));
        tau(1,0)(3,1) = tau(0,1)(2,0) = tau(0,1)(2,2) = tau(1,0)(1,1) =
                (beta * c2 * sin_beta_2 / 16.0) * (c3 * (8.0 - beta2 - 8.0 * cos2_beta_2 - 2.0 * b_sin_cos_beta_2) +
                                                 beta * c4 * (beta + 2.0 * beta * cos2_beta_2 - 6.0 * sin_beta_2 * cos_beta_2));
        tau(1,0)(3,2) = tau(1,0)(2,2) = tau(0,1)(1,0) = tau(0,1)(1,1) =
                (c2 / 96.0) * (c3 * (beta * (15.0 - 2.0 * beta2 - 84.0 * cos2_beta_2) * sin_beta_2 + 3.0 * (58.0 - 7.0 * beta2 - (58.0 - 4.0 * beta2) *
                             cos2_beta_2) * cos_beta_2) + c4 * (beta * (15.0 + beta2) * sin_beta_2 - 6.0 * beta * (25.0 - 2.0 * beta2) *
                             sin_beta_2 * cos2_beta_2 + (30.0 - 12.0 * beta2 - beta4 - (30.0 - 72.0 * beta2) * cos2_beta_2) * cos_beta_2));
        tau(1,0)(3,3) = tau(0,1)(0,0) = 0.0;
        tau(0,1)(2,1) = tau(1,0)(2,1) =
                (1.0 / 96.0) * (-3.0 * c3_2* (4.0 - beta2 - 4.0 * cos_beta - beta * sin_beta) - 6.0 * c3 * c4 * (2.0 * beta2 - 3.0 * beta * sin_beta +
                              beta2 * cos_beta) + beta * c4_2 * (beta3 - 3.0 * beta2 * sin_beta - 6.0 * beta * cos_beta + 6.0 * sin_beta));


        tau(0,2)(3,0) = tau(2,0)(3,0) = tau(0,2)(3,3) = tau(2,0)(0,0) =
                (-c1 / 48.0) * (3 * c3 * (3 * beta * sin_beta_2 - 2 * cos_beta_2 - beta2  * cos_beta_2 + 2 * cos3_beta_2) +
                     c4 * (-beta * (36 +beta2) * sin_beta_2 + 3 * (8 + 3 * beta2 - 8 * cos2_beta_2 + 2 * b_sin_cos_beta_2) * cos_beta_2));
        tau(1,1)(3,0) = (c1 / 48.0) *
                (-3 * c3 * (beta * sin_beta_2 + (2 - beta2 - 2 * cos2_beta_2) * cos_beta_2) + c4 * (beta * (6 + beta2) * sin_beta_2 + 3 *
                                                                  (4 - 3 * beta2 - 4 * cos2_beta_2 + 2 * b_sin_cos_beta_2) * cos_beta_2));
        tau(0,2)(3,1) = tau(0,2)(3,2) = tau(2,0)(2,0) = tau(2,0)(1,0) =
                (c2 / 96.0) * (c3 * (beta * (51 - 4 * beta2 + 12 * cos2_beta_2) * sin_beta_2 - 7 * (6 + 3 * beta2 - 6 * cos2_beta_2) *
                       cos_beta_2) - c4 * (3 * beta * (69 - 5 * beta2 + 2 * cos2_beta_2) * sin_beta_2 - (162 + 78 * beta2 + beta4 -
                                           6 * (27 + 2 * beta2) * cos2_beta_2) * cos_beta_2));
        tau(1,1)(3,1) = tau(1,1)(2,0) =
                (c2 / 96.0) * (c3 * (beta * (15 - 4 * beta2) * sin_beta_2 - 3 * (6 - beta2 + 2 * (2 * beta * sin_beta_2 - 3 * cos_beta_2) *
                             cos_beta_2) * cos_beta_2) - c4 * (beta * (63 - 19 * beta2) * sin_beta_2 - (66 - 18 * beta2 + beta4 + 6 *
                                                    (7 * beta * sin_beta_2 - (11 - 2 * beta2) * cos_beta_2) * cos_beta_2) * cos_beta_2));
        tau(2,0)(3,1) = tau(0,2)(2,0) = tau(0,2)(2,2) = tau(2,0)(1,1) =
                (c2 / 16.0) * (c3 * (-beta3 * sin_beta_2 - (16 - 4 * beta2 - 2 * (2 * beta * sin_beta_2 + (8 - beta2) * cos_beta_2) *
                             cos_beta_2) * cos_beta_2) + c4 * (5 * beta3 * sin_beta_2 + 2 * beta * (-4 * beta + ((6 - beta2) * sin_beta_2 +
                                                               beta * cos_beta_2) * cos_beta_2) * cos_beta_2));
        tau(1,1)(3,2) = tau(1,1)(1,0) =
                (c2 / 48.0) * (3 * c3 * (beta* sin_beta_2 + (2 - beta2 - 2 * cos2_beta_2) * cos_beta_2) - c4 * (beta * (6 + beta2) *
                             sin_beta_2 + 3 * (4 - 3 * beta2 + 2 * (beta * sin_beta_2 - 2 * cos_beta_2) * cos_beta_2) * cos_beta_2));
        tau(2,0)(3,2) = tau(2,0)(2,2) = tau(0,2)(1,0) = tau(0,2)(1,1) =
                (c2 / 96.0) * (-c3 * (beta * (15 + 2 * beta2) * sin_beta_2 - 3 * (38 - beta2 - 2 * (2 * beta * sin_beta_2 + (19 + 2 * beta2) *
                         cos_beta_2) * cos_beta_2) * cos_beta_2) + c4 * (beta * (63  + 11 * beta2) * sin_beta_2 - (66 + 6 * beta2 + beta4 +
                         6 * (beta * (19 + 2 * beta2) * sin_beta_2 - (11 + 8 * beta2) * cos_beta_2) * cos_beta_2) * cos_beta_2));
        tau(0,2)(2,1) = tau(1,1)(2,1) = tau(2,0)(2,1) =
                (beta / 96.0) * (3 * c3_2 * (beta - sin_beta) - 6 * c3 * c4 * (4 * beta - 3 * sin_beta + beta * cos_beta) +
                               c4_2 * (48 * beta + beta3 - 30 * sin_beta + 3 * beta2 * sin_beta - 18 * beta * cos_beta));
        tau(1,1)(2,2) = tau(1,1)(1,1) =
                (c2 / 96.0) * (-c3 * (beta * (15 + 2 * beta2) * sin_beta_2 - (18 + 9 * beta2 + 6 * (2 * beta * sin_beta_2 - (2 * beta2 + 3) *
                         cos_beta_2) * cos_beta_2) * cos_beta_2) + c4 * (beta * (63 + 11 * beta2) * sin_beta_2 - (66 + 18 * beta2 + beta4 +
                         (6 * beta * (7 + 2 * beta2) * sin_beta_2 - 6 * (11 + 4 * beta2) * cos_beta_2) * cos_beta_2) * cos_beta_2));
        tau(1,1)(3,3) = tau(2,0)(3,3) = tau(0,2)(0,0) = tau(1,1)(0,0) = 0;

        break;
    }

    GLboolean is_parallel = GL_TRUE;

    if (is_parallel)
    {
        GLint num_of_faces = _face.size();

        #pragma omp parallel for
        for (GLint i = 0; i < num_of_faces; i++)
        {
            Face* fit = &_face[i];
            LinearCombination3 *lc[3] = {nullptr, nullptr, nullptr};

            pair<GLint, GLint> index_pair((*fit)[0], (*fit)[1]);
            lc[0] = _edges[index_pair]->boundary;
            lc[1] = _edges[index_pair]->next->boundary;
            lc[2] = _edges[index_pair]->next->next->boundary;

            switch(type)
            {
            case CUBIC_POLYNOMIAL:
                fit->C0 = new CubicPolinomialSurface();
                break;
            case QUARTIC_POLYNOMIAL:
                fit->C0 = new QuarticPolinomialSurface();
                break;
            case SECOND_ORDER_TRIGONOMETRIC:
                fit->C0 = new SecondOrderTrigonometricSurface(beta);
                break;
            case SECOND_ORDER_HYPERBOLIC:
                fit->C0 = new SecondOrderHyperbolicSurface(beta);
                break;
            case FIRST_ORDER_ALGEBRAIC_TRIGONOMETRIC:
                fit->C0 = new FirstOrderAlgebraicTrigonometricSurface(beta);
                break;
            }

            fit->C0->SetData(0, 0, (*lc[0])[0]);
            fit->C0->SetData(1, 0, (*lc[0])[1]);
            fit->C0->SetData(2, 0, (*lc[0])[2]);
            fit->C0->SetData(3, 0, (*lc[1])[0]);
            fit->C0->SetData(3, 1, (*lc[1])[1]);
            fit->C0->SetData(3, 2, (*lc[1])[2]);
            fit->C0->SetData(3, 3, (*lc[2])[0]);
            fit->C0->SetData(2, 2, (*lc[2])[1]);
            fit->C0->SetData(1, 1, (*lc[2])[2]);

            DCoordinate3 &p_111 = (*fit->C0)(2, 1);

            GLdouble denom = 0;
            GLuint gamma = 3;

            for (GLuint g = 1; g <= gamma; g ++)
            {
                GLdouble sg = 0;
                for(GLuint z = 0; z <= g; z++)
                {
                    sg += combination(g, z) * tau(z, g - z)(2, 1);
                }
               denom += sg * epsilon[g];
            }

            denom = -denom;

            for(GLuint g = 1; g <= gamma; g++)
            {
                DCoordinate3 sp_111(0,0,0);
                for (GLuint z = 0; z <= g; z++)
                {
                    DCoordinate3 ssp_111(0,0,0);
                    for (GLuint r = 0; r <= 3; r++)
                    {
                        if (r != 2)
                        {
                            for (GLuint s = 0; s <= r; s++)
                            {
                                ssp_111 += tau(z, g - z)(r, s) * (*fit->C0)(r, s);
                            }
                        }
                    }

                    for (GLuint s = 0; s <= 2; s++) {
                        if( s != 1)
                        {
                            ssp_111 += tau(z, g - z)(2, s) * (*fit->C0)(2, s);
                        }
                    }
                    sp_111 += ssp_111 * combination(g, z);
                }
                p_111 += sp_111 * epsilon[g];
            }

            //_data(2, 1) = p / denom;
            p_111 /= denom;

            /*fit->img_C0 = fit->C0->GenerateImage(div_point_count, usage_flag);

            if (!fit->img_C0)
            {
                cout<<"ERROR! Could not generate image";
            }

            if (!fit->img_C0->UpdateVertexBufferObjects(usage_flag)){
                cout<<"Could not update vertex buffer objects"<<endl;
            }*/
        }

        for ( vector<Face>::iterator fit = _face.begin(); fit != _face.end(); fit++)
        {
            fit->img_C0 = fit->C0->GenerateImage(div_point_count, usage_flag);

            if (!fit->img_C0)
            {
                cout<<"ERROR! Could not generate image";
            }

            if (!fit->img_C0->UpdateVertexBufferObjects(usage_flag)){
                cout<<"Could not update vertex buffer objects"<<endl;
            }
        }
    }
    else
    {
        for ( vector<Face>::iterator fit = _face.begin(); fit != _face.end(); fit++)
        {
            LinearCombination3 *lc[3] = {nullptr, nullptr, nullptr};

            pair<GLint, GLint> index_pair((*fit)[0], (*fit)[1]);
            lc[0] = _edges[index_pair]->boundary;
            lc[1] = _edges[index_pair]->next->boundary;
            lc[2] = _edges[index_pair]->next->next->boundary;

            switch(type)
            {
            case CUBIC_POLYNOMIAL:
                fit->C0 = new CubicPolinomialSurface();
                break;
            case QUARTIC_POLYNOMIAL:
                fit->C0 = new QuarticPolinomialSurface();
                break;
            case SECOND_ORDER_TRIGONOMETRIC:
                fit->C0 = new SecondOrderTrigonometricSurface(beta);
                break;
            case SECOND_ORDER_HYPERBOLIC:
                fit->C0 = new SecondOrderHyperbolicSurface(beta);
                break;
            case FIRST_ORDER_ALGEBRAIC_TRIGONOMETRIC:
                fit->C0 = new FirstOrderAlgebraicTrigonometricSurface(beta);
                break;
            }

            fit->C0->SetData(0, 0, (*lc[0])[0]);
            fit->C0->SetData(1, 0, (*lc[0])[1]);
            fit->C0->SetData(2, 0, (*lc[0])[2]);
            fit->C0->SetData(3, 0, (*lc[1])[0]);
            fit->C0->SetData(3, 1, (*lc[1])[1]);
            fit->C0->SetData(3, 2, (*lc[1])[2]);
            fit->C0->SetData(3, 3, (*lc[2])[0]);
            fit->C0->SetData(2, 2, (*lc[2])[1]);
            fit->C0->SetData(1, 1, (*lc[2])[2]);

            DCoordinate3 &p_111 = (*fit->C0)(2, 1);

            GLdouble denom = 0;
            GLuint gamma = 3;

            for (GLuint g = 1; g <= gamma; g ++)
            {
                GLdouble sg = 0;
                for(GLuint z = 0; z <= g; z++)
                {
                    sg += combination(g, z) * tau(z, g - z)(2, 1);
                }
               denom += sg * epsilon[g];
            }

            denom = -denom;

            for(GLuint g = 1; g <= gamma; g++)
            {
                DCoordinate3 sp_111(0,0,0);
                for (GLuint z = 0; z <= g; z++)
                {
                    DCoordinate3 ssp_111(0,0,0);
                    for (GLuint r = 0; r <= 3; r++)
                    {
                        if (r != 2)
                        {
                            for (GLuint s = 0; s <= r; s++)
                            {
                                ssp_111 += tau(z, g - z)(r, s) * (*fit->C0)(r, s);
                            }
                        }
                    }

                    for (GLuint s = 0; s <= 2; s++) {
                        if( s != 1)
                        {
                            ssp_111 += tau(z, g - z)(2, s) * (*fit->C0)(2, s);
                        }
                    }
                    sp_111 += ssp_111 * combination(g, z);
                }
                p_111 += sp_111 * epsilon[g];
            }

            //_data(2, 1) = p / denom;
            p_111 /= denom;

            fit->img_C0 = fit->C0->GenerateImage(div_point_count, usage_flag);

            if (!fit->img_C0)
            {
                cout<<"ERROR! Could not generate image";
            }

            if (!fit->img_C0->UpdateVertexBufferObjects(usage_flag)){
                cout<<"Could not update vertex buffer objects"<<endl;
            }
        }
    }

}


void TriangulatedHalfEdgeDataStructure3::GenerateG1TriangularSurfaces(BoundaryType type, GLdouble beta,
                                                                      GLdouble theta_1, GLdouble theta_2, GLdouble theta_3,
                                                                      GLuint div_point_count, GLenum usage_flag)
{
    // phi lookup tables
    GLuint rho = 3;

    RowMatrix<GLdouble> theta(rho + 1);

    theta[0] = 0;
    theta[1] = theta_1;
    theta[2] = theta_2;
    theta[3] = theta_3;

    RowMatrix< TriangularMatrix<GLdouble> > phi(rho + 1); // phi[0](k, l), phi[1](k, l), phi[2](k, l), ..., phi[rho](k, l)
    TriangularMatrix<GLdouble> PHI(4);

    for (GLuint r = 0; r <= rho; r++) {
        phi[r].ResizeRows(4);
    }

    RowMatrix<GLdouble> c, ch;
    c.ResizeColumns(5);
    ch.ResizeColumns(5);

    GLdouble sin_beta = sin(beta), sin_beta2 = sin_beta * sin_beta, sin_beta3 = sin_beta * sin_beta2,
            //sin_beta4 = sin_beta2 * sin_beta2,
            sin_2beta = sin(2.0 * beta), sin_3beta = sin(3.0 * beta), sin_4beta = sin(4.0 * beta),
            cos_beta = cos(beta), cos_beta2 = cos_beta * cos_beta, cos_beta3 = cos_beta * cos_beta2,
            cos_beta4 = cos_beta2 * cos_beta2,
            cos_2beta = cos(2.0 * beta);

    GLdouble sin_beta_2 = sin(beta / 2.0), sin_beta_2_2 = sin_beta_2 * sin_beta_2, sin_beta_2_3 = sin_beta_2 * sin_beta_2_2,
            sin_beta_2_4 = sin_beta_2_2 * sin_beta_2_2, sin_beta_2_8 = sin_beta_2_4 * sin_beta_2_4,
            sin_3beta_2 = sin(3.0 * beta / 2.0), sin_5beta_2 = sin(5.0 * beta / 2.0), sin_7beta_2 = sin(7.0 * beta / 2.0);
    GLdouble cos_beta_2 = cos(beta / 2.0), cos_beta_2_2 = cos_beta_2 * cos_beta_2,
            cos_3beta_2 = cos(3.0 * beta / 2.0);
    GLdouble denom = 96.0 * sin_beta_2_8;
    GLdouble denom2 = 6.0 * (1.0 - 4.0 * cos_beta + 6.0 * cos_beta2 - 4.0 * cos_beta3 + cos_beta4);
    GLdouble denom3 = denom2 / 2.0;
    GLdouble d1 = beta - sin_beta, d2 = d1 * (2.0 * sin_beta - beta - beta * cos_beta);

    c[0] = 1.0 / sin_beta_2_4;
    c[1] = 1.0 / sin_beta_2_4 * 4.0 * cos_beta_2;
    c[2] = 1.0 / sin_beta_2_4 * (4.0 * cos_beta_2_2 + 2.0);
    c[3] = c[1];
    c[4] = c[0];


    GLdouble cosh1 = cosh(beta), cosh2 = cosh1 * cosh1, cosh3 = cosh2 * cosh1,
            cosh4 = cosh2 * cosh2;
    GLdouble sinh1 = sinh(beta), sinh2 = sinh(2.0 * beta), sinh3 = sinh(3.0 * beta),
            sinh4 = sinh(4.0 * beta);
    GLdouble sinh_beta = sinh(beta/2), sinh_beta2 = sinh_beta * sinh_beta, sinh_beta4 = sinh_beta2 * sinh_beta2,
            sinh_beta8 = sinh_beta4 * sinh_beta4;
    GLdouble denomh = 1536.0 * sinh_beta8;
    GLdouble denomh2 = 1.0 - 4.0 * cosh1 + 6.0 * cosh2 - 4.0 * cosh3 + cosh4;

    GLdouble a = 2.0 * sin_beta - beta - beta * cos_beta,
            b = beta - sin_beta, b2 = b * b;
    GLdouble denom_a = a * b;

    switch (type)
    {
    case CUBIC_POLYNOMIAL:
        //phi[1](0,0) = phi[1](3,0) = phi[1](3,3) = 0.0;
        phi[1](1,0) = phi[1](3,2) = -9.0 / 10.0;
        phi[1](2,0) = phi[1](3,1) = -3.0 / 5.0;
        phi[1](1,1) = phi[1](2,2) = 6.0 / 5.0;
        phi[1](2,1) = 3.0 / 10.0;

        phi[2](1,0) = phi[2](3,2) = phi[2](2,1) = -18.0;
        phi[2](2,0) = phi[2](3,1) = 0.0;
        phi[2](1,1) = phi[2](2,2) = 36.0;
        //phi[2](3,3) = phi[2](3,0) = phi[2](0,0)

        phi[3](0,0) = phi[3](3,3) =  36.0;
        phi[3](1,0) = phi[3](3,2) = -108.0;
        phi[3](1,1) = phi[3](2,2) = 324.0;
        phi[3](2,0) = phi[3](3,1) = 108.0;
        phi[3](2,1) = -324.0;
        phi[3](3,0) = -36.0;
        break;

    case QUARTIC_POLYNOMIAL:
        //phi[1](0,0) = phi[1](3,0) = phi[1](3,3) = 0.0;
        phi[1](1,0) = phi[1](3,2) = -52.0 / 35.0;
        phi[1](2,0) = phi[1](3,1) = -24.0 / 35.0;
        phi[1](1,1) = phi[1](2,2) = 66.0 / 35.0;
        phi[1](2,1) = 2.0 / 7.0;

        phi[2](1,0) = phi[2](3,2) = -204.0 / 5.0;
        phi[2](2,0) = phi[2](3,1) = 36.0 / 5.0;
        phi[2](1,1) = phi[2](2,2) = 324.0 / 5.0;
        phi[2](2,1) = -156.0 / 5.0;
        //phi[2](3,3) = phi[2](3,0) = phi[2](0,0)

        phi[3](0,0) = phi[3](3,3) = 192.0;
        phi[3](1,0) = phi[3](3,2) = -336.0;
        phi[3](1,1) = phi[3](2,2) = 624.0;
        phi[3](2,0) = phi[3](3,1) = 240.0;
        phi[3](2,1) = - 528.0;
        phi[3](3,0) = - 96.0;
        break;

    case SECOND_ORDER_TRIGONOMETRIC:
        //phi[1](0,0) = phi[1](3,0) = phi[1](3,3) = 0.0;
        phi[1](1,0) = phi[1](3,2) = ((24.0 + 4.0 * cos_beta - 18.0 * cos_beta2 + 5.0 * cos_beta3) * sin_beta -
                3.0 * beta * (6.0 + 2.0 * cos_beta - 3.0 * cos_beta2)) / denom;
        phi[1](2,0) = phi[1](3,1) = ((-12.0 + 24.0 * cos_beta + 2.0 * cos_beta2 + cos_beta3) * sin_beta +
                3.0 * beta * (2.0 - 2.0 * cos_beta - 5.0 * cos_beta2)) / denom;
        phi[1](1,1) = phi[1](2,2) = ((-32.0 - 12.0 * cos_beta + 34.0 * cos_beta2 - 5.0 * cos_beta3) * sin_beta +
                                     3.0 * beta * (8.0 + 4.0 * cos_beta - 5.0 * cos_beta2 - 2.0 * cos_beta3)) / denom2;
        phi[1](2,1) = ((20.0 - 16.0 * cos_beta - 18 * cos_beta2 - cos_beta3) * sin_beta -
                       3.0 * beta * (4.0  - 7.0 * cos_beta2 - 2.0 * cos_beta3)) / denom2;

        phi[2](1,0) = phi[2](3,2) = ((12.0 + 5.0 * cos_beta + 9.0 * cos_beta2 - 14.0 * cos_beta3) * sin_beta -
                                     3.0 * beta * (6.0 + cos_beta - 3.0 * cos_beta2)) / denom3;
        phi[2](2,0) = phi[2](3,1) = ((-18.0 + 21.0 * cos_beta + 7.0 * cos_beta2 + 2.0 * cos_beta3) * sin_beta +
                                     3.0 * beta * (4.0 - cos_beta - 7.0 * cos_beta2)) / denom3;
        phi[2](1,1) = phi[2](2,2) = ((-16.0 -24.0 * cos_beta + 11.0 * cos_beta2 + 17.0 * cos_beta3) * sin_beta +
                                     3.0 * beta * (10.0 - 7.0 * cos_beta2 + 2.0 * cos_beta - cos_beta3)) / denom3;
        phi[2](2,1) = ((22.0 - 2.0 * cos_beta - 27.0 * cos_beta2 - 5.0 * cos_beta3) * sin_beta -
                       3.0 * beta * (8.0 - 11.0 * cos_beta2 - cos_beta3)) / denom3;
        //phi[2](3,3) = phi[2](3,0) = phi[2](0,0)

        /*phi[3](0,0) = phi[3](3,3) = c[0] * c[0] * 1.0 / 24.0 * (15 * beta - 16.0 * sin_beta3 -
                                                                3.0 * sin_4beta - 3.0 * sin_beta * cos_beta);
        phi[3](3,0) = c[0] * c[4] * 1.0 / 24.0 * (3.0 * beta * cos_beta + 12.0 * beta * cos_2beta + sin_beta * (13.0 - 28.0 * cos_beta));
        phi[3](1,0) = phi[3](3,2) = 1.0 / 4.0 * c[4] * (c[1] * 5.0 / 24.0  * (18.0 * sin_beta - 9.0 * sin_2beta + 4.0 * (sin_3beta - 3.0 * beta) * cos_beta_2)+
                c[2] * 1.0 / 48.0 * (6.0 * beta - 30.0 * sin_beta - 3.0 * sin_2beta - 8.0 * sin_3beta + 54.0 * beta * cos_beta));
        phi[3](3,1) = 1.0 / 4.0 * c[4] * (c[1] * 1.0 / 4.0 * (-7.0 * sin_beta_2 +
                                     17.0 / 3.0 * sin_3beta_2 + 2.0 * sin_5beta_2 - beta * cos_beta_2 * (17.0 * cos_beta - 7.0)) +
                c[2] * 1.0 / 48.0 * (6.0 * beta - 30.0 * sin_beta - 3.0 * sin_2beta - 8.0 * sin_3beta + 54.0 * beta * cos_beta));
        phi[3](2,0) = 1.0 / 4.0 * c[0] * (c[3] * 1.0 / 24.0 * (-42.0 * sin_beta_2 + 34.0 * sin_3beta_2 + 12.0 * sin_5beta_2 -
                                                               9.0 * beta * cos_beta_2 - 51.0 * beta * cos_3beta_2) +
                c[2] * 1.0 / 48.0 * (6.0 * beta - 30.0 * sin_beta - 3.0 * sin_2beta - 8.0 * sin_3beta + 54.0 * beta * cos_beta));
        phi[3](1,1) = phi[3](2,2) = c[3] * c[3] * 1.0 / 16.0 * 1.0 / 48.0 * (444.0 * beta - 435 * sin_beta + 78 * sin_2beta - 67 * sin_3beta + 36 * beta * cos_beta) +
                c[2] * c[3] * 1.0 / 8.0 * (-3.0/4.0 * sin_beta_2_3 * (9.0 * cos_beta + 11.0) + 5.0 / 96.0 * (82.0 * sin_beta_2 + 3.0 * sin_3beta_2 + sin_5beta_2
                                                                                                             - 48.0 * beta * cos_beta) +
                                           1.0 / 4.0 * sin_beta_2_2 * (2.0 * (6.0 * sin_beta_2 + sin_3beta_2) - 3* beta * cos_beta_2) -
                                           5.0 / 48.0 * (2.0 * (sin_beta_2 - 1.0 * sin_3beta_2 + sin_5beta_2) + 21.0 * beta * cos_beta_2 +
                                                         3.0 * beta * cos_3beta_2)) +
                c[2] * c[2] * 1.0 / 16.0 * 1.0 / 12.0 * (27.0 * beta - 19.0 * sin_beta + 3.0 * beta * cos_beta - 11.0 * sin_beta * cos_beta);

        phi[3](2,1) = c[2] * c[2] * 1.0 / 16.0 * 1.0 / 12.0 * (27.0 * beta - 19.0 * sin_beta + 3.0 * beta * cos_beta - 11.0 * sin_beta * cos_beta) +
                c[2] * c[3] * 1.0 / 4.0 * (-3.0/4.0 * sin_beta_2_3 * (9.0 * cos_beta + 11.0) + 5.0 / 96.0 * (82.0 * sin_beta_2 + 3.0 * sin_3beta_2 + sin_5beta_2
                                                                                                             - 48.0 * beta * cos_beta) +
                                           1.0 / 4.0 * sin_beta_2_2 * (2.0 * (6.0 * sin_beta_2 + sin_3beta_2) - 3* beta * cos_beta_2) -
                                           5.0 / 48.0 * (2.0 * (sin_beta_2 - 1.0 * sin_3beta_2 + sin_5beta_2) + 21.0 * beta * cos_beta_2 +
                                                         3.0 * beta * cos_3beta_2)) +
                c[1] * c[3] *1.0 / 16.0 * (5.0 / 32.0 * (9.0 * beta - 16.0 * sin_beta - 16.0 * sin_2beta + 36.0 * beta * cos_beta + 3.0 * beta * cos_2beta) +
                                           5.0 / 96.0 * (3.0 * beta * sin_beta2 - 3.0 * beta * (cos_beta - 9.0) * (cos_beta + 1) + 4 * sin_beta * (cos_beta-13.0)) -
                                           9.0 / 4.0 * sin_beta_2_3 * (beta * sin_beta_2 + 4.0 * cos_beta_2) -
                                           3.0 / 8.0 * sin_beta_2_2 * (5.0 * beta - 26.0 * sin_beta + beta * cos_beta));*/
        break;

    case SECOND_ORDER_HYPERBOLIC:
        //phi[1](0,0) = phi[1](3,0) = phi[1](3,3) = 0.0;
        phi[1](1,0) = phi[1](3,2) = (48.0 * beta * (6.0 + 2.0 * cosh1 - 3.0 * cosh2)  - 312.0 * sinh1 -
                                     52.0 * sinh2 + 72.0 * sinh3 - 10.0 * sinh4) / denomh;
        phi[1](2,0) = phi[1](3,1) = (-48.0 * beta * (2.0 - 2.0 * cosh1 - 5.0 * cosh2) + 184.0 * sinh1 -
                                     196.0 * sinh2 - 8.0 * sinh3 - 2.0 * sinh4) / denomh;
        phi[1](1,1) = phi[1](2,2) = (-48.0 * beta * (8.0 + 4.0 * cosh1 - 5.0 * cosh2 - 2.0 * cosh3) +
                                     376.0 * sinh1 + 116.0 * sinh2 - 136.0 * sinh3 + 10.0 * sinh4) / ( 96.0 * denomh2);
        phi[1](2,1) = (48.0 * beta * (4.0 - 7.0 * cosh2 - 2.0 * cosh3) - 248.0 * sinh1 + 132.0 * sinh2 +
                       72.0 * sinh3 + 2.0 * sinh4) / ( 96.0 * denomh2);


        phi[2](1,0) = phi[2](3,2) = (-24.0 * beta * (6.0 + cosh1 - 3.0 * cosh2) + 114.0 * sinh1 -
                                     8.0 * sinh2 + 18.0 * sinh3 - 14.0 * sinh4) / ( 24.0 * denomh2);
        phi[2](2,0) = phi[2](3,1) = (24.0 * beta * (4.0 - cosh1 - 7.0 * cosh2) -130.0 * sinh1 + 88.0 * sinh2 +
                                     14.0 * sinh3 + 2.0 * sinh4) / ( 24.0 * denomh2);
        phi[2](1,1) = phi[2](2,2) = (48.0 * beta * (10.0 + 2.0 * cosh1 - 7.0 * cosh2 - cosh3) -212.0 * sinh1 -
                                     124.0 * sinh2 + 44.0 * sinh3 + 34.0 * sinh4) / ( 48.0 * denomh2);
        phi[2](2,1) = (-48.0 * beta * (8.0 - 11.0 * cosh2 - cosh3) + 244.0 * sinh1 - 36.0 * sinh2 - 108.0 * sinh3 -
                       10.0 * sinh4) / ( 48.0 * denomh2);
        //phi[2](3,3) = phi[2](3,0) = phi[2](0,0)

        break;

    case FIRST_ORDER_ALGEBRAIC_TRIGONOMETRIC:
        //phi[1](0,0) = phi[1](3,0) = phi[1](3,3) = 0.0;
        GLdouble beta2 = beta * beta, beta3 = beta * beta2;
        phi[1](1,0) = phi[1](3,2) = (sin_beta * (-3.0 * beta + 6.0 * sin_beta + 2.0 * beta * cos_beta -
                                             3.0 * sin_2beta + beta * cos_2beta)) / (4.0 * a * b2);
        phi[1](2,0) = phi[1](3,1) = (sin_beta * (- beta - sin_beta + beta2 * sin_beta +
                                             beta * cos_beta + cos_beta * sin_beta)) / (2.0 * a * b2);
        phi[1](1,1) = phi[1](2,2) = (sin_beta2 * (2.0 * beta3 - 4.0 * sin_beta - 4.0 * beta2 * sin_beta +
                                             4.0 * beta * cos_beta + 2.0 * sin_2beta - beta2 * sin_2beta -
                                             4.0 * beta * cos_2beta)) / (4.0 * denom_a * denom_a);
        phi[1](2,1) = (sin_beta2 * (6.0 * beta - 2.0 * sin_beta - 3.0 * beta2 * sin_beta -
                               6.0 * beta * cos_beta + beta3 * cos_beta + sin_2beta)) /(2.0 * denom_a * denom_a);

        phi[2](1,0) = phi[2](3,2) = (sin_beta * (- beta - 2.0 * sin_beta + 2.0 * beta * cos_beta +
                                             sin_2beta - beta * cos_2beta)) / (4.0 * a * b2);
        phi[2](2,0) = phi[2](3,1) = (sin_beta * (- beta - sin_beta + beta2 * sin_beta + beta * cos_beta +
                                             cos_beta * sin_beta)) / (2.0 * a * b2);
        phi[2](1,1) = phi[2](2,2) = (sin_beta2 * (2.0 * beta + 2.0 * beta3 + 4.0 * sin_beta -
                                             4.0 * beta2 * sin_beta - 4.0 * beta * cos_beta - 2.0 * sin_2beta +
                                             beta2 * sin_2beta + 2.0 * beta * cos_2beta)) / (4.0 * denom_a * denom_a);
        phi[2](2,1) = (sin_beta2 * (beta + 2.0 * sin_beta - beta2 * sin_beta - 2.0 * beta * cos_beta +
                               beta3 * cos_beta - sin_2beta + beta * cos_2beta)) / (2.0 * denom_a * denom_a);
        //phi[2](3,3) = phi[2](3,0) = phi[2](0,0)


        phi[3](0,0) = phi[3](3,3) = 1.0 / (2.0 * d1 * d1) * (beta + sin_beta * cos_beta);
        phi[3](1,0) = phi[3](3,2) = sin_beta / (d1 * d2) * sin_beta_2_2 * (sin_beta - beta * (cos_beta + 2.0));
        phi[3](1,1) = phi[3](2,2) = sin_beta2 / (4.0 * d2 * d2) * (2.0 * (beta3 + 3.0 * beta - 2.0 * sin_beta + sin_2beta) -
                                                     beta * (beta * (4.0 * sin_beta + sin_2beta) + 4.0 * cos_beta + 2.0 * cos_2beta));
        phi[3](2,0) = phi[3](3,1) = sin_beta / (2.0 * d1 * d2) * ((beta2 + 1) * sin_beta - beta + (beta - sin_beta) * cos_beta);
        phi[3](2,1) = sin_beta2 / (d2 * d2) * (-1.0 / 2.0 * (3.0 * beta2 + 2.0) * sin_beta +
                                               (1.0 / 2.0 * beta * (beta2 - 2) + sin_beta) * cos_beta + beta + beta * sin_beta2);
        phi[3](3,0) = -1.0 / (2.0 * d1 * d1) * (sin_beta + beta * cos_beta);

        break;
    }

    for (GLuint r = 1 ; r <= rho; r++) {
        for (GLuint k = 0; k <= 3; k++) {
            for (GLuint l = 0; l <= k; l++)
                PHI(k, l) += theta[r] * phi[r](k, l);
        }
    }

    for (vector<Face>::iterator fit = _face.begin(); fit != _face.end(); fit++)
    {
        // calculating number of vertices, unit normal vectors and texture coordinates
        GLuint vertex_count = div_point_count * (div_point_count + 1) / 2;

        // calculating number of triangular faces
        GLuint face_count = (div_point_count - 1) * (div_point_count - 1);

        TriangulatedMesh3 *result = nullptr;

        result = new TriangulatedMesh3(vertex_count, face_count);

        if (!result)
        {
            return;
        }

        //for face indexing
        GLuint current_face = 0;

        pair<GLint, GLint> index_pair((*fit)[0], (*fit)[1]);

        TriangulatedHalfEdgeDataStructure3::HalfEdge *he_ij = _edges[index_pair];
        TriangulatedHalfEdgeDataStructure3::HalfEdge *he_jk = he_ij->next;
        TriangulatedHalfEdgeDataStructure3::HalfEdge *he_ki = he_ij->next->next;

        DCoordinate3 p_i = *he_ij->vertex;//_vertex[(*fit)[0]];//(*fit->C0)(0,0);
        DCoordinate3 p_j = *he_jk->vertex;//_vertex[(*fit)[1]];//(*fit->C0)(3,0);
        DCoordinate3 p_k = *he_ki->vertex;//_vertex[(*fit)[2]];//(*fit->C0)(3,3);

        LinearCombination3* c_jk = he_jk->boundary;
        LinearCombination3* c_ki = he_ki->boundary;
        LinearCombination3* c_ij = he_ij->boundary;

        LinearCombination3* g_i_jk = nullptr;
        LinearCombination3* g_j_ki = nullptr;
        LinearCombination3* g_k_ij = nullptr;

        //TriangularSurface3 *s_ijk = fit->C0;
        TriangularSurface3 *s_i_jk = fit->C0;
        TriangularSurface3 *s_j_ki = nullptr;
        TriangularSurface3 *s_k_ij = nullptr;

        TriangularSurface3 *s_o_ij = (he_ij->opposite ? he_ij->opposite->face->C0 : nullptr);
        TriangularSurface3 *s_o_jk = (he_jk->opposite ? he_jk->opposite->face->C0 : nullptr);
        TriangularSurface3 *s_o_ki = (he_ki->opposite ? he_ki->opposite->face->C0 : nullptr);

        TriangularSurface3 *s_x_kj = nullptr;
        TriangularSurface3 *s_x_ik = nullptr;
        TriangularSurface3 *s_x_ji = nullptr;


        switch (type) {
        case CUBIC_POLYNOMIAL:
            s_x_ik = new CubicPolinomialSurface();
            s_x_ji = new CubicPolinomialSurface();
            s_x_kj = new CubicPolinomialSurface();

            s_j_ki = new CubicPolinomialSurface();
            s_k_ij = new CubicPolinomialSurface();

            g_i_jk = new CubicPolinomialCurve();
            g_j_ki = new CubicPolinomialCurve();
            g_k_ij = new CubicPolinomialCurve();
            break;
        case QUARTIC_POLYNOMIAL:
            s_x_ik = new QuarticPolinomialSurface();
            s_x_ji = new QuarticPolinomialSurface();
            s_x_kj = new QuarticPolinomialSurface();

            s_j_ki = new QuarticPolinomialSurface();
            s_k_ij = new QuarticPolinomialSurface();

            g_i_jk = new QuarticPolinomialCurve();
            g_j_ki = new QuarticPolinomialCurve();
            g_k_ij = new QuarticPolinomialCurve();
            break;
        case SECOND_ORDER_TRIGONOMETRIC:
            s_x_ik = new SecondOrderTrigonometricSurface(beta);
            s_x_ji = new SecondOrderTrigonometricSurface(beta);
            s_x_kj = new SecondOrderTrigonometricSurface(beta);

            s_j_ki = new SecondOrderTrigonometricSurface(beta);
            s_k_ij = new SecondOrderTrigonometricSurface(beta);

            g_i_jk = new SecondOrderTrigonometricCurve(beta);
            g_j_ki = new SecondOrderTrigonometricCurve(beta);
            g_k_ij = new SecondOrderTrigonometricCurve(beta);
            break;
        case SECOND_ORDER_HYPERBOLIC:
            s_x_ik = new SecondOrderHyperbolicSurface(beta);
            s_x_ji = new SecondOrderHyperbolicSurface(beta);
            s_x_kj = new SecondOrderHyperbolicSurface(beta);

            s_j_ki = new SecondOrderHyperbolicSurface(beta);
            s_k_ij = new SecondOrderHyperbolicSurface(beta);

            g_i_jk = new SecondOrderHyperbolicCurve(beta);
            g_j_ki = new SecondOrderHyperbolicCurve(beta);
            g_k_ij = new SecondOrderHyperbolicCurve(beta);
            break;
        case FIRST_ORDER_ALGEBRAIC_TRIGONOMETRIC:
            s_x_ik = new FirstOrderAlgebraicTrigonometricSurface(beta);
            s_x_ji = new FirstOrderAlgebraicTrigonometricSurface(beta);
            s_x_kj = new FirstOrderAlgebraicTrigonometricSurface(beta);

            s_j_ki = new FirstOrderAlgebraicTrigonometricSurface(beta);
            s_k_ij = new FirstOrderAlgebraicTrigonometricSurface(beta);

            g_i_jk = new FirstOrderAlgebraicTrigonometricCurve(beta);
            g_j_ki = new FirstOrderAlgebraicTrigonometricCurve(beta);
            g_k_ij = new FirstOrderAlgebraicTrigonometricCurve(beta);
            break;

        }

        if(!g_i_jk)
        {
            cout<<"no g_i_jk boundary"<<endl;
        }

        if(s_x_kj && s_o_jk)
        {
            LinearCombination3 *lc_x_i[3] = {nullptr, nullptr, nullptr};
            lc_x_i[0] = he_jk->opposite->next->next->boundary;
            lc_x_i[1] = he_jk->opposite->boundary;
            lc_x_i[2] = he_jk->opposite->next->boundary;
            //TriangularSurface3* s_x_i = he_jk->opposite->face->C0;

            s_x_kj->SetData(0, 0, (*lc_x_i[0])[0]);
            s_x_kj->SetData(1, 0, (*lc_x_i[0])[1]);
            s_x_kj->SetData(2, 0, (*lc_x_i[0])[2]);
            s_x_kj->SetData(3, 0, (*lc_x_i[1])[0]);
            s_x_kj->SetData(3, 1, (*lc_x_i[1])[1]);
            s_x_kj->SetData(3, 2, (*lc_x_i[1])[2]);
            s_x_kj->SetData(3, 3, (*lc_x_i[2])[0]);
            s_x_kj->SetData(2, 2, (*lc_x_i[2])[1]);
            s_x_kj->SetData(1, 1, (*lc_x_i[2])[2]);

            s_x_kj->SetData(2, 1, (*s_o_jk)(2,1));

        }

        if(s_x_ik && s_o_ki)
        {
            LinearCombination3 *lc_x_j[3] = {nullptr, nullptr, nullptr};
            lc_x_j[0] = he_ki->opposite->next->next->boundary;
            lc_x_j[1] = he_ki->opposite->boundary;
            lc_x_j[2] = he_ki->opposite->next->boundary;

            s_x_ik->SetData(0, 0, (*lc_x_j[0])[0]);
            s_x_ik->SetData(1, 0, (*lc_x_j[0])[1]);
            s_x_ik->SetData(2, 0, (*lc_x_j[0])[2]);
            s_x_ik->SetData(3, 0, (*lc_x_j[1])[0]);
            s_x_ik->SetData(3, 1, (*lc_x_j[1])[1]);
            s_x_ik->SetData(3, 2, (*lc_x_j[1])[2]);
            s_x_ik->SetData(3, 3, (*lc_x_j[2])[0]);
            s_x_ik->SetData(2, 2, (*lc_x_j[2])[1]);
            s_x_ik->SetData(1, 1, (*lc_x_j[2])[2]);

            s_x_ik->SetData(2, 1, (*s_o_ki)(2,1));
        }

        if(s_x_ji && s_o_ij)
        {
            LinearCombination3 *lc_x_k[3] = {nullptr, nullptr, nullptr};
            lc_x_k[0] = he_ij->opposite->next->next->boundary;
            lc_x_k[1] = he_ij->opposite->boundary;
            lc_x_k[2] = he_ij->opposite->next->boundary;

            s_x_ji->SetData(0, 0, (*lc_x_k[0])[0]);
            s_x_ji->SetData(1, 0, (*lc_x_k[0])[1]);
            s_x_ji->SetData(2, 0, (*lc_x_k[0])[2]);
            s_x_ji->SetData(3, 0, (*lc_x_k[1])[0]);
            s_x_ji->SetData(3, 1, (*lc_x_k[1])[1]);
            s_x_ji->SetData(3, 2, (*lc_x_k[1])[2]);
            s_x_ji->SetData(3, 3, (*lc_x_k[2])[0]);
            s_x_ji->SetData(2, 2, (*lc_x_k[2])[1]);
            s_x_ji->SetData(1, 1, (*lc_x_k[2])[2]);

            s_x_ji->SetData(2, 1, (*s_o_ij)(2,1));

        }

        s_j_ki->SetData(0, 0, (*c_jk)[0]);
        s_j_ki->SetData(1, 0, (*c_jk)[1]);
        s_j_ki->SetData(2, 0, (*c_jk)[2]);
        s_j_ki->SetData(3, 0, (*c_ki)[0]);
        s_j_ki->SetData(3, 1, (*c_ki)[1]);
        s_j_ki->SetData(3, 2, (*c_ki)[2]);
        s_j_ki->SetData(3, 3, (*c_ij)[0]);
        s_j_ki->SetData(2, 2, (*c_ij)[1]);
        s_j_ki->SetData(1, 1, (*c_ij)[2]);

        s_j_ki->SetData(2, 1, (*s_i_jk)(2, 1));

        s_k_ij->SetData(0, 0, (*c_ki)[0]);
        s_k_ij->SetData(1, 0, (*c_ki)[1]);
        s_k_ij->SetData(2, 0, (*c_ki)[2]);
        s_k_ij->SetData(3, 0, (*c_ij)[0]);
        s_k_ij->SetData(3, 1, (*c_ij)[1]);
        s_k_ij->SetData(3, 2, (*c_ij)[2]);
        s_k_ij->SetData(3, 3, (*c_jk)[0]);
        s_k_ij->SetData(2, 2, (*c_jk)[1]);
        s_k_ij->SetData(1, 1, (*c_jk)[2]);

        s_k_ij->SetData(2, 1, (*s_i_jk)(2, 1));

        DCoordinate3 n_i = *he_ij->normal;
        DCoordinate3 f_i = -n_i;
        DCoordinate3 n_j = *he_jk->normal;
        DCoordinate3 f_j = -n_j;
        DCoordinate3 n_k = *he_ki->normal;
        DCoordinate3 f_k = -n_k;

        DCoordinate3 p(1.0, 0, 0), q(0, 1.0, 0), r(0, 0, 1.0);
        GLdouble ds = (1.0 - 2.0 * EPS) / (GLdouble)(div_point_count - 1);
        GLfloat sd = 1.0f / (div_point_count - 1);

        GLboolean is_parallel = GL_FALSE;

        if(is_parallel)
        {
            #pragma omp parallel for
            for (GLint i_j = 0; i_j < (GLint)vertex_count; i_j++)
            {
                GLint i = (GLint)sqrt(i_j * 2);
                if (i * (i + 1) / 2 > i_j)
                {
                    i = i - 1;
                }

                GLint j = i_j - i * (i + 1) / 2;

                GLfloat ss = min(i * sd, 1.0f);
                GLdouble s = min(EPS + i * ds, 1.0 - EPS);

                DCoordinate3 x = p + (q - p) * s; // a [p,q] szakaszon felvett egyenközű osztópont
                DCoordinate3 y = p + (r - p) * s; // a [p,r] szakaszon felvett egyenközű osztópont

                GLdouble dt = i ? (1.0 - 2.0 * EPS) / (GLdouble)i : 0.0;

                GLfloat tt = min(j * sd, 1.0f);
                GLdouble t = min(EPS + j * dt, 1.0 - EPS);

                DCoordinate3 b = x + (y - x) * t;

                GLdouble u = min(beta * b[2] / (1.0 - b[0]), beta);
                GLdouble v = min(beta * b[0] / (1.0 - b[1]), beta);
                GLdouble w = min(beta * b[1] / (1.0 - b[2]), beta);

                /***********************************for p_i and p_jk**********************************/

                TriangularSurface3::PartialDerivatives pd_s_i_jk, pd_s_x_kj;

                if(s_i_jk->CalculatePartialDerivatives(1, 0, beta - u, pd_s_i_jk) == GL_FALSE)
                {
                    cout<<"Could not calculate partial derivatives"<<endl;
                    //return;
                }

                LinearCombination3::Derivatives d_c_jk;

                if(c_jk->CalculateDerivatives(0, u, d_c_jk) == GL_FALSE)
                {
                    cout<<"Could not calculate derivatives"<<endl;
                }

                DCoordinate3 p_jk = d_c_jk(0);//pd_s_i_j_k(0,0);//d_c_j_k(0);
                DCoordinate3 n_jk = pd_s_i_jk(1,0);
                n_jk ^= pd_s_i_jk(1,1);
                n_jk.normalize();

                //if hed_j_k->opposite
                if(s_x_kj->CalculatePartialDerivatives(1, 0, u, pd_s_x_kj) == GL_FALSE)
                {
                    cout<<"Could not calculate partial derivatives"<<endl;
                }

                DCoordinate3 n_kj = pd_s_x_kj(1,0);
                n_kj ^= pd_s_x_kj(1,1);
                n_kj.normalize();

                if (he_jk->opposite)
                {
                    n_jk += n_kj;
                    n_jk.normalize();
                }

                DCoordinate3 b_i = (p_jk - p_i) ^ f_i;
                b_i.normalize();
                DCoordinate3 t_i = f_i ^ b_i;

                DCoordinate3 f_jk = -n_jk;
                DCoordinate3 b_jk = (p_i - p_jk) ^ f_jk;
                b_jk.normalize();
                DCoordinate3 t_jk = f_jk ^ b_jk;

                GLdouble dot_t_i_t_jk = t_i * t_jk;

                GLdouble delta_i = PHI(1, 1) * PHI(2, 2) - dot_t_i_t_jk * dot_t_i_t_jk * PHI(2, 1) * PHI(2, 1);

                RealSquareMatrix inverse_i(2);

                inverse_i(0, 0) = PHI(2, 2);
                inverse_i(0, 1) = -dot_t_i_t_jk * PHI(2, 1);
                inverse_i(1, 0) = inverse_i(0, 1);
                inverse_i(1, 1) = PHI(1, 1);

                for (GLuint m = 0; m <= 1; m++)
                {
                    for (GLuint n = 0; n <= 1; n++)
                    {
                        inverse_i(m, n) /= -delta_i;
                    }
                }

                ColumnMatrix<GLdouble> B_i(2);

                B_i[0] = (p_i * (PHI(1, 0) + PHI(1, 1)) + p_jk * (PHI(2, 1) + PHI(3, 1))) * t_i;
                B_i[1] = (p_i * (PHI(2, 0) + PHI(2, 1)) + p_jk * (PHI(2, 2) + PHI(3, 2))) * t_jk;

                ColumnMatrix<GLdouble> lambda_i(2);

                for (GLuint n = 0; n <= 1; n++)
                {
                    for (GLuint m = 0; m <= 1; m++)
                    {
                        lambda_i[n] += inverse_i(n, m) * B_i[m];
                    }
                }

                (*g_i_jk)[0] = p_i;
                (*g_i_jk)[1] = p_i + lambda_i[0] * t_i;
                (*g_i_jk)[2] = p_jk + lambda_i[1] * t_jk;
                (*g_i_jk)[3] = p_jk;

                /***********************************for p_j and p_ki**********************************/

                TriangularSurface3::PartialDerivatives pd_s_j_ki, pd_s_x_ik;

                if(s_j_ki->CalculatePartialDerivatives(1, 0, beta - v, pd_s_j_ki) == GL_FALSE)
                {
                    cout<<"Could not calculate partial derivatives"<<endl;
                    //return;
                }

                LinearCombination3::Derivatives d_c_ki;

                if(c_ki->CalculateDerivatives(0, v, d_c_ki) == GL_FALSE)
                {
                    cout<<"Could not calculate derivatives"<<endl;
                }

                DCoordinate3 p_ki = d_c_ki(0);//pd_s_j_k_i(0,0);//d_c_k_i(0);
                DCoordinate3 n_ki = pd_s_j_ki(1,0);
                n_ki ^= pd_s_j_ki(1,1);
                n_ki.normalize();

                //if hed_j_k->opposite

                if(s_x_ik->CalculatePartialDerivatives(1, 0, v, pd_s_x_ik) == GL_FALSE)
                {
                    cout<<"Could not calculate partial derivatives"<<endl;
                }

                DCoordinate3 n_ik = pd_s_x_ik(1,0);
                n_ik ^= pd_s_x_ik(1,1);
                n_ik.normalize();

                if(he_ki->opposite)
                {
                    n_ki += n_ik;
                    n_ki.normalize();
                }

                DCoordinate3 b_j = (p_ki - p_j) ^ f_j;
                b_j.normalize();
                DCoordinate3 t_j = f_j ^ b_j;

                DCoordinate3 f_ki = - n_ki;
                DCoordinate3 b_ki = (p_j - p_ki) ^ f_ki;
                b_ki.normalize();
                DCoordinate3 t_ki = f_ki ^ b_ki;

                GLdouble dot_t_j_t_ki = t_j * t_ki;
                GLdouble delta_j = PHI(1, 1) * PHI(2, 2) - dot_t_j_t_ki * dot_t_j_t_ki * PHI(2, 1) * PHI(2, 1);

                RealSquareMatrix inverse_j(2);

                inverse_j(0, 0) = PHI(2, 2);
                inverse_j(0, 1) = -dot_t_j_t_ki * PHI(2, 1);
                inverse_j(1, 0) = inverse_j(0, 1);
                inverse_j(1, 1) = PHI(1, 1);

                for (GLuint m = 0; m <= 1; m++)
                {
                    for (GLuint n = 0; n <= 1; n++)
                    {
                        inverse_j(m, n) /= -delta_j;
                    }
                }

                ColumnMatrix<GLdouble> B_j(2);

                B_j[0] = (p_j * (PHI(1, 0) + PHI(1, 1)) + p_ki * (PHI(2, 1) + PHI(3, 1))) * t_j;
                B_j[1] = (p_j * (PHI(2, 0) + PHI(2, 1)) + p_ki * (PHI(2, 2) + PHI(3, 2))) * t_ki;

                ColumnMatrix<GLdouble> lambda_j(2);

                for (GLuint n = 0; n <= 1; n++)
                {
                    for (GLuint m = 0; m <= 1; m++)
                    {
                        lambda_j[n] += inverse_j(n, m) * B_j[m];
                    }
                }

                (*g_j_ki)[0] = p_j;
                (*g_j_ki)[1] = p_j + lambda_j[0] * t_j;
                (*g_j_ki)[2] = p_ki + lambda_j[1] * t_ki;
                (*g_j_ki)[3] = p_ki;

                /***********************************for p_k and p_ij**********************************/

                TriangularSurface3::PartialDerivatives pd_s_k_ij, pd_s_x_ji;

                if(s_k_ij->CalculatePartialDerivatives(1, 0, beta - w, pd_s_k_ij) == GL_FALSE)
                {
                    cout<<"Could not calculate partial derivatives"<<endl;
                    //return;
                }


                LinearCombination3::Derivatives d_c_ij;

                if(c_ij->CalculateDerivatives(0, w, d_c_ij) == GL_FALSE)
                {
                    cout<<"Could not calculate derivatives"<<endl;
                }

                DCoordinate3 p_ij = d_c_ij(0);//pd_s_k_i_j(0,0);//d_c_i_j(0);
                DCoordinate3 n_ij = pd_s_k_ij(1,0);
                n_ij ^= pd_s_k_ij(1,1);
                n_ij.normalize();

                //if hed_j_k->opposite
                if(s_x_ji->CalculatePartialDerivatives(1, 0, w, pd_s_x_ji) == GL_FALSE)
                {
                    cout<<"Could not calculate partial derivatives"<<endl;
                }

                DCoordinate3 n_ji = pd_s_x_ji(1,0);
                n_ji ^= pd_s_x_ji(1,1);
                n_ji.normalize();

                if(he_ij->opposite)
                {
                    n_ij += n_ji;
                    n_ij.normalize();
                }

                DCoordinate3 b_k = (p_ij - p_k) ^ f_k;
                b_k.normalize();
                DCoordinate3 t_k = f_k ^ b_k;

                DCoordinate3 f_ij = - n_ij;
                DCoordinate3 b_ij = (p_k - p_ij) ^ f_ij;
                b_ij.normalize();
                DCoordinate3 t_ij = f_ij ^ b_ij;

                GLdouble dot_t_k_t_ij = t_k * t_ij;
                GLdouble delta_k = PHI(1, 1) * PHI(2, 2) - dot_t_k_t_ij * dot_t_k_t_ij * PHI(2, 1) * PHI(2, 1);

                RealSquareMatrix inverse_k(2);

                inverse_k(0, 0) = PHI(2, 2);
                inverse_k(0, 1) = -dot_t_k_t_ij * PHI(2, 1);
                inverse_k(1, 0) = inverse_k(0, 1);
                inverse_k(1, 1) = PHI(1, 1);

                for (GLuint m = 0; m <= 1; m++)
                {
                    for (GLuint n = 0; n <= 1; n++)
                    {
                        inverse_k(m, n) /= -delta_k;
                    }
                }

                ColumnMatrix<GLdouble> B_k(2);

                B_k[0] = (p_k * (PHI(1, 0) + PHI(1, 1)) + p_ij * (PHI(2, 1) + PHI(3, 1))) * t_k;
                B_k[1] = (p_k * (PHI(2, 0) + PHI(2, 1)) + p_ij * (PHI(2, 2) + PHI(3, 2))) * t_ij;

                ColumnMatrix<GLdouble> lambda_k(2);

                for (GLuint n = 0; n <= 1; n++)
                {
                    for (GLuint m = 0; m <= 1; m++)
                    {
                        lambda_k[n] += inverse_k(n, m) * B_k[m];
                    }
                }

                (*g_k_ij)[0] = p_k;
                (*g_k_ij)[1] = p_k + lambda_k[0] * t_k;
                (*g_k_ij)[2] = p_ij + lambda_k[1] * t_ij;
                (*g_k_ij)[3] = p_ij;

                /**************************************calculate averaged surface point********************************/
                LinearCombination3::Derivatives d_g_i_jk, d_g_j_ki, d_g_k_ij;

                if(g_i_jk->CalculateDerivatives(0, min(max(0.0,beta * (1.0 - b[0])), beta), d_g_i_jk) == GL_FALSE)
                {
                    cout<<"Could not calculate derivatives, g_i_jk"<<endl;
                }

                DCoordinate3 s_i = d_g_i_jk[0];

                if(g_j_ki->CalculateDerivatives(0, min(max(0.0,beta * (1.0 - b[1])), beta), d_g_j_ki) == GL_FALSE)
                {
                    cout<<"Could not calculate derivatives, g_j_ki"<<endl;
                }

                DCoordinate3 s_j = d_g_j_ki[0];

                if(g_k_ij->CalculateDerivatives(0, min(max(0.0,beta * (1.0 - b[2])), beta), d_g_k_ij) == GL_FALSE)
                {
                    cout<<"Could not calculate derivatives, g_k_ij"<<endl;
                }

                DCoordinate3 s_k = d_g_k_ij[0];

                GLdouble denomiator = b[0] * b[1] + b[0] * b[2] + b[1] * b[2];
                GLdouble w_i = b[1] * b[2] / denomiator;
                GLdouble w_j = b[0] * b[2] / denomiator;
                GLdouble w_k = b[0] * b[1] / denomiator;

                DCoordinate3 p_s_i_j_k = w_i * s_i + w_j * s_j + w_k * s_k;

                GLuint index[4];

                index[0] = i * (i + 1) / 2 + j;
                index[1] = index[0] + i + 1;
                index[2] = index[1] + 1;
                index[3] = index[0] + 1;

                (*result)._vertex[index[0]] = p_s_i_j_k;

                // texture coordinates
                (*result)._tex[index[0]].s() = ss;
                (*result)._tex[index[0]].t() = tt;

                // faces
                #pragma omp critical
                if (i < (GLint)div_point_count - 1 && j < (GLint)div_point_count - 1)
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

                //normals
                if (j == 0)
                {
                    (*result)._normal[index[0]] = n_ij;
                }

                if (i == div_point_count - 1)
                {
                    (*result)._normal[index[0]] = n_jk;
                }

                if (i == j)
                {
                    (*result)._normal[index[0]] = n_ki;
                }

            }

        }else {

            for(GLuint i = 0; i < div_point_count; i++)
            {
                GLfloat ss = min(i * sd, 1.0f);
                GLdouble s = min(EPS + i * ds, 1.0 - EPS);

                DCoordinate3 x = p + (q - p) * s; // a [p,q] szakaszon felvett egyenközű osztópont
                DCoordinate3 y = p + (r - p) * s; // a [p,r] szakaszon felvett egyenközű osztópont

                GLdouble dt = i ? (1.0 - 2.0 * EPS) / (GLdouble)i : 0.0;

                for(GLuint j = 0; j <= i; j++)
                {
                    GLfloat tt = min(j * sd, 1.0f);
                    GLdouble t = min(EPS + j * dt, 1.0 - EPS);

                    DCoordinate3 b = x + (y - x) * t;

                    GLdouble u = min(beta * b[2] / (1.0 - b[0]), beta);
                    GLdouble v = min(beta * b[0] / (1.0 - b[1]), beta);
                    GLdouble w = min(beta * b[1] / (1.0 - b[2]), beta);

                    /***********************************for p_i and p_jk**********************************/

                    TriangularSurface3::PartialDerivatives pd_s_i_jk, pd_s_x_kj;

                    if(s_i_jk->CalculatePartialDerivatives(1, 0, beta - u, pd_s_i_jk) == GL_FALSE)
                    {
                        cout<<"Could not calculate partial derivatives"<<endl;
                        //return;
                    }

                    LinearCombination3::Derivatives d_c_jk;

                    if(c_jk->CalculateDerivatives(0, u, d_c_jk) == GL_FALSE)
                    {
                        cout<<"Could not calculate derivatives"<<endl;
                    }

                    DCoordinate3 p_jk = d_c_jk(0);//pd_s_i_j_k(0,0);//d_c_j_k(0);
                    DCoordinate3 n_jk = pd_s_i_jk(1,0);
                    n_jk ^= pd_s_i_jk(1,1);
                    n_jk.normalize();

                    //if hed_j_k->opposite
                    if(s_x_kj->CalculatePartialDerivatives(1, 0, u, pd_s_x_kj) == GL_FALSE)
                    {
                        cout<<"Could not calculate partial derivatives"<<endl;
                    }

                    DCoordinate3 n_kj = pd_s_x_kj(1,0);
                    n_kj ^= pd_s_x_kj(1,1);
                    n_kj.normalize();

                    if (he_jk->opposite)
                    {
                        n_jk += n_kj;
                        n_jk.normalize();
                    }

                    DCoordinate3 b_i = (p_jk - p_i) ^ f_i;
                    b_i.normalize();
                    DCoordinate3 t_i = f_i ^ b_i;

                    DCoordinate3 f_jk = -n_jk;
                    DCoordinate3 b_jk = (p_i - p_jk) ^ f_jk;
                    b_jk.normalize();
                    DCoordinate3 t_jk = f_jk ^ b_jk;

                    GLdouble dot_t_i_t_jk = t_i * t_jk;

                    GLdouble delta_i = PHI(1, 1) * PHI(2, 2) - dot_t_i_t_jk * dot_t_i_t_jk * PHI(2, 1) * PHI(2, 1);

                    RealSquareMatrix inverse_i(2);

                    inverse_i(0, 0) = PHI(2, 2);
                    inverse_i(0, 1) = -dot_t_i_t_jk * PHI(2, 1);
                    inverse_i(1, 0) = inverse_i(0, 1);
                    inverse_i(1, 1) = PHI(1, 1);

                    for (GLuint m = 0; m <= 1; m++)
                    {
                        for (GLuint n = 0; n <= 1; n++)
                        {
                            inverse_i(m, n) /= -delta_i;
                        }
                    }

                    ColumnMatrix<GLdouble> B_i(2);

                    B_i[0] = (p_i * (PHI(1, 0) + PHI(1, 1)) + p_jk * (PHI(2, 1) + PHI(3, 1))) * t_i;
                    B_i[1] = (p_i * (PHI(2, 0) + PHI(2, 1)) + p_jk * (PHI(2, 2) + PHI(3, 2))) * t_jk;

                    ColumnMatrix<GLdouble> lambda_i(2);

                    for (GLuint n = 0; n <= 1; n++)
                    {
                        for (GLuint m = 0; m <= 1; m++)
                        {
                            lambda_i[n] += inverse_i(n, m) * B_i[m];
                        }
                    }

                    (*g_i_jk)[0] = p_i;
                    (*g_i_jk)[1] = p_i + lambda_i[0] * t_i;
                    (*g_i_jk)[2] = p_jk + lambda_i[1] * t_jk;
                    (*g_i_jk)[3] = p_jk;

                    /***********************************for p_j and p_ki**********************************/

                    TriangularSurface3::PartialDerivatives pd_s_j_ki, pd_s_x_ik;

                    if(s_j_ki->CalculatePartialDerivatives(1, 0, beta - v, pd_s_j_ki) == GL_FALSE)
                    {
                        cout<<"Could not calculate partial derivatives"<<endl;
                        //return;
                    }

                    LinearCombination3::Derivatives d_c_ki;

                    if(c_ki->CalculateDerivatives(0, v, d_c_ki) == GL_FALSE)
                    {
                        cout<<"Could not calculate derivatives"<<endl;
                    }

                    DCoordinate3 p_ki = d_c_ki(0);//pd_s_j_k_i(0,0);//d_c_k_i(0);
                    DCoordinate3 n_ki = pd_s_j_ki(1,0);
                    n_ki ^= pd_s_j_ki(1,1);
                    n_ki.normalize();

                    //if hed_j_k->opposite

                    if(s_x_ik->CalculatePartialDerivatives(1, 0, v, pd_s_x_ik) == GL_FALSE)
                    {
                        cout<<"Could not calculate partial derivatives"<<endl;
                    }

                    DCoordinate3 n_ik = pd_s_x_ik(1,0);
                    n_ik ^= pd_s_x_ik(1,1);
                    n_ik.normalize();

                    if(he_ki->opposite)
                    {
                        n_ki += n_ik;
                        n_ki.normalize();
                    }

                    DCoordinate3 b_j = (p_ki - p_j) ^ f_j;
                    b_j.normalize();
                    DCoordinate3 t_j = f_j ^ b_j;

                    DCoordinate3 f_ki = - n_ki;
                    DCoordinate3 b_ki = (p_j - p_ki) ^ f_ki;
                    b_ki.normalize();
                    DCoordinate3 t_ki = f_ki ^ b_ki;

                    GLdouble dot_t_j_t_ki = t_j * t_ki;
                    GLdouble delta_j = PHI(1, 1) * PHI(2, 2) - dot_t_j_t_ki * dot_t_j_t_ki * PHI(2, 1) * PHI(2, 1);

                    RealSquareMatrix inverse_j(2);

                    inverse_j(0, 0) = PHI(2, 2);
                    inverse_j(0, 1) = -dot_t_j_t_ki * PHI(2, 1);
                    inverse_j(1, 0) = inverse_j(0, 1);
                    inverse_j(1, 1) = PHI(1, 1);

                    for (GLuint m = 0; m <= 1; m++)
                    {
                        for (GLuint n = 0; n <= 1; n++)
                        {
                            inverse_j(m, n) /= -delta_j;
                        }
                    }

                    ColumnMatrix<GLdouble> B_j(2);

                    B_j[0] = (p_j * (PHI(1, 0) + PHI(1, 1)) + p_ki * (PHI(2, 1) + PHI(3, 1))) * t_j;
                    B_j[1] = (p_j * (PHI(2, 0) + PHI(2, 1)) + p_ki * (PHI(2, 2) + PHI(3, 2))) * t_ki;

                    ColumnMatrix<GLdouble> lambda_j(2);

                    for (GLuint n = 0; n <= 1; n++)
                    {
                        for (GLuint m = 0; m <= 1; m++)
                        {
                            lambda_j[n] += inverse_j(n, m) * B_j[m];
                        }
                    }

                    (*g_j_ki)[0] = p_j;
                    (*g_j_ki)[1] = p_j + lambda_j[0] * t_j;
                    (*g_j_ki)[2] = p_ki + lambda_j[1] * t_ki;
                    (*g_j_ki)[3] = p_ki;

                    /***********************************for p_k and p_ij**********************************/

                    TriangularSurface3::PartialDerivatives pd_s_k_ij, pd_s_x_ji;

                    if(s_k_ij->CalculatePartialDerivatives(1, 0, beta - w, pd_s_k_ij) == GL_FALSE)
                    {
                        cout<<"Could not calculate partial derivatives"<<endl;
                        //return;
                    }


                    LinearCombination3::Derivatives d_c_ij;

                    if(c_ij->CalculateDerivatives(0, w, d_c_ij) == GL_FALSE)
                    {
                        cout<<"Could not calculate derivatives"<<endl;
                    }

                    DCoordinate3 p_ij = d_c_ij(0);//pd_s_k_i_j(0,0);//d_c_i_j(0);
                    DCoordinate3 n_ij = pd_s_k_ij(1,0);
                    n_ij ^= pd_s_k_ij(1,1);
                    n_ij.normalize();

                    //if hed_j_k->opposite
                    if(s_x_ji->CalculatePartialDerivatives(1, 0, w, pd_s_x_ji) == GL_FALSE)
                    {
                        cout<<"Could not calculate partial derivatives"<<endl;
                    }

                    DCoordinate3 n_ji = pd_s_x_ji(1,0);
                    n_ji ^= pd_s_x_ji(1,1);
                    n_ji.normalize();

                    if(he_ij->opposite)
                    {
                        n_ij += n_ji;
                        n_ij.normalize();
                    }

                    DCoordinate3 b_k = (p_ij - p_k) ^ f_k;
                    b_k.normalize();
                    DCoordinate3 t_k = f_k ^ b_k;

                    DCoordinate3 f_ij = - n_ij;
                    DCoordinate3 b_ij = (p_k - p_ij) ^ f_ij;
                    b_ij.normalize();
                    DCoordinate3 t_ij = f_ij ^ b_ij;

                    GLdouble dot_t_k_t_ij = t_k * t_ij;
                    GLdouble delta_k = PHI(1, 1) * PHI(2, 2) - dot_t_k_t_ij * dot_t_k_t_ij * PHI(2, 1) * PHI(2, 1);

                    RealSquareMatrix inverse_k(2);

                    inverse_k(0, 0) = PHI(2, 2);
                    inverse_k(0, 1) = -dot_t_k_t_ij * PHI(2, 1);
                    inverse_k(1, 0) = inverse_k(0, 1);
                    inverse_k(1, 1) = PHI(1, 1);

                    for (GLuint m = 0; m <= 1; m++)
                    {
                        for (GLuint n = 0; n <= 1; n++)
                        {
                            inverse_k(m, n) /= -delta_k;
                        }
                    }

                    ColumnMatrix<GLdouble> B_k(2);

                    B_k[0] = (p_k * (PHI(1, 0) + PHI(1, 1)) + p_ij * (PHI(2, 1) + PHI(3, 1))) * t_k;
                    B_k[1] = (p_k * (PHI(2, 0) + PHI(2, 1)) + p_ij * (PHI(2, 2) + PHI(3, 2))) * t_ij;

                    ColumnMatrix<GLdouble> lambda_k(2);

                    for (GLuint n = 0; n <= 1; n++)
                    {
                        for (GLuint m = 0; m <= 1; m++)
                        {
                            lambda_k[n] += inverse_k(n, m) * B_k[m];
                        }
                    }

                    (*g_k_ij)[0] = p_k;
                    (*g_k_ij)[1] = p_k + lambda_k[0] * t_k;
                    (*g_k_ij)[2] = p_ij + lambda_k[1] * t_ij;
                    (*g_k_ij)[3] = p_ij;

                    /**************************************calculate averaged surface point********************************/
                    LinearCombination3::Derivatives d_g_i_jk, d_g_j_ki, d_g_k_ij;

                    if(g_i_jk->CalculateDerivatives(0, min(max(0.0,beta * (1.0 - b[0])), beta), d_g_i_jk) == GL_FALSE)
                    {
                        cout<<"Could not calculate derivatives, g_i_jk"<<endl;
                    }

                    DCoordinate3 s_i = d_g_i_jk[0];

                    if(g_j_ki->CalculateDerivatives(0, min(max(0.0,beta * (1.0 - b[1])), beta), d_g_j_ki) == GL_FALSE)
                    {
                        cout<<"Could not calculate derivatives, g_j_ki"<<endl;
                    }

                    DCoordinate3 s_j = d_g_j_ki[0];

                    if(g_k_ij->CalculateDerivatives(0, min(max(0.0,beta * (1.0 - b[2])), beta), d_g_k_ij) == GL_FALSE)
                    {
                        cout<<"Could not calculate derivatives, g_k_ij"<<endl;
                    }

                    DCoordinate3 s_k = d_g_k_ij[0];

                    GLdouble denomiator = b[0] * b[1] + b[0] * b[2] + b[1] * b[2];
                    GLdouble w_i = b[1] * b[2] / denomiator;
                    GLdouble w_j = b[0] * b[2] / denomiator;
                    GLdouble w_k = b[0] * b[1] / denomiator;

                    DCoordinate3 p_s_i_j_k = w_i * s_i + w_j * s_j + w_k * s_k;

                    GLuint index[4];

                    index[0] = i * (i + 1) / 2 + j;
                    index[1] = index[0] + i + 1;
                    index[2] = index[1] + 1;
                    index[3] = index[0] + 1;

                    (*result)._vertex[index[0]] = p_s_i_j_k;

                    // texture coordinates
                    (*result)._tex[index[0]].s() = ss;
                    (*result)._tex[index[0]].t() = tt;

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

                    //normals
                    if (j == 0)
                    {
                        (*result)._normal[index[0]] = n_ij;
                    }

                    if (i == div_point_count - 1)
                    {
                        (*result)._normal[index[0]] = n_jk;
                    }

                    if (i == j)
                    {
                        (*result)._normal[index[0]] = n_ki;
                    }

                }
            }
        }



        //normals
        for (GLuint i = 2; i < div_point_count - 1; i++)
        {
            for ( GLuint j = 1; j < i; j++)
            {
                GLuint index[4];

                index[0] = i * (i + 1) / 2 + j;
                index[1] = index[0] + i + 1;
                index[2] = index[1] + 1;
                index[3] = index[0] + 1;

                DCoordinate3 n = (*result)._vertex[index[1]];
                n -= (*result)._vertex[index[0]];

                DCoordinate3 p = (*result)._vertex[index[2]];
                p -= (*result)._vertex[index[0]];

                n ^= p;

                DCoordinate3 n2 = (*result)._vertex[index[2]];
                n2 -= (*result)._vertex[index[0]];

                DCoordinate3 p2 = (*result)._vertex[index[3]];
                p2 -= (*result)._vertex[index[0]];

                n2 ^= p2;

                for (GLint node = 0; node < 3; ++node){
                    (*result)._normal[index[node]] += n;
                    if (node == 0 || node == 2)
                    {
                        (*result)._normal[index[node]] += n2;
                    }
                }

                (*result)._normal[index[3]] += n2;
            }
        }

        for (vector<DCoordinate3>::iterator nit = (*result)._normal.begin(); nit != (*result)._normal.end(); ++nit)
            nit->normalize();

        fit->G1 = result;
        fit->G1->UpdateVertexBufferObjects(usage_flag);

    }

}


GLboolean TriangulatedHalfEdgeDataStructure3::RenderCurves(GLuint order, GLenum render_mode)
{
    for (map<pair<GLint, GLint>, HalfEdge*>::iterator it = _edges.begin(); it != _edges.end(); it++)
    {
        //if (it->second->is_img_boundary_renderable)
        {
            glColor3f(1.0f, 0.0f, 0.0f);
            if (!it->second->img_boundary->RenderDerivatives(order, render_mode))
            {
                cout<<"ERROR! Cannot render derivatives."<<endl;
                return GL_FALSE;
            }
        }
    }

    return GL_TRUE;
}

GLboolean TriangulatedHalfEdgeDataStructure3::RenderControlPolygons(GLenum render_mode)
{

    for (map<pair<GLint, GLint>, HalfEdge*>::iterator it = _edges.begin(); it != _edges.end(); it++)
    {
        glColor3f(0.1f, 0.5f, 0.4f);
        if (!it->second->boundary || !it->second->boundary->RenderData(render_mode))
        {
            cout<<"cannot render boundary"<<endl;
            return GL_FALSE;
        }

        glPointSize(15);
        if (!it->second->boundary || !it->second->boundary->RenderData(GL_POINTS))
        {
            cout<<"cannot render boundary control points"<<endl;
            return GL_FALSE;
        }
    }

    return GL_TRUE;
}

GLboolean TriangulatedHalfEdgeDataStructure3::RenderHalfEdges()
{

    if (_edges.empty())
    {
        return GL_FALSE;
    }

    glBegin(GL_LINES);
        for (map< pair<GLint, GLint>, HalfEdge* >::const_iterator it = _edges.begin(); it != _edges.end(); it++)
        {
            glColor3f(0.0f,0.0f,0.0f);
            //glColor3f(0.4f, 0.5f, 0.6f);
            if (it->second->vertex && it->second->next->vertex)
            {
                glVertex3dv(&(*(it->second->vertex))[0]);
                glVertex3dv(&(*(it->second->next->vertex))[0]);
            }
        }
    glEnd();

    return GL_TRUE;
}

GLboolean TriangulatedHalfEdgeDataStructure3::RenderFaces()
{
    if (_face.empty())
    {
        return GL_FALSE;
    }

    glBegin(GL_TRIANGLES);
        for (vector<Face>::const_iterator fit = _face.begin(); fit != _face.end(); ++fit)
        {
                glVertex3dv(&(_vertex[(*fit)[0]])[0]);
                glVertex3dv(&(_vertex[(*fit)[1]])[0]);
                glVertex3dv(&(_vertex[(*fit)[2]])[0]);

        }
    glEnd();

    return GL_TRUE;
}

GLboolean TriangulatedHalfEdgeDataStructure3::RenderC0Meshes(GLenum render_mode)
{
    for (vector<Face>::const_iterator fit = _face.begin(); fit != _face.end(); ++fit)
    {
        if (!fit->img_C0->Render(render_mode))
        {
            cout<<"ERROR! Could not render image"<<endl;
            return GL_FALSE;
        }

    }
    return GL_TRUE;
}


GLboolean TriangulatedHalfEdgeDataStructure3::RenderG1Meshes(GLenum render_mode)
{
    for (vector<Face>::const_iterator fit = _face.begin(); fit != _face.end(); ++fit)
    {
        if (!fit->G1->Render(render_mode))
        {
            cout<<"ERROR! Could not render image"<<endl;
            return GL_FALSE;
        }
    }
    return GL_TRUE;
}


GLboolean TriangulatedHalfEdgeDataStructure3::RenderC0ControlNets()
{
    GLuint i = 0;
    for (vector<Face>::const_iterator fit = _face.begin(); fit != _face.end(); ++fit, i++)
    {
        if (!fit->C0->UpdateVertexBufferObjectsOfData(GL_STATIC_DRAW))
        {
            cout<<"Error! could not update vertex objects"<<endl;
        }

        glColor3f(0.1f, 0.5f, 0.9f);
        if (!fit->C0->RenderData(GL_LINE_STRIP)) {
            cout<<"ERROR!, Could not render data"<<endl;
            return GL_FALSE;
        }
    }
    return GL_TRUE;
}

GLboolean TriangulatedHalfEdgeDataStructure3::RenderC0BoundaryNormals(GLdouble scale)
{
    if (_face.empty())
    {
        return GL_FALSE;
    }

    for (vector<Face>::const_iterator fit = _face.begin(); fit != _face.end(); fit++)
    {
        if (!fit->img_C0)
        {
            return GL_FALSE;
        }

        GLuint div_point_count = (GLuint)sqrt(fit->img_C0->FaceCount()) + 1;
        for (GLuint i  = 0; i < div_point_count; i++)
        {
            {
                glColor3d(1.0, 0.0, 0.0);
                glBegin(GL_LINES);
                    GLuint index = i * (i + 1) / 2 + 0;
                    glVertex3dv(&fit->img_C0->_vertex[index][0]);
                    DCoordinate3 sum = fit->img_C0->_vertex[index];
                    sum += scale * fit->img_C0->_normal[index];
                    glVertex3dv(&sum[0]);
                glEnd();
            }

            {
                glColor3d(0.1, 0.5, 1.0);
                glBegin(GL_LINES);
                    GLuint index = i * (i + 1) / 2 + i;
                    glVertex3dv(&fit->img_C0->_vertex[index][0]);
                    DCoordinate3 sum = fit->img_C0->_vertex[index];
                    sum += scale * fit->img_C0->_normal[index];
                    glVertex3dv(&sum[0]);
                glEnd();
            }

            {
                glColor3d(0.1, 1.0, 0.5);
                glBegin(GL_LINES);
                    GLuint index = (div_point_count - 1) * div_point_count / 2 + i;
                    glVertex3dv(&fit->img_C0->_vertex[index][0]);
                    DCoordinate3 sum = fit->img_C0->_vertex[index];
                    sum += scale * fit->img_C0->_normal[index];
                    glVertex3dv(&sum[0]);
                glEnd();
            }
        }
    }
    return GL_TRUE;
}

GLboolean TriangulatedHalfEdgeDataStructure3::RenderG1BoundaryNormals(GLdouble scale)
{
    if (_face.empty())
    {
        return GL_FALSE;
    }

    for (vector<Face>::const_iterator fit = _face.begin(); fit != _face.end(); fit++)
    {
        if (!fit->G1)
        {
            return GL_FALSE;
        }

        GLuint div_point_count = (GLuint)sqrt(fit->G1->FaceCount()) + 1;
        for (GLuint i  = 0; i < div_point_count; i++)
        {
            {
                glColor3d(1.0, 0.0, 0.0);
                glBegin(GL_LINES);
                    GLuint index = i * (i + 1) / 2 + 0;
                    glVertex3dv(&fit->G1->_vertex[index][0]);
                    DCoordinate3 sum = fit->G1->_vertex[index];
                    sum += scale * fit->G1->_normal[index];
                    glVertex3dv(&sum[0]);
                glEnd();
            }

            {
                glColor3d(0.1, 0.5, 1.0);
                glBegin(GL_LINES);
                    GLuint index = i * (i + 1) / 2 + i;
                    glVertex3dv(&fit->G1->_vertex[index][0]);
                    DCoordinate3 sum = fit->G1->_vertex[index];
                    sum += scale * fit->G1->_normal[index];
                    glVertex3dv(&sum[0]);
                glEnd();
            }

            {
                glColor3d(0.1, 1.0, 0.5);
                glBegin(GL_LINES);
                    GLuint index = (div_point_count - 1) * div_point_count / 2 + i;
                    glVertex3dv(&fit->G1->_vertex[index][0]);
                    DCoordinate3 sum = fit->G1->_vertex[index];
                    sum += scale * fit->G1->_normal[index];
                    glVertex3dv(&sum[0]);
                glEnd();
            }
        }
    }
    return GL_TRUE;
}

GLboolean TriangulatedHalfEdgeDataStructure3::RenderApproximatedUnitTangents(GLdouble scale)
{
    if(_edges.empty())
    {
        return GL_FALSE;
    }

    glColor3f(0.1f, 0.8f, 0.3f);
    glBegin(GL_LINES);
    for (map< pair<GLint, GLint>, HalfEdge* >::const_iterator it = _edges.begin(); it != _edges.end(); it++)
    {
        glVertex3dv(&(*it->second->vertex)[0]);
        DCoordinate3 sum = *it->second->vertex;
        sum += (it->second->tangent) * scale;
        glVertex3dv(&sum[0]);
    }
    glEnd();

    return GL_TRUE;
}

GLboolean TriangulatedHalfEdgeDataStructure3::RenderVertexNormals(GLdouble scale)
{
    if (_vertex.empty())
    {
        return GL_FALSE;
    }

    glColor3f(0.6f, 0.1f, 0.8f);
    glBegin(GL_LINES);
    for (map< pair<GLint, GLint>, HalfEdge* >::const_iterator it = _edges.begin(); it != _edges.end(); it++)
    {
        glVertex3dv(&(*(it->second->vertex))[0]);
        DCoordinate3 sum = *it->second->vertex;
        sum += (*it->second->normal) * scale;
        glVertex3dv(&sum[0]);
    }
    glEnd();

    return GL_TRUE;
}

GLboolean TriangulatedHalfEdgeDataStructure3::RenderVertices()
{
    if (_vertex.empty())
    {
        return GL_FALSE;
    }

    glPointSize(20);
    glColor3f(0.5f, 0.6f, 0.5f);
    glBegin(GL_POINTS);
    for (vector<DCoordinate3>::const_iterator vit = _vertex.begin(); vit != _vertex.end(); vit++)
    {
        glVertex3f((*vit).x(), (*vit).y(), (*vit).z());
    }
    glEnd();

    return GL_TRUE;
}

GLvoid TriangulatedHalfEdgeDataStructure3::DeleteBoundaryCurves()
{
    for(map< pair<GLint, GLint>, HalfEdge* >::const_iterator it = _edges.begin(); it != _edges.end(); it++)
    {
        if (it->second->boundary)
        {
            delete it->second->boundary; it->second->boundary = nullptr;
        }

        if (it->second->img_boundary)
        {
            delete it->second->img_boundary; it->second->img_boundary = nullptr;
        }
    }

}

GLvoid TriangulatedHalfEdgeDataStructure3::DeleteC0Surfaces()
{
    for (vector<Face>::iterator fit = _face.begin(); fit != _face.end(); ++fit)
    {
        if (fit->C0)
        {
            delete fit->C0; fit->C0 = nullptr;
        }

        if (fit->img_C0)
        {
            delete fit->img_C0; fit->img_C0 = nullptr;
        }
    }
}

GLvoid TriangulatedHalfEdgeDataStructure3::DeleteG1Surfaces()
{
    for (vector<Face>::iterator fit = _face.begin(); fit != _face.end(); ++fit)
    {
        if (fit->G1)
        {
            delete fit->G1; fit->G1 = nullptr;
        }
    }
}

GLvoid TriangulatedHalfEdgeDataStructure3::CleanAll()
{
    DeleteG1Surfaces();
    DeleteC0Surfaces();
    DeleteBoundaryCurves();
    _vertex.clear();
    _normal.clear();
    _face.clear();
    _edges.clear();
}

TriangulatedHalfEdgeDataStructure3::~TriangulatedHalfEdgeDataStructure3()
{
    _vertex.clear();
    _normal.clear();
    _face.clear();
    _edges.clear();
}
