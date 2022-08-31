#include "Bicu_bic_B_Spline_Patch3.h"

#include <Core/Materials.h>
#include "ShaderPrograms.h"
#include "DCoordinates3.h"
#include "Matrices.h"
//#include <algorithm>
#include <math.h>

namespace cagd
{
    class BicubicBSplineCompositeSurface3
    {
    public:
        enum Direction{N, NW, W, SW, S, SE, E, NE};

        class PatchAttributes
        {
        public:
            BicubicBSplinePatch3    *patch;
            TriangulatedMesh3       *image;
            Material                *material;
            ShaderProgram           *shader;
            RowMatrix<GenericCurve3*>   *u_isoparametric_lines;
            RowMatrix<GenericCurve3*>   *v_isoparametric_lines;

            RowMatrix<PatchAttributes*> neighbours;

            PatchAttributes(Material &mat=MatFBBrass);
            PatchAttributes(const PatchAttributes &pa);

            PatchAttributes& operator = (const PatchAttributes &pa);

            BicubicBSplinePatch3* getPatch() const;
            GLvoid setPatch(BicubicBSplinePatch3* p);

            TriangulatedMesh3* getImage()const;
            GLvoid setImage(TriangulatedMesh3* img);

            RowMatrix<GenericCurve3*>* getUIsoparametricLines() const;
            RowMatrix<GenericCurve3*>* getVIsoparametricLines() const;

            //GLboolean setNeighbour(const PatchAttributes &pa);

            ~PatchAttributes();
            // to do: constructor, copy constructor, assignment operator, destructor
                    //        setters/getters

       };

    protected:
        //Matrix<DCoordinate3>            _big_control_net; // 9es
        //GLdouble                        _u_alpha, _v_alpha;
        //std::vector<PatchAttributes>    _attributes;
        Matrix<DCoordinate3>       _control_net; // ebből kell kinyerni az összes 4×4-es kontrollpont-blokkot
        Matrix<PatchAttributes>    _attributes;  //feluleti foltok
        GLuint _row_count, _column_count; //feluleti foltokhoz

        GLdouble                        _scale_iso_lines;
        GLboolean                       _show_u_iso_lines, _show_v_iso_lines;


    public:
        BicubicBSplineCompositeSurface3(GLuint row_count, GLuint column_count); // a _control_net méretét állítja
        DCoordinate3  operator ()(GLuint row, GLuint column) const; // kontrollpont érték szerinti lekérdezése
        DCoordinate3& operator ()( GLuint row, GLuint column);      // kontrollpont referencia szerinti lekérdezése

        GLvoid render();  //osszetett folt kirajzolasa
        GLvoid renderControlPoints();
        GLvoid renderIsoLines();

        GLvoid setScaleIsoLines(double value);
        GLvoid setShowUIsoLines();
        GLvoid setShowVIsoLines();
        GLvoid updateIsoLines();
        GLuint getRowCount();

        GLvoid SetSelectedPointU(GLint u);
        GLvoid SetSelectedPointV(GLint v);
        GLint GetSelectedPointU();
        GLint GetSelectedPointV();


        GLvoid updateControllPoint(int row, int col, DCoordinate3 point);
        GLvoid setControlPoint(int row, int col, DCoordinate3 point);

        GLvoid UpdateCorner(GLint i_begin, GLint i_end, GLint j_begin, GLint j_end, GLint row, GLint column);

        GLvoid generateImage();

    };
};
