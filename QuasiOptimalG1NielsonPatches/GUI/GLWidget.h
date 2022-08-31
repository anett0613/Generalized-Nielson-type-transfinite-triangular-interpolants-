#pragma once

#include <GL/glew.h>
#include <QGLWidget>
#include <QGLFormat>
#include "../Core/DCoordinates3.h"
//#include "../Parametric/ParametricCurves3.h"
#include "../Core/Materials.h"
#include "../Core/TriangulatedMeshes3.h"
//#include "../Parametric/ParametricSurfaces3.h"
#include "../Core/ShaderPrograms.h"
//#include "Test/TestFunctions.h"
#include <QTimer>
//#include "../Core/ToroidalCompositeBSplinePatch3.h"
//#include "../BSpline/SpecializedBSplineCurves3.h"

//#include "../Cyclic/CyclicCurves3.h"
#include "../QuasiOptimalNielsonPatches/TriangulatedHalfEdgeDataStructures.h"

namespace cagd
{
    class GLWidget: public QGLWidget
    {
        Q_OBJECT

    private:

        // variables defining the projection matrix
        double       _aspect;            // aspect ratio of the rendering window
        double       _fovy;              // field of view in direction y
        double       _z_near, _z_far;    // distance of near and far clipping planes

        // variables defining the model-view matrix
        double       _eye[3], _center[3], _up[3];

        // variables needed by transformations
        int         _angle_x, _angle_y, _angle_z;
        double      _zoom;
        double      _trans_x, _trans_y, _trans_z;

        // your other declarations

        int _type_index;
		int _combo_index;
		int _color_index;
        int _normals_surface;
        double _first_derivative;
        double _second_derivative;
        double _scale_value=1.0;

        GLuint _ps_selected;
        GLuint _pc_count;
		RowMatrix<GLdouble> _scale;

        GLuint _ps_count;

		GLfloat _angle = 0;
		QTimer *_timer;

        GLuint _n;

        //ShaderProgram _shader;
        RowMatrix<ShaderProgram>                _shader;
        int                                     _shader_index;
        bool                                    _shader_enabled;
        double                                  _r, _g, _b, _a;

        GLfloat _shading_parameter;
        GLfloat _scale_factor_parameter;
        GLfloat _smoothing_parameter;

        //Nielson-type
        int                                                 _nielson_index = 0;
        int                                                 _model_index = 0;
        GLdouble                                            _beta = 1.0;
        bool                                                _beta_enabled = false;
        TriangulatedHalfEdgeDataStructure3                  _heds; //elephant
        TriangulatedHalfEdgeDataStructure3::BoundaryType    _boundary_type = TriangulatedHalfEdgeDataStructure3::CUBIC_POLYNOMIAL;
        int                                                 _div_point_count = 20;
        int                                                 _order = 0;
        GLdouble                                            _theta1 = 1.0;
        GLdouble                                            _theta2 = 0.0;
        GLdouble                                            _theta3 = 0.0;
        GLdouble                                            _epsilon1 = 1.0;
        GLdouble                                            _epsilon2 = 0.0;
        GLdouble                                            _epsilon3 = 0.0;
        GLdouble                                            _scale_parameter = 1.0;

        bool                                                _show_control_net = false;
        bool                                                _show_control_polygons = false;
        bool                                                _show_boundaries = false;
        bool                                                _show_halfEdges = false;
        bool                                                _show_vertices = false;
        bool                                                _show_vertex_normals = false;
        bool                                                _show_averaged_normals = false;
        bool                                                _show_C0_normals = false;
        bool                                                _show_tangents = false;

        bool                                                _dark_mode = false;


    public:
        // special and default constructor
        // the format specifies the properties of the rendering window
        GLWidget(QWidget* parent = 0, const QGLFormat& format = QGL::Rgba | QGL::DepthBuffer | QGL::DoubleBuffer);

        // redeclared virtual functions
        GLvoid initializeGL();
        void paintGL();
        void resizeGL(int w, int h);
        virtual ~GLWidget();

        void switch_parametric_surface(GLuint index);
        void switch_parametric_curve(GLuint index);

	private slots:
        //void animate();

    public slots:
        // public event handling methods/slots
        void set_angle_x(int value);
        void set_angle_y(int value);
        void set_angle_z(int value);

        void set_zoom_factor(double value);

        void set_trans_x(double value);
        void set_trans_y(double value);
        void set_trans_z(double value);

		void set_combo_index(int index);
        void set_type_index(int index);

		void set_color(int index);
		void set_color_index(int index);
        void set_scale(double value);

        void setShaderNumber(int);
        void setShading(double);

        void set_shader(int index);

        void set_shader_enabled();
        void set_shader_scale(double value);
        void set_shader_smoothing(double value);
        void set_shader_shading(double value);

        void set_shader_color();
        void set_shader_red(double value);
        void set_shader_green(double value);
        void set_shader_blue(double value);
        void set_shader_alpha(double value);

        //Nielson type triangular surfaces
        void show_Nielson_type();
        void show_faces();
        void show_C0_surfaces();
        void show_G1_surfaces();
        void show_boundary_curves();
        void show_halfEdges();

        void show_control_net();
        void show_C0_normals();
        void show_control_polygons();
        void show_G1_normals();
        void show_vertex_normals();
        void show_vertices();
        void show_boundaries();
        void render_halfEdges();
        void show_tangents();

        void set_boundary_type(int value);
        void set_div_point_count(int value);
        void set_beta(double value);
        void set_order(int value);
        void set_theta1(double value);
        void set_theta2(double value);
        void set_theta3(double value);
        void set_epsilon1(double value);
        void set_epsilon2(double value);
        void set_epsilon3(double value);
        void set_Nielson_index(int value);
        void set_scale_parameter(double value);
        void set_model_index(int value);

        void set_dark_mode();

    signals:
            // B-Spline curve signals
        void x_changed(double);
        void y_changed(double);
        void z_changed(double);

        //Nielson-type surfaces signals
        void boundary_type_changed(bool);

	private:

        GLvoid initShader();
        GLvoid paintPoint(DCoordinate3);
        GLvoid init_Nielson_type_surfaces();
    };
}
