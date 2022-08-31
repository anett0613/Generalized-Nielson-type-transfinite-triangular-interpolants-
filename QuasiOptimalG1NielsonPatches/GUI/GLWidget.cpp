#include "GLWidget.h"

#if !defined(__APPLE__)
#include <GL/glu.h>
#endif

#include <iostream>
#include "../Core/Matrices.h"
#include "../Test/TestFunctions.h"
#include "../Core/Constants.h"
#include <omp.h>
using namespace std;
using namespace cagd;


#include <Core/Exceptions.h>

namespace cagd
{
    //--------------------------------
    // special and default constructor
    //--------------------------------
    GLWidget::GLWidget(QWidget *parent, const QGLFormat &format): QGLWidget(format, parent)
    {
		_timer = new QTimer(this);
		_timer->setInterval(0);
        //connect(_timer, SIGNAL(timeout()), this, SLOT(animate()));
    }

    //--------------------------------------------------------------------------------------
    // this virtual function is called once before the first call to paintGL() or resizeGL()
    //--------------------------------------------------------------------------------------

    GLWidget::~GLWidget(){

    }

    void GLWidget::initializeGL()
    {
        // creating a perspective projection matrix
        glMatrixMode(GL_PROJECTION);

        glLoadIdentity();

        _aspect = (double)width() / (double)height();
        _z_near = 1.0;
        _z_far  = 1000.0;
        _fovy   = 45.0;

        gluPerspective(_fovy, _aspect, _z_near, _z_far);

        // setting the model view matrix
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        _eye[0] = _eye[1] = 0.0; _eye[2] = 6.0;
        _center[0] = _center[1] = _center[2] = 0.0;
        _up[0] = _up[2] = 0.0; _up[1] = 1.0;

        gluLookAt(_eye[0], _eye[1], _eye[2], _center[0], _center[1], _center[2], _up[0], _up[1], _up[2]);

        // enabling the depth test
        glEnable(GL_DEPTH_TEST);

        // setting the background color
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        //set_dark_mode();

        // initial values of transformation parameters
        _angle_x = _angle_y = _angle_z = 0.0;
        _trans_x = _trans_y = _trans_z = 0.0;
        _zoom = 1.0;

        // ...

        try
        {
            // initializing the OpenGL Extension Wrangler library
            GLenum error = glewInit();

            if (error != GLEW_OK)
            {
                throw Exception("Could not initialize the OpenGL Extension Wrangler Library!");
            }

            if (!glewIsSupported("GL_VERSION_2_0"))
            {
                throw Exception("Your graphics card is not compatible with OpenGL 2.0+! "
                                "Try to update your driver or buy a new graphics adapter!");
            }

            // create and store your geometry in display lists or vertex buffer objects
            // ...
            /*if (!_heds.LoadFromOFF("Models/elephant_v162_f320.off", GL_TRUE))
            {
                throw Exception("Could not load the model file!");
            }
            */

            //_heds.GenerateBoundaryCurves(TriangulatedHalfEdgeDataStructure3::SECOND_ORDER_TRIGONOMETRIC, _beta, 1, 0, 0, 10, GL_STATIC_DRAW);
            //_heds.GenerateC0TriangularSurfaces(TriangulatedHalfEdgeDataStructure3::SECOND_ORDER_TRIGONOMETRIC, _beta, 1, 0, 0, 10, GL_STATIC_DRAW);
            //_heds.GenerateG1TriangularSurfaces(TriangulatedHalfEdgeDataStructure3::CUBIC_POLYNOMIAL, _beta, 1, 0, 0, 10, GL_STATIC_DRAW);
            init_Nielson_type_surfaces();
            initShader();
        }
        catch (Exception &e)
        {
            cout << e << endl;
        }
        glEnable(GL_POINT_SMOOTH );
        glHint(GL_POINT_SMOOTH_HINT,GL_NICEST);
        glEnable(GL_LINE_SMOOTH );
        glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
        glEnable(GL_POLYGON_SMOOTH );
        glHint(GL_POLYGON_SMOOTH_HINT,GL_NICEST);

        glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
        glEnable(GL_DEPTH_TEST);

        glewInit();

        _type_index = 0;
		_combo_index = 0;
        _pc_count = 7;


    }

    //-----------------------
    // the rendering function
    //-----------------------
    void GLWidget::paintGL()
    {
        // clears the color and depth buffers
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // stores/duplicates the original model view matrix
        glPushMatrix();

            // applying transformations
            glRotatef(_angle_x, 1.0, 0.0, 0.0);
            glRotatef(_angle_y, 0.0, 1.0, 0.0);
            glRotatef(_angle_z, 0.0, 0.0, 1.0);
            glTranslated(_trans_x, _trans_y, _trans_z);
            glScaled(_zoom, _zoom, _zoom);

            // render your geometry (this is oldest OpenGL rendering technique, later we will use some advanced methods)

          /*  glColor3f(1.0f, 1.0f, 1.0f);
            glBegin(GL_LINES);
                glVertex3f(0.0f, 0.0f, 0.0f);
                glVertex3f(1.1f, 0.0f, 0.0f);

                glVertex3f(0.0f, 0.0f, 0.0f);
                glVertex3f(0.0f, 1.1f, 0.0f);

                glVertex3f(0.0f, 0.0f, 0.0f);
                glVertex3f(0.0f, 0.0f, 1.1f);
            glEnd();

            glBegin(GL_TRIANGLES);
                // attributes
                glColor3f(1.0f, 0.0f, 0.0f);
                // associated with position
                glVertex3f(1.0f, 0.0f, 0.0f);

                // attributes
                glColor3f(0.0, 1.0, 0.0);
                // associated with position
                glVertex3f(0.0, 1.0, 0.0);

                // attributes
                glColor3f(0.0f, 0.0f, 1.0f);
                // associated with position
                glVertex3f(0.0f, 0.0f, 1.0f);
            glEnd();
            if(_image_of_pc){
                glColor3f(1.0f,0.0f,0.0f);
                _image_of_pc->RenderDerivatives(0,GL_LINE_STRIP);

                glPointSize(5.0f);
                    glColor3f(0.0f, 0.5f, 0.0f);
                    _image_of_pc->RenderDerivatives(1, GL_LINES);
                    _image_of_pc->RenderDerivatives(1,GL_POINTS);

                    glColor3f(0.1f, 0.5f, 0.9f);
                    _image_of_pc->RenderDerivatives(2, GL_LINES);
                    _image_of_pc->RenderDerivatives(2,GL_POINTS);
                glPointSize(1.0f);
            }*/
            /*_shader[0].Enable(GL_TRUE);
            glEnable(GL_LIGHTING);
            glEnable(GL_LIGHT0);
            glEnable(GL_NORMALIZE);

            MatFBBrass.Apply();
            _heds.RenderC0Meshes();
            //_heds.RenderG1Meshes();
            _shader[0].Disable();
            glDisable(GL_LIGHT0);
            glDisable(GL_LIGHTING);
            glDisable(GL_NORMALIZE);*/


            //_heds.RenderCurves(0, GL_LINE_STRIP);
            //_heds.RenderCurves(1, GL_LINES);
            //_heds.RenderCurves(2, GL_LINES);

            //MatFBBrass.Apply();

            //_heds.RenderFaces();




            //glColor3f(0.0f, 1.0f, 0.0f);
            //_heds.RenderVertexNormals(0.5);
            //_heds.RenderC0ControlNets();
            //_heds.RenderControlPolygons();


            //glColor3f(0.0f, 0.0f, 0.0f);
            //_heds.RenderHalfEdges();
            //_heds.RenderC0BoundaryNormals(0.1);
            //_heds.RenderG1BoundaryNormals(0.1);
            //

            show_Nielson_type();

            //MatFBBrass.Apply();
            //_heds.RenderC0Meshes();
        // pops the current matrix stack, replacing the current matrix with the one below it on the stack,
        // i.e., the original model view matrix is restored
        glPopMatrix();
    }

    //----------------------------------------------------------------------------
    // when the main window is resized one needs to redefine the projection matrix
    //----------------------------------------------------------------------------
    void GLWidget::resizeGL(int w, int h)
    {
        // setting the new size of the rendering context
        glViewport(0, 0, w, h);

        // redefining the projection matrix
        glMatrixMode(GL_PROJECTION);

        glLoadIdentity();

        _aspect = (double)w / (double)h;

        gluPerspective(_fovy, _aspect, _z_near, _z_far);

        // switching back to the model view matrix
        glMatrixMode(GL_MODELVIEW);

        updateGL();
    }

    //-----------------------------------
    // implementation of the public slots
    //-----------------------------------

    void GLWidget::set_angle_x(int value)
    {
        if (_angle_x != value)
        {
            _angle_x = value;
            updateGL();
        }
    }

    void GLWidget::set_angle_y(int value)
    {
        if (_angle_y != value)
        {
            _angle_y = value;
            updateGL();
        }
    }

    void GLWidget::set_angle_z(int value)
    {
        if (_angle_z != value)
        {
            _angle_z = value;
            updateGL();
        }
    }

    void GLWidget::set_zoom_factor(double value)
    {
        if (_zoom != value)
        {
            _zoom = value;
            updateGL();
        }
    }

    void GLWidget::set_trans_x(double value)
    {
        if (_trans_x != value)
        {
            _trans_x = value;
            updateGL();
        }
    }

    void GLWidget::set_trans_y(double value)
    {
        if (_trans_y != value)
        {
            _trans_y = value;
            updateGL();
        }
    }

    void GLWidget::set_trans_z(double value)
    {
        if (_trans_z != value)
        {
            _trans_z = value;
            updateGL();
        }
    }

    void GLWidget::set_color(int index)
    {
		switch (index) {
		case 0: {
			MatFBBrass.Apply();
			break;
		}
		case 1: {
			MatFBEmerald.Apply();
			break;
		}
		case 2: {
			MatFBGold.Apply();
			break;
		}
		case 3: {
			MatFBPearl.Apply();
			break;
		}
		case 4: {
			MatFBRuby.Apply();
			break;
		}
		case 5: {
			MatFBSilver.Apply();
			break;
		}
		case 6: {
			MatFBTurquoise.Apply();
			break;
		}
		}
	}

    /*void GLWidget::animate() {
		GLfloat *vertex = _mouse.MapVertexBuffer(GL_READ_WRITE);
		GLfloat *normal = _mouse.MapNormalBuffer(GL_READ_ONLY);

		_angle += DEG_TO_RADIAN;
		if (_angle >= TWO_PI) {
			_angle -= TWO_PI;
		}
		GLfloat scale = sin(_angle) / 3000.0;
		for (GLuint i = 0; i < _mouse.VertexCount(); ++i) {

			for (GLuint coordinate = 0; coordinate < 3; ++coordinate, ++vertex, ++normal)
				*vertex += scale * (*normal);
		}
		_mouse.UnmapVertexBuffer();
		_mouse.UnmapNormalBuffer();
		updateGL();
    }*/

	void GLWidget::set_type_index(int index){
		_type_index = index;
		_combo_index = 0;
		updateGL();
	}
	void GLWidget::set_combo_index(int index){
		if (index > -1) {
			_combo_index = index;
			updateGL();
		}
	}
	void GLWidget::set_color_index(int index) {
		if (index > -1) {
			_color_index = index;
			updateGL();
		}
	}

    void GLWidget::set_shader_enabled()
    {
        _shader_enabled = !_shader_enabled;
        updateGL();
    }

    void GLWidget::set_shader_scale(double value)
    {
        _shader[1].Enable();
        _shader[1].SetUniformVariable1f("scale_factor", value);
        _shader[1].Disable();
        updateGL();
    }

    void GLWidget::set_shader_smoothing(double value)
    {
        _shader[1].Enable();
        _shader[1].SetUniformVariable1f("smoothing", value);
        _shader[1].Disable();
        updateGL();
    }

    void GLWidget::set_shader_shading(double value)
    {
        _shader[1].Enable();
        _shader[1].SetUniformVariable1f("shading", value);
        _shader[1].Disable();
        updateGL();
    }

    void GLWidget::set_shader_red(double value)
    {
        _r = value;
        set_shader_color();
        updateGL();
    }

    void GLWidget::set_shader_green(double value)
    {
        _g = value;
        set_shader_color();
        updateGL();
    }

    void GLWidget::set_shader_blue(double value)
    {
        _b = value;
        set_shader_color();
        updateGL();
    }

    void GLWidget::set_shader_alpha(double value)
    {
        _a = value;
        set_shader_color();
        updateGL();
    }

    void GLWidget::set_shader_color()
    {
        _shader[2].Enable();
        _shader[2].SetUniformVariable4f("default_outline_color", _r, _g, _b, _a);
        _shader[2].Disable();
    }


    void GLWidget::set_scale(double value){
        if(_scale_value!=value){
            _scale_value=value;
            updateGL();
        }

    }

    void GLWidget::initShader(){
      _shader_index = 0;
      _shader.ResizeColumns(4);

      try {
        if (!_shader[0].InstallShaders("../Shaders/directional_light.vert",
                                        "../Shaders/directional_light.frag")) {
          throw Exception("Directional light failed to load.");
        }

        if (!_shader[1].InstallShaders("../Shaders/reflection_lines.vert",
                                        "../Shaders/reflection_lines.frag")) {
          throw Exception("Reflection lines failed to load.");
        } else {
          _shader[1].Enable();
          _shader[1].SetUniformVariable1f("scale_factor", 4.0f);
          _shader[1].SetUniformVariable1f("smoothing", 2.0f);
          _shader[1].SetUniformVariable1f("shading", 1.0f);
          _shader[1].Disable();
        }

        if (!_shader[2].InstallShaders("../Shaders/toon.vert",
                                        "../Shaders/toon.frag")) {
          throw Exception("Toon failed to load.");
        }
        {
          _shader[2].Enable();
          _shader[2].SetUniformVariable4f(
            "default_outline_color", 1.0f, 1.0f, 1.0f, 1.0f);
          _shader[2].Disable();
        }

        if (!_shader[3].InstallShaders("../Shaders/two_sided_lighting.vert",
                                        "../Shaders/two_sided_lighting.frag")) {
          throw Exception("Two sided lighting failed to load.");
        }
      } catch (Exception& e) {
        cerr << e << endl;
      }
    }

    void GLWidget::setShaderNumber(int index){
      _shader_index = index;
      updateGL();
    }

    void GLWidget::setShading(double shading){
      this->_shading_parameter = (float)shading;
      updateGL();
    }

    void GLWidget::set_shader(int index)
    {
        if (_shader_index != index)
        {
            _shader_index = index;
            updateGL();
        }
    }

    GLvoid GLWidget::paintPoint(DCoordinate3 point)
    {
      glPointSize(12);
      glBegin(GL_POINTS);
      glColor3f(0, 1, 1);
      glVertex3f(point.x(), point.y(), point.z());
      glEnd();
      glPointSize(2);
    }

    //Nielson-type triangular surfaces

    GLvoid GLWidget::init_Nielson_type_surfaces()
    {

        switch(_model_index)
        {
        case 0:
            if(!_heds.LoadFromOFF("Models/cube.off", GL_TRUE))
            {
                throw Exception("Could not load the model file!");
            }
            break;
        case 1:
            if(!_heds.LoadFromOFF("Models/icosahedron.off", GL_TRUE))
            {
                throw Exception("Could not load the model file!");
            }
            break;
        case 2:
            if(!_heds.LoadFromOFF("Models/elephant_v162_f320.off", GL_TRUE))
            {
                throw Exception("Could not load the model file!");
            }
            break;
        case 3:
            if(!_heds.LoadFromOFF("Models/cone.off", GL_TRUE))
            {
                throw Exception("Could not load the model file!");
            }
            break;
        /*case 4:
            if(!_heds.LoadFromOFF("Models/bunny.off", GL_TRUE))
            {
                throw Exception("Could not load the model file!");
            }
            break;*/
        case 4:
            if(!_heds.LoadFromOFF("Models/star.off", GL_TRUE))
            {
                throw Exception("Could not load the model file!");
            }
            break;
        /*case 6:
            if(!_heds.LoadFromOFF("Models/sphere.off", GL_TRUE))
            {
                throw Exception("Could not load the model file!");
            }
            break;
        case 7:
            if(!_heds.LoadFromOFF("Models/goblet.off", GL_TRUE))
            {
                throw Exception("Could not load the model file!");
            }
            break;*/
        default:
            if(!_heds.LoadFromOFF("Models/elephant_v162_f320.off", GL_TRUE))
            {
                throw Exception("Could not load the model file!");
            }
            break;
        }

        GLdouble time1 = omp_get_wtime();
        _heds.GenerateBoundaryCurves(_boundary_type, _beta, _theta1, _theta2, _theta3, _div_point_count, GL_STATIC_DRAW);
        cout<<"GenerateBoundary "<< omp_get_wtime() - time1<<endl;
        GLdouble time2 = omp_get_wtime();
        _heds.GenerateC0TriangularSurfaces(_boundary_type, _beta, _epsilon1, _epsilon2, _epsilon3, _div_point_count, GL_STATIC_DRAW);
        cout<<"GenerateC0 "<< omp_get_wtime() - time2<<endl;
        GLdouble time3 = omp_get_wtime();
        _heds.GenerateG1TriangularSurfaces(_boundary_type, _beta, _theta1, _theta2, _theta3, _div_point_count, GL_STATIC_DRAW);

        time3 = omp_get_wtime() - time3;
        cout<<time3<<endl;
    }

    void GLWidget::show_Nielson_type()
    {
        switch(_nielson_index)
        {
        case 0:
            show_halfEdges();
            break;
        case 1:
            show_boundary_curves();
            break;
        case 2:
            glEnable(GL_LIGHTING);
            glEnable(GL_NORMALIZE);
            glEnable(GL_LIGHT0);
            _shader[_shader_index].Enable(GL_TRUE);
            set_color(_color_index);

            show_faces();

            _shader[_shader_index].Disable();
            glDisable(GL_LIGHT0);
            glDisable(GL_LIGHTING);
            glDisable(GL_NORMALIZE);

            break;
        case 3:
            glEnable(GL_LIGHTING);
            glEnable(GL_NORMALIZE);
            glEnable(GL_LIGHT0);
            _shader[_shader_index].Enable(GL_TRUE);
            set_color(_color_index);

            show_C0_surfaces();

            _shader[_shader_index].Disable();
            glDisable(GL_LIGHT0);
            glDisable(GL_LIGHTING);
            glDisable(GL_NORMALIZE);

            break;
        case 4:
            glEnable(GL_LIGHTING);
            glEnable(GL_NORMALIZE);
            glEnable(GL_LIGHT0);
            _shader[_shader_index].Enable(GL_TRUE);
            set_color(_color_index);

            show_G1_surfaces();

            _shader[_shader_index].Disable();
            glDisable(GL_LIGHT0);
            glDisable(GL_LIGHTING);
            glDisable(GL_NORMALIZE);

            break;
        }

        if(_show_control_net)
        {
            _heds.RenderC0ControlNets();
        }

        if(_show_control_polygons)
        {
            _heds.RenderControlPolygons();
        }

        if(_show_vertex_normals)
        {
            _heds.RenderVertexNormals(0.1);
        }

        if(_show_vertices)
        {
            _heds.RenderVertices();
        }

        if(_show_halfEdges)
        {
            _heds.RenderHalfEdges();
        }

        if(_show_boundaries)
        {
            _heds.RenderCurves(_order);
        }

        if(_show_C0_normals)
        {
            _heds.RenderC0BoundaryNormals(_scale_parameter);
        }

        if(_show_averaged_normals)
        {
            _heds.RenderG1BoundaryNormals(_scale_parameter);
        }

        if(_show_tangents)
        {
            _heds.RenderApproximatedUnitTangents(_scale_parameter);
        }
    }

    void GLWidget::show_control_net()
    {
        //_heds.RenderC0ControlNets();
        _show_control_net = !_show_control_net;
        updateGL();
    }

    void GLWidget::show_C0_normals()
    {
        //_heds.RenderC0BoundaryNormals();
        _show_C0_normals = !_show_C0_normals;
        updateGL();
    }

    void GLWidget::show_C0_surfaces()
    {
        _shader[_shader_index].Enable(GL_TRUE);
        _heds.RenderC0Meshes();
        _shader[_shader_index].Disable();
        glDisable(GL_LIGHT0);
        glDisable(GL_LIGHTING);
        glDisable(GL_NORMALIZE);
        //updateGL();
    }

    void GLWidget::show_G1_surfaces()
    {
        _shader[_shader_index].Enable(GL_TRUE);
        _heds.RenderG1Meshes();
        _shader[_shader_index].Disable();
        glDisable(GL_LIGHT0);
        glDisable(GL_LIGHTING);
        glDisable(GL_NORMALIZE);
        //updateGL();
    }

    void GLWidget::show_halfEdges()
    {
        _heds.RenderHalfEdges();
        //updateGL();
    }

    void GLWidget::show_boundary_curves()
    {
        if (!_order)
        {
            _heds.RenderCurves(_order);
        }
        else
        {
            _heds.RenderCurves(_order, GL_LINES);
            _heds.RenderCurves(_order, GL_POINTS);

        }
        //updateGL();
    }

    void GLWidget::show_faces()
    {
        _shader[_shader_index].Enable(GL_TRUE);

        _heds.RenderFaces();

        _shader[_shader_index].Disable();
        glDisable(GL_LIGHT0);
        glDisable(GL_LIGHTING);
        glDisable(GL_NORMALIZE);
        //updateGL();
    }

    void GLWidget::show_G1_normals()
    {
        //_heds.RenderG1BoundaryNormals(_scale_parameter);
        _show_averaged_normals = !_show_averaged_normals;
        updateGL();
    }

    void GLWidget::show_control_polygons()
    {
        _show_control_polygons = !_show_control_polygons;
        updateGL();
    }

    void GLWidget::show_vertices()
    {
        _show_vertices = !_show_vertices;
        updateGL();
    }

    void GLWidget::show_vertex_normals()
    {
        _show_vertex_normals = !_show_vertex_normals;
        updateGL();
    }

    void GLWidget::show_boundaries()
    {
        _show_boundaries = !_show_boundaries;
        updateGL();
    }

    void GLWidget::render_halfEdges()
    {
        _show_halfEdges = !_show_halfEdges;
        updateGL();
    }

    void GLWidget::show_tangents()
    {
        _show_tangents = !_show_tangents;
        updateGL();
    }

    void GLWidget::set_boundary_type(int boundary_type)
    {
        switch(boundary_type)
        {
        case 0:
            _boundary_type = TriangulatedHalfEdgeDataStructure3::CUBIC_POLYNOMIAL;

            if (_beta_enabled)
            {
                _beta_enabled = !_beta_enabled;
                emit boundary_type_changed(_beta_enabled);
                //cout<<"boundary type changed"<<endl;
            }
            break;
        case 1:
            _boundary_type = TriangulatedHalfEdgeDataStructure3::QUARTIC_POLYNOMIAL;

            if (_beta_enabled)
            {
                _beta_enabled = !_beta_enabled;
                emit boundary_type_changed(_beta_enabled);
            }
            break;
        case 2:
            _boundary_type = TriangulatedHalfEdgeDataStructure3::SECOND_ORDER_TRIGONOMETRIC;

            if (!_beta_enabled)
            {
                _beta_enabled = !_beta_enabled;
                emit boundary_type_changed(_beta_enabled);
                //cout<<"boundary type changed"<<endl;
            }
            break;
        case 3:
            _boundary_type = TriangulatedHalfEdgeDataStructure3::SECOND_ORDER_HYPERBOLIC;

            if (!_beta_enabled)
            {
                _beta_enabled = !_beta_enabled;
                emit boundary_type_changed(_beta_enabled);
            }
            break;
        case 4:
            _boundary_type = TriangulatedHalfEdgeDataStructure3::FIRST_ORDER_ALGEBRAIC_TRIGONOMETRIC;

            if (!_beta_enabled)
            {
                _beta_enabled = !_beta_enabled;
                emit boundary_type_changed(_beta_enabled);
            }
            break;
        }

        init_Nielson_type_surfaces();
        //show_Nielson_type(_nielson_index);
        updateGL();
    }

    void GLWidget::set_beta(double value)
    {
        if (_beta != value)
        {
            _beta = value;
            init_Nielson_type_surfaces();
            updateGL();
        }        
    }

    void GLWidget::set_div_point_count(int value)
    {
        if (_div_point_count != value)
        {
            _div_point_count = value;
            updateGL();
        }
    }

    void GLWidget::set_order(int value)
    {
        if (_order != value)
        {
            _order = value;
            updateGL();
        }        
    }

    void GLWidget::set_theta1(double value)
    {
        if (_theta1 != value)
        {
            _theta1 = value;
            init_Nielson_type_surfaces();
            updateGL();
        }        
    }

    void GLWidget::set_theta2(double value)
    {
        if (_theta2 != value)
        {
            _theta2 = value;
            init_Nielson_type_surfaces();
            updateGL();
        }        
    }

    void GLWidget::set_theta3(double value)
    {
        if (_theta3 != value)
        {
            _theta3 = value;
            init_Nielson_type_surfaces();
            updateGL();
        }        
    }

    void GLWidget::set_epsilon1(double value)
    {
        if (_epsilon1 != value)
        {
            _epsilon1 = value;
            init_Nielson_type_surfaces();
            updateGL();
        }
    }

    void GLWidget::set_epsilon2(double value)
    {
        if (_epsilon2 != value)
        {
            _epsilon2 = value;
            init_Nielson_type_surfaces();
            updateGL();
        }
    }

    void GLWidget::set_epsilon3(double value)
    {
        if (_epsilon3 != value)
        {
            _epsilon3 = value;
            init_Nielson_type_surfaces();
            updateGL();
        }
    }

    void GLWidget::set_scale_parameter(double value)
    {
        if (_scale_parameter != value)
        {
            _scale_parameter = value;
            updateGL();
        }
    }

    void GLWidget::set_Nielson_index(int value)
    {
        if (_nielson_index != value)
        {
            _nielson_index = value;
            show_Nielson_type();
            updateGL();
        }
    }

    void GLWidget::set_model_index(int value)
    {
        if (_model_index != value)
        {
            _model_index = value;
            _heds.CleanAll();
            init_Nielson_type_surfaces();
            show_Nielson_type();
            updateGL();
        }
    }


    void GLWidget::set_dark_mode()
    {
        _dark_mode = !_dark_mode;
        // setting the background color
        if (_dark_mode)
        {
            glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
        } else {
            glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        }

        updateGL();
    }
}







