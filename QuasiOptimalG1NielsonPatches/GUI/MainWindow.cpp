#include "MainWindow.h"

namespace cagd
{
    MainWindow::MainWindow(QWidget *parent): QMainWindow(parent)
    {
        setupUi(this);

    /*

      the structure of the main window's central widget

     *---------------------------------------------------*
     |                 central widget                    |
     |                                                   |
     |  *---------------------------*-----------------*  |
     |  |     rendering context     |   scroll area   |  |
     |  |       OpenGL widget       | *-------------* |  |
     |  |                           | | side widget | |  |
     |  |                           | |             | |  |
     |  |                           | |             | |  |
     |  |                           | *-------------* |  |
     |  *---------------------------*-----------------*  |
     |                                                   |
     *---------------------------------------------------*

    */
        _side_widget = new SideWidget(this);

        _scroll_area = new QScrollArea(this);
        _scroll_area->setWidget(_side_widget);
        _scroll_area->setSizePolicy(_side_widget->sizePolicy());
        _scroll_area->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);

        _gl_widget = new GLWidget(this);

        centralWidget()->setLayout(new QHBoxLayout());
        centralWidget()->layout()->addWidget(_gl_widget);
        centralWidget()->layout()->addWidget(_scroll_area);

        // creating a signal slot mechanism between the rendering context and the side widget
        connect(_side_widget->rotate_x_slider, SIGNAL(valueChanged(int)), _gl_widget, SLOT(set_angle_x(int)));
        connect(_side_widget->rotate_y_slider, SIGNAL(valueChanged(int)), _gl_widget, SLOT(set_angle_y(int)));
        connect(_side_widget->rotate_z_slider, SIGNAL(valueChanged(int)), _gl_widget, SLOT(set_angle_z(int)));

        connect(_side_widget->zoom_factor_spin_box, SIGNAL(valueChanged(double)), _gl_widget, SLOT(set_zoom_factor(double)));

        connect(_side_widget->trans_x_spin_box, SIGNAL(valueChanged(double)), _gl_widget, SLOT(set_trans_x(double)));
        connect(_side_widget->trans_y_spin_box, SIGNAL(valueChanged(double)), _gl_widget, SLOT(set_trans_y(double)));
        connect(_side_widget->trans_z_spin_box, SIGNAL(valueChanged(double)), _gl_widget, SLOT(set_trans_z(double)));

        connect(_side_widget->model_comboBox, SIGNAL(currentIndexChanged(int)), _gl_widget, SLOT(set_model_index(int)));

        connect(_side_widget->material_comboBox, SIGNAL(currentIndexChanged(int)), _gl_widget, SLOT(set_color_index(int)));
        connect(_side_widget->scale_parameter_spinBox, SIGNAL(valueChanged(double)),_gl_widget,SLOT(set_scale_parameter(double)));

        connect(_side_widget->comboBox_Shader, SIGNAL(currentIndexChanged(int)),_gl_widget,SLOT(set_shader(int)));
        connect(_side_widget->shader_enabled_checkbox, SIGNAL(stateChanged(int)), _gl_widget, SLOT(set_shader_enabled()));
        connect(_side_widget->scale_spinbox, SIGNAL(valueChanged(double)), _gl_widget, SLOT(set_shader_scale(double)));
        connect(_side_widget->smoothing_spinbox, SIGNAL(valueChanged(double)), _gl_widget, SLOT(set_shader_smoothing(double)));
        connect(_side_widget->shading_spinbox, SIGNAL(valueChanged(double)), _gl_widget, SLOT(set_shader_shading(double)));

        connect(_side_widget->red_spinbox, SIGNAL(valueChanged(double)), _gl_widget, SLOT(set_shader_red(double)));
        connect(_side_widget->green_spinbox, SIGNAL(valueChanged(double)), _gl_widget, SLOT(set_shader_green(double)));
        connect(_side_widget->blue_spinbox, SIGNAL(valueChanged(double)), _gl_widget, SLOT(set_shader_blue(double)));
        connect(_side_widget->alpha_spinbox, SIGNAL(valueChanged(double)), _gl_widget, SLOT(set_shader_alpha(double)));

        //Nielson-type triangular surfaces
        connect(_side_widget->dark_mode, SIGNAL(stateChanged(int)), _gl_widget, SLOT(set_dark_mode()));

        connect(_side_widget->theta1_spinBox, SIGNAL(valueChanged(double)), _gl_widget, SLOT(set_theta1(double)));
        connect(_side_widget->theta2_spinBox, SIGNAL(valueChanged(double)), _gl_widget, SLOT(set_theta2(double)));
        connect(_side_widget->theta3_spinBox, SIGNAL(valueChanged(double)), _gl_widget, SLOT(set_theta3(double)));

        connect(_side_widget->beta_spinBox, SIGNAL(valueChanged(double)), _gl_widget, SLOT(set_beta(double)));
        connect(_side_widget->div_point_spinBox, SIGNAL(valueChanged(int)), _gl_widget, SLOT(set_div_point_count(int)));

        connect(_side_widget->boundaryCurve_checkBox, SIGNAL(stateChanged(int)), _gl_widget, SLOT(show_boundaries()));
        connect(_side_widget->halfEdges_checkBox, SIGNAL(stateChanged(int)), _gl_widget, SLOT(render_halfEdges()));
        connect(_side_widget->controlNet_checkBox, SIGNAL(stateChanged(int)), _gl_widget, SLOT(show_control_net()));
        connect(_side_widget->controlPolygon_checkBox, SIGNAL(stateChanged(int)), _gl_widget, SLOT(show_control_polygons()));
        connect(_side_widget->vertex_checkBox, SIGNAL(stateChanged(int)), _gl_widget, SLOT(show_vertices()));
        connect(_side_widget->vertexNormals_checkBox, SIGNAL(stateChanged(int)), _gl_widget, SLOT(show_vertex_normals()));
        connect(_side_widget->averagedNormals_checkBox, SIGNAL(stateChanged(int)), _gl_widget, SLOT(show_G1_normals()));
        connect(_side_widget->boundaryNormals_checkBox, SIGNAL(stateChanged(int)), _gl_widget, SLOT(show_C0_normals()));
        connect(_side_widget->tangent_checkBox, SIGNAL(stateChanged(int)), _gl_widget, SLOT(show_tangents()));

        connect(_side_widget->nielson_comboBox, SIGNAL(currentIndexChanged(int)), _gl_widget, SLOT(set_Nielson_index(int)));

        connect(_side_widget->basis_function_comboBox, SIGNAL(currentIndexChanged(int)), _gl_widget, SLOT(set_boundary_type(int)));

        connect(_gl_widget, SIGNAL(boundary_type_changed(bool)), _side_widget->beta_spinBox, SLOT(setEnabled(bool)));

    }

    //--------------------------------
    // implementation of private slots
    //--------------------------------
    void MainWindow::on_action_Quit_triggered()
    {
        qApp->exit(0);
    }
}
