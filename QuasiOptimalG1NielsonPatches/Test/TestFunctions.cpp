#include <cmath>
#include "TestFunctions.h"
#include "../Core/Constants.h"

using namespace cagd;
using namespace std;

GLdouble spiral_on_cone::u_min = -TWO_PI;
GLdouble spiral_on_cone::u_max = +TWO_PI;

DCoordinate3 spiral_on_cone::d0(GLdouble u)
{
    return DCoordinate3(u * cos(u), u * sin(u), u);
}

DCoordinate3 spiral_on_cone::d1(GLdouble u)
{
    GLdouble c = cos(u), s = sin(u);
    return DCoordinate3(c - u * s, s + u * c, 1.0);
}

DCoordinate3 spiral_on_cone::d2(GLdouble u)
{
    GLdouble c = cos(u), s = sin(u);
    return DCoordinate3(-2.0 * s - u * c, 2.0 * c - u * s, 0);
}

//Torus

GLdouble torus::u_min=0;
GLdouble torus::u_max=6*PI;

DCoordinate3 torus::d0(GLdouble u) {
    GLdouble c = cos(u), c2=cos(2.0*u/3), s=sin(u), s2=sin(2.0*u/3);
    return DCoordinate3((2.0+c2)*c,(2.0+c2)*s, s2);
}

DCoordinate3 torus::d1(GLdouble u) {
    GLdouble c = cos(u), c2=cos(2.0*u/3.0), s=sin(u), s2=sin(2.0*u/3.0);
    return DCoordinate3(-(2.0+c2)*s-2.0/3.0*s2*c,(2.0+c2)*c-2.0/3.0*s2*s,2.0/3.0*c2);
}
DCoordinate3 torus::d2(GLdouble u) {
	GLdouble c = cos(u), c2 = cos(2.0*u / 3.0), s = sin(u), s2 = sin(2.0*u / 3.0);
    return DCoordinate3(-(2.0 + c2)*c + 2.0 / 3.0*s2*s - 2.0 / 3.0*(c2*c - s2 * s), -(2.0 + c2)*s - 2.0 / 3.0*s2*c - 2.0 / 3.0*(c2*s + s2 * c),
		-4.0 / 9.0*s2);
}

//Lissajous

GLdouble lissajous::u_min = -1.0;
GLdouble lissajous::u_max = 1.0;

DCoordinate3 lissajous::d0(GLdouble u) {
	return DCoordinate3(sin(5.0*u + PI / 2.0), sin(4.0*u), u);
}

DCoordinate3 lissajous::d1(GLdouble u) {
	return DCoordinate3(5.0*cos(5.0*u + PI / 2.0), 4.0*cos(4.0*u), 1.0);
}

DCoordinate3 lissajous::d2(GLdouble u) {
    return DCoordinate3(-25.0*sin(5.0*u + PI / 2.0), -16.0*sin(4.0*u), 0.0);
}

//Hypo

GLdouble hypo::u_min = -3.0;
GLdouble hypo::u_max = 3.0;

DCoordinate3 hypo::d0(GLdouble u) {
	return DCoordinate3(5.0*cos(u) + cos(5.0*u), 5.0*sin(u) - sin(5.0*u),u);
}

DCoordinate3 hypo::d1(GLdouble u) {
	return DCoordinate3(-5.0*sin(u) - 5.0*sin(5.0*u), 5.0*cos(u) - 5.0*cos(5.0*u), 1.0);
}

DCoordinate3 hypo::d2(GLdouble u) {
	return DCoordinate3(-5.0*cos(u) - 25.0*cos(5.0*u), -5.0*sin(u) + 25.0*sin(5.0*u), 0.0);
}

//Cyclo

GLdouble cyclo::u_min = 0;
GLdouble cyclo::u_max = 2.0*PI;

DCoordinate3 cyclo::d0(GLdouble u) {
	return DCoordinate3(2.0*(u - sin(u)), 2.0*(1.0 - cos(u)), u);
}

DCoordinate3 cyclo::d1(GLdouble u) {
	return DCoordinate3(2.0*(1.0 - cos(u)), 2.0*sin(u), 1.0);
}

DCoordinate3 cyclo::d2(GLdouble u) {
    return DCoordinate3(2.0*sin(u), 2.0*cos(u), 0.0);
}

//Ellipse

GLdouble ellipse::u_min = -6.0;
GLdouble ellipse::u_max = 6.0;

DCoordinate3 ellipse::d0(GLdouble u) {
	return DCoordinate3(6.0*cos(u), 4.0*sin(u), u);
}
DCoordinate3 ellipse::d1(GLdouble u){
    return DCoordinate3(-6.0*sin(u),4.0*cos(u),1.0);
}
DCoordinate3 ellipse::d2(GLdouble u){
    return DCoordinate3(-6.0*cos(u),-4.0*sin(u),0.0);
}

//Rose

GLdouble rose::u_min = -PI;
GLdouble rose::u_max = PI;

DCoordinate3 rose::d0(GLdouble u) {
	return DCoordinate3(1.0/2.0+cos(3.0*u)*cos(u), 1.0 / 2.0 + cos(3.0*u)*sin(u), u);
}
DCoordinate3 rose::d1(GLdouble u) {
	return DCoordinate3(-3.0*sin(3.0*u)*cos(u) - (1.0 / 2.0 + cos(3.0*u))*sin(u),
		(1.0 / 2.0 + cos(3.0*u))*cos(u) - 3.0*sin(3.0*u)*sin(u), 1.0);

}
DCoordinate3 rose::d2(GLdouble u) {
	GLdouble x = -3 * (3 * cos(3 * u) * cos(u) - sin(u) * sin(3 * u)) - 
		1 * cos(u) * (1 / 2 + cos(3 * u)) + 3 * sin(u) * sin(u);
	GLdouble y = -3 * (3 * cos(3 * u) * sin(u)) + sin(3 * u) * cos(u) -1 * sin(u) * (1 / 2 + cos(3 * u)) -
		3 * cos(u) * sin(3 * u);
	return DCoordinate3(x, y, 0.0);
}

//SURFACES

//Torus

GLdouble first::u_min = 0;
GLdouble first::u_max = PI;
GLdouble first::v_min = 0;
GLdouble first::v_max = PI;

DCoordinate3 first::d00(GLdouble u, GLdouble v){

	return DCoordinate3(sin(2 * u) * cos(v) * cos(v),sin(u) * sin(2 * v) / 2,cos(2 * u) * cos(v) * cos(v));
}

DCoordinate3 first::d10(GLdouble u, GLdouble v){

	return DCoordinate3(cos(v) * cos(v) * 2 * cos(2 * u),sin(2 * v) / 2 * cos(u),cos(v) * cos(v) * -2 * sin(2 * u));
}

DCoordinate3 first::d01(GLdouble u, GLdouble v){

	return DCoordinate3(sin(2 * u) * -2 * cos(v) * sin(v),sin(u) * cos(2 * v),cos(2 * u) * -2 * cos(v) * sin(v));
}

//Rose

GLdouble rose_surface::u_min = 0;
GLdouble rose_surface::u_max =  PI;
GLdouble rose_surface::v_min = 0;
GLdouble rose_surface::v_max =  TWO_PI;


DCoordinate3 rose_surface::d00(GLdouble u, GLdouble v) {
    GLdouble x = cos(4.0*u)*cos(u)*pow(cos(v) , 9);
    GLdouble y= cos(4.0*u)*sin(u)*pow(cos(v) , 9);
    GLdouble z = cos(4.0*u)*sin(v)*pow(cos(v) , 8);
	return DCoordinate3(x,y,z);
}
DCoordinate3 rose_surface::d01(GLdouble u, GLdouble v) {
    GLdouble x = (pow(cos(v) , 9))*(-4.0*sin(4.0*u)*cos(u) - cos(4.0*u)*sin(u));
    GLdouble y = pow(cos(v) , 9)*(cos(4.0*u)*cos(u) - 4.0*sin(4.0*u)*sin(u));
    GLdouble z = -4.0*sin(4.0*u)*sin(v)*pow(cos(v) , 8);
	return DCoordinate3(x, y, z);
}
DCoordinate3 rose_surface::d10(GLdouble u, GLdouble v) {
    GLdouble x = -9.0*cos(4.0*u)*cos(u)*pow(cos(v) , 8)*sin(v);
    GLdouble y = -9.0*cos(4.0*u)*sin(u)*pow(cos(v) , 8)*sin(v);
    GLdouble z = cos(4.0*u)*(cos(v)*pow(cos(v) , 8)-8.0*sin(v)*sin(v)*pow(cos(v) , 7));
	return DCoordinate3(x, y, z);
}

//Third

GLdouble third::u_min = -1.0;
GLdouble third::u_max = 1.0;
GLdouble third::v_min = -1.0;
GLdouble third::v_max = 1.0;

DCoordinate3 third::d00(GLdouble u, GLdouble v){

	return DCoordinate3(u, v, u * v);
}

DCoordinate3 third::d10(GLdouble u, GLdouble v){

	return DCoordinate3(1.0, 0.0, v);
}

DCoordinate3 third::d01(GLdouble u, GLdouble v){

	return DCoordinate3(0.0, -2.0 * sin(v), 2.0 * cos(v));
}

//Fourth

GLdouble fourth::u_min = -0.0;
GLdouble fourth::u_max = 3.0;
GLdouble fourth::v_min = 0.0;
GLdouble fourth::v_max = TWO_PI;

DCoordinate3 fourth::d00(GLdouble u, GLdouble v)
{
	return DCoordinate3(2.0 * sin(v), 2.0 * cos(v), u);
}

DCoordinate3 fourth::d10(GLdouble u, GLdouble v)
{
	return DCoordinate3(0.0, 0.0, 1.0);
}

DCoordinate3 fourth::d01(GLdouble u, GLdouble v)
{
	return DCoordinate3(2.0 * cos(v), -2.0 * sin(v), 0.0);
}

//Fifth

GLdouble fifth::u_min = -PI;
GLdouble fifth::u_max = PI;
GLdouble fifth::v_min = -PI;
GLdouble fifth::v_max = PI;

DCoordinate3 fifth::d00(GLdouble u, GLdouble v) {
	return DCoordinate3(cos(u) * cos(v), u, cos(u) * sin(v));
}

DCoordinate3 fifth::d10(GLdouble u, GLdouble v) {
	return DCoordinate3(cos(v) * -1 * sin(u), 1, -sin(u) * sin(v));
}

DCoordinate3 fifth::d01(GLdouble u, GLdouble v) {
	return DCoordinate3(cos(u) * -1 * sin(v), 0, cos(u) * cos(v));
}

//Sixth

GLdouble sixth::u_min = 0;
GLdouble sixth::u_max = PI;
GLdouble sixth::v_min = 0;
GLdouble sixth::v_max = TWO_PI;

DCoordinate3 sixth::d00(GLdouble u, GLdouble v){
	return DCoordinate3(sin(u) * cos(v), sin(u) * sin(v), cos(u));
}

DCoordinate3 sixth::d10(GLdouble u, GLdouble v){
	return DCoordinate3(cos(v) * cos(u), cos(u) * sin(v), -1 * sin(u));
}

DCoordinate3 sixth::d01(GLdouble u, GLdouble v){
	return DCoordinate3(-1 * sin(v) * sin(u), sin(u) * cos(v), 0);
}

//Dini
GLdouble dini::u_min = -TWO_PI;
GLdouble dini::u_max =  TWO_PI;
GLdouble dini::v_min = 0;
GLdouble dini::v_max =  TWO_PI;


DCoordinate3 dini::d00(GLdouble u, GLdouble v) {
    GLdouble x = cos(u)*sin(v);
    GLdouble y= sin(u)*sin(v);
    GLdouble z = cos(v)+log(tan(v/2.0))+u;
    return DCoordinate3(x,y,z);
}
DCoordinate3 dini::d10(GLdouble u, GLdouble v) {
    GLdouble x = -sin(u)*sin(v);
    GLdouble y = cos(u)*sin(v);
    GLdouble z = 1.0;
    return DCoordinate3(x, y, z);
}
DCoordinate3 dini::d01(GLdouble u, GLdouble v) {
    GLdouble x = cos(u)*cos(v);
    GLdouble y = sin(u)*cos(v);
    GLdouble z = -sin(v)+1.0/2.0*tan(v/2.0)+1.0/2.0;
    return DCoordinate3(x, y, z);
}

//Torus
GLdouble torus_surf::u_min = 0;
GLdouble torus_surf::u_max =  TWO_PI;
GLdouble torus_surf::v_min = 0;
GLdouble torus_surf::v_max =  TWO_PI;


DCoordinate3 torus_surf::d00(GLdouble u, GLdouble v) {
    GLdouble x = (2.0+cos(u))*cos(v);
    GLdouble y= (2.0+cos(u))*sin(v);
    GLdouble z = 2.0+sin(u);
    return DCoordinate3(x,y,z);
}
DCoordinate3 torus_surf::d01(GLdouble u, GLdouble v) {
    GLdouble x = -sin(u)*cos(v);
    GLdouble y = -sin(u)*sin(v);
    GLdouble z = cos(u);
    return DCoordinate3(x, y, z);
}
DCoordinate3 torus_surf::d10(GLdouble u, GLdouble v) {
    GLdouble x = -(2.0+cos(u))*sin(v);
    GLdouble y = (2.0+cos(u))*cos(v);
    GLdouble z = 0;
    return DCoordinate3(x, y, z);
}
