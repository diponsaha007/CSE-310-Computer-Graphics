#include <bits/stdc++.h>
#include <GL/glut.h>

using namespace std;
const double pi = 2 * acos(0.0);
const int total_size = 36;
int cur_sq = 17;

struct Point {
    double x, y, z;

    Point() {}

    Point(double x, double y, double z) : x(x), y(y), z(z) {}

    Point(const Point &p) {
        x = p.x;
        y = p.y;
        z = p.z;
    }

    void normalize() {
        double len = sqrt(x * x + y * y + z * z);
        x /= len;
        y /= len;
        z /= len;
    }

    Point &operator+=(const Point& t) {
        x += t.x;
        y += t.y;
        z += t.z;
        return *this;
    }

    Point &operator-=(const Point& t) {
        x -= t.x;
        y -= t.y;
        z -= t.z;
        return *this;
    }

    Point &operator*=(double t) {
        x *= t;
        y *= t;
        z *= t;
        return *this;
    }

    Point &operator/=(double t) {
        x /= t;
        y /= t;
        z /= t;
        return *this;
    }

    Point operator+(const Point& t) const {
        return Point(*this) += t;
    }

    Point operator-(const Point& t) const {
        return Point(*this) -= t;
    }

    Point operator*(double t) const {
        return Point(*this) *= t;
    }

    Point operator/(double t) const {
        return Point(*this) /= t;
    }
};


class Camera {
public:
    Point pos, u, r, l;
    const int MOVE_CONSTANT = 2;
    const double ROTATION_CONSTANT = 3.00 * pi / 180.00;

    Camera() {
        double val = 1.00 / sqrt(2);
        pos = Point(100, 100, 0);
        u = Point(0, 0, 1);
        r = Point(-val, val, 0);
        l = Point(-val, -val, 0);
    }

    void move(const Point& p, bool pos_dir) {
        pos = pos + p * (pos_dir ? 1 : -1) * MOVE_CONSTANT;
    }

    void rotate(Point &p1, Point &p2, bool pos_rotation) const {
        double angle = (pos_rotation ? 1 : -1) * ROTATION_CONSTANT;
        Point p3, p4;
        p3 = p1 * cos(angle) - p2 * sin(angle);
        p4 = p1 * sin(angle) + p2 * cos(angle);
        p1 = p3;
        p2 = p4;
        p1.normalize();
        p2.normalize();
    }
} cam;


void drawAxes() {
    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINES);
    {
        glVertex3f(100, 0, 0);
        glVertex3f(-100, 0, 0);

        glVertex3f(0, -100, 0);
        glVertex3f(0, 100, 0);

        glVertex3f(0, 0, 100);
        glVertex3f(0, 0, -100);
    }
    glEnd();
}


void drawSquare(double a) {
    glBegin(GL_QUADS);
    {
        glVertex3f(a, a, 0);
        glVertex3f(a, -a, 0);
        glVertex3f(-a, -a, 0);
        glVertex3f(-a, a, 0);
    }
    glEnd();
}

void drawSphereOneEighth(double radius, int slices = 24, int stacks = 24) {
    vector<vector<Point> > points(stacks + 1, vector<Point>(slices + 1));
    int i, j;
    double h, r;
    for (i = 0; i <= stacks; i++) {
        h = radius * sin(((double) i / (double) stacks) * (pi / 2));
        r = radius * cos(((double) i / (double) stacks) * (pi / 2));
        for (j = 0; j <= slices; j++) {
            points[i][j].x = r * cos(((double) j / (double) slices) * (pi / 2));
            points[i][j].y = r * sin(((double) j / (double) slices) * (pi / 2));
            points[i][j].z = h;
        }
    }

    for (i = 0; i < stacks; i++) {
        for (j = 0; j < slices; j++) {
            glBegin(GL_QUADS);
            {
                glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
                glVertex3f(points[i][j + 1].x, points[i][j + 1].y, points[i][j + 1].z);
                glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i + 1][j + 1].z);
                glVertex3f(points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z);
            }
            glEnd();
        }
    }
}

void drawCylinderOneFourth(double radius, int height, int slices = 24, int stacks = 24) {
    vector<vector<Point> > points(stacks + 1, vector<Point>(slices + 1));
    int i, j;
    double h, r;
    //generate points
    for (i = 0; i <= stacks; i++) {
        h = height * sin(((double) i / (double) stacks) * (pi / 2));
        r = radius;
        for (j = 0; j <= slices; j++) {
            points[i][j].x = r * cos(((double) j / (double) slices) * (pi / 2));
            points[i][j].y = r * sin(((double) j / (double) slices) * (pi / 2));
            points[i][j].z = h;
        }
    }
    //draw quads using generated points
    for (i = 0; i < stacks; i++) {
        //glColor3f((double) i / (double) stacks, (double) i / (double) stacks, (double) i / (double) stacks);
        for (j = 0; j < slices; j++) {
            glBegin(GL_QUADS);
            {
                //upper hemisphere
                glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
                glVertex3f(points[i][j + 1].x, points[i][j + 1].y, points[i][j + 1].z);
                glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i + 1][j + 1].z);
                glVertex3f(points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z);
                glVertex3f(points[i][j + 1].x, points[i][j + 1].y, -points[i][j + 1].z);
                glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, -points[i + 1][j + 1].z);
                glVertex3f(points[i + 1][j].x, points[i + 1][j].y, -points[i + 1][j].z);
            }
            glEnd();
        }
    }
}


void drawBoxSquares() {
    glColor3f(1, 1, 1);
    glPushMatrix();
    {
        glTranslated(0, 0, total_size);
        drawSquare(cur_sq);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslated(0, 0, -total_size);
        drawSquare(cur_sq);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslated(0, total_size, 0);
        glRotated(90, 1, 0, 0);
        drawSquare(cur_sq);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslated(0, -total_size, 0);
        glRotated(-90, 1, 0, 0);
        drawSquare(cur_sq);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslated(total_size, 0, 0);
        glRotated(90, 0, 1, 0);
        drawSquare(cur_sq);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslated(-total_size, 0, 0);
        glRotated(-90, 0, 1, 0);
        drawSquare(cur_sq);
    }
    glPopMatrix();
}

void drawBoxSpheres() {
    glColor3f(1, 0, 0);
    int val = total_size - cur_sq;
    int val2 = cur_sq;
    //upper side
    glPushMatrix();
    {
        glTranslated(val2, val2, val2);
        drawSphereOneEighth(val);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslated(-val2, val2, val2);
        glRotated(90, 0, 0, 1);
        drawSphereOneEighth(val);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslated(-val2, -val2, val2);
        glRotated(180, 0, 0, 1);
        drawSphereOneEighth(val);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslated(val2, -val2, val2);
        glRotated(270, 0, 0, 1);
        drawSphereOneEighth(val);
    }
    glPopMatrix();

    //lower side
    glPushMatrix();
    {
        glTranslated(val2, val2, -val2);
        glRotated(270, 1, 0, 0);
        drawSphereOneEighth(val);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslated(-val2, val2, -val2);
        glRotated(90, 0, 0, 1);
        glRotated(270, 1, 0, 0);
        drawSphereOneEighth(val);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslated(-val2, -val2, -val2);
        glRotated(180, 0, 0, 1);
        glRotated(270, 1, 0, 0);
        drawSphereOneEighth(val);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslated(val2, -val2, -val2);
        glRotated(270, 0, 0, 1);
        glRotated(270, 1, 0, 0);
        drawSphereOneEighth(val);
    }
    glPopMatrix();
}

void drawBoxCylinders() {
    glColor3f(0, 1, 0);
    int val = total_size - cur_sq;
    int val2 = cur_sq;
    // 4 corner
    glPushMatrix();
    {
        glTranslated(val2, val2, 0);
        drawCylinderOneFourth(val, val2);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslated(-val2, val2, 0);
        glRotated(90, 0, 0, 1);
        drawCylinderOneFourth(val, val2);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslated(-val2, -val2, 0);
        glRotated(180, 0, 0, 1);
        drawCylinderOneFourth(val, val2);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslated(val2, -val2, 0);
        glRotated(270, 0, 0, 1);
        drawCylinderOneFourth(val, val2);
    }
    glPopMatrix();

    //upper 4
    glPushMatrix();
    {
        glTranslated(0, val2, val2);
        glRotated(90, 0, 0, 1);
        glRotated(90, 1, 0, 0);
        drawCylinderOneFourth(val, val2);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslated(-val2, 0, val2);
        glRotated(90 * 2, 0, 0, 1);
        glRotated(90, 1, 0, 0);
        drawCylinderOneFourth(val, val2);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslated(0, -val2, val2);
        glRotated(90 * 3, 0, 0, 1);
        glRotated(90, 1, 0, 0);
        drawCylinderOneFourth(val, val2);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslated(val2, 0, val2);
        glRotated(90 * 4, 0, 0, 1);
        glRotated(90, 1, 0, 0);
        drawCylinderOneFourth(val, val2);
    }
    glPopMatrix();

    //lower 4
    glPushMatrix();
    {
        glTranslated(0, val2, -val2);
        glRotated(90, 0, 0, 1);
        glRotated(270, 1, 0, 0);
        drawCylinderOneFourth(val, val2);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslated(-val2, 0, -val2);
        glRotated(90 * 2, 0, 0, 1);
        glRotated(270, 1, 0, 0);
        drawCylinderOneFourth(val, val2);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslated(0, -val2, -val2);
        glRotated(90 * 3, 0, 0, 1);
        glRotated(270, 1, 0, 0);
        drawCylinderOneFourth(val, val2);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslated(val2, 0, -val2);
        glRotated(90 * 4, 0, 0, 1);
        glRotated(270, 1, 0, 0);
        drawCylinderOneFourth(val, val2);
    }
    glPopMatrix();
}

void drawBox() {
    drawBoxSquares();
    drawBoxSpheres();
    drawBoxCylinders();
}

void keyboardListener(unsigned char key, int x, int y) {
    switch (key) {
        case '1':
            cam.rotate(cam.r, cam.l, false);
            break;
        case '2':
            cam.rotate(cam.r, cam.l, true);
            break;
        case '3':
            cam.rotate(cam.l, cam.u, false);
            break;
        case '4':
            cam.rotate(cam.l, cam.u, true);
            break;
        case '5':
            cam.rotate(cam.u, cam.r, false);
            break;
        case '6':
            cam.rotate(cam.u, cam.r, true);
            break;
        default:
            break;
    }
}


void specialKeyListener(int key, int x, int y) {
    switch (key) {
        case GLUT_KEY_DOWN:
            cam.move(cam.l, false);
            break;
        case GLUT_KEY_UP:
            cam.move(cam.l, true);
            break;
        case GLUT_KEY_RIGHT:
            cam.move(cam.r, true);
            break;
        case GLUT_KEY_LEFT:
            cam.move(cam.r, false);
            break;
        case GLUT_KEY_PAGE_UP:
            cam.move(cam.u, true);
            break;
        case GLUT_KEY_PAGE_DOWN:
            cam.move(cam.u, false);
            break;
        case GLUT_KEY_HOME:
            if (cur_sq > 0)
                cur_sq--;
            break;
        case GLUT_KEY_END:
            if (cur_sq < total_size)
                cur_sq++;
            break;
        default:
            break;
    }
}


void display() {

    //clear the display
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0, 0, 0, 0);    //color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    /********************
    / set-up camera here
    ********************/
    //load the correct matrix -- MODEL-VIEW matrix
    glMatrixMode(GL_MODELVIEW);

    //initialize the matrix
    glLoadIdentity();

    //now give three info
    //1. where is the camera (viewer)?
    //2. where is the camera looking?
    //3. Which direction is the camera's UP direction?
//    gluLookAt(100, 100, 100, 0, 0, 0, 0, 0, 1);
    gluLookAt(cam.pos.x, cam.pos.y, cam.pos.z, cam.pos.x + cam.l.x, cam.pos.y + cam.l.y, cam.pos.z + cam.l.z,
              cam.u.x, cam.u.y, cam.u.z);

    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);


    /****************************
    / Add your objects from here
    ****************************/
    //add objects

    drawAxes();
    drawBox();

    //ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    glutSwapBuffers();
}


void animate() {
    //codes for any changes in Models, Camera
    glutPostRedisplay();
}

void init() {
    //codes for initialization
    //clear the screen
    glClearColor(0, 0, 0, 0);

    /************************
    / set-up projection here
    ************************/
    //load the PROJECTION matrix
    glMatrixMode(GL_PROJECTION);

    //initialize the matrix
    glLoadIdentity();

    //give PERSPECTIVE parameters
    gluPerspective(80, 1, 1, 1000.0);
    //field of view in the Y (vertically)
    //aspect ratio that determines the field of view in the X direction (horizontally)
    //near distance
    //far distance
}

int main(int argc, char **argv) {
    glutInit(&argc, argv);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);    //Depth, Double buffer, RGB color

    glutCreateWindow("Fully Controllable Camera & Sphere to/from Cube");

    init();

    glEnable(GL_DEPTH_TEST);    //enable Depth Testing

    glutDisplayFunc(display);    //display callback function
    glutIdleFunc(animate);        //what you want to do in the idle time (when no drawing is occuring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMainLoop();        //The main loop of OpenGL

    return 0;
}
