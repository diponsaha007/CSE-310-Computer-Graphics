#include <bits/stdc++.h>
#include <GL/glut.h>

using namespace std;
const double pi = 2 * acos(0.0);

double cameraHeight = 140.0;
double cameraAngle = 1.0;

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

    Point &operator+=(const Point &t) {
        x += t.x;
        y += t.y;
        z += t.z;
        return *this;
    }

    Point &operator-=(const Point &t) {
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

    Point operator+(const Point &t) const {
        return Point(*this) += t;
    }

    Point operator-(const Point &t) const {
        return Point(*this) -= t;
    }

    Point operator*(double t) const {
        return Point(*this) *= t;
    }

    Point operator/(double t) const {
        return Point(*this) /= t;
    }
};


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

void drawGrid() {
    int i;
    glColor3f(0.4, 0.4, 0.4);    //grey
    int lim = 18;
    glBegin(GL_LINES);
    {
        for (i = -lim; i <= lim; i++) {
            //lines parallel to Y-axis
            glVertex3f(i * 15, -lim * 15, 0);
            glVertex3f(i * 15, lim * 15, 0);
            //lines parallel to X-axis
            glVertex3f(-lim * 15, i * 15, 0);
            glVertex3f(lim * 15, i * 15, 0);
        }
    }
    glEnd();

}

void drawSquare(double a) {
    glBegin(GL_QUADS);
    {
        glVertex3f(a, a, 2);
        glVertex3f(a, -a, 2);
        glVertex3f(-a, -a, 2);
        glVertex3f(-a, a, 2);
    }
    glEnd();
}

class WheelMovement {
public:
    Point pos;
    Point dir;
    int angle;
    int angle2;
    const int MOVEMENT_UNIT = 3;
    const int ROTATION_UNIT_DEG = 2;
    const double ROTATION_UNIT_RAD = ROTATION_UNIT_DEG * pi / 180.00;


    WheelMovement() {
        pos = Point(0, 0, 0);
        dir = Point(1, 0, 0);
        angle = 90;
        angle2 = 0;
    }

    void move(bool positive) {
        pos = pos + dir * (positive ? 1 : -1) * MOVEMENT_UNIT;
        angle2 += (positive ? 1 : -1) * ROTATION_UNIT_DEG * 2;
        angle2 %= 360;
    }

    void rotate(bool anticlockwise) {
        int multiply = (anticlockwise ? 1 : -1);
        angle += multiply * ROTATION_UNIT_DEG;
        angle %= 360;
        Point tmp(0, 0, 0);
        tmp.x = dir.x * cos(ROTATION_UNIT_RAD * multiply) - dir.y * sin(ROTATION_UNIT_RAD * multiply);
        tmp.y = dir.x * sin(ROTATION_UNIT_RAD * multiply) + dir.y * cos(ROTATION_UNIT_RAD * multiply);
        dir = tmp;
        dir.normalize();
    }
} wheel;

void drawWheel(double radius, int height, int segments = 50) {
    vector<Point> points(segments + 1);
    for (int i = 0; i <= segments; i++) {
        points[i].x = radius * cos(((double) i / (double) segments) * 2 * pi);
        points[i].y = radius * sin(((double) i / (double) segments) * 2 * pi);
    }
    for (int i = 0; i < segments; i++) {
        double shade;
        if (i < segments / 2)shade = 2 * (double) i / (double) segments;
        else shade = 2 * (1.0 - (double) i / (double) segments);
        shade *= 0.8;
        glColor3f(shade, shade, shade);
        glBegin(GL_QUADS);
        {
            glVertex3f(points[i].x, points[i].y, 0);
            glVertex3f(points[i].x, points[i].y, height);
            glVertex3f(points[i + 1].x, points[i + 1].y, height);
            glVertex3f(points[i + 1].x, points[i + 1].y, 0);
        }
        glEnd();
    }
    glColor3f(0.55, 0.55, 0.55);
    glBegin(GL_QUADS);
    {
        glVertex3f(0, radius, 0);
        glVertex3f(0, radius, height);
        glVertex3f(0, -radius, height);
        glVertex3f(0, -radius, 0);
    }
    glEnd();
    glBegin(GL_QUADS);
    {
        glVertex3f(radius, 0, 0);
        glVertex3f(radius, 0, height);
        glVertex3f(-radius, 0, height);
        glVertex3f(-radius, 0, 0);
    }
    glEnd();
}

void drawScene() {
    //drawAxes();
    int radius = 32;
    int height = 14;
    drawGrid();
    glPushMatrix();
    {
        glTranslated(wheel.pos.x, wheel.pos.y, 0);
        glRotated(wheel.angle, 0, 0, 1);
        glTranslated(-height / 2.00, 0, radius);
        glRotated(90, 0, 1, 0);
        glRotated(wheel.angle2, 0, 0, 1);
        drawWheel(radius, height);
    }
    glPopMatrix();

}

void keyboardListener(unsigned char key, int x, int y) {
    switch (key) {
        case 'w':
            wheel.move(0);
            break;
        case 's':
            wheel.move(1);
            break;
        case 'a':
            wheel.rotate(1);
            break;
        case 'd':
            wheel.rotate(0);
            break;
        default:
            break;
    }
}


void specialKeyListener(int key, int x, int y) {
    switch (key) {
        case GLUT_KEY_DOWN:        //down arrow key
            cameraHeight -= 3.0;
            break;
        case GLUT_KEY_UP:        // up arrow key
            cameraHeight += 3.0;
            break;

        case GLUT_KEY_RIGHT:
            cameraAngle += 0.03;
            break;
        case GLUT_KEY_LEFT:
            cameraAngle -= 0.03;
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

    //gluLookAt(100,100,100,	0,0,0,	0,0,1);
    gluLookAt(200 * cos(cameraAngle), 200 * sin(cameraAngle), cameraHeight, 0, 0, 0, 0, 0, 1);
    //gluLookAt(0,0,200,	0,0,0,	1,0,0);

    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);
    /****************************
    / Add your objects from here
    ****************************/
    //add objects
    drawScene();
    //ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    glutSwapBuffers();
}


void animate() {
    //codes for any changes in Models, Camera
    glutPostRedisplay();
}

void init() {

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
    glutCreateWindow("Wheel");
    init();
    glEnable(GL_DEPTH_TEST);    //enable Depth Testing
    glutDisplayFunc(display);    //display callback function
    glutIdleFunc(animate);        //what you want to do in the idle time (when no drawing is occuring)
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMainLoop();        //The main loop of OpenGL

    return 0;
}
