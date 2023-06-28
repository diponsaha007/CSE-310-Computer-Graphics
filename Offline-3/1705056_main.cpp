#include <bits/stdc++.h>
#include <GL/glut.h>
#include "1705056_classes.h"
#include "bitmap_image.hpp"

using namespace std;

int recursion_level;
int dimension;
int imageID = 11;

void drawAxes() {
    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINES);
    {
        glVertex3f(500, 0, 0);
        glVertex3f(-500, 0, 0);

        glVertex3f(0, -500, 0);
        glVertex3f(0, 500, 0);

        glVertex3f(0, 0, 500);
        glVertex3f(0, 0, -500);
    }
    glEnd();
}

int getBitmapColor(double color) {
    if (color < 0)
        return 0;
    else if (color >= 1)
        return 255;
    else
        return round(color * 255);
}

void capture() {
    string filename = "output_" + to_string(imageID) + ".bmp";
    imageID++;

    bitmap_image bitmapImage(dimension, dimension);

    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            bitmapImage.set_pixel(i, j, 0, 0, 0);
        }
    }
    double windowHeight = 500;
    double windowWidth = 500;
    double planeDistance = windowHeight / (2.0 * tan(80.00 / 2.0 * M_PI / 180.0));
    Point topleft = camera.pos + camera.l * planeDistance - camera.r * windowWidth / 2.00 +
                    camera.u * windowHeight / 2.00;
    double du = windowWidth * 1.00 / dimension;
    double dv = windowHeight * 1.00 / dimension;
    topleft = topleft + camera.r * (0.5 * du) - camera.u * (0.5 * dv);
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            Point p = topleft + camera.r * (i * du) - camera.u * (j * dv);
            Ray ray = Ray(camera.pos, p - camera.pos);
            double tMin = inf;
            int nearest = -1;
            for (int k = 0; k < objects.size(); k++) {
                double now = objects[k]->intersect(ray);
                if (now > 0 && now < inf) {
                    if (now < tMin) {
                        tMin = now;
                        nearest = k;
                    }
                }
            }
            if (nearest == -1)
                continue;
            Color color = objects[nearest]->recursiveIntersection(nearest, recursion_level, ray);
            bitmapImage.set_pixel(i, j, getBitmapColor(color.rgb[0]), getBitmapColor(color.rgb[1]),
                                  getBitmapColor(color.rgb[2]));
        }
    }

    bitmapImage.save_image(filename);
    cout << "Bitmap image saved. Filename = " << filename << "\n";
}

void keyboardListener(unsigned char key, int x, int y) {
    switch (key) {
        case '0':
            capture();
            break;
        case '1':
            camera.rotate(camera.r, camera.l, false);
            break;
        case '2':
            camera.rotate(camera.r, camera.l, true);
            break;
        case '3':
            camera.rotate(camera.l, camera.u, false);
            break;
        case '4':
            camera.rotate(camera.l, camera.u, true);
            break;
        case '5':
            camera.rotate(camera.u, camera.r, false);
            break;
        case '6':
            camera.rotate(camera.u, camera.r, true);
            break;
        default:
            break;
    }
}


void specialKeyListener(int key, int x, int y) {
    switch (key) {
        case GLUT_KEY_DOWN:
            camera.move(camera.l, false);
            break;
        case GLUT_KEY_UP:
            camera.move(camera.l, true);
            break;
        case GLUT_KEY_RIGHT:
            camera.move(camera.r, true);
            break;
        case GLUT_KEY_LEFT:
            camera.move(camera.r, false);
            break;
        case GLUT_KEY_PAGE_UP:
            camera.move(camera.u, true);
            break;
        case GLUT_KEY_PAGE_DOWN:
            camera.move(camera.u, false);
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
    gluLookAt(camera.pos.x, camera.pos.y, camera.pos.z, camera.pos.x + camera.l.x, camera.pos.y + camera.l.y,
              camera.pos.z + camera.l.z,
              camera.u.x, camera.u.y, camera.u.z);

    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);


    /****************************
    / Add your objects from here
    ****************************/
    //add objects

//    drawAxes();
    for (int i = 0; i < objects.size(); i++) {
        objects[i]->draw();
    }

    for (int i = 0; i < plights.size(); i++) {
        plights[i].draw();
    }

    for (int i = 0; i < slights.size(); i++) {
        slights[i].draw();
    }
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

void loadData() {
    InsertFloor();
    ifstream input("input.txt");
    input >> recursion_level >> dimension;

    int n;
    input >> n;

    for (int i = 0; i < n; i++) {
        string s;
        input >> s;
        if (s == "sphere") {

            Sphere s;
            input >> s;
            Object *now = new Sphere(s.center, s.radius, s.color, s.ambient, s.diffuse, s.specular, s.shine,
                                     s.recursive);
            objects.push_back(now);
        } else if (s == "triangle") {
            Triangle t;
            input >> t;
            Object *now = new Triangle(t.p1, t.p2, t.p3, t.color, t.ambient, t.diffuse, t.specular, t.shine,
                                       t.recursive);
            objects.push_back(now);
        } else if (s == "general") {
            GeneralQuadricSurfaces g;
            input >> g;
            Object *now = new GeneralQuadricSurfaces(g.a, g.b, g.c, g.d, g.e, g.f, g.g, g.h, g.i, g.j,
                                                     g.cube_reference_point, g.length, g.width, g.height, g.color,
                                                     g.ambient, g.diffuse, g.specular, g.shine, g.recursive);
            objects.push_back(now);
        }
    }

    input >> n;
    for (int i = 0; i < n; i++) {
        PointLight pl;
        input >> pl;
        plights.push_back(pl);
    }
    input >> n;
    for (int i = 0; i < n; i++) {
        SpotLight sl;
        input >> sl;
        slights.push_back(sl);
    }
}

void ClearVectors() {
    for (int i = 0; i < objects.size(); i++) {
        delete objects[i];
    }
    objects.clear();
    plights.clear();
    slights.clear();
}

int main(int argc, char **argv) {
    glutInit(&argc, argv);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);    //Depth, Double buffer, RGB color

    glutCreateWindow("Ray Tracing");

    init();

    glEnable(GL_DEPTH_TEST);    //enable Depth Testing

    glutDisplayFunc(display);    //display callback function
    glutIdleFunc(animate);        //what you want to do in the idle time (when no drawing is occuring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);

    loadData();
    atexit(ClearVectors);
    glutMainLoop();        //The main loop of OpenGL

    return 0;
}
