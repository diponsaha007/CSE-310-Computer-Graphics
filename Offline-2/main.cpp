#include <bits/stdc++.h>
#include "bitmap_image.hpp"

using namespace std;

const double pi = 2 * acos(0.0);
const double eps = 1e-9;
const int N = 4;
mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

int getrand(int a, int b) {
    int x = uniform_int_distribution<int>(a, b)(rng);
    return x;
}

int sign(double x) { return (x > eps) - (x < -eps); }

struct Point {
    double x, y, z;

    Point() = default;

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

    bool operator==(const Point &a) const { return sign(a.x - x) == 0 && sign(a.y - y) == 0 && sign(a.z - z) == 0; }

    bool operator!=(const Point &a) const { return !(*this == a); }

    bool operator<(const Point &a) const {
        if (sign(a.x - x) == 0 && sign(a.y - y) == 0)
            return z < a.z;
        else if (sign(a.x - x) == 0)
            return y < a.y;
        return x < a.x;
    }

    bool operator>(const Point &a) const {
        if (sign(a.x - x) == 0 && sign(a.y - y) == 0)
            return z > a.z;
        else if (sign(a.x - x) == 0)
            return y > a.y;
        return x > a.x;
    }
};

istream &operator>>(istream &is, Point &p) {
    is >> p.x >> p.y >> p.z;
    return is;
}

ostream &operator<<(ostream &os, Point &p) {
    os << p.x << " " << p.y << " " << p.z;
    return os;
}

double dot(const Point& a, const Point& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Point cross(const Point& a, const Point& b) {
    return Point(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

// returns true if  point p is on line segment ab
bool is_point_on_seg(const Point& a, const Point& b, const Point& p) {
    if (fabs(cross(p - b, a - b).z) < eps) {
        if (p.x < min(a.x, b.x) || p.x > max(a.x, b.x)) return false;
        if (p.y < min(a.y, b.y) || p.y > max(a.y, b.y)) return false;
        return true;
    }
    return false;
}

// intersection point between segment ab and segment cd assuming unique intersection exists
bool seg_seg_intersection(const Point& a, const Point& b, const Point& c, const Point& d, Point &ans) {
    double oa = cross(d - c, a - c).z, ob = cross(d - c, b - c).z;
    double oc = cross(b - a, c - a).z, od = cross(b - a, d - a).z;
    if (oa * ob < 0 && oc * od < 0) {
        ans = (a * ob - b * oa) / (ob - oa);
        return true;
    } else return false;
}

// intersection point between segment ab and segment cd assuming unique intersection may not exists
// se.size()==0 means no intersection
// se.size()==1 means one intersection
// se.size()==2 means range intersection
set<Point> seg_seg_intersection_inside(const Point& a, const Point& b, const Point& c, const Point& d) {
    Point ans;
    if (seg_seg_intersection(a, b, c, d, ans)) return {ans};
    set<Point> se;
    if (is_point_on_seg(c, d, a)) se.insert(a);
    if (is_point_on_seg(c, d, b)) se.insert(b);
    if (is_point_on_seg(a, b, c)) se.insert(c);
    if (is_point_on_seg(a, b, d)) se.insert(d);
    return se;
}

double to_rad(double ang) {
    return ang * pi / 180.00;
}

struct Matrix {
    vector<vector<double>> a;

    Matrix() = default;

    Matrix(vector<vector<double>> &b) {
        a = b;
    }

    Matrix multiply(Matrix m1) {
        vector<vector<double>> &b = m1.a;
        assert(a[0].size() == b.size());
        int n = a.size(), m = a[0].size(), l = b[0].size();
        vector<vector<double> > ret(n, vector<double>(l));
        for (int i = 0; i < n; i++) {
            for (int k = 0; k < m; k++) {
                for (int j = 0; j < l; j++) {
                    ret[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        return {ret};
    }

    void normalize() {
        //Normalize a 4*1 Point
        assert(a[0].size() == 1);
        for (int i = 0; i < N; i++) {
            a[i][0] /= a[N - 1][0];
        }
    }
};

Matrix getIdentityMatrix() {
    vector<vector<double>> a(N, vector<double>(N, 0));
    for (int i = 0; i < N; i++) {
        a[i][i] = 1;
    }
    return {a};
}

Matrix getTranslationMatrix(const Point& t) {
    vector<vector<double>> a(N, vector<double>(N, 0));
    for (int i = 0; i < N; i++) {
        a[i][i] = 1;
    }
    a[0][N - 1] = t.x;
    a[1][N - 1] = t.y;
    a[2][N - 1] = t.z;
    return {a};
}

Matrix getScalingMatrix(const Point& s) {
    vector<vector<double>> a(N, vector<double>(N, 0));
    a[0][0] = s.x;
    a[1][1] = s.y;
    a[2][2] = s.z;
    a[3][3] = 1;
    return {a};
}

Point Rodrigues(const Point& x, const Point& a, double angle) {
    return x * cos(to_rad(angle)) + a * (dot(a, x)) * (1 - cos(to_rad(angle))) + cross(a, x) * sin(to_rad(angle));
}

Matrix getRotationMatrix(double angle, Point a) {
    a.normalize();
    Point c1 = Rodrigues(Point(1.0, 0.0, 0.0), a, angle);
    Point c2 = Rodrigues(Point(0.0, 1.0, 0.0), a, angle);
    Point c3 = Rodrigues(Point(0.0, 0.0, 1.0), a, angle);

    vector<vector<double>> matrix(N, vector<double>(N, 0));
    matrix[0][0] = c1.x;
    matrix[1][0] = c1.y;
    matrix[2][0] = c1.z;

    matrix[0][1] = c2.x;
    matrix[1][1] = c2.y;
    matrix[2][1] = c2.z;

    matrix[0][2] = c3.x;
    matrix[1][2] = c3.y;
    matrix[2][2] = c3.z;

    matrix[3][3] = 1;
    return {matrix};
}

Matrix getViewTransformationMatrix(const Point& look, const Point& eye, const Point& up) {
    Point l = look - eye;
    l.normalize();
    Point r = cross(l, up);
    r.normalize();
    Point u = cross(r, l);
    u.normalize();

    Matrix T = getTranslationMatrix(eye * -1);
    vector<vector<double>> rot(N, vector<double>(N, 0));
    rot[0][0] = r.x;
    rot[0][1] = r.y;
    rot[0][2] = r.z;

    rot[1][0] = u.x;
    rot[1][1] = u.y;
    rot[1][2] = u.z;

    rot[2][0] = -l.x;
    rot[2][1] = -l.y;
    rot[2][2] = -l.z;

    rot[3][3] = 1;
    Matrix R(rot);
    Matrix V = R.multiply(T);
    return V;
}

Matrix getProjectionMatrix(double fovY, double aspectRatio, double near, double far) {
    double fovX = fovY * aspectRatio;
    double t = near * tan(to_rad(fovY / 2.0));
    double r = near * tan(to_rad(fovX / 2.0));
    vector<vector<double>> a(N, vector<double>(N, 0));
    a[0][0] = near / r;
    a[1][1] = near / t;
    a[2][2] = -(far + near) / (far - near);
    a[2][3] = -(2 * far * near) / (far - near);
    a[3][2] = -1;
    return {a};
}

struct Stack {
    stack<Matrix> st;
    Matrix cur_transformation;

    Stack() {
        cur_transformation = getIdentityMatrix();
    }

    void push() {
        st.push(cur_transformation);
    }

    void pop() {
        if (st.empty()) {
            cout << "Error! Stack is empty!";
            return;
        }
        cur_transformation = st.top();
        st.pop();
    }

    void applyTranslation(const Point& t) {
        Matrix now = getTranslationMatrix(t);
        cur_transformation = cur_transformation.multiply(now);
    }

    void applyScale(const Point& s) {
        Matrix now = getScalingMatrix(s);
        cur_transformation = cur_transformation.multiply(now);
    }

    void applyRotation(double angle, const Point& a) {
        Matrix now = getRotationMatrix(angle, a);
        cur_transformation = cur_transformation.multiply(now);
    }
};

Point applyTransformation(Matrix cur_transformation, const Point& p) {
    vector<vector<double>> a(N, vector<double>(1, 0));
    a[0][0] = p.x;
    a[1][0] = p.y;
    a[2][0] = p.z;
    a[3][0] = 1;
    Matrix m = cur_transformation.multiply({a});
    m.normalize();
    return {m.a[0][0], m.a[1][0], m.a[2][0]};
}

Point eye, look, up;
double fovY, aspectRatio, near, far;
ifstream scene;
Stack st;
vector<vector<Point>> triangles;

void Stage1() {
    ofstream stage1("stage1.txt");
    stage1 << fixed << setprecision(7);

    while (true) {
        string command;
        scene >> command;
        if (command == "triangle") {
            vector<Point> p(3);
            scene >> p[0] >> p[1] >> p[2];
            for (int i = 0; i < 3; i++) {
                p[i] = applyTransformation(st.cur_transformation, p[i]);
                stage1 << p[i] << "\n";
            }
            stage1 << "\n";
            triangles.push_back(p);
        } else if (command == "translate") {
            Point t;
            scene >> t;
            st.applyTranslation(t);
        } else if (command == "scale") {
            Point s;
            scene >> s;
            st.applyScale(s);
        } else if (command == "rotate") {
            double angle;
            Point a;
            scene >> angle >> a;
            st.applyRotation(angle, a);
        } else if (command == "push") {
            st.push();
        } else if (command == "pop") {
            st.pop();
        } else if (command == "end") {
            break;
        }
    }
    stage1.close();
}

void Stage2() {
    ofstream stage2("stage2.txt");
    stage2 << fixed << setprecision(7);
    Matrix V = getViewTransformationMatrix(look, eye, up);
    for (auto &triangle: triangles) {
        for (auto &j: triangle) {
            j = applyTransformation(V, j);
            stage2 << j << "\n";
        }
        stage2 << "\n";
    }
    stage2.close();
}

void Stage3() {
    ofstream stage3("stage3.txt");
    stage3 << fixed << setprecision(7);
    Matrix P = getProjectionMatrix(fovY, aspectRatio, near, far);
    for (auto &triangle: triangles) {
        for (auto &j: triangle) {
            j = applyTransformation(P, j);
            stage3 << j << "\n";
        }
        stage3 << "\n";
    }
    stage3.close();
}

struct Color {
    int r, g, b;
};

Color getColor() {
    vector<int> color(3);
    for (int i = 0; i < 3; i++) {
        color[i] = getrand(0, 255);
    }
    if (color[0] == 0 && color[1] == 0 && color[2] == 0) {
        color[2] = getrand(1, 255);
    }
    return {color[0], color[1], color[2]};
}

struct Triangle {
    vector<Point> points;
    Color color;

    Triangle() {
        points.resize(3);
        color = getColor();
    }
};


void Stage4() {
    // Read data
    ifstream config("config.txt");
    int width, height;
    double maxx, minx, maxy, miny, maxz, minz;
    config >> width >> height >> minx >> miny >> minz >> maxz;
    maxx = -minx;
    maxy = -miny;
    vector<Triangle> a(triangles.size());
    for (int i = 0; i < triangles.size(); i++) {
        a[i].points = triangles[i];
    }
    config.close();

    //Initialize Z-buffer and Frame buffer
    double dx, dy, top_y, left_x;
    dx = (maxx - minx) / width;
    dy = (maxy - miny) / height;
    top_y = maxy - dy / 2;
    left_x = minx + dx / 2;
    vector<vector<double>> zbuffer(height, vector<double>(width, maxz));
    vector<vector<Color>> framebuffer(height, vector<Color>(width, {0, 0, 0}));

    //Apply procedure
    for (const Triangle &it: a) {
        Triangle now = it;
        double low_y = min({now.points[0].y, now.points[1].y, now.points[2].y});
        double high_y = max({now.points[0].y, now.points[1].y, now.points[2].y});
        int top_scanline = max(0, (int) round((top_y - high_y) / dy));
        int bottom_scanline = min(height - 1, (int) round((top_y - low_y) / dy));

        for (int row = top_scanline; row <= bottom_scanline; row++) {
            double low_x = 1e9, high_x = -1e9;
            set<Point> st;
            //scanline segment
            Point p1(-1e9, top_y - row * dy, 0);
            Point p2(1e9, top_y - row * dy, 0);

            for (int i = 0; i < 3; i++) {
                Point p3(now.points[i].x, now.points[i].y, 0);
                Point p4(now.points[(i + 1) % 3].x, now.points[(i + 1) % 3].y, 0);
                set<Point> cur = seg_seg_intersection_inside(p1, p2, p3, p4);
                for (const Point &pt: cur)
                    st.insert(pt);
            }
            if (st.empty())
                continue;
            for (const Point& pt: st) {
                low_x = min(low_x, pt.x);
                high_x = max(high_x, pt.x);
            }
            int col1 = max(0, (int) round((low_x - left_x) / dx));
            int col2 = min(width - 1, (int) round((high_x - left_x) / dx));
            if (col1 > col2)
                continue;
            double za = -1e9, zb = -1e9;
            bool f1 = false, f2 = false;
            for (int i = 0; i < 3; i++) {
                Point p3(now.points[i].x, now.points[i].y, 0);
                Point p4(now.points[(i + 1) % 3].x, now.points[(i + 1) % 3].y, 0);
                set<Point> cur = seg_seg_intersection_inside(p1, p2, p3, p4);
                if (cur.size() != 1)
                    continue;
                Point point = *cur.begin();
                p3.z = now.points[i].z;
                p4.z = now.points[(i + 1) % 3].z;
                if (p3.y < p4.y)
                    swap(p3, p4);
                if (sign(point.x - low_x) == 0) {
                    za = p3.z - (p3.z - p4.z) * ((p3.y - p1.y) / (p3.y - p4.y));
                    f1 = true;
                }
                if (sign(point.x - high_x) == 0) {
                    zb = p3.z - (p3.z - p4.z) * ((p3.y - p1.y) / (p3.y - p4.y));
                    f2 = true;
                }
            }
            if (!f1 || !f2) {
                continue;
            }
            for (int col = col1; col <= col2; col++) {
                double xp = left_x + col * dx;
                double zp = zb - (zb - za) * ((high_x - xp) / (high_x - low_x));
                if (zbuffer[row][col] > zp) {
                    zbuffer[row][col] = zp;
                    framebuffer[row][col] = now.color;
                }
            }

        }
    }

    //Printing Z buffer
    ofstream z_buffer("z_buffer.txt");
    z_buffer << fixed << setprecision(6);

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (zbuffer[i][j] < maxz) {
                z_buffer << zbuffer[i][j] << "\t";
            }
        }
        z_buffer << "\n";
    }
    z_buffer.close();
    //saving image
    bitmap_image bitmapImage(width, height);

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            bitmapImage.set_pixel(j, i, framebuffer[i][j].r, framebuffer[i][j].g, framebuffer[i][j].b);
        }
    }
    bitmapImage.save_image("out.bmp");

}

int main() {
    scene.open("scene.txt");
    scene >> eye >> look >> up >> fovY >> aspectRatio >> near >> far;
    Stage1();
    scene.close();
    Stage2();
    Stage3();
    Stage4();
    return 0;
}
