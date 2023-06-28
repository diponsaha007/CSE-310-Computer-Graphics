#include <utility>

using namespace std;
const double pi = 2 * acos(0.0);
const double eps = 1e-8;
const double inf = 1e18;
const int floor_width = 500;
const int tile_width = 10;
const double floor_ambient = 0.35;
const double floor_diffuse = 0.35;
const double floor_specular = 0.35;
const int floor_shine = 30;
const double floor_recursive = 0.35;

int sign(double x) { return (x > eps) - (x < -eps); }

struct Point {
    double x, y, z;

    Point() {
        x = 0;
        y = 0;
        z = 0;
    }

    Point(double x, double y, double z) : x(x), y(y), z(z) {}

    Point(const Point &p) {
        x = p.x;
        y = p.y;
        z = p.z;
    }

    double len() const {
        return sqrt(x * x + y * y + z * z);
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

double dot(const Point &a, const Point &b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Point cross(const Point &a, const Point &b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

double to_deg(double ang) {
    return ang * 180.00 / pi;
}

double dist(const Point &a, const Point &b) {
    return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

Point translate(const Point &a, const Point &b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

double angle(const Point &a, const Point &b) // in radian
{
    if (a.len() == 0 || b.len() == 0)
        return inf;
    double ang = dot(a, b) / a.len() / b.len();
    ang = acos(ang);
    return ang;
}

typedef Point Vector;

struct Camera {
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

    void move(const Point &p, bool pos_dir) {
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
};

struct Color {
    array<double, 3> rgb;

    Color() { rgb = {0, 0, 0}; }

    Color(double r, double g, double b) {
        rgb = {r, g, b};
    }

    void roundValues() {
        for (int i = 0; i < 3; i++) {
            if (rgb[i] < 0)
                rgb[i] = 0;
            if (rgb[i] > 1)
                rgb[i] = 1;
        }
    }

    friend istream &operator>>(istream &is, Color &c) {
        is >> c.rgb[0] >> c.rgb[1] >> c.rgb[2];
        return is;
    }
};


void drawSphere(double radius, int slices, int stacks) {
    vector<vector<Point>> points(stacks + 1, vector<Point>(slices + 1));
    for (int i = 0; i <= stacks; i++) {
        double h = radius * sin(((double) i / (double) stacks) * (pi / 2));
        double r = radius * cos(((double) i / (double) stacks) * (pi / 2));
        for (int j = 0; j <= slices; j++) {
            points[i][j].x = r * cos(((double) j / (double) slices) * 2 * pi);
            points[i][j].y = r * sin(((double) j / (double) slices) * 2 * pi);
            points[i][j].z = h;
        }
    }

    for (int i = 0; i < stacks; i++) {
        for (int j = 0; j < slices; j++) {
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

struct PointLight {
    Point light_pos;
    Color color;

    PointLight() = default;

    PointLight(const Point &lightPos, const Color &color) : light_pos(lightPos), color(color) {}

    void draw() {
        glPushMatrix();
        glColor3f(color.rgb[0], color.rgb[1], color.rgb[2]);
        glTranslatef(light_pos.x, light_pos.y, light_pos.z);
        drawSphere(1, 20, 20);
        glPopMatrix();
    }

    friend istream &operator>>(istream &is, PointLight &pl) {
        is >> pl.light_pos >> pl.color;
        return is;
    }
};

struct SpotLight {
    PointLight point_light;
    Vector light_direction;
    double cutoff_angle;

    SpotLight() = default;

    SpotLight(PointLight pointLight, const Vector &lightDirection, double cutoffAngle) : point_light(
            std::move(pointLight)),
                                                                                         light_direction(
                                                                                                 lightDirection),
                                                                                         cutoff_angle(
                                                                                                 cutoffAngle) {}

    void draw() {
        glPushMatrix();
        glColor3f(point_light.color.rgb[0], point_light.color.rgb[1], point_light.color.rgb[2]);
        glTranslatef(point_light.light_pos.x, point_light.light_pos.y, point_light.light_pos.z);
        drawSphere(3, 20, 20);
        glPopMatrix();
    }

    friend istream &operator>>(istream &is, SpotLight &sl) {
        is >> sl.point_light >> sl.light_direction >> sl.cutoff_angle;
        return is;
    }
};

struct Ray {
    Point p;
    Vector v;

    Ray() = default;

    Ray(const Point &p, const Vector &v) : p(p), v(v) {
        this->v.normalize();
    }
};

struct Object {
    Color color;
    double ambient, diffuse, specular, recursive;
    int shine;

    virtual void draw() {
    }

    //if it doesn't intersect, return -1
    virtual double intersect(Ray ray) {
        return -1;
    }

    virtual Vector getNormal(Point intersectionPoint, Ray ray) {
        return {0, 0, 0};
    }

    virtual Color getColor(Point intersectionPoint) {
        return color;
    }

    //Takes intersection point , ray and object color currently at intersection point
    // and returns the color of the intersection point
    Color getIntersectionColor(Point intersectionPoint, Ray ray, Color color1);

    Color recursiveIntersection(int index, int level, Ray ray);
};

vector<Object *> objects;
vector<PointLight> plights;
vector<SpotLight> slights;
Camera camera;


Color Object::getIntersectionColor(Point intersectionPoint, Ray ray, Color color1) {
    //ambient
    Color ret = color1;
    for (int i = 0; i < 3; i++) {
        ret.rgb[i] *= ambient;
    }
    //diffuse and specular
    for (int i = 0; i < plights.size(); i++) {
        Ray r1 = Ray(plights[i].light_pos, intersectionPoint - plights[i].light_pos);
        double t_min = inf;
        //check if r1 is obscured by any object
        for (int j = 0; j < objects.size(); j++) {
            double t = objects[j]->intersect(r1);
            if (t > 0 && t < t_min)
                t_min = t;
        }
        if (t_min <= 0 || t_min >= inf) {
            continue;
        }
        Point intersectionPoint1 = r1.p + r1.v * t_min;
        if (intersectionPoint1 == intersectionPoint) {
            Vector normal = getNormal(intersectionPoint, ray);
            Vector light_dir = plights[i].light_pos - intersectionPoint;
            light_dir.normalize();
            double diffuse_intensity = dot(light_dir, normal);
            if (diffuse_intensity > 0) {
                for (int j = 0; j < 3; j++) {
                    ret.rgb[j] += plights[i].color.rgb[j] * diffuse * diffuse_intensity * color1.rgb[j];
                }
            }
            Vector r2 = r1.v - normal * 2 * dot(r1.v, normal);
            r2.normalize();
            Vector r3 = Point(0, 0, 0) - ray.v;
            r3.normalize();
            double specular_intensity = pow(dot(r2, r3), shine);
            if (specular_intensity > 0) {
                for (int j = 0; j < 3; j++) {
                    ret.rgb[j] += plights[i].color.rgb[j] * specular * specular_intensity * color1.rgb[j];
                }
            }
        }
    }

    for (int i = 0; i < slights.size(); i++) {
        Ray r1 = Ray(slights[i].point_light.light_pos, intersectionPoint - slights[i].point_light.light_pos);
        double ang = to_deg(angle(slights[i].light_direction, r1.v));
        if (ang > slights[i].cutoff_angle) {
            continue;
        }
        double t_min = inf;
        //check if r1 is obscured by any object
        for (int j = 0; j < objects.size(); j++) {
            double t = objects[j]->intersect(r1);
            if (t > 0 && t < t_min)
                t_min = t;
        }
        if (t_min <= 0 || t_min >= inf) {
            continue;
        }

        Point intersectionPoint1 = r1.p + r1.v * t_min;
        if (intersectionPoint1 == intersectionPoint) {

            Vector normal = getNormal(intersectionPoint, ray);
            Vector light_dir = slights[i].point_light.light_pos - intersectionPoint;
            light_dir.normalize();
            double diffuse_intensity = dot(light_dir, normal);
            if (diffuse_intensity > 0) {
                for (int j = 0; j < 3; j++) {
                    ret.rgb[j] += slights[i].point_light.color.rgb[j] * diffuse * diffuse_intensity * color1.rgb[j];
                }
            }
            Vector r2 = r1.v - normal * 2 * dot(r1.v, normal);
            r2.normalize();
            Vector r3 = Point(0, 0, 0) - ray.v;
            r3.normalize();
            double specular_intensity = pow(dot(r2, r3), shine);
            if (specular_intensity > 0) {
                for (int j = 0; j < 3; j++) {
                    ret.rgb[j] += slights[i].point_light.color.rgb[j] * specular * specular_intensity * color1.rgb[j];
                }
            }
        }
    }
    return ret;

}


Color Object::recursiveIntersection(int index, int level, Ray ray) {
    double t = objects[index]->intersect(ray);
    Color ret = Color(0, 0, 0);
    Point intersectionPoint;
    intersectionPoint = ray.p + ray.v * t;
    ret = objects[index]->getColor(intersectionPoint);
    ret = objects[index]->getIntersectionColor(intersectionPoint, ray, ret);
    ret.roundValues();


    if (level <= 0) {
        return ret;
    }
    //get reflected ray
    Vector normal = objects[index]->getNormal(intersectionPoint, ray);
    Vector r = ray.v - normal * 2 * dot(ray.v, normal);
    r.normalize();
    Ray r1 = Ray(intersectionPoint + r * 0.1, r);

    double tMin = inf;
    int nearest = -1;
    for (int i = 0; i < objects.size(); i++) {
        double now = objects[i]->intersect(r1);
        if (now > 0 && now < inf) {
            if (now < tMin) {
                tMin = now;
                nearest = i;
            }
        }
    }
    if (nearest == -1) {
        return ret;
    }
    Color color2 = recursiveIntersection(nearest, level - 1, r1);
    for (int i = 0; i < 3; i++) {
        ret.rgb[i] += color2.rgb[i] * recursive;
    }
    ret.roundValues();
    return ret;
}


struct Sphere : public Object {
    Point center;
    double radius;

    Sphere() = default;

    Sphere(const Point &center, double radius, const Color &color, double ambient, double diffuse, double specular,
           int shine, double recursive) : center(center), radius(radius) {
        this->color = color;
        this->ambient = ambient;
        this->diffuse = diffuse;
        this->specular = specular;
        this->shine = shine;
        this->recursive = recursive;
    }

    void draw() override {
        glPushMatrix();
        glColor3f(color.rgb[0], color.rgb[1], color.rgb[2]);
        glTranslatef(center.x, center.y, center.z);
        drawSphere(radius, 20, 20);
        glPopMatrix();
    }


    double intersect(Ray ray) override {
        Ray r = ray;
        r.p = translate(r.p, Point(0, 0, 0) - center);
        double a = dot(r.v, r.v);
        double b = 2 * dot(r.v, r.p);
        double c = dot(r.p, r.p) - radius * radius;
        double d = b * b - 4 * a * c;
        if (d < 0) return -1;
        double t1 = (-b + sqrt(d)) / (2 * a);
        double t2 = (-b - sqrt(d)) / (2 * a);
        double t = min(t1, t2);
        if (t < 0)
            t = max(t1, t2);
        return t;
    }

    Vector getNormal(Point intersectionPoint, Ray ray) override {
        Vector normal = intersectionPoint - center;
        normal.normalize();
        //camera is inside the sphere
        if (sign(dist(camera.pos, center) - radius) <= 0)
            normal = Point(0, 0, 0) - normal;
        return normal;
    }

    friend istream &operator>>(istream &in, Sphere &s) {
        in >> s.center >> s.radius >> s.color >> s.ambient >> s.diffuse >> s.specular >> s.recursive >> s.shine;
        return in;
    }
};

struct Floor : public Object {
    int floor_width;
    int tile_width;

    Floor() = default;

    Floor(int floor_width, int tile_width, const Color &color, double ambient, double diffuse, double specular,
          int shine, double recursive) : floor_width(floor_width), tile_width(tile_width) {
        this->color = color;
        this->ambient = ambient;
        this->diffuse = diffuse;
        this->specular = specular;
        this->shine = shine;
        this->recursive = recursive;
    }

    void draw() override {
        int lim = floor_width / tile_width;
        for (int row = 0; row < lim; row++) {
            for (int col = 0; col < lim; col++) {
                Color now = color;
                if ((row + col) % 2)
                    now = Color(1 - now.rgb[0], 1 - now.rgb[1], 1 - now.rgb[2]);
                glColor3f(now.rgb[0], now.rgb[1], now.rgb[2]);
                Point p = Point(-floor_width / 2.00 + row * tile_width, -floor_width / 2.00 + col * tile_width, 0);
                glBegin(GL_QUADS);
                glVertex3f(p.x, p.y, p.z);
                glVertex3f(p.x + tile_width, p.y, p.z);
                glVertex3f(p.x + tile_width, p.y + tile_width, p.z);
                glVertex3f(p.x, p.y + tile_width, p.z);
                glEnd();
            }
        }
    }

    bool checkInside(Point p) const {
        return (p.x >= -floor_width / 2.00 && p.x <= floor_width / 2.00 && p.y >= -floor_width / 2.00 &&
                p.y <= floor_width / 2.00);
    }

    double intersect(Ray ray) override {
        if (ray.v.z == 0)
            return -1;
        Point normal = getNormal(Point(0, 0, 0), Ray());
        double t = -(dot(normal, ray.p)) / dot(normal, ray.v);
        if (t < 0 || !checkInside(ray.p + ray.v * t))
            return -1;
        return t;
    }

    Vector getNormal(Point intersectionPoint, Ray ray) override {
        Vector v1 = camera.pos - intersectionPoint;
        if (dot(v1, Vector(0, 0, 1)) > 0)
            return {0, 0, 1};
        else
            return {0, 0, -1};
    }

    Color getColor(Point intersectionPoint) override {
        int i = (intersectionPoint.x + floor_width / 2.00) / tile_width;
        int j = (intersectionPoint.y + floor_width / 2.00) / tile_width;
        Color now = color;
        if ((i + j) % 2)
            now = Color(1 - now.rgb[0], 1 - now.rgb[1], 1 - now.rgb[2]);
        return now;
    }
};

void InsertFloor() {
    Object *floor = new Floor(floor_width, tile_width, Color(0, 0, 0), floor_ambient, floor_diffuse, floor_specular,
                              floor_shine, floor_recursive);
    objects.push_back(floor);
}

struct Triangle : public Object {
    Point p1, p2, p3;

    Triangle() = default;

    Triangle(const Point &p1, const Point &p2, const Point &p3, const Color &color, double ambient, double diffuse,
             double specular, int shine, double recursive) : p1(p1), p2(p2), p3(p3) {
        this->color = color;
        this->ambient = ambient;
        this->diffuse = diffuse;
        this->specular = specular;
        this->shine = shine;
        this->recursive = recursive;
    }

    void draw() override {
        glBegin(GL_TRIANGLES);
        glColor3f(color.rgb[0], color.rgb[1], color.rgb[2]);
        glVertex3f(p1.x, p1.y, p1.z);
        glVertex3f(p2.x, p2.y, p2.z);
        glVertex3f(p3.x, p3.y, p3.z);
        glEnd();
    }

    double intersect(Ray ray) override {
        //Ray triangle intersection
        Vector e1 = p2 - p1;
        Vector e2 = p3 - p1;
        Vector h = cross(ray.v, e2);
        double a = dot(e1, h);
        if (a > -eps && a < eps)
            return -1;
        double f = 1.0 / a;
        Vector s = ray.p - p1;
        double u = f * dot(s, h);
        if (u < 0.0 || u > 1.0)
            return -1;
        Vector q = cross(s, e1);
        double v = f * dot(ray.v, q);
        if (v < 0.0 || u + v > 1.0)
            return -1;
        double t = f * dot(e2, q);
        if (t > eps)
            return t;
        else
            return -1;

    }

    Vector getNormal(Point intersectionPoint, Ray ray) override {
        Vector e1 = p2 - p1;
        Vector e2 = p3 - p1;
        Vector normal = cross(e1, e2);
        normal.normalize();
        if (dot(normal, ray.v) > 0)
            normal = Point(0, 0, 0) - normal;
        return normal;
    }

    //overload the >> operator for input stream
    friend istream &operator>>(istream &in, Triangle &s) {
        in >> s.p1 >> s.p2 >> s.p3 >> s.color >> s.ambient >> s.diffuse >> s.specular >> s.recursive >> s.shine;
        return in;
    }
};

struct GeneralQuadricSurfaces : public Object {
    double a, b, c, d, e, f, g, h, i, j;
    Point cube_reference_point;
    double length, width, height;

    GeneralQuadricSurfaces() = default;

    GeneralQuadricSurfaces(double a, double b, double c, double d, double e, double f, double g, double h, double i,
                           double j, const Point &cube_reference_point, double length, double width, double height,
                           const Color &color, double ambient, double diffuse, double specular, int shine,
                           double recursive) : a(a), b(b), c(c), d(d), e(e), f(f), g(g), h(h), i(i), j(j),
                                               cube_reference_point(cube_reference_point), length(length), width(width),
                                               height(height) {
        this->color = color;
        this->ambient = ambient;
        this->diffuse = diffuse;
        this->specular = specular;
        this->shine = shine;
        this->recursive = recursive;
    }

    bool safe(const Point &p) const {
        bool flag = true;
        if (!(length == 0 || (p.x >= cube_reference_point.x && p.x <= cube_reference_point.x + length))) {
            flag = false;
        }
        if (!(width == 0 || (p.y >= cube_reference_point.y && p.y <= cube_reference_point.y + width))) {
            flag = false;
        }
        if (!(height == 0 || (p.z >= cube_reference_point.z && p.z <= cube_reference_point.z + height))) {
            flag = false;
        }
        return flag;
    }

    double intersect(Ray ray) override {
        double aq = a * ray.v.x * ray.v.x + b * ray.v.y * ray.v.y + c * ray.v.z * ray.v.z + d * ray.v.x * ray.v.y +
                    e * ray.v.x * ray.v.z + f * ray.v.y * ray.v.z;
        double bq = 2 * (a * ray.v.x * ray.p.x + b * ray.v.y * ray.p.y + c * ray.v.z * ray.p.z)
                    + d * (ray.v.x * ray.p.y + ray.v.y * ray.p.x) + e * (ray.v.x * ray.p.z + ray.v.z * ray.p.x) +
                    f * (ray.v.y * ray.p.z + ray.v.z * ray.p.y) + g * ray.v.x + h * ray.v.y + i * ray.v.z;
        double cq = a * ray.p.x * ray.p.x + b * ray.p.y * ray.p.y + c * ray.p.z * ray.p.z + d * ray.p.x * ray.p.y +
                    e * ray.p.x * ray.p.z + f * ray.p.y * ray.p.z + g * ray.p.x + h * ray.p.y + i * ray.p.z + j;

        double t;
        if (aq == 0 && bq == 0)
            return -1;
        if (aq == 0) {
            t = -cq / bq;
            if (t <= 0 || !safe(ray.p + ray.v * t)) {
                t = -1;
            }
        } else {
            double discriminant = bq * bq - 4 * aq * cq;
            if (discriminant < 0)
                return -1;

            t = (-bq - sqrt(discriminant)) / (2 * aq);
            if (t <= 0 || !safe(ray.p + ray.v * t)) {
                t = (-bq + sqrt(discriminant)) / (2 * aq);

                if (t <= 0 || !safe(ray.p + ray.v * t)) {
                    t = -1;
                }
            }
        }
        if (t <= 0)
            return -1;

        return t;
    }

    Vector getNormal(Point intersectionPoint, Ray ray) override {
        double xn = 2 * a * intersectionPoint.x + d * intersectionPoint.y + e * intersectionPoint.z + g;
        double yn = 2 * b * intersectionPoint.y + d * intersectionPoint.x + f * intersectionPoint.z + h;
        double zn = 2 * c * intersectionPoint.z + e * intersectionPoint.x + f * intersectionPoint.y + i;
        Vector normal = Vector(xn, yn, zn);
        normal.normalize();
        if (dot(ray.v * -1, normal) <= 0)
            normal = normal * (-1.0);
        return normal;
    }

    friend istream &operator>>(istream &in, GeneralQuadricSurfaces &s) {
        in >> s.a >> s.b >> s.c >> s.d >> s.e >> s.f >> s.g >> s.h >> s.i >> s.j >> s.cube_reference_point >> s.length
           >>
           s.width >> s.height >> s.color >> s.ambient >> s.diffuse >> s.specular >> s.recursive >> s.shine;
        return in;
    }
};
