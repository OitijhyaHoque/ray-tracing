#ifndef OBJECTS_H
#define OBJECTS_H

#ifdef __linux__
#include <GL/glut.h> // For Linux systems
#elif defined(_WIN32) || defined(WIN32)
#include <windows.h>
#include <GL/glut.h> // For Windows systems
#elif defined(__APPLE__)
#include <GLUT/glut.h> // For macOS systems
#else
#include <GL/glut.h> // Default fallback
#endif

#include "bitmap_image.hpp" // For bitmap image handling

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

#define EPSILON 1e-6

// Optimized Vector struct using fixed-size arrays
struct Vector
{
    double vec[3];
    
    Vector() { vec[0] = vec[1] = vec[2] = 0.0; }
    Vector(double x, double y, double z) { vec[0] = x; vec[1] = y; vec[2] = z; }
    
    inline void normalize()
    {
        double normValue = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
        if (normValue > EPSILON) {
            double inv = 1.0 / normValue;
            vec[0] *= inv;
            vec[1] *= inv;
            vec[2] *= inv;
        }
    }
    
    inline double norm() const
    {
        return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    }
    
    inline double normSquared() const
    {
        return vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
    }
};

// Inline optimized vector operations
inline double dotProduct(const Vector &a, const Vector &b)
{
    return a.vec[0] * b.vec[0] + a.vec[1] * b.vec[1] + a.vec[2] * b.vec[2];
}

inline Vector crossProduct(const Vector &a, const Vector &b)
{
    return Vector(
        a.vec[1] * b.vec[2] - a.vec[2] * b.vec[1],
        a.vec[2] * b.vec[0] - a.vec[0] * b.vec[2],
        a.vec[0] * b.vec[1] - a.vec[1] * b.vec[0]
    );
}

inline Vector add(const Vector &a, const Vector &b)
{
    return Vector(a.vec[0] + b.vec[0], a.vec[1] + b.vec[1], a.vec[2] + b.vec[2]);
}

inline Vector subtract(const Vector &a, const Vector &b)
{
    return Vector(a.vec[0] - b.vec[0], a.vec[1] - b.vec[1], a.vec[2] - b.vec[2]);
}

inline Vector multiply(double scalar, const Vector &v)
{
    return Vector(scalar * v.vec[0], scalar * v.vec[1], scalar * v.vec[2]);
}

inline Vector multiply(const Vector &v, double scalar)
{
    return Vector(scalar * v.vec[0], scalar * v.vec[1], scalar * v.vec[2]);
}

Vector rotation(const Vector &target, const Vector &axis, double angle);

class Ray
{
public:
    Vector start; // Starting point of the ray
    Vector dir;   // Direction of the ray

    Ray(const Vector &start, const Vector &dir);
    // void normalize();
    // double intersect(Object *obj, double *current_color, int level);
    // Vector getNormal(Object *obj, const Vector &point);
};

class Object
{
public:
    Vector reference_point;       // Reference point for the object
    double height, width, length; // Dimensions of the object
    Vector color;
    double ambient, diffuse, specular, reflection;
    int shine;

    Object();
    virtual ~Object() = default;

    virtual void draw() = 0;

    void setColor(double c[3]);
    void setShine(int s);
    void setCoefficients(double ambient, double diffuse, double specular, double reflection);

    double intersect(const Ray &ray, double *current_color, int level);
    virtual double solveIntersection(const Ray &ray) = 0;
    virtual Vector getNormal(const Vector &point, const Vector &direction = Vector()) = 0;
    virtual Vector getColorAt(const Vector &point);
};

class Sphere : public Object
{
public:
    Sphere(const Vector &center, double radius);
    void draw() override;
    double solveIntersection(const Ray &ray) override;
    Vector getNormal(const Vector &point, const Vector &direction = Vector()) override;
};

class Triangle : public Object
{
public:
    Vector p1, p2, p3;

    Triangle(const Vector &p1, const Vector &p2, const Vector &p3);
    void draw() override;
    double solveIntersection(const Ray &ray) override;
    Vector getNormal(const Vector &point, const Vector &direction = Vector()) override;
    Vector getColorAt(const Vector &point) override;
};

class GeneralQuadratic : public Object
{
public:
    double coefficients[10];
    GeneralQuadratic(double coeffs[10]);
    void draw() override;
    double solveIntersection(const Ray &ray) override;
    Vector getNormal(const Vector &point, const Vector &direction = Vector()) override;
    Vector getColorAt(const Vector &point) override;
    bool isPointInBoundingBox(const Vector &point) const;
};

class Floor : public Object
{
public:
    double floorWidth;
    double tileWidth;

    double primary[3];
    double secondary[3];

    std::string texture;

    unsigned char *textureData;
    int textureWidth;
    int textureHeight;
    int textureChannels;
    bool textureLoaded;

    double leftX, leftY;

    Floor(double floorWidth, double tileWidth);
    ~Floor();
    void draw() override;
    double solveIntersection(const Ray &ray) override;
    Vector getNormal(const Vector &point, const Vector &direction = Vector()) override;
    Vector getColorAt(const Vector &point) override;

    bool loadTexture(const std::string &filename);
    void freeTexture();
};

class PointLight
{
public:
    Vector position; // Position of the point light
    Vector color;    // Color of the light

    PointLight(const Vector &pos, const double col[3]);
    void draw();
};

class SpotLight
{
public:
    PointLight pointLight; // Point light properties
    Vector direction;      // Direction of the spotlight
    double angle;          // Cutoff angle

    SpotLight(const PointLight &pl, const Vector &dir, double ang);
    void draw();
};

class Camera
{
public:
    Vector eye;    // Camera position
    Vector center; // Look-at point
    Vector up;     // Up vector

    double movementSensitivity; // Movement sensitivity
    double rotationSensitivity; // Rotation sensitivity

    Camera(const Vector &eye_pos, const Vector &center_pos, const Vector &up_vec);

    void print();
    Vector getLookDirection();
    void moveFB(int dir);      // Move forward/backward
    void moveLR(int dir);      // Move left/right
    void moveUD(int dir);      // Move up/down
    void rotateYaw(int dir);   // Rotate around up vector
    void rotatePitch(int dir); // Rotate around right vector
    void rotateRoll(int dir);  // Rotate around look vector
};

struct Environment
{
    int recursionDepth;
    int pixelCount;
    std::vector<Object *> objects;
    std::vector<PointLight *> pointLights;
    std::vector<SpotLight *> spotLights;

    Environment();
    ~Environment();
};

void capture();

void toggleTexture();

void cleanupTextures();

#endif // OBJECTS_H
