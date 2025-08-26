#include "2005049_obj.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

extern Camera camera;
extern Environment environment;

bool applyTexture = false;
int textureId = 0;
int textureCount = 2;

Vector::Vector(int s) : size(s)
{
    vec = std::vector<double>(s, 0);
}

void Vector::normalize()
{
    double normValue = norm();

    if (normValue <= EPSILON)
        throw std::runtime_error("Cannot normalize a zero vector.");

    for (double &val : vec)
        val /= normValue;
}

double Vector::norm()
{
    double norm = 0;
    for (double val : vec)
        norm += val * val;
    return sqrt(norm);
}

double dotProduct(const Vector &a, const Vector &b)
{
    if (a.size != b.size)
        throw std::invalid_argument("Vectors must be of the same size.");

    double result = 0;
    for (int i = 0; i < a.size; i++)
    {
        result += a.vec[i] * b.vec[i];
    }
    return result;
}

Vector crossProduct(const Vector &a, const Vector &b)
{
    if (a.size != 3 || b.size != 3)
        throw std::invalid_argument("Cross product is only defined for 3D vectors.");

    Vector result(3);
    result.vec[0] = a.vec[1] * b.vec[2] - a.vec[2] * b.vec[1];
    result.vec[1] = a.vec[2] * b.vec[0] - a.vec[0] * b.vec[2];
    result.vec[2] = a.vec[0] * b.vec[1] - a.vec[1] * b.vec[0];

    return result;
}

Vector add(const Vector &a, const Vector &b)
{
    if (a.size != b.size)
        throw std::invalid_argument("Vectors must be of the same size.");

    Vector result(a.size);
    for (int i = 0; i < a.size; i++)
    {
        result.vec[i] = a.vec[i] + b.vec[i];
    }
    return result;
}

Vector multiply(double scalar, const Vector &v)
{
    Vector result(v.size);
    for (int i = 0; i < v.size; i++)
    {
        result.vec[i] = scalar * v.vec[i];
    }
    return result;
}

Vector rotation(const Vector &target, const Vector &axis, double angle)
{
    if (axis.size != 3 || target.size != 3)
        throw std::invalid_argument("Axis and target must be 3D");

    // assume normalized axis

    // Rodrigues' rotation formula
    // v_rotated = v * cos(angle) + (k x v) * sin(angle) + k * (k . v) * (1 - cos(angle))

    // radians
    double cosAngle = cos(angle);
    double sinAngle = sin(angle);

    // k x v
    Vector cross = crossProduct(axis, target);

    // k . v
    double dot = dotProduct(axis, target);

    Vector term1 = multiply(cosAngle, target);
    Vector term2 = multiply(sinAngle, cross);
    Vector term3 = multiply(dot * (1 - cosAngle), axis);
    Vector result = add(term1, term2);
    result = add(result, term3);
    return result;
}

Object::Object() : reference_point(3), color(3)
{
    height = width = length = 0.0;
    color.vec[0] = color.vec[1] = color.vec[2] = 0.0;
    ambient = diffuse = specular = reflection = 0.0;
    shine = 0;
}

void Object::setColor(double c[3])
{
    color.vec[0] = c[0];
    color.vec[1] = c[1];
    color.vec[2] = c[2];
}

void Object::setShine(int s)
{
    shine = s;
}

void Object::setCoefficients(double ambient, double diffuse, double specular, double reflection)
{
    this->ambient = ambient;
    this->diffuse = diffuse;
    this->specular = specular;
    this->reflection = reflection;
}

double Object::intersect(const Ray &ray, double *current_color, int level)
{
    double t = solveIntersection(ray);

    if (t < 0)
        return -1.0;

    if (level == 0)
        return t;

    // color computation
    Vector intersection_point = add(ray.start, multiply(t - EPSILON, ray.dir));
    Vector normal = getNormal(intersection_point, ray.dir);

    Vector resultColor(3);
    resultColor.vec[0] = resultColor.vec[1] = resultColor.vec[2] = 0.0;

    // add ambient color
    Vector colorAtPoint = getColorAt(intersection_point);
    resultColor = add(resultColor, multiply(ambient, colorAtPoint));

    for (int i = 0; i < environment.pointLights.size(); i++)
    {
        PointLight *light = environment.pointLights[i];

        // intersection point to light
        Vector lightDir = add(light->position, multiply(-1, intersection_point));
        double lightDistance = lightDir.norm();
        lightDir.normalize();

        // ray from intersection point to light
        // start point = intersection_point + EPSILON * lightDir to avoid self-shadowing
        Ray shadowRay(add(intersection_point, multiply(EPSILON, lightDir)), lightDir);

        bool inShadow = false;
        for (Object *obj : environment.objects)
        {
            double shadowT = obj->intersect(shadowRay, NULL, 0);
            if (shadowT > EPSILON && shadowT < lightDistance - EPSILON)
            {
                inShadow = true;
                break;
            }
        }

        if (inShadow)
            continue;

        // Diffuse component
        // Id = kd * IL * max(0, N.L)
        double dotNL = dotProduct(normal, lightDir);
        if (dotNL > 0)
        {
            Vector diffuseColor = multiply(diffuse * dotNL, light->color);
            resultColor = add(resultColor, diffuseColor);
        }

        // Specular component
        // Is = ks * IL * max(0, R.V)^shine
        Vector viewDir = multiply(-1, ray.dir);
        viewDir.normalize();
        // R = 2 * (L.N) * N - L
        Vector reflectDir = add(multiply(2 * dotProduct(lightDir, normal), normal), multiply(-1, lightDir));
        reflectDir.normalize();
        double spec = pow(std::max(dotProduct(viewDir, reflectDir), 0.0), shine);
        Vector specularColor = multiply(specular * spec, light->color);
        resultColor = add(resultColor, specularColor);
    }

    for (int i = 0; i < environment.spotLights.size(); i++)
    {
        SpotLight *sl = environment.spotLights[i];

        // intersection point to light
        Vector lightDir = add(sl->pointLight.position, multiply(-1, intersection_point));
        double lightDistance = lightDir.norm();
        lightDir.normalize();

        // check cutoff angle
        double angleCos = dotProduct(multiply(-1, sl->direction), lightDir);
        if (angleCos < cos(sl->angle * M_PI / 180.0))
            continue;

        // spotFactor * attenuation
        double spotFactor = pow(angleCos, 10) / (EPSILON + lightDistance*lightDistance);

        Ray shadowRay(add(intersection_point, multiply(EPSILON, lightDir)), lightDir);

        bool inShadow = false;
        for (Object *obj : environment.objects)
        {
            double shadowT = obj->intersect(shadowRay, NULL, 0);
            if (shadowT > 0 && shadowT < lightDistance)
            {
                inShadow = true;
                break;
            }
        }
        if (inShadow)
            continue;

        // Diffuse component - correctly check if light is in front of surface
        double dotNL = dotProduct(normal, lightDir);
        if (dotNL > 0)
        {
            Vector diffuseColor = multiply(diffuse * dotNL * spotFactor, sl->pointLight.color);
            resultColor = add(resultColor, diffuseColor);
        }

        // Specular component
        // Is = ks * IL * max(0, R.V)^shine
        Vector viewDir = multiply(-1, ray.dir);
        viewDir.normalize();
        // R = 2 * (L.N) * N - L
        Vector reflectDir = add(multiply(2 * dotProduct(lightDir, normal), normal), multiply(-1, lightDir));
        reflectDir.normalize();
        double spec = pow(std::max(dotProduct(viewDir, reflectDir), 0.0), shine);
        Vector specularColor = multiply(specular * spec * spotFactor, sl->pointLight.color);
        resultColor = add(resultColor, specularColor);
    }

    if (level >= environment.recursionDepth)
    {
        current_color[0] = resultColor.vec[0];
        current_color[1] = resultColor.vec[1];
        current_color[2] = resultColor.vec[2];
        return t;
    }

    // reflected ray
    Vector reflected_ray_dir = add(ray.dir, multiply(-2.0 * dotProduct(ray.dir, normal), normal));
    reflected_ray_dir.normalize();
    Vector reflected_ray_start = add(intersection_point, multiply(EPSILON, reflected_ray_dir));
    Ray reflected_ray(reflected_ray_start, reflected_ray_dir);

    // find nearest object for reflection
    int nearest_object_index = -1;
    double t_min = std::numeric_limits<double>::max();

    for (size_t k = 0; k < environment.objects.size(); ++k)
    {
        double dummy_color[3];
        double t_reflected = environment.objects[k]->intersect(reflected_ray, dummy_color, 0);
        if (t_reflected > EPSILON && t_reflected < t_min + EPSILON)
        {
            t_min = t_reflected;
            nearest_object_index = k;
        }
    }

    if (nearest_object_index != -1)
    {
        double reflected_color_arr[3] = {0.0, 0.0, 0.0};
        environment.objects[nearest_object_index]->intersect(reflected_ray, reflected_color_arr, level + 1);

        Vector reflected_color(3);
        reflected_color.vec[0] = reflected_color_arr[0];
        reflected_color.vec[1] = reflected_color_arr[1];
        reflected_color.vec[2] = reflected_color_arr[2];

        resultColor = add(resultColor, multiply(reflection, reflected_color));
    }

    current_color[0] = resultColor.vec[0];
    current_color[1] = resultColor.vec[1];
    current_color[2] = resultColor.vec[2];

    return t;
}

Vector Object::getColorAt(const Vector &point)
{
    return color;
}

Sphere::Sphere(const Vector &center, double radius)
{
    reference_point = center;
    length = radius; // Using length to represent the radius
}

void Sphere::draw()
{
    glColor3f(color.vec[0], color.vec[1], color.vec[2]);
    glPushMatrix();
    glTranslatef(reference_point.vec[0], reference_point.vec[1], reference_point.vec[2]);
    glutSolidSphere(length, 200, 200); // Using length as radius
    glPopMatrix();
}

double Sphere::solveIntersection(const Ray &ray)
{
    Vector o = add(ray.start, multiply(-1, reference_point)); // o = ray.start - sphere.center
    double b = 2 * dotProduct(o, ray.dir);
    double c = dotProduct(o, o) - length * length; // length is radius

    double discriminant = b * b - 4 * c;

    if (discriminant < 0)
        return -1.0; // No intersection

    double t1 = (-b - sqrt(discriminant)) / 2.0;
    double t2 = (-b + sqrt(discriminant)) / 2.0;

    if (t1 > EPSILON)
        return t1;
    else if (t2 > EPSILON)
        return t2;

    return -1.0; // Both intersections are behind the ray
}

Vector Sphere::getNormal(const Vector &point, const Vector &direction)
{
    Vector normal = add(point, multiply(-1, reference_point));
    normal.normalize();
    return normal;
}

Triangle::Triangle(const Vector &p1, const Vector &p2, const Vector &p3) : p1(p1), p2(p2), p3(p3) {}

void Triangle::draw()
{
    glColor3f(color.vec[0], color.vec[1], color.vec[2]);
    glBegin(GL_TRIANGLES);
    glVertex3f(p1.vec[0], p1.vec[1], p1.vec[2]);
    glVertex3f(p2.vec[0], p2.vec[1], p2.vec[2]);
    glVertex3f(p3.vec[0], p3.vec[1], p3.vec[2]);
    glEnd();
}

double Triangle::solveIntersection(const Ray &ray)
{
    // Implement Moller-Trumbore ray-triangle intersection algorithm
    Vector edge1 = add(p2, multiply(-1, p1));
    Vector edge2 = add(p3, multiply(-1, p1));

    Vector h = crossProduct(ray.dir, edge2);
    double a = dotProduct(edge1, h);

    // Ray is parallel to the triangle
    if (a > -EPSILON && a < EPSILON)
        return -1.0;

    double f = 1.0 / a;
    Vector s = add(ray.start, multiply(-1, p1));
    double u = f * dotProduct(s, h);

    // Ray doesn't intersect the triangle
    if (u < 0.0 || u > 1.0)
        return -1.0;

    Vector q = crossProduct(s, edge1);
    double v = f * dotProduct(ray.dir, q);

    // Ray doesn't intersect the triangle
    if (v < 0.0 || u + v > 1.0)
        return -1.0;

    // Calculate t, ray intersection
    double t = f * dotProduct(edge2, q);

    if (t > EPSILON)
        return t;

    return -1.0; // No intersection or behind the ray
}

Vector Triangle::getNormal(const Vector &point, const Vector &direction)
{
    Vector edge1 = add(p2, multiply(-1, p1));
    Vector edge2 = add(p3, multiply(-1, p1));
    Vector normal = crossProduct(edge1, edge2);
    normal.normalize();

    if (dotProduct(normal, direction) > 0)
    {
        normal = multiply(-1, normal);
    }

    return normal;
}

Vector Triangle::getColorAt(const Vector &point)
{
    return color;
}

GeneralQuadratic::GeneralQuadratic(double coeffs[10])
{
    for (int i = 0; i < 10; ++i)
    {
        coefficients[i] = coeffs[i];
    }
    // set from outside
    reference_point = Vector(3);
    height = width = length = 1.0;
}

void GeneralQuadratic::draw()
{
}

bool GeneralQuadratic::isPointInBoundingBox(const Vector &point) const
{
    // Calculate half dimensions of bounding box
    // double halfLength = length / 2;
    // double halfWidth = width / 2;
    // double halfHeight = height / 2;

    double halfLength = length;
    double halfWidth = width;
    double halfHeight = height;

    // If any dimension is 0, no clipping along that dimension
    bool withinX = (fabs(length) < EPSILON) ||
                   (point.vec[0] >= reference_point.vec[0] - halfLength &&
                    point.vec[0] <= reference_point.vec[0] + halfLength);

    bool withinY = (fabs(width) < EPSILON) ||
                   (point.vec[1] >= reference_point.vec[1] - halfWidth &&
                    point.vec[1] <= reference_point.vec[1] + halfWidth);

    bool withinZ = (fabs(height) < EPSILON) ||
                   (point.vec[2] >= reference_point.vec[2] - halfHeight &&
                    point.vec[2] <= reference_point.vec[2] + halfHeight);

    return withinX && withinY && withinZ;
}

double GeneralQuadratic::solveIntersection(const Ray &ray)
{
    // General quadratic equation: Ax^2 + By^2 + Cz^2 + Dxy + Exz + Fyz + Gx + Hy + Iz + J = 0
    // coefficients[0] = A, coefficients[1] = B, ..., coefficients[9] = J

    double A = coefficients[0];
    double B = coefficients[1];
    double C = coefficients[2];
    double D = coefficients[3];
    double E = coefficients[4];
    double F = coefficients[5];
    double G = coefficients[6];
    double H = coefficients[7];
    double I = coefficients[8];
    double J = coefficients[9];

    double x0 = ray.start.vec[0];
    double y0 = ray.start.vec[1];
    double z0 = ray.start.vec[2];
    double dx = ray.dir.vec[0];
    double dy = ray.dir.vec[1];
    double dz = ray.dir.vec[2];

    // Compute the coefficients of the quadratic equation at² + bt + c = 0
    double a = A * dx * dx + B * dy * dy + C * dz * dz +
               D * dx * dy + E * dx * dz + F * dy * dz;

    double b = 2 * A * x0 * dx + 2 * B * y0 * dy + 2 * C * z0 * dz +
               D * (x0 * dy + y0 * dx) + E * (x0 * dz + z0 * dx) +
               F * (y0 * dz + z0 * dy) + G * dx + H * dy + I * dz;

    double c = A * x0 * x0 + B * y0 * y0 + C * z0 * z0 +
               D * x0 * y0 + E * x0 * z0 + F * y0 * z0 +
               G * x0 + H * y0 + I * z0 + J;

    // Solve the quadratic equation
    if (fabs(a) < EPSILON)
    {
        // If a is approximately zero, we have a linear equation
        if (fabs(b) < EPSILON)
        {
            return -1.0; // No solution
        }
        double t = -c / b;
        if (t > EPSILON)
        {
            Vector point = add(ray.start, multiply(t, ray.dir));
            if (isPointInBoundingBox(point))
            {
                return t;
            }
        }
        return -1.0;
    }

    double discriminant = b * b - 4 * a * c;

    if (discriminant < 0)
    {
        return -1.0; // No real solutions
    }

    double t1 = (-b - sqrt(discriminant)) / (2 * a);
    double t2 = (-b + sqrt(discriminant)) / (2 * a);

    // Check if intersection points are within the bounding box
    bool t1Valid = false;
    bool t2Valid = false;

    if (t1 > EPSILON)
    {
        Vector point1 = add(ray.start, multiply(t1, ray.dir));
        t1Valid = isPointInBoundingBox(point1);
    }

    if (t2 > EPSILON)
    {
        Vector point2 = add(ray.start, multiply(t2, ray.dir));
        t2Valid = isPointInBoundingBox(point2);
    }

    if (t1Valid && t2Valid)
    {
        return std::min(t1, t2); // Return the closer intersection
    }
    else if (t1Valid)
    {
        return t1;
    }
    else if (t2Valid)
    {
        return t2;
    }

    return -1.0; // Both solutions are either negative or outside the bounding box
}

Vector GeneralQuadratic::getNormal(const Vector &point, const Vector &direction)
{
    // The gradient of the quadratic function gives the normal vector
    double A = coefficients[0];
    double B = coefficients[1];
    double C = coefficients[2];
    double D = coefficients[3];
    double E = coefficients[4];
    double F = coefficients[5];
    double G = coefficients[6];
    double H = coefficients[7];
    double I = coefficients[8];

    double x = point.vec[0];
    double y = point.vec[1];
    double z = point.vec[2];

    Vector normal(3);
    normal.vec[0] = 2 * A * x + D * y + E * z + G;
    normal.vec[1] = 2 * B * y + D * x + F * z + H;
    normal.vec[2] = 2 * C * z + E * x + F * y + I;

    normal.normalize();

    if (dotProduct(normal, direction) > 0)
    {
        normal = multiply(-1, normal); // Ensure the normal points outward
    }

    return normal;
}

Vector GeneralQuadratic::getColorAt(const Vector &point)
{
    return color;
}

Floor::Floor(double floorWidth, double tileWidth)
{
    this->floorWidth = floorWidth;
    this->tileWidth = tileWidth;

    this->leftX = -floorWidth / 2;
    this->leftY = -floorWidth / 2;

    primary[0] = 1.0;
    primary[1] = 1.0;
    primary[2] = 1.0;

    secondary[0] = 0.0;
    secondary[1] = 0.0;
    secondary[2] = 0.0;

    texture = "";

    textureData = nullptr;
    textureWidth = 0;
    textureHeight = 0;
    textureChannels = 0;
    textureLoaded = false;
}

Floor::~Floor()
{
    freeTexture();
}

bool Floor::loadTexture(const std::string &filename)
{
    freeTexture();

    FILE *f = fopen(filename.c_str(), "rb");

    if (!f)
    {
        std::cerr << "Failed to open texture file: " << filename << std::endl;
        return false;
    }

    textureData = stbi_load_from_file(f, &textureWidth, &textureHeight, &textureChannels, 0);

    if (textureData)
    {
        textureLoaded = true;
        std::cout << "Texture loaded: " << filename << " (" << textureWidth << "x" << textureHeight
                  << ", " << textureChannels << " channels)" << std::endl;
        return true;
    }
    else
    {
        std::cerr << "Failed to load texture: " << filename << std::endl;
        return false;
    }
}

void Floor::freeTexture()
{
    if (textureData)
    {
        stbi_image_free(textureData);
        textureData = NULL;
        textureLoaded = false;
    }
}

double Floor::solveIntersection(const Ray &ray)
{
    // The floor is a plane at z=0, so we need to find the intersection of the ray with this plane
    double dz = ray.dir.vec[2];

    // parallel ray - won't intersect
    if (fabs(dz) < EPSILON)
    {
        return -1.0;
    }

    double t = -ray.start.vec[2] / dz;

    // check if within boundary
    if (t > EPSILON)
    {
        double x = ray.start.vec[0] + t * ray.dir.vec[0];
        double y = ray.start.vec[1] + t * ray.dir.vec[1];

        // Check if the point is within the floor's boundaries
        if (x >= leftX && x <= leftX + floorWidth &&
            y >= leftY && y <= leftY + floorWidth)
        {
            return t;
        }
    }

    return -1.0;
}

Vector Floor::getNormal(const Vector &point, const Vector &direction)
{
    Vector normal(3);
    normal.vec[0] = 0;
    normal.vec[1] = 0;
    normal.vec[2] = 1;
    return normal;
}

Vector Floor::getColorAt(const Vector &point)
{
    int tileX = static_cast<int>((point.vec[0] - leftX) / tileWidth);
    int tileY = static_cast<int>((point.vec[1] - leftY) / tileWidth);

    if (point.vec[0] < leftX)
        tileX--;
    if (point.vec[1] < leftY)
        tileY--;

    Vector tileColor(3);

    if (applyTexture && textureLoaded)
    {
        double localX = (point.vec[0] - (leftX + tileX * tileWidth)) / tileWidth;
        double localY = (point.vec[1] - (leftY + tileY * tileWidth)) / tileWidth;

        int texX = static_cast<int>(localX * textureWidth);
        int texY = static_cast<int>(localY * textureHeight);

        texX = std::max(0, std::min(texX, textureWidth - 1));
        texY = std::max(0, std::min(texY, textureHeight - 1));

        int index = (texY * textureWidth + texX) * textureChannels;
        tileColor.vec[0] = textureData[index] / 255.0;
        tileColor.vec[1] = textureData[index + 1] / 255.0;
        tileColor.vec[2] = textureData[index + 2] / 255.0;
    }
    else
    {
        bool isPrimary = ((tileX + tileY) % 2 == 0);

        if (isPrimary)
        {
            tileColor.vec[0] = primary[0];
            tileColor.vec[1] = primary[1];
            tileColor.vec[2] = primary[2];
        }
        else
        {
            tileColor.vec[0] = secondary[0];
            tileColor.vec[1] = secondary[1];
            tileColor.vec[2] = secondary[2];
        }
    }

    return tileColor;
}

void Floor::draw()
{
    int numTiles = static_cast<int>(floorWidth / tileWidth);

    for (int i = 0; i < numTiles; ++i)
    {
        for (int j = 0; j < numTiles; ++j)
        {
            if ((i + j) % 2 == 0)
            {
                glColor3f(primary[0], primary[1], primary[2]); // Primary color
            }
            else
            {
                glColor3f(secondary[0], secondary[1], secondary[2]); // Secondary color
            }

            glBegin(GL_QUADS);
            glVertex3f(leftX + i * tileWidth, leftY + j * tileWidth, 0);
            glVertex3f(leftX + (i + 1) * tileWidth, leftY + j * tileWidth, 0);
            glVertex3f(leftX + (i + 1) * tileWidth, leftY + (j + 1) * tileWidth, 0);
            glVertex3f(leftX + i * tileWidth, leftY + (j + 1) * tileWidth, 0);
            glEnd();
        }
    }
}

PointLight::PointLight(const Vector &pos, const double col[3]) : position(pos), color(3)
{
    color.vec[0] = col[0];
    color.vec[1] = col[1];
    color.vec[2] = col[2];
}

void PointLight::draw()
{
    glColor3f(color.vec[0], color.vec[1], color.vec[2]);
    glPushMatrix();
    glTranslatef(position.vec[0], position.vec[1], position.vec[2]);
    glutSolidSphere(1, 20, 20);
    glPopMatrix();
}

SpotLight::SpotLight(const PointLight &pl, const Vector &dir, double ang)
    : pointLight(pl), direction(dir), angle(ang)
{
    if (angle < 0 || angle > 180)
        throw std::invalid_argument("Angle must be between 0 and 180 degrees.");
    direction.normalize();
}

void SpotLight::draw()
{
    pointLight.draw();
    glColor3f(pointLight.color.vec[0], pointLight.color.vec[1], pointLight.color.vec[2]);
    glPushMatrix();
    glTranslatef(pointLight.position.vec[0], pointLight.position.vec[1], pointLight.position.vec[2]);
    glRotatef(-angle / 2, direction.vec[0], direction.vec[1], direction.vec[2]);
    glutSolidCone(1, 1.0, 20, 20);
    glPopMatrix();
}

Camera::Camera(const Vector &eye_pos, const Vector &center_pos, const Vector &up_vec)
    : eye(eye_pos), center(center_pos), up(up_vec)
{
    movementSensitivity = 0.4;
    rotationSensitivity = 0.01 * M_PI;
}

void Camera::print()
{
    std::cout << "--------------------------------" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "(" << eye.vec[0] << ", " << eye.vec[1] << ", " << eye.vec[2] << ")\tCamera Position" << std::endl;
    std::cout << "(" << center.vec[0] << ", " << center.vec[1] << ", " << center.vec[2] << ")\tLook-at Point" << std::endl;
    std::cout << "(" << up.vec[0] << ", " << up.vec[1] << ", " << up.vec[2] << ")\tUp Vector" << std::endl;
}

Vector Camera::getLookDirection()
{
    Vector look(3);
    look.vec[0] = center.vec[0] - eye.vec[0];
    look.vec[1] = center.vec[1] - eye.vec[1];
    look.vec[2] = center.vec[2] - eye.vec[2];
    return look;
}

void Camera::moveFB(int dir)
{
    Vector look = getLookDirection();
    look.normalize();

    Vector movement = multiply(dir * movementSensitivity, look);
    eye = add(eye, movement);
    center = add(center, movement);
}

void Camera::moveLR(int dir)
{
    Vector look = getLookDirection();
    Vector right = crossProduct(look, up);
    right.normalize();

    Vector movement = multiply(dir * movementSensitivity, right);
    eye = add(eye, movement);
    center = add(center, movement);
}

void Camera::moveUD(int dir)
{
    Vector up_normalized = up;
    up_normalized.normalize();

    Vector movement = multiply(dir * movementSensitivity, up_normalized);
    eye = add(eye, movement);
    center = add(center, movement);
}

void Camera::rotateYaw(int dir)
{
    Vector look = getLookDirection();
    Vector up_normalized = up;
    up_normalized.normalize();

    Vector new_look = rotation(look, up_normalized, dir * rotationSensitivity);
    center.vec[0] = eye.vec[0] + new_look.vec[0];
    center.vec[1] = eye.vec[1] + new_look.vec[1];
    center.vec[2] = eye.vec[2] + new_look.vec[2];
}

void Camera::rotatePitch(int dir)
{
    Vector look = getLookDirection();
    Vector right = crossProduct(look, up);
    right.normalize();

    Vector new_look = rotation(look, right, dir * rotationSensitivity);
    Vector new_up = rotation(up, right, dir * rotationSensitivity);

    center.vec[0] = eye.vec[0] + new_look.vec[0];
    center.vec[1] = eye.vec[1] + new_look.vec[1];
    center.vec[2] = eye.vec[2] + new_look.vec[2];

    up = new_up;
}

void Camera::rotateRoll(int dir)
{
    Vector look = getLookDirection();
    look.normalize();

    up = rotation(up, look, dir * rotationSensitivity);
}

Ray::Ray(const Vector &start, const Vector &dir) : start(start), dir(dir)
{
    // Normalize the direction vector
    this->dir.normalize();
}

Environment::Environment() : recursionDepth(0), pixelCount(0) {}

Environment::~Environment()
{
    for (Object *obj : objects)
        delete obj;
    for (PointLight *pl : pointLights)
        delete pl;
    for (SpotLight *sl : spotLights)
        delete sl;
}

void capture()
{
    std::cout << "Capturing image..." << std::endl;

    if (applyTexture)
    {
        for (Object *obj : environment.objects)
        {
            Floor *floor = dynamic_cast<Floor *>(obj);
            if (floor)
            {

                std::string textureFile = "./textures/tex0" + std::to_string(textureId + 1) + ".jpg";
                floor->loadTexture(textureFile);
                break;
            }
        }
    }

    // initialize image
    bitmap_image image(environment.pixelCount, environment.pixelCount);
    for (int i = 0; i < environment.pixelCount; i++)
    {
        for (int j = 0; j < environment.pixelCount; j++)
        {
            image.set_pixel(i, j, 0, 0, 0);
        }
    }

    // double viewAngle = 80.0;
    double viewAngle = 45.0;
    double planeDistance = (environment.pixelCount / 2.0) / tan((viewAngle * M_PI / 180.0) / 2.0);

    Vector l = camera.getLookDirection();
    l.normalize();
    Vector r = crossProduct(l, camera.up);
    r.normalize();
    // up - making sure the 3 vectors are orthogonal
    Vector u = crossProduct(r, l);
    u.normalize();

    double windowHeight = environment.pixelCount;
    double windowWidth = environment.pixelCount;

    int imageWidth = environment.pixelCount;
    int imageHeight = environment.pixelCount;

    double du = (double)windowWidth / imageWidth;
    double dv = (double)windowHeight / imageHeight;

    // topleft = eye + l * planeDistance - r * windowWidth /2 + u * windowHeight /2
    Vector topleft = add(camera.eye, multiply(planeDistance, l));
    topleft = add(topleft, multiply(-windowWidth / 2.0, r));
    topleft = add(topleft, multiply(windowHeight / 2.0, u));

    topleft = add(topleft, multiply(0.5 * du, r));
    topleft = add(topleft, multiply(-0.5 * dv, u));

    for (int i = 0; i < imageWidth; ++i)
    {
        for (int j = 0; j < imageHeight; ++j)
        {
            // calculate curpixel
            Vector curPixel = add(topleft, multiply(i * du, r));
            curPixel = add(curPixel, multiply(-j * dv, u));

            // cast ray
            Vector rayDir = add(curPixel, multiply(-1, camera.eye));
            Ray ray(camera.eye, rayDir);

            int nearest_object_index = -1;
            double t_min = std::numeric_limits<double>::max();

            for (size_t k = 0; k < environment.objects.size(); ++k)
            {
                double dummy_color[3];
                double t = environment.objects[k]->intersect(ray, dummy_color, 0);
                if (t > 0 && t < t_min)
                {
                    t_min = t;
                    nearest_object_index = k;
                }
            }

            if (nearest_object_index != -1)
            {
                double color[3] = {0.0, 0.0, 0.0};
                environment.objects[nearest_object_index]->intersect(ray, color, 1);

                // Clamp color values to be within [0, 1]
                for (int c_idx = 0; c_idx < 3; ++c_idx)
                {
                    color[c_idx] = std::max(0.0, std::min(1.0, color[c_idx]));
                }

                // Update the bitmap image with the calculated color
                image.set_pixel(i, j, (unsigned char)(color[0] * 255), (unsigned char)(color[1] * 255), (unsigned char)(color[2] * 255));
            }
        }
    }

    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y%m%d_%H%M%S");
    std::string filename = "output_" + oss.str() + ".bmp";
    image.save_image(filename);
    std::cout << "Image saved as " << filename << std::endl;

    // Free textures after capture is done
    cleanupTextures();
}

void cleanupTextures()
{
    std::cout << "Cleaning up textures..." << std::endl;
    for (Object *obj : environment.objects)
    {
        Floor *floor = dynamic_cast<Floor *>(obj);
        if (floor)
        {
            floor->freeTexture();
        }
    }
}

void toggleTexture()
{
    applyTexture = !applyTexture;
    std::cout << "Texture mode " << (applyTexture ? "enabled" : "disabled") << std::endl;
}
