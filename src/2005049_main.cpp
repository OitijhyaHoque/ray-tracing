#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <stdexcept>
#include <string>

#include "2005049_obj.h"

using namespace std;

extern int textureId;
extern int textureCount;

Vector eye_pos(3), center_pos(3), up_vec(3);
Camera camera(eye_pos, center_pos, up_vec);

Environment environment;

vector<string> readFileLines(const string &filename)
{
    vector<string> lines;
    ifstream infile(filename);

    if (!infile)
        throw runtime_error("Could not open file: " + filename);

    string line;
    while (getline(infile, line))
    {
        bool hasNonWhitespace = false;
        for (char c : line)
        {
            if (!isspace(c))
            {
                hasNonWhitespace = true;
                break;
            }
        }

        if (hasNonWhitespace)
        {
            lines.push_back(line);
        }
    }
    return lines;
}

vector<string> processLine(const string &line)
{
    vector<string> tokens;

    istringstream iss(line);
    string token;
    while (iss >> token)
    {
        tokens.push_back(token);
    }

    return tokens;
}

void loadData(string filename)
{
    vector<string> lines = readFileLines(filename);

    int i = 0, n = lines.size();
    bool parsedObjects = false, parsedPointlights = false, parsedSpotlights = false;

    while (i < n)
    {
        vector<string> tokens = processLine(lines[i]);

        if (tokens.size() == 2 && tokens[0] == "Input" && tokens[1] == "explanation")
        {
            break;
        }

        if (i == 0)
        {
            environment.recursionDepth = stoi(tokens[0]);
            cout << "Recursion Depth: " << environment.recursionDepth << endl;

            i++;
            continue;
        }

        if (i == 1)
        {
            environment.pixelCount = stoi(tokens[0]);
            cout << "Pixel Count: " << environment.pixelCount << endl;

            i++;
            continue;
        }

        if (i == 2)
        {
            int objectCount = stoi(tokens[0]);
            cout << "Object Count: " << objectCount << endl;
            parsedObjects = true;

            i++;

            for (int j = 0; j < objectCount; j++)
            {
                vector<string> tokens = processLine(lines[i]);

                if (tokens[0] == "sphere")
                {
                    vector<string> l1 = processLine(lines[i + 1]);
                    vector<string> l2 = processLine(lines[i + 2]);
                    vector<string> l3 = processLine(lines[i + 3]);
                    vector<string> l4 = processLine(lines[i + 4]);
                    vector<string> l5 = processLine(lines[i + 5]);

                    Vector center(3);
                    center.vec[0] = stod(l1[0]);
                    center.vec[1] = stod(l1[1]);
                    center.vec[2] = stod(l1[2]);

                    double radius = stod(l2[0]);
                    double color[3];
                    color[0] = stod(l3[0]);
                    color[1] = stod(l3[1]);
                    color[2] = stod(l3[2]);

                    double ambient = stod(l4[0]);
                    double diffuse = stod(l4[1]);
                    double specular = stod(l4[2]);
                    double reflection = stod(l4[3]);

                    int shine = stoi(l5[0]);

                    Sphere *sphere = new Sphere(center, radius);
                    sphere->setColor(color);
                    sphere->setCoefficients(ambient, diffuse, specular, reflection);
                    sphere->setShine(shine);

                    environment.objects.push_back(sphere);

                    i += 6;
                }
                else if (tokens[0] == "triangle")
                {
                    vector<string> l1 = processLine(lines[i + 1]);
                    vector<string> l2 = processLine(lines[i + 2]);
                    vector<string> l3 = processLine(lines[i + 3]);
                    vector<string> l4 = processLine(lines[i + 4]);
                    vector<string> l5 = processLine(lines[i + 5]);
                    vector<string> l6 = processLine(lines[i + 6]);

                    Vector p1(3), p2(3), p3(3);
                    p1.vec[0] = stod(l1[0]);
                    p1.vec[1] = stod(l1[1]);
                    p1.vec[2] = stod(l1[2]);
                    p2.vec[0] = stod(l2[0]);
                    p2.vec[1] = stod(l2[1]);
                    p2.vec[2] = stod(l2[2]);
                    p3.vec[0] = stod(l3[0]);
                    p3.vec[1] = stod(l3[1]);
                    p3.vec[2] = stod(l3[2]);

                    double color[3];
                    color[0] = stod(l4[0]);
                    color[1] = stod(l4[1]);
                    color[2] = stod(l4[2]);

                    double ambient = stod(l5[0]);
                    double diffuse = stod(l5[1]);
                    double specular = stod(l5[2]);
                    double reflection = stod(l5[3]);

                    int shine = stoi(l6[0]);

                    Triangle *triangle = new Triangle(p1, p2, p3);
                    triangle->setColor(color);
                    triangle->setCoefficients(ambient, diffuse, specular, reflection);
                    triangle->setShine(shine);

                    environment.objects.push_back(triangle);

                    i += 7;
                }
                else if (tokens[0] == "general")
                {
                    vector<string> l1 = processLine(lines[i + 1]);
                    vector<string> l2 = processLine(lines[i + 2]);
                    vector<string> l3 = processLine(lines[i + 3]);
                    vector<string> l4 = processLine(lines[i + 4]);
                    vector<string> l5 = processLine(lines[i + 5]);

                    double coeffs[10];
                    for (int k = 0; k < 10; k++)
                    {
                        coeffs[k] = stod(l1[k]);
                    }

                    Vector ref_pt(3);
                    ref_pt.vec[0] = stod(l2[0]);
                    ref_pt.vec[1] = stod(l2[1]);
                    ref_pt.vec[2] = stod(l2[2]);
                    double length = stod(l2[3]);
                    double width = stod(l2[4]);
                    double height = stod(l2[5]);

                    double color[3];
                    color[0] = stod(l3[0]);
                    color[1] = stod(l3[1]);
                    color[2] = stod(l3[2]);

                    double ambient = stod(l4[0]);
                    double diffuse = stod(l4[1]);
                    double specular = stod(l4[2]);
                    double reflection = stod(l4[3]);

                    int shine = stoi(l5[0]);

                    GeneralQuadratic *general = new GeneralQuadratic(coeffs);
                    general->reference_point = ref_pt;
                    general->length = length;
                    general->width = width;
                    general->height = height;
                    general->setColor(color);
                    general->setCoefficients(ambient, diffuse, specular, reflection);
                    general->setShine(shine);

                    environment.objects.push_back(general);

                    i += 6;
                }
                else
                {
                    throw runtime_error("Unknown object type: " + tokens[0]);
                    return;
                }
            }

            continue;
        }

        if (i > 2 && parsedObjects && !parsedPointlights)
        {
            int pointLightCount = stoi(tokens[0]);
            cout << "Point Light Count: " << pointLightCount << endl;

            parsedPointlights = true;

            i++;

            for (int j = 0; j < pointLightCount; j++)
            {
                vector<string> l1 = processLine(lines[i]);
                vector<string> l2 = processLine(lines[i + 1]);

                Vector position(3);
                position.vec[0] = stod(l1[0]);
                position.vec[1] = stod(l1[1]);
                position.vec[2] = stod(l1[2]);

                double color[3];
                color[0] = stod(l2[0]);
                color[1] = stod(l2[1]);
                color[2] = stod(l2[2]);

                PointLight *pointLight = new PointLight(position, color);
                environment.pointLights.push_back(pointLight);

                i += 2;
            }

            continue;
        }

        if (i > 2 && parsedObjects && parsedPointlights)
        {
            int spotLightCount = stoi(tokens[0]);
            cout << "Spot Light Count: " << spotLightCount << endl;

            parsedSpotlights = true;

            i++;

            for (int j = 0; j < spotLightCount; j++)
            {
                vector<string> l1 = processLine(lines[i]);
                vector<string> l2 = processLine(lines[i + 1]);
                vector<string> l3 = processLine(lines[i + 2]);
                vector<string> l4 = processLine(lines[i + 3]);

                Vector position(3);
                position.vec[0] = stod(l1[0]);
                position.vec[1] = stod(l1[1]);
                position.vec[2] = stod(l1[2]);

                double color[3];
                color[0] = stod(l2[0]);
                color[1] = stod(l2[1]);
                color[2] = stod(l2[2]);

                Vector direction(3);
                direction.vec[0] = stod(l3[0]);
                direction.vec[1] = stod(l3[1]);
                direction.vec[2] = stod(l3[2]);

                double angle = stod(l4[0]);

                SpotLight *spotLight = new SpotLight(PointLight(position, color), direction, angle);
                environment.spotLights.push_back(spotLight);

                i += 4;
            }

            continue;
        }
    }

    Floor *floor = new Floor(1000.0, 20.0);

    floor->setCoefficients(0.4, 0.2, 0.2, 0.2);
    floor->setShine(5);

    environment.objects.push_back(floor); // Add a default floor object
}

void initGL()
{
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black background
    glEnable(GL_DEPTH_TEST);              // Enable depth testing for z-culling
}

void display()
{
    // Clear color and depth buffers
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Set up the model-view matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // Position camera using the eye, center and up vectors
    gluLookAt(camera.eye.vec[0], camera.eye.vec[1], camera.eye.vec[2],          // Camera position
              camera.center.vec[0], camera.center.vec[1], camera.center.vec[2], // Look-at point
              camera.up.vec[0], camera.up.vec[1], camera.up.vec[2]);            // Up vector

    // Draw objects based on visibility flags
    for (auto &obj : environment.objects)
    {
        obj->draw();
    }

    for (auto &pl : environment.pointLights)
    {
        pl->draw();
    }

    for (auto &sl : environment.spotLights)
    {
        sl->draw();
    }

    // Swap buffers (double buffering)
    glutSwapBuffers();
}

void reshapeListener(GLsizei width, GLsizei height)
{
    // Prevent division by zero
    if (height == 0)
        height = 1;

    // Calculate aspect ratio
    GLfloat aspect = (GLfloat)width / (GLfloat)height;

    // Set viewport to cover entire window
    glViewport(0, 0, width, height);

    // Set up perspective projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    // 45-degree field of view, aspect ratio, near and far clipping planes
    gluPerspective(45.0f, aspect, 0.1f, 10000.0f);
}

void keyboardListener(unsigned char key, int x, int y)
{
    switch (key)
    {
    // --- Camera Rotation Controls ---
    case '1':
        camera.rotateYaw(1);
        break;
    case '2':
        camera.rotateYaw(-1);
        break;
    case '3':
        camera.rotatePitch(1);
        break;
    case '4':
        camera.rotatePitch(-1);
        break;
    case '5':
        camera.rotateRoll(1);
        break;
    case '6':
        camera.rotateRoll(-1);
        break;
    case '0':
        capture();
        break;
    case 'T':
        toggleTexture();
        break;
    case 't':
        toggleTexture();
        break;
    case 'U':
        textureId = (textureId + 1) % textureCount;
        break;
    case 'u':
        textureId = (textureId + 1) % textureCount;
        break;

    // --- Program Control ---
    case 27: // ESC key: exit program
        exit(0);
        break;
    }

    // camera.print();      // Print camera state for debugging
    glutPostRedisplay(); // Request a screen refresh
}

void specialKeyListener(int key, int x, int y)
{
    switch (key)
    {
    case GLUT_KEY_UP: // Move forward
        camera.moveFB(1);
        break;
    case GLUT_KEY_DOWN: // Move backward
        camera.moveFB(-1);
        break;
    case GLUT_KEY_LEFT: // Move left
        camera.moveLR(-1);
        break;
    case GLUT_KEY_RIGHT: // Move right
        camera.moveLR(1);
        break;
    case GLUT_KEY_PAGE_UP: // Move upward
        camera.moveUD(1);
        break;
    case GLUT_KEY_PAGE_DOWN: // Move downward
        camera.moveUD(-1);
        break;
    }

    // camera.print();      // Print camera state for debugging
    glutPostRedisplay(); // Request a screen refresh
}

int main(int argc, char *argv[])
{
    loadData("./input/scene.txt");

    // Initialize camera vectors
    eye_pos.vec[0] = 0.0;
    eye_pos.vec[1] = -50.0;
    eye_pos.vec[2] = 0.0;
    center_pos.vec[0] = 0.0;
    center_pos.vec[1] = 0.0;
    center_pos.vec[2] = 0.0;
    up_vec.vec[0] = 0.0;
    up_vec.vec[1] = 0.0;
    up_vec.vec[2] = 1.0;

    camera = Camera(eye_pos, center_pos, up_vec);

    // Initialize GLUT
    glutInit(&argc, argv);

    // Configure display mode and window
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(environment.pixelCount, environment.pixelCount);
    glutInitWindowPosition(50, 50);

    glutCreateWindow("Ray tracing");

    // Register callback functions
    glutDisplayFunc(display);
    glutReshapeFunc(reshapeListener);
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);

    // Initialize OpenGL settings
    initGL();

    // Enter the GLUT event loop
    glutMainLoop();

    return 0;
}
