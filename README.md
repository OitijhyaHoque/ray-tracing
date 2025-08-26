# Ray Tracing

This project is an implementation of ray tracing. It renders 3D scenes described in a text file. The application uses OpenGL for real-time camera control and visualization, and it can capture the ray-traced scene to a BMP image file.

![Ray Traced Scene](./outputs/output_20250711_205553.bmp)

## Features

  * **Ray-Object Intersection**: Implements ray intersection logic for multiple geometric primitives:
      * Spheres
      * Triangles
      * General Quadric Surfaces
      * A checkerboard floor
  * **Illumination Model**: Utilizes the Phong Lighting Model to calculate lighting effects, including:
      * Ambient, diffuse, and specular components.
      * Support for multiple light sources (Point Lights and Spotlights).
      * Shadow calculation by casting shadow rays towards light sources.
  * **Recursive Reflection**: Simulates realistic reflections by recursively tracing reflected rays up to a user-defined depth.
  * **Texture Mapping**: Supports applying image textures to the floor object, with a mechanism to toggle between the texture and the default checkerboard pattern.
  * **Fully Controllable Camera**: A 3D camera with functionalities for movement and rotation:
      * Move: Forward, backward, left, right, up, and down.
      * Rotate: Yaw (look left/right), Pitch (look up/down), and Roll (tilt).
  * **Scene Capture**: Renders the current camera view using ray tracing and saves it as a high-quality BMP image file.

## Run

This project requires **OpenGL** and **GLUT** to run. Use `run.sh` to compile and run the project on linux. Make changes to `input/scene.txt` to modify the scene. The scene file also includes how it should be formatted. Please note that the general quadratic curve is not visible while running the program, it is only visible in the captured image.

Compiling on **Windows**: `g++ ./src/2005049_main.cpp ./src/2005049_obj.cpp -o r.exe -lfreeglut -lglew32 -lopengl32 -lglu32`

## Controls

### Camera Movement

| Key | Action |
| :--- | :--- |
| **Up Arrow** | Move Forward |
| **Down Arrow** | Move Backward |
| **Left Arrow** | Move Left |
| **Right Arrow**| Move Right |
| **Page Up** | Move Up |
| **Page Down** | Move Down |

### Camera Rotation

| Key | Action |
| :--- | :--- |
| **1** | Rotate Left (Yaw) |
| **2** | Rotate Right (Yaw) |
| **3** | Look Up (Pitch) |
| **4** | Look Down (Pitch) |
| **5** | Tilt Clockwise (Roll) |
| **6** | Tilt Counter-clockwise (Roll) |

### Application Controls

| Key | Action |
| :--- | :--- |
| **0** | **Capture Scene**: Renders the current view and saves it as `output_YYYYMMDD_HHMMSS.bmp`. |
| **T / t** | **Toggle Texture**: Switches the floor rendering between the checkerboard and the loaded texture. |
| **U / u** | **Update Texture**: Cycles to the next available texture file (e.g., `tex01.jpg`, `tex02.jpg`). |
| **ESC** | **Exit** the application. |

## Sample Outputs

![Ray Traced Scene](./outputs/output_20250711_205553.bmp)
![Texture 1](./outputs/output_20250711_205928.bmp)
![Texture 2](./outputs/output_20250711_210743.bmp)
