//Standard C/C++
#ifdef _WIN32
#include <cstdlib>
#define M_PI 3.14159265358979323846
#endif
#include <fstream>
using namespace std;

//GL
#ifdef __linux__
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#elif _WIN32
#include <glut.h>
#endif
#define ESC 27

//MathGeoLib
#include <MathGeoLib.h>

//SOIL
#include <SOIL.h>

//Project
#include "particleSystem.h"

struct CelestialObject
{
    // Orbital characteristics
    double eccentricity; // []
    double semiMajorAxis; // [Mm]
    double inclination; // [degrees]
    double longitude; // [degrees]
    double argument; // [degrees]

    // Physical characteristics
    double meanRadius; // [Mm]
    double axialTilt; // [degrees]
    double rotationSpeed; // [rad/s]
    double meanMotion; // [rad/s]

    // Appearance
    GLfloat r, g, b, a;
    const char *name;
    const char *texName;

    // State
    double meanAnomaly; // [rad]
    double rotation; // [rad]
    float4x4 orbitMat;
    GLuint texId;
    GLUquadricObj *obj;
};

//Objects
static const GLfloat lightPosition[4] = {0.0, 0.0, 0.0, 1.0};
static CelestialObject Sun =     {   0.0,        0.0,   0.0,      0.0,      0.0, 696.3420,     0.0,          0.0,         0.0,  1.0, 0.68, 0.35, 1.0, "Sun",         "textures/sun.jpg"};
static CelestialObject Mercury = {0.2056,   57909.05,  3.38,   48.331,   29.124,   2.4397,   0.034,  1.24002e-06, 8.26626e-07, 0.41, 0.41, 0.41, 1.0, "Mercury", "textures/mercury.jpg"};
static CelestialObject Venus =   {0.0067,  108208.00,  3.86,   76.678,   55.186,   6.0518,  177.36, -2.99247e-07, 3.23624e-07, 0.98, 0.98, 0.82, 1.0, "Venus",     "textures/venus.jpg"};
static CelestialObject Earth =   {0.0167,  149513.00, 7.155, -11.2606, 102.9471,   6.3710, 23.4393,  7.29212e-05, 1.99256e-07,  0.0,  0.0,  0.5, 1.0, "Earth",     "textures/earth.jpg"};
static CelestialObject Mars =    {0.0935,  227939.10,  5.65,   49.562,  286.537,   3.3895,   25.19,  7.08821e-05, 1.05852e-07,  1.0,  0.5, 0.31, 1.0, "Mars",       "textures/mars.jpg"};
static CelestialObject Jupiter = {0.0488,  778547.20,  6.09,  100.492,  275.066,  69.9110,    3.13,  1.75853e-04, 1.67768e-08, 0.87, 0.72, 0.53, 1.0, "Jupiter", "textures/jupiter.jpg"};
static CelestialObject Saturn =  {0.0557, 1433449.37,  5.51, 113.6428, 336.0139,  58.2320,   26.73,  1.63788e-04, 6.71301e-09, 0.93, 0.91, 0.66, 1.0, "Saturn",   "textures/saturn.jpg"};
static CelestialObject Uranus =  {0.0472, 2870671.40,  6.48,  73.9993,  96.9989,  25.3620,   97.77, -1.01238e-04, 2.36844e-09, 0.67, 0.85, 0.90, 1.0, "Uranus",   "textures/uranus.jpg"};
static CelestialObject Neptune = {0.0087, 4498542.60,  6.43, 131.7830, 273.2194,  24.6220,   28.32,  1.08338e-04, 1.20735e-09, 0.27, 0.51, 0.70, 1.0, "Neptune", "textures/neptune.jpg"};
static CelestialObject Pluto =   {0.2488, 5874000.00, 11.88, 110.2868, 113.7635,   1.1840, 119.591, -1.13856e-05, 8.09148e-10,  0.4, 0.31, 0.21, 1.0, "Pluto",     "textures/pluto.jpg"};
static CelestialObject *celestialObjects[] = {&Sun, &Mercury, &Venus, &Earth, &Mars, &Jupiter, &Saturn, &Uranus, &Neptune, &Pluto, 0};
static GLuint particleTexId;
static ParticleSystem particles;
static VecArray KBOs;
static GLuint saturnRingTexId;
static GLUquadricObj *saturnRingDisk;

//Drawing options
static bool axesOn = false;
static bool camFollow = false;
static bool HUDOn = true;
static bool KBOsOn = true;
static bool orbitsOn = true;
static bool particlesOn = true;
static bool paused = false;
//Quanta
static const GLdouble motionQuantum = 0.0075;
static const GLdouble rotateQuantum = 0.03;
static const GLdouble timeQuantum = 1.0e-03; // GLUT time resolution: 1 ms
static const GLdouble translateQuantum = 1000.0;
static const GLdouble zoomQuantum = 1000.0;
//View options
static const GLdouble defaultCamPos = 3.0e+05;
static const GLfloat planetScale = 1.0e+03;
static const GLfloat sunScale = 2.5e+01;
static const CelestialObject *camTarget = 0;
static GLdouble camX = defaultCamPos;
static GLdouble camY = defaultCamPos;
static GLdouble camZ = defaultCamPos;
static bool fullScreen = false;
static int motionX;
static int motionY;
static GLdouble multiplier = 1.0e+02;
static GLdouble tick = 8.64e+01; // 1 day/second
static GLdouble transX = 0.0;
static GLdouble transY = 0.0;
static GLdouble transZ = 0.0;
//Window properties
static const int defaultWinWidth = 800;
static const int defaultWinHeight = 600;
static int lastWinWidth;
static int lastWinHeight;
static int winWidth;
static int winHeight;
static int winId;

//Helper functions
void init();
void cleanup();
bool loadKBOs(const char *fileName);
GLuint loadTexture(const char *fileName);
GLdouble camRho();
GLdouble camTheta();
GLdouble camPhi();
//Compute functions
double eccentricAnomaly(double meanAnomaly, double eccentricity);
float3 heliocentricCoordinates(const CelestialObject &co);
float4x4 orbitalMatrix(const CelestialObject &co);
float4x4 rotationMatrix(const CelestialObject &co);
//View functions
void resetCamPos();
void updateProjection();
void updateCamPos();
void updateCartesianCamPos(GLdouble x, GLdouble y, GLdouble z);
void updateSphericalCamPos(GLdouble rho, GLdouble theta, GLdouble phi);
//Plot functions
void drawAxes(GLfloat length);
void drawEllipse(GLdouble a, GLdouble e, GLuint segmentCount);
void drawHUD();
void drawKBOs();
void drawParticleSystem();
void drawPlanet(const CelestialObject &co);
void drawPlanets();
void drawSaturnRings();
void drawString(const char *str);
void drawSun();
//Callback functions
void displayCallback();
void idleCallback();
void keyboardCallback(unsigned char key, int x, int y);
void menuCallback(int value);
void motionCallback(int x, int y);
void mouseCallback(int button, int state, int x, int y);
void reshapeCallback(int width, int height);
void specialCallback(int key, int x, int y);

//Helper functions
void init()
{
    //setup menus
    glutCreateMenu(menuCallback);
    glutAttachMenu(GLUT_RIGHT_BUTTON);
    //setup planets
    for (int i = 0; celestialObjects[i]; ++i) {
        glutAddMenuEntry(celestialObjects[i]->name, i);
        if (i > 0) {
            celestialObjects[i]->orbitMat = orbitalMatrix(*celestialObjects[i]);
        }
        celestialObjects[i]->texId = loadTexture(celestialObjects[i]->texName);
        celestialObjects[i]->obj = gluNewQuadric();
        gluQuadricTexture(celestialObjects[i]->obj, GLU_TRUE);
    }
    //setup misc textures and objects
    particleTexId = loadTexture("textures/particle.bmp");
    saturnRingTexId = loadTexture("textures/saturn_ring.png");
    saturnRingDisk = gluNewQuadric();
    gluQuadricTexture(saturnRingDisk, GLU_TRUE);
    //setup lights
    GLfloat sceneAmbient[] = {0.0, 0.0, 0.0, 1.0};
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, sceneAmbient);
    GLfloat sunAmbient[] = {0.1, 0.1, 0.1, 1.0};
    GLfloat sunDiffuse[] = {1.0, 1.0, 1.0, 1.0};
    GLfloat sunSpecular[] = {1.0, 1.0, 1.0, 1.0};
    glLightfv(GL_LIGHT0, GL_AMBIENT, sunAmbient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, sunDiffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, sunSpecular);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);
    //setup materials
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    //enable depth buffer
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_DEPTH_TEST);
    //enable culling
    glFrontFace(GL_CCW);
    glEnable(GL_CULL_FACE);
    //black canvas
    glClearColor(0.0, 0.0, 0.0, 1.0);
    //set our perspective
    updateProjection();
}
void cleanup()
{
    for (int i = 0; celestialObjects[i]; ++i) {
        gluDeleteQuadric(celestialObjects[i]->obj);
        celestialObjects[i]->obj = 0;
    }
    gluDeleteQuadric(saturnRingDisk);
    saturnRingDisk = 0;
}
bool loadKBOs(const char *fileName)
{
    ifstream file(fileName);
    if (!file.is_open()) {
        cerr << "Error loading Kepler Belt Objects from file: " << fileName << endl;
        return false;
    }

    KBOs.clear();

    string line;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }

        CelestialObject KBO;
        istringstream str(line);
        str >> KBO.semiMajorAxis >> KBO.eccentricity >> KBO.inclination >> KBO.longitude >> KBO.argument >> KBO.meanAnomaly;
        KBO.semiMajorAxis *= 149597.871; // Astronomical Units (AU) -> Megameters (Mm)
        KBO.meanAnomaly *= M_PI / 180.f; // Degrees -> Radians
        KBO.orbitMat = orbitalMatrix(KBO);
        KBOs.push_back(heliocentricCoordinates(KBO));
    }
    return true;
}
GLuint loadTexture(const char *fileName)
{
    GLuint texId = SOIL_load_OGL_texture(fileName, SOIL_LOAD_AUTO, SOIL_CREATE_NEW_ID, SOIL_FLAG_INVERT_Y | SOIL_FLAG_MIPMAPS | SOIL_FLAG_TEXTURE_REPEATS);
    if (!texId) {
        cerr << "Error loading texture: " << fileName << endl;
        cerr << "SOIL error: " << SOIL_last_result() << endl;
    }
    return texId;
}
GLdouble camRho()
{
    return sqrt(pow(camX, 2) + pow(camY, 2) + pow(camZ, 2));
}
GLdouble camTheta()
{
    GLdouble rho = camRho();
    return acos(camZ / (rho != 0.0 ? rho : DBL_MIN));
}
GLdouble camPhi()
{
    return atan2(camY, camX);
}
//Compute functions
double eccentricAnomaly(double meanAnomaly, double eccentricity)
{
    return atan2(sin(meanAnomaly), cos(meanAnomaly) - eccentricity);
}
float3 heliocentricCoordinates(const CelestialObject &co)
{
    // Eccentricity (e)
    const GLdouble &e = co.eccentricity;
    // Semi-major axis length (α)
    const GLdouble &a = co.semiMajorAxis;
    // Mean anomaly (M)
    const GLdouble &M = co.meanAnomaly;
    // Eccentric anomaly (E)
    const GLdouble E = eccentricAnomaly(M, e);
    // Semi-minor axis length (β)
    const GLdouble b = a * sqrt(1.0 - e * e);

    return co.orbitMat.TransformPos(a * cos(E), b * sin(E), 0.0);
}
float4x4 orbitalMatrix(const CelestialObject &co)
{
    // Eccentricity (e)
    const GLdouble &e = co.eccentricity;
    // Semi-major axis length (α)
    const GLdouble &a = co.semiMajorAxis;
    // Mean anomaly (M)
    const GLdouble &M = co.meanAnomaly;
    // Eccentric anomaly (E)
    const GLdouble E = eccentricAnomaly(M, e);
    // Semi-minor axis length (β)
    const GLdouble b = a * sqrt(1.0 - e * e);

    // Quaternion for argument of the perihelion (ω) rotation
    Quat q1 = Quat::RotateZ(DegToRad(co.argument));
    // Quaternion for orbital plane inclination (i) rotation
    Quat q2 = Quat::RotateX(DegToRad(co.inclination));
    // Quaternion for longitude of the ascending node (Ω) rotation
    Quat q3 = Quat::RotateZ(DegToRad(co.longitude));
    // Quaternion containing all the necessary orbital rotations:
    // 1. Argument of the perihelion (ω)
    // 2. Orbital plane inclination (i)
    // 3. Longitude of the ascending node (Ω)
    Quat qOrbit = (q3 * q2 * q1).Normalized();

    // Matrix containing all the necessary orbital transformations:
    // 1. Translate to focal point (Sun), since f = α * e
    // 2. Perform orbital rotations (around focal point)
    return qOrbit.ToFloat4x4() * float4x4::Translate(-a * e, 0.0, 0.0);
}
float4x4 rotationMatrix(const CelestialObject &co)
{
    // Quaternion for planet self-rotation
    Quat q1 = Quat::RotateZ(co.rotation);
    // Quaternion for axial tilt rotation
    Quat q2 = Quat::RotateX(-DegToRad(co.axialTilt));
    // Quaternion containing all the necessary planetary rotations:
    // 1. Planet spin around its rotation axis
    // 2. Axial tilt of the planet
    Quat qPlanet = (q2 * q1).Normalized();
    // Matrix containing all the necessary planetary rotations:
    return qPlanet.ToFloat4x4();
}
//View functions
void resetCamPos()
{
    camX = camY = camZ = (camTarget ? 5 * camTarget->meanRadius * planetScale : defaultCamPos);
    transX = transY = transZ = 0.0;
}
void updateProjection()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glViewport(0, 0, winWidth, winHeight);
    gluPerspective(45.0, ((GLdouble)winWidth) / winHeight, 1000.0, 100000000.0);
}
void updateCamPos()
{
    float3 eye(camX, camY, camZ), center(0.0, 0.0, 0.0), up(0.0, 0.0, 1.0);
    if (camTarget) {
        center = heliocentricCoordinates(*camTarget);
        if (camFollow) {
            float4x4 rotMat = rotationMatrix(*camTarget);
            eye = rotMat.TransformPos(eye) + center;
            up = rotMat.TransformPos(up);
        } else {
            eye = eye + center;
        }
    }

    //set our viewpoint
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(eye.x, eye.y, eye.z, center.x, center.y, center.z, up.x, up.y, up.z);
}
void updateCartesianCamPos(GLdouble x, GLdouble y, GLdouble z)
{
    camX = x;
    camY = y;
    camZ = z;
    updateCamPos();
}
void updateSphericalCamPos(GLdouble rho, GLdouble theta, GLdouble phi)
{
    if (rho < 0.0 || theta < 0.0 || theta > M_PI) {
        return;
    }
    updateCartesianCamPos(rho * sin(theta) * cos(phi), rho * sin(theta) * sin(phi), rho * cos(theta));
}
//Plot functions
void drawAxes(GLfloat length)
{
    GLboolean lightState = glIsEnabled(GL_LIGHTING);
    if (lightState) {
        glDisable(GL_LIGHTING);
    }

    glBegin(GL_LINES);
    //[X]
    glColor3f(1.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(length, 0.0, 0.0);
    //[Y]
    glColor3f(0.0, 1.0, 0.0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, length, 0.0);
    //[Z]
    glColor3f(0.0, 0.0, 1.0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, length);
    glEnd();

    if (lightState) {
        glEnable(GL_LIGHTING);
    }
}
void drawEllipse(GLdouble a, GLdouble e, GLuint segmentCount)
{
    GLboolean lightState = glIsEnabled(GL_LIGHTING);
    if (lightState) {
        glDisable(GL_LIGHTING);
    }

    // Semi-minor axis length (β)
    const GLdouble b = a * sqrt(1.0 - e * e);

    glBegin(GL_LINE_LOOP);
    for (GLfloat t = 0.0f; t < 2 * M_PI; t += 2 * M_PI / segmentCount) {
        glVertex2f(a * cos(t), b * sin(t));
    }
    glEnd();

    if (lightState) {
        glEnable(GL_LIGHTING);
    }
}
void drawHUD()
{
    GLboolean lightState = glIsEnabled(GL_LIGHTING);
    if (lightState) {
        glDisable(GL_LIGHTING);
    }

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();

    gluOrtho2D(0, winWidth, winHeight, 0);
    glColor3f(0.0, 1.0, 0.0);
    char str[64];

    sprintf(str, "Planet Scale: %gx", planetScale);
    glRasterPos2i(4, 15);
    drawString(str);

    sprintf(str, "Sun Scale: %gx", sunScale);
    glRasterPos2i(4, 30);
    drawString(str);

    if (paused) {
        sprintf(str, "Time Flow: PAUSED");
    } else {
        sprintf(str, "Time Flow: %gx", tick / timeQuantum);
    }
    glRasterPos2i(4, 45);
    drawString(str);

    sprintf(str, "Multiplier: %gx", multiplier);
    glRasterPos2i(4, 60);
    drawString(str);

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    if (lightState) {
        glEnable(GL_LIGHTING);
    }
}
void drawKBOs()
{
    GLboolean lightState = glIsEnabled(GL_LIGHTING);
    if (lightState) {
        glDisable(GL_LIGHTING);
    }

    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_POINTS);
    for (VecArray::const_iterator it = KBOs.begin(); it != KBOs.end(); ++it) {
        glVertex3fv(it->ptr());
    }
    glEnd();

    if (lightState) {
        glEnable(GL_LIGHTING);
    }
}
void drawParticleSystem()
{
    //enable blending
    glBlendFunc(GL_SRC_COLOR, GL_ONE);
    glEnable(GL_BLEND);
    //disable depth
    glDepthMask(GL_FALSE);
    //enable texturing
    glColor3f(1.0, 1.0, 1.0);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, particleTexId);

    const GLfloat f = 150.0 * sunScale;
    for (int i = 0; i < PARTICLE_COUNT; ++i) {
        const Particle &p = particles.nextParticle();
        glPushMatrix();
            glRotatef(p.phi, 0.0, 1.0, 0.0);
            glRotatef(p.theta, 0.0, 0.0, 1.0);
            glTranslated((Sun.meanRadius + p.rho) * sunScale, 0.0, 0.0);
            glRotatef(90.0, 0.0, 1.0, 0.0);

            glBegin(GL_TRIANGLE_STRIP);
                glTexCoord2i(1, 1);
                glVertex2f(f, f); // Top Right
                glTexCoord2i(0, 1);
                glVertex2f(-f, f); // Top Left
                glTexCoord2i(1, 0);
                glVertex2f(f, -f); // Bottom Right
                glTexCoord2i(0, 0);
                glVertex2f(-f, -f); // Bottom Left
            glEnd();

            glBegin(GL_TRIANGLE_STRIP);
                glTexCoord2i(1, 1);
                glVertex2f(-f, f); // Top Right
                glTexCoord2i(0, 1);
                glVertex2f(f, f); // Top Left
                glTexCoord2i(1, 0);
                glVertex2f(-f, -f); // Bottom Right
                glTexCoord2i(0, 0);
                glVertex2f(f, -f); // Bottom Left
            glEnd();
        glPopMatrix();
    }

    //disable texturing
    glDisable(GL_TEXTURE_2D);
    //enable depth
    glDepthMask(GL_TRUE);
    //disable blending
    glDisable(GL_BLEND);
}
void drawPlanet(const CelestialObject &co)
{
    if (orbitsOn) {
        glPushMatrix();
            glMultMatrixf(co.orbitMat.Transposed().ptr());
            glColor3fv(&co.r);
            //Pluto's orbit is so enormous that we need greater
            //resolution to allow him to stay on his trajectory!
            drawEllipse(co.semiMajorAxis, co.eccentricity, strcmp(co.name, "Pluto") == 0 ? 300 : 120);
        glPopMatrix();
    }

    // Matrix containing all the necessary planetary transformations:
    // 1. Apply all planetary rotations
    // 2. Translate to planet's current position
    float4x4 mPlanet = float4x4::Translate(heliocentricCoordinates(co)) * rotationMatrix(co);

    glPushMatrix();
        glMultMatrixf(mPlanet.Transposed().ptr());
        if (axesOn) {
            drawAxes(3 * co.meanRadius * planetScale);
        }
        glColor3f(1.0, 1.0, 1.0);
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, co.texId);
        gluSphere(co.obj, co.meanRadius * planetScale, 60, 60);
        glDisable(GL_TEXTURE_2D);
        if (strcmp(co.name, "Saturn") == 0) {
            drawSaturnRings();
        }
    glPopMatrix();
}
void drawSaturnRings()
{
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
    glDisable(GL_CULL_FACE);

    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, saturnRingTexId);
    gluDisk(saturnRingDisk, 0, (80 + Saturn.meanRadius) * planetScale, 60, 60);
    glDisable(GL_TEXTURE_2D);

    glEnable(GL_CULL_FACE);
    glDisable(GL_BLEND);
}
void drawPlanets()
{
    for (int i = 1; celestialObjects[i]; ++i) {
        drawPlanet(*celestialObjects[i]);
    }
}
void drawString(const char *str)
{
    if (!str) {
        return;
    }

    while (*str) {
        glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *(str++));
    }
}
void drawSun()
{
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);

    if (axesOn) {
        drawAxes(3 * Sun.meanRadius * sunScale);
    }
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, Sun.texId);
    gluSphere(Sun.obj, Sun.meanRadius * sunScale, 60, 60);
    glDisable(GL_TEXTURE_2D);
}
//Callback functions
void displayCallback()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    updateCamPos();

    glPushMatrix();
        glTranslated(transX, transY, transZ);
        if (KBOsOn) {
            drawKBOs();
        }
        drawSun();
        drawPlanets();
        if (particlesOn) {
            drawParticleSystem();
        }
        if (HUDOn) {
            drawHUD();
        }
    glPopMatrix();

    glutSwapBuffers();
}
void idleCallback()
{
    static unsigned t0 = 0;
    const unsigned dt = glutGet(GLUT_ELAPSED_TIME) - t0;
    t0 += dt;

    if (paused) {
        return;
    }

    if (particlesOn) {
        particles.update();
    }

    for (int i = 1; celestialObjects[i]; ++i) {
        celestialObjects[i]->meanAnomaly += celestialObjects[i]->meanMotion * dt * tick;
        celestialObjects[i]->rotation += celestialObjects[i]->rotationSpeed * dt * tick;
    }
    glutPostRedisplay();
}
void keyboardCallback(unsigned char key, int x, int y)
{
    switch (key) {
    case ESC:
    case 'q':
    case 'Q':
        glutDestroyWindow(winId);
        exit(EXIT_SUCCESS);
    case 'f':
    case 'F':
        if (!fullScreen) {
            lastWinWidth = winWidth;
            lastWinHeight = winHeight;
            glutFullScreen();
        } else {
            glutReshapeWindow(lastWinWidth, lastWinHeight);
        }
        fullScreen = !fullScreen;
        break;
    case ' ':
        paused = !paused;
        break;
    case '0':
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9':
        multiplier = pow(10.f, key - '0');
        break;
    case 'a':
    case 'A':
        axesOn = !axesOn;
        break;
    case 'c':
    case 'C':
        camFollow = !camFollow;
        resetCamPos();
        break;
    case 'h':
    case 'H':
        HUDOn = !HUDOn;
        break;
    case 'k':
    case 'K':
        KBOsOn = !KBOsOn;
        break;
    case 'o':
    case 'O':
        orbitsOn = !orbitsOn;
        break;
    case 'p':
    case 'P':
        particlesOn = !particlesOn;
        if (!particlesOn) {
            particles.reset();
        }
        break;
    case 'x':
        transX += translateQuantum;
        break;
    case 'X':
        transX -= translateQuantum;
        break;
    case 'y':
        transY += translateQuantum;
        break;
    case 'Y':
        transY -= translateQuantum;
        break;
    case 'z':
        transZ += translateQuantum;
        break;
    case 'Z':
        transZ -= translateQuantum;
        break;
    case '+':
        updateSphericalCamPos(camRho() - zoomQuantum * multiplier, camTheta(), camPhi());
        break;
    case '-':
        updateSphericalCamPos(camRho() + zoomQuantum * multiplier, camTheta(), camPhi());
        break;
    case '*':
        tick += timeQuantum * multiplier;
        return;
    case '/':
        tick -= timeQuantum * multiplier;
        if (tick < timeQuantum) {
            tick = timeQuantum;
        }
        return;
    case '=':
        resetCamPos();
        break;
    default:
        return;
    }
    glutPostRedisplay();
}
void menuCallback(int value)
{
    camTarget = (value > 0 ? celestialObjects[value] : 0);
    resetCamPos();
    glutPostRedisplay();
}
void motionCallback(int x, int y)
{
    updateSphericalCamPos(camRho(), camTheta() + (motionY - y) * motionQuantum, camPhi() + (motionX - x) * motionQuantum);
    motionX = x;
    motionY = y;
    glutPostRedisplay();
}
void mouseCallback(int button, int state, int x, int y)
{
    if (state != GLUT_DOWN) {
        return;
    }

    switch (button) {
    case GLUT_LEFT_BUTTON:
    case GLUT_MIDDLE_BUTTON:
    case GLUT_RIGHT_BUTTON:
        motionX = x;
        motionY = y;
        break;
    case 0x0003: // wheel scroll up - non-standard behaviour
        updateSphericalCamPos(camRho() - zoomQuantum * multiplier, camTheta(), camPhi());
        glutPostRedisplay();
        break;
    case 0x0004: // wheel scroll down - non-standard behaviour
        updateSphericalCamPos(camRho() + zoomQuantum * multiplier, camTheta(), camPhi());
        glutPostRedisplay();
        break;
    }
}
void reshapeCallback(int width, int height)
{
    winWidth = width;
    winHeight = (height != 0 ? height : 1);
    updateProjection();
}
void specialCallback(int key, int x, int y)
{
    switch (key) {
    case GLUT_KEY_UP:
        updateSphericalCamPos(camRho(), camTheta() - rotateQuantum, camPhi());
        break;
    case GLUT_KEY_DOWN:
        updateSphericalCamPos(camRho(), camTheta() + rotateQuantum, camPhi());
        break;
    case GLUT_KEY_LEFT:
        updateSphericalCamPos(camRho(), camTheta(), camPhi() - rotateQuantum);
        break;
    case GLUT_KEY_RIGHT:
        updateSphericalCamPos(camRho(), camTheta(), camPhi() + rotateQuantum);
        break;
    default:
        return;
    }
    glutPostRedisplay();
}

int main(int argc, char **argv)
{
    loadKBOs("L7SyntheticModel-v09.txt");

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize(defaultWinWidth, defaultWinHeight);
    winId = glutCreateWindow("Solar System");
    glutDisplayFunc(displayCallback);
    glutIdleFunc(idleCallback);
    glutKeyboardFunc(keyboardCallback);
    glutMotionFunc(motionCallback);
    glutMouseFunc(mouseCallback);
    glutReshapeFunc(reshapeCallback);
    glutSpecialFunc(specialCallback);
    init();
    atexit(cleanup);
    glutMainLoop();
}
