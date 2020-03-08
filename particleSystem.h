//Standard C/C++
#include <cstdlib>
#include <ctime>

#define PARTICLE_COUNT 5000

struct Particle
{
    int lifetime;
    float rho, theta, phi;
    float dRho, dTheta, dPhi;
};

class ParticleSystem
{
public:
    ParticleSystem();

    Particle &nextParticle();
    void reset();
    void update();
private:
    float randomNumber(float min, float max) const;
    Particle newParticle() const;

    Particle m_particles[PARTICLE_COUNT];
    int m_curParticle;
};
