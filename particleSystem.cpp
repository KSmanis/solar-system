//Own
#include "particleSystem.h"

//Particle System Contants
static const int lifetimeMin = 10;
static const int lifetimeMax = 20;
static const float distMin = 0;
static const float distMax = 30;
static const float angleMin = 0;
static const float angleMax = 360;
static const float dDistMin = 1;
static const float dDistMax = 15;
static const float dAngleMin = 0.5;
static const float dAngleMax = 3.0;

ParticleSystem::ParticleSystem()
{
    reset();
}

Particle &ParticleSystem::nextParticle()
{
    m_curParticle %= PARTICLE_COUNT;
    return m_particles[m_curParticle++];
}
void ParticleSystem::reset()
{
    for (int i = 0; i < PARTICLE_COUNT; ++i) {
        m_particles[i] = newParticle();
    }
    m_curParticle = 0;
}
void ParticleSystem::update()
{
    for (int i = 0; i < PARTICLE_COUNT; ++i) {
        if (m_particles[i].lifetime > 0) {
            m_particles[i].rho += m_particles[i].dRho;
            m_particles[i].theta += m_particles[i].dTheta;
            m_particles[i].phi += m_particles[i].dPhi;
            m_particles[i].lifetime--;
        } else {
            m_particles[i] = newParticle();
        }
    }
}

float ParticleSystem::randomNumber(float min, float max) const
{
    static bool firstTime = true;
    if (firstTime) {
        srand(time(0));
        firstTime = false;
    }

    return min + rand() * (max - min) / RAND_MAX;
}
Particle ParticleSystem::newParticle() const
{
    Particle p;
    p.lifetime = (int)randomNumber(lifetimeMin, lifetimeMax);
    p.rho = randomNumber(distMin, distMax);
    p.theta = randomNumber(angleMin, angleMax);
    p.phi = randomNumber(angleMin, angleMax);
    p.dRho = randomNumber(dDistMin, dDistMax);
    p.dTheta = randomNumber(dAngleMin, dAngleMax);
    p.dPhi = randomNumber(dAngleMin, dAngleMax);
    return p;
}
