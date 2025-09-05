#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <ctime>
#include <string>
#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
using namespace std;

const double G = 5.0;
const double dt = 0.005;
const int steps = 400000;
const double box_size = 300.0;
const double radius = 1.5;

struct Vec2 {
    double x, y;

    Vec2 operator+(const Vec2& o) const { return { x + o.x, y + o.y }; }
    Vec2 operator-(const Vec2& o) const { return { x - o.x, y - o.y }; }
    Vec2 operator*(double s) const { return { x * s, y * s }; }
    Vec2& operator+=(const Vec2& o) { x += o.x; y += o.y; return *this; }
    Vec2& operator-=(const Vec2& o) { x -= o.x; y -= o.y; return *this; }
};

double norm(const Vec2& v) {
    return sqrt(v.x * v.x + v.y * v.y) + 1e-5;
}

struct Particle {
    Vec2 pos;
    Vec2 vel;
    double mass;
    bool active = true;
};

void compute_gravity(vector<Particle>& p, vector<Vec2>& forces) {
    int N = p.size();
    forces.assign(N, { 0.0, 0.0 });
    for (int i = 0; i < N; ++i) {
        if (!p[i].active) continue;
        for (int j = i + 1; j < N; ++j) {
            if (!p[j].active) continue;
            Vec2 r = p[j].pos - p[i].pos;
            double d = norm(r);
            Vec2 dir = r * (1.0 / d);

            Vec2 f = { 0.0, 0.0 };

            if (d > radius && d < 300.0) {
                double Fg = G / (d * d);
                f += dir * Fg;
            }

            forces[i] += f;
            forces[j] -= f;
        }
    }
}

void bounces(vector<Particle>& p) {
    int N = p.size();
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            Vec2 r = p[j].pos - p[i].pos;
            double d = norm(r);
            if (d > 1.0 * radius && d < 2 * radius) {
                Vec2 n = r * (1.0 / d);
                double v1n = p[i].vel.x * n.x + p[i].vel.y * n.y;
                double v2n = p[j].vel.x * n.x + p[j].vel.y * n.y;
                if (v1n - v2n < 0) continue;
                double a = 0.0;
                if (v1n - v2n > 10) a = 0.1;
                double m1 = p[i].mass;
                double m2 = p[j].mass;
                double v1n_new = (v1n * (m1 - m2) + 2 * m2 * v2n) / (m1 + m2);
                double v2n_new = (v2n * (m2 - m1) + 2 * m1 * v1n) / (m1 + m2);
                Vec2 dv1 = n * (v1n_new - v1n) * (1 - a);
                Vec2 dv2 = n * (v2n_new - v2n) * (1 - a);
                p[i].vel += dv1;
                p[j].vel += dv2;
            }
        }
    }
}

void PollEvents(sf::RenderWindow& window, bool& paused, int& sps) {
    while (const optional event = window.pollEvent()) {
        if (event->is<sf::Event::Closed>()) {
            window.close();
        }
        if (const auto* keyPressed = event->getIf<sf::Event::KeyPressed>()) {
            if (keyPressed->scancode == sf::Keyboard::Scancode::Space) {
                paused = !paused;
            }
        }
        if (const auto* keyPressed = event->getIf<sf::Event::KeyPressed>()) {
            if (keyPressed->scancode == sf::Keyboard::Scancode::NumpadPlus) {
                sps += 100;
            }
        }
        if (const auto* keyPressed = event->getIf<sf::Event::KeyPressed>()) {
            if (keyPressed->scancode == sf::Keyboard::Scancode::NumpadMinus) {
                sps += 100;
            }
        }
    }
}

int main() {
    const int N_start = 300;
    int sps = 1000; // steps per second
	unsigned int step = 0;
	bool paused = false;
	sf::Vector2i lastMousePos;
    vector<Particle> particles;
    mt19937 rng(time(nullptr));
    uniform_real_distribution<> pos_dist(-box_size, box_size);
    uniform_real_distribution<> vel_dist(0.0, 0.4);

    sf::Clock fpsClock;
    int frameCount = 0;

    for (int i = 0; i < N_start; ++i) {
        particles.push_back({
            {pos_dist(rng), pos_dist(rng)},
            {vel_dist(rng), vel_dist(rng)},
            1.0
            });
    }

    vector<Vec2> forces;

    unsigned int width = 1600, height = 800;
	unsigned int width2 = 800, height2 = 800;
    sf::RenderWindow* window = new sf::RenderWindow(sf::VideoMode({ width, height }), "Particle Simulation");
    window->setFramerateLimit(sps);

    sf::RenderWindow* Energy = new sf::RenderWindow(sf::VideoMode({ width2, height2 }), "Graph of energy");
    Energy->setFramerateLimit(sps);


    sf::CircleShape circle(4.0f);
    circle.setOrigin(circle.getGeometricCenter());
    circle.setPosition({ width / 2.0f, height / 2.0f });
    circle.setFillColor(sf::Color::Blue);

	sf::CircleShape Ek(0.3f);
	Ek.setOrigin(Ek.getGeometricCenter());
	Ek.setFillColor(sf::Color::Blue);

    sf::CircleShape Ep(0.3f);
    Ep.setOrigin(Ep.getGeometricCenter());
    Ep.setFillColor(sf::Color::Yellow);

    sf::CircleShape Ec(0.3f);
    Ec.setOrigin(Ec.getGeometricCenter());
    Ec.setFillColor(sf::Color::Green);

    while (window->isOpen()) {
        PollEvents(*window, paused, sps);

        step++;

        int x = sf::Mouse::getPosition(*window).x;
        int y = sf::Mouse::getPosition(*window).y;

        if (!paused) {
            compute_gravity(particles, forces);
            bounces(particles);

            for (size_t i = 0; i < particles.size(); ++i) {
                Particle& pt = particles[i];
                pt.vel += forces[i] * (dt / pt.mass);
                pt.pos += pt.vel * dt;
            }
        }

        double kinetic = 0.0;
        for (const auto& pt : particles) {
            kinetic += 0.5 * pt.mass * (pt.vel.x * pt.vel.x + pt.vel.y * pt.vel.y);
        }

        double potential = 0.0;
        for (size_t i = 0; i < particles.size(); ++i) {
            for (size_t j = i + 1; j < particles.size(); ++j) {
                Vec2 r = particles[j].pos - particles[i].pos;
                double d = norm(r);
                potential -= G * particles[j].mass * particles[i].mass / d;
            }
        }

        Ek.setPosition({
            static_cast<float>(step * 0.01), static_cast<float>(-kinetic * 0.01 + height2 / 2)
            });
        Ep.setPosition({
            static_cast<float>(step * 0.01), static_cast<float>(-potential * 0.01 + height2 / 2)
            });
        Ec.setPosition({
            static_cast<float>(step * 0.01), static_cast<float>(-(kinetic + potential) * 0.01 + height2 / 2)
            });

        window->clear(sf::Color::Black);

        for (size_t i = 0; i < particles.size(); ++i) {
            circle.setPosition({
                static_cast<float>(particles[i].pos.x * (width / (2.0 * box_size)) + x),
                static_cast<float>(particles[i].pos.y * (height / (2.0 * box_size)) + y)
                });
            window->draw(circle);
        }

        if (!paused) {
            Energy->draw(Ek);
            Energy->draw(Ep);
            Energy->draw(Ec);
        }

        frameCount++;
        if (fpsClock.getElapsedTime().asSeconds() >= 1.0f) {
            cout << "FPS: " << frameCount << ", x: " << x << ", y: " << y << " Step: " << step << " Energy: " << kinetic << ", " << potential << ", " << kinetic + potential << endl;
            frameCount = 0;
            fpsClock.restart();
        }
        window->display();
		Energy->display();

    }


    delete window;
	delete Energy;

    return 0;
}
