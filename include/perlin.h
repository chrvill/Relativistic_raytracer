#include <cmath>
#include <vector>
#include <random>
#include <algorithm>

class PerlinNoise {
public:
    PerlinNoise(unsigned int seed = 0) {
        reseed(seed);
    }

    void reseed(unsigned int seed) {
        std::mt19937 generator(seed);
        std::uniform_real_distribution<double> distribution(0.0, 1.0);

        for (unsigned int i = 0; i < 256; ++i) {
            p.push_back(i);
            g1.push_back(distribution(generator) * 2.0 - 1.0);
        }

        // Shuffle the permutation vector
        std::shuffle(p.begin(), p.end(), generator);

        // Duplicate the permutation vector
        for (unsigned int i = 0; i < 256; ++i) {
            p.push_back(p[i]);
            g1.push_back(g1[i]);
        }
    }

    double noise(double x, double y, double z) const {
        int X = (int)std::floor(x) & 255;
        int Y = (int)std::floor(y) & 255;
        int Z = (int)std::floor(z) & 255;

        x -= std::floor(x);
        y -= std::floor(y);
        z -= std::floor(z);

        double u = fade(x);
        double v = fade(y);
        double w = fade(z);

        int A = p[X] + Y;
        int AA = p[A] + Z;
        int AB = p[A + 1] + Z;
        int B = p[X + 1] + Y;
        int BA = p[B] + Z;
        int BB = p[B + 1] + Z;

        return lerp(w, lerp(v, lerp(u, grad(p[AA], x, y, z),
                                      grad(p[BA], x - 1, y, z)),
                              lerp(u, grad(p[AB], x, y - 1, z),
                                      grad(p[BB], x - 1, y - 1, z))),
                       lerp(v, lerp(u, grad(p[AA + 1], x, y, z - 1),
                                      grad(p[BA + 1], x - 1, y, z - 1)),
                              lerp(u, grad(p[AB + 1], x, y - 1, z - 1),
                                      grad(p[BB + 1], x - 1, y - 1, z - 1))));
    }

private:
    std::vector<int> p;
    std::vector<double> g1;

    static double fade(double t) {
        return t * t * t * (t * (t * 6 - 15) + 10);
    }

    static double lerp(double t, double a, double b) {
        return a + t * (b - a);
    }

    double grad(int hash, double x, double y, double z) const {
        int h = hash & 15;
        double u = h < 8 ? x : y;
        double v = h < 4 ? y : h == 12 || h == 14 ? x : z;
        return ((h & 1) ? -u : u) + ((h & 2) ? -v : v);
    }
};
