#pragma once

#include "array2.h"

// Color
class Color8
{
public:
    unsigned char r, g, b, a;

    explicit inline Color8()
    {
    }
    explicit inline Color8(unsigned char v)
    {
        r = g = b = v;
        a = 255;
    }
    explicit inline Color8(unsigned char r, unsigned char g, unsigned char b, unsigned char a) : r(r), g(g), b(b), a(a)
    {
    }
    static Color8 Lerp(const double& t, const Color8& a, const Color8& b);
};

inline Color8 Color8::Lerp(const double& t, const Color8& a, const Color8& b) {
    return Color8((1.0 - t) * a.r + t * b.r, (1.0 - t) * a.g + t * b.g, (1.0 - t) * a.b + t * b.b, (1.0 - t) * a.a + t * b.a);
}

// CPU Texture
class Texture2D : public Array2
{
private:
    std::vector<Color8> colors;

public:
    inline Texture2D() { }
    Texture2D(int w, int h);
    Texture2D(const std::vector<Color8>& cols, int w, int h);
    inline ~Texture2D() { }

    void Fill(const Color8& c);
    const Color8* Data() const;
};
