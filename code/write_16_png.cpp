#include "write_16_png.h"
#include "lodepng.h"
#include <iostream>

void write_16_png(const char* filename, const std::vector<unsigned char>& img, unsigned int width, unsigned int height) {

    if (img.size() != width * height * 2) {
        std::cerr << "img size (" << img.size() << ") != width (" << width << ") * height (" << height << ") * 2" << std::endl;
        return;
    }

    std::vector<unsigned char> image;
    image.resize(width * height * 2); // 2 bytes per pixel
    for (unsigned y = 0; y < height; y++)
        for (unsigned x = 0; x < width; x++)
        {
            image[2 * width * y + 2 * x + 0] = y; // most significant
            image[2 * width * y + 2 * x + 1] = x; // least significant
        }
    lodepng::encode(filename, img, width, height, LCT_GREY, 16);
    //lodepng::save_file(img, "SomeImage.png");
}