#include "tools/mipmapGenerator.h"

#include <vector>
#include "utils/console.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"


int mipmapLevels(int width, int height)
{
	if ((width == 0) || (width & (width - 1))) {
		Console::err << std::endl
			<< "Heightfield must be of width power of two."
			<< std::endl << std::endl;
		return 0;
	}
	if ((height == 0) || (height & (height - 1))) {
		Console::err << std::endl
			<< "Heightfield must be of height power of two."
			<< std::endl << std::endl;
		return 0;
	}
	return std::min(log2(width), log2(height));
}

void decreaseLevel(const std::vector<uint8_t>& pixels, int width, int height, std::vector<uint8_t>& newPixels)
{
	for (int y = 0; y < height - 1; y += 2)
	{
		for (int x = 0; x < width - 1; x += 2)
		{
			int index_top_left = y * width + x;
			int index_top_right = y * width + (x + 1);
			int index_bottom_left = (y + 1) * width + x;
			int index_bottom_right = (y + 1) * width + (x + 1);

			uint8_t g1 = pixels[index_top_left];
			uint8_t g2 = pixels[index_top_right];
			uint8_t g3 = pixels[index_bottom_left];
			uint8_t g4 = pixels[index_bottom_right];

			uint8_t g = (g1 + g2 + g3 + g4) / 4;
			newPixels.push_back(g);
		}
	}
}

void MipmapGenerator::generate(const std::string& filename, const std::string& inputFolder, const std::string& outputFolder)
{
	if (! CreateDirectoryA(outputFolder.c_str(), NULL) && ERROR_ALREADY_EXISTS != GetLastError())
	{
		Console::err << "Can't create folder " << outputFolder << std::endl;
		return;
	}

	int width, height, channels;
	unsigned char *data = stbi_load((inputFolder + filename).c_str(), &width, &height, &channels, 1);
	stbi_write_png((outputFolder + "0_" + filename).c_str(),
		width, height, channels, data, width * sizeof(char));

	std::vector<uint8_t> pixels;
	for (int i = 0; i < width * height; ++i) {
		pixels.push_back((uint8_t)data[i]);
	}
	stbi_image_free(data);

	int levels = mipmapLevels(width, height);
	for (int l = 1; l <= levels; ++l)
	{
		std::vector<uint8_t> grey;
		decreaseLevel(pixels, width, height, grey);
		width /= 2;
		height /= 2;
		stbi_write_png((outputFolder + std::to_string(l) + "_" + filename).c_str(),
					   width, height, channels, grey.data(), width * sizeof(uint8_t));
		pixels = grey;
	}

	Console::succ << std::endl
		<< "Image saved to " << outputFolder << filename
		<< std::endl << std::endl;
}