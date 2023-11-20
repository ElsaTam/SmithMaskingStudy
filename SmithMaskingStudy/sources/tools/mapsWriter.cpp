#include "tools/mapsWriter.h"

#include <omp.h>
#include "utils/console.h"
#include "utils/params.h"

using namespace gdt;

struct Frame {
    gdt::vec3sc t, bt, n;
};

static gdt::vec3sc toWorld(const Frame &f, const gdt::vec3sc &v) {
    return f.t * v.x + f.bt * v.y + f.n * v.z;
}

// forward-declaration
bool createVisibilityMap(const std::string& filename,
                         const gdt::vec3sc &v,
                         const TriangleMesh* const shape,
                         bool fromDir, int nSamples);




bool MapsWriter::normalMap(const std::string& filename,
                           const TriangleMesh* const shape) {
    Console::print(OutLevel::ERR, "No library found for writing png files.");
    return false;
    
    /*
    int width = shape->bounds.size().x;
    int height = shape->bounds.size().y;

    FILE* fp;
    fopen_s(&fp, filename.c_str(), "wb");
    if(!fp) return false;

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) return false;

    png_infop info = png_create_info_struct(png);
    if (!info) return false;

    if (setjmp(png_jmpbuf(png))) return false;

    png_init_io(png, fp);

    // Output is 8bit depth, RGBA format.
    png_set_IHDR(
        png,
        info,
        width, height,
        8,
        PNG_COLOR_TYPE_RGBA,
        PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_DEFAULT,
        PNG_FILTER_TYPE_DEFAULT
    );
    png_write_info(png, info);

    png_bytep *row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * height);
    for(int y = 0; y < height; ++y) {
        row_pointers[y] = (png_byte*)malloc(png_get_rowbytes(png, info));
    }

    std::vector<vec3sc> data(width * height, vec3sc(0));
    const std::vector<vec3sc>& positions = shape->meshes[0]->vertex;
    const std::vector<vec3sc>& normals = shape->meshes[0]->normal;
    for(size_t i = 0; i < positions.size(); ++i) {
        int x = floor(positions[i].x - shape->bounds.lower.x);
        int y = floor(positions[i].y - shape->bounds.lower.y);
        if (x >= 0 && x < width && y >= 0 && y < height) {
            const vec3sc& n = normals[i];
            data[y * width + x] += n;
        }
    }

    int x, y;
#ifdef WITH_OPENMP
    #pragma omp parallel for collapse(2)
#endif
    for (y = 0; y < height; ++y) {
        for (x = 0; x < width; ++x) {
            vec3sc n = data[y * width + x];
            (&(row_pointers[y][x * 4]))[0] = 127 + n.x * 127; // r
            (&(row_pointers[y][x * 4]))[1] = 127 + n.y * 127; // g
            (&(row_pointers[y][x * 4]))[2] = 127 + n.z * 127; // b
            (&(row_pointers[y][x * 4]))[3] = 255;             // alpha : no transparency
        }
    }

    png_write_image(png, row_pointers);
    png_write_end(png, NULL);

    fclose(fp);

    png_destroy_write_struct(&png, &info);
    for(int y = 0; y < height; y++) {
        free(row_pointers[y]);
    }
    free(row_pointers);

    return true;
    */
}

bool MapsWriter::aoMap(const std::string& filename,
                       const TriangleMesh* shape,
                       int nSamples) {
    Console::print(OutLevel::ERR, "No library found for writing png files.");
    return false;

    /*
    vec2i size(shape->bounds.size().x, shape->bounds.size().y);
    float *img_data = new float[size.x * size.y];

    float z = shape->bounds.upper.z + 1.;
    float xMin = shape->bounds.lower.x;
    float yMin = shape->bounds.lower.y;

    int x, y;
#ifdef WITH_OPENMP
    #pragma omp parallel for collapse(2)
#endif
    for(y = 0; y < size.y; ++y) {
        for(x = 0; x < size.x; ++x) {
            img_data[y * size.x + x] = 0;
            for (int i = 0; i < nSamples; ++i) {
                vec2f d = math::randomSample<vec2f>();
                Ray ray(vec3sc(xMin + (float)x + d.x, yMin + (float)y + d.y, (float)z),
                        vec3sc(0, 0, -1));
                img_data[y * size.x + x] += Li(ray, shape);
            }
            img_data[y * size.x + x] /= static_cast<float>(nSamples);
        }
    }

    FILE* fp;
    fopen_s(&fp, filename.c_str(), "wb");
    if(!fp) return false;

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) return false;

    png_infop info = png_create_info_struct(png);
    if (!info) return false;

    if (setjmp(png_jmpbuf(png))) return false;

    png_init_io(png, fp);

    // Output is 8bit depth, RGBA format.
    png_set_IHDR(
        png,
        info,
        size.x, size.y,
        8,
        PNG_COLOR_TYPE_RGBA,
        PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_DEFAULT,
        PNG_FILTER_TYPE_DEFAULT
    );
    png_write_info(png, info);

    png_bytep *row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * size.y);
    for(int y = 0; y < size.y; ++y) {
        row_pointers[y] = (png_byte*)malloc(png_get_rowbytes(png, info));
    }

    for (int y = 0; y < size.y; ++y) {
        for (int x = 0; x < size.x; ++x) {
            float g = img_data[y * size.x + x];

            (&(row_pointers[y][x * 4]))[0] = 127 + g * 127; // r
            (&(row_pointers[y][x * 4]))[1] = 127 + g * 127; // g
            (&(row_pointers[y][x * 4]))[2] = 127 + g * 127; // b
            (&(row_pointers[y][x * 4]))[3] = 255;           // alpha : no transparency
        }
    }

    png_write_image(png, row_pointers);
    png_write_end(png, NULL);

    fclose(fp);

    png_destroy_write_struct(&png, &info);
    for(int y = 0; y < size.y; y++) {
        free(row_pointers[y]);
    }
    free(row_pointers);

    delete[] img_data;

    return true;
    */
}

bool MapsWriter::visibilityMapFromPoint(const std::string& filename,
                                        const gdt::vec3sc &pos,
                                        const TriangleMesh* const shape,
                                        int nSamples) {
    return createVisibilityMap(filename, pos, shape, false, nSamples);
}
bool MapsWriter::visibilityMapFromDir(const std::string& filename,
                                      const gdt::vec3sc &dir,
                                      const TriangleMesh* const shape,
                                      int nSamples) {
    return createVisibilityMap(filename, dir, shape, true, nSamples);
}


// "private" methods

bool createVisibilityMap(const std::string& filename,
                         const gdt::vec3sc &v,
                         const TriangleMesh* const shape,
                         bool fromDir, int nSamples) {

    Console::print(OutLevel::ERR, "No library found for writing png files.");
    return false;

    /*
    int width = shape->bounds.size().x;
    int height = shape->bounds.size().y;

    FILE* fp;
    fopen_s(&fp, filename.c_str(), "wb");
    if(!fp) return false;

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) return false;

    png_infop info = png_create_info_struct(png);
    if (!info) return false;

    if (setjmp(png_jmpbuf(png))) return false;

    png_init_io(png, fp);

    // Output is 8bit depth, RGBA format.
    png_set_IHDR(
        png,
        info,
        width, height,
        8,
        PNG_COLOR_TYPE_RGBA,
        PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_DEFAULT,
        PNG_FILTER_TYPE_DEFAULT
    );
    png_write_info(png, info);

    png_bytep *row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * height);
    for(int y = 0; y < height; ++y) {
        row_pointers[y] = (png_byte*)malloc(png_get_rowbytes(png, info));
    }

    ScalarRng rng(0, 1);
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            float V = 0;
            for (int i = 0; i < nSamples; ++i) {
                vec2f sample(x + rng(), y + rng());
                V += shape->isVisible(sample, v, fromDir);
            }
            V /= (float)nSamples;
            
            (&(row_pointers[y][x * 4]))[0] = V * 255; // r
            (&(row_pointers[y][x * 4]))[1] = V * 255; // g
            (&(row_pointers[y][x * 4]))[2] = V * 255; // b
            (&(row_pointers[y][x * 4]))[3] = 255;           // alpha : no transparency
        }
    }

    png_write_image(png, row_pointers);
    png_write_end(png, NULL);

    fclose(fp);

    png_destroy_write_struct(&png, &info);
    for(int y = 0; y < height; y++) {
        free(row_pointers[y]);
    }
    free(row_pointers);

    return true;
    */
}
