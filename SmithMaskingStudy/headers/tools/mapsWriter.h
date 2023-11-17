#pragma once

#include <iostream>
#include "shapes/TriangleMesh.h"

/**
 * @brief Namespace containing methods used to write maps and textures.
 * 
 */
namespace MapsWriter {

    /**
     * @brief Write the normal map in a .png file
     * 
     * @param filename name of the .png file
     * @param shape heightfield associated
     * @return true if writing went alright
     * @return false if a problem occured
     */
    bool normalMap(const std::string& filename,
                   const TriangleMesh* const shape);

    /**
     * @brief 
     * 
     * @param filename name of the .png file
     * @param shape heightfield associated
     * @param nSamples number of samples used by texel
     * @return true if writing went alright
     * @return false if a problem occured
     */
	bool aoMap(const std::string& filename,
               const TriangleMesh*shape,
               int nSamples = 1);

    /**
     * @brief 
     * 
     * @param filename name of the .png file
     * @param pos position of the observation point (or ponctual light source)
     * @param shape heightfield associated
     * @param nSamples number of samples used by texel
     * @return true if writing went alright
     * @return false if a problem occured
     */
    bool visibilityMapFromPoint(const std::string& filename,
                                const gdt::vec3sc &pos,
                                const TriangleMesh* const shape,
                                int nSamples = 1);

    /**
     * @brief 
     * 
     * @param filename name of the .png file
     * @param dir view direction (or light direction)
     * @param shape heightfield associated
     * @param nSamples number of samples used by texel
     * @return true if writing went alright
     * @return false if a problem occured
     */
    bool visibilityMapFromDir(const std::string& filename,
                              const gdt::vec3sc &dir,
                              const TriangleMesh* const shape,
                              int nSamples = 1);

};