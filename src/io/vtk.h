// Copyright (c) 2024-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later
#ifndef SFCGAL_IO_VTK_H_
#define SFCGAL_IO_VTK_H_

#include "SFCGAL/Geometry.h"
#include <istream>
#include <memory>
#include <ostream>
#include <string>

namespace SFCGAL::io::VTK {

/**
 * @brief Saves a geometry to a legacy VTK format stream.
 *
 * @param[in] geom The geometry to save
 * @param[out] out The output stream
 * @throws SFCGAL::Exception If the geometry is invalid or unsupported
 */
SFCGAL_API void
save(const Geometry &geom, std::ostream &out);

/**
 * @brief Saves a geometry to a legacy VTK file.
 *
 * @param[in] geom The geometry to save
 * @param[in] filename The name of the file to save to
 * @throws SFCGAL::Exception If the file cannot be opened or the geometry is
 * invalid
 */
SFCGAL_API void
save(const Geometry &geom, const std::string &filename);

/**
 * @brief Saves a geometry to a legacy VTK format string.
 *
 * @param[in] geom The geometry to save
 * @return The legacy VTK format string
 * @throws SFCGAL::Exception If the geometry is invalid or unsupported
 */
SFCGAL_API auto
saveToString(const Geometry &geom) -> std::string;

/**
 * @brief Saves a geometry to a legacy VTK format buffer (C API).
 *
 * @param[in] geom The geometry to save
 * @param[out] buffer The buffer to write to
 * @param[in,out] size On input, the size of the buffer. On output, the number
 * of bytes written (or required if buffer is null)
 * @throws SFCGAL::Exception If the geometry is invalid or unsupported
 */
SFCGAL_API void
saveToBuffer(const Geometry &geom, char *buffer, size_t *size);

/**
 * @brief Loads a geometry from a legacy VTK format stream.
 *
 * @param[in] in The input stream
 * @return The loaded geometry
 * @throws SFCGAL::Exception If the stream is invalid or malformed
 */
SFCGAL_API auto
load(std::istream &in) -> std::unique_ptr<Geometry>;

/**
 * @brief Loads a geometry from a legacy VTK format string.
 *
 * @param[in] vtk The VTK format string
 * @return The loaded geometry
 * @throws SFCGAL::Exception If the string is invalid or malformed
 */
SFCGAL_API auto
load(const std::string &vtk) -> std::unique_ptr<Geometry>;

/**
 * @brief Loads a geometry from a legacy VTK file.
 *
 * @param[in] filename The name of the file to load from
 * @return The loaded geometry
 * @throws SFCGAL::Exception If the file cannot be opened or is invalid
 */
SFCGAL_API auto
loadFromFile(const std::string &filename) -> std::unique_ptr<Geometry>;

} // namespace SFCGAL::io::VTK

#endif // SFCGAL_IO_VTK_H_
