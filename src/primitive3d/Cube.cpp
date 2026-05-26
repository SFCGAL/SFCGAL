// Copyright (c) 2025-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/primitive3d/Cube.h"

namespace SFCGAL {

Cube::Cube(const Kernel::FT &size) { m_parameters["size"] = size; }

auto
Cube::primitiveType() const -> std::string
{
  return "Cube";
}

auto
Cube::primitiveTypeId() const -> PrimitiveType
{
  return PrimitiveType::TYPE_CUBE;
}

void
Cube::invalidateCache()
{
  Primitive::invalidateCache();
  m_box.reset();
}

void
Cube::setSize(const Kernel::FT &size)
{
  validateAndSetParameter("size", size);
}

void
Cube::validateParameters(
    std::unordered_map<std::string, PrimitiveParameter> const &tempParameters)
    const
{
  const double size =
      CGAL::to_double(std::get<Kernel::FT>(tempParameters.at("size")));

  if (size < 0.) {
    BOOST_THROW_EXCEPTION(Exception("Cube size cannot be negative."));
  }
}

auto
Cube::generateBox() const -> Box
{
  if (!m_box) {
    m_box.emplace(size(), size(), size());
    m_box->setTransformation(m_transform);
  }

  return *m_box;
}

auto
Cube::generatePolyhedralSurface() const -> PolyhedralSurface
{
  return generateBox().generatePolyhedralSurface();
}

auto
Cube::area3D(bool withDiscretization) const -> double
{
  return generateBox().area3D(withDiscretization);
}

auto
Cube::volume(bool withDiscretization) const -> double
{
  return generateBox().volume(withDiscretization);
}

auto
Cube::size() const -> const Kernel::FT &
{
  return std::get<Kernel::FT>(m_parameters.at("size"));
}

auto
Cube::toString() const -> std::string
{
  std::ostringstream stringStream;
  stringStream << "[Primitive type: " << primitiveType() << ", size: " << size()
               << "]";

  return stringStream.str();
}

} // namespace SFCGAL
