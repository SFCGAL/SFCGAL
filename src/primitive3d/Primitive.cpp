// Copyright (c) 2025-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/primitive3d/Primitive.h"
#include "SFCGAL/Exception.h"

namespace SFCGAL {

void
Primitive::setParameter(const std::string        &name,
                        const PrimitiveParameter &parameter)
{
  auto parameterIt = m_parameters.find(name);

  // check that parameter name exists
  if (parameterIt == m_parameters.end()) {
    BOOST_THROW_EXCEPTION(
        Exception((boost::format("%s does not have a parameter named %s") %
                   primitiveType() % name)
                      .str()));
  }

  // check the type of the parameter
  const bool typeMatches = std::visit(
      [&](auto &&currentParameter) {
        using expectedType = std::decay_t<decltype(currentParameter)>;
        return std::holds_alternative<expectedType>(parameter);
      },
      parameterIt->second);

  if (!typeMatches) {
    BOOST_THROW_EXCEPTION(Exception(
        (boost::format("Wrong type for parameter %s.") % name).str()));
    return;
  }

  // Validate the new value
  validateAndSetParameter(name, parameter);
}

void
Primitive::validateAndSetParameter(const std::string        &name,
                                   const PrimitiveParameter &parameter)
{
  // check if valid
  auto tempParameters     = m_parameters;
  tempParameters.at(name) = parameter;
  validateParameters(tempParameters);

  // assign parameter
  m_parameters.at(name) = parameter;

  invalidateCache();
}

auto
Primitive::parameter(const std::string &name) const -> PrimitiveParameter
{
  auto parameterIt = m_parameters.find(name);
  // check that parameter name exists
  if (parameterIt == m_parameters.end()) {
    BOOST_THROW_EXCEPTION(
        Exception((boost::format("%s does not have a parameter named %s") %
                   primitiveType() % name)
                      .str()));
  }

  return parameterIt->second;
}

auto
Primitive::parameters() const
    -> std::unordered_map<std::string, PrimitiveParameter>
{
  return m_parameters;
}

void
Primitive::setTransformation(const Kernel::Aff_transformation_3 &transform)
{
  m_transform = transform;
  invalidateCache();
}

auto
Primitive::transformation() const -> Kernel::Aff_transformation_3
{
  return m_transform;
}

auto
Primitive::almostEqual(const Primitive &other, double epsilon) const -> bool
{
  using FT = Kernel::FT;

  if (primitiveTypeId() != other.primitiveTypeId()) {
    return false;
  }

  if (epsilon <= 0.0) {
    return (*this) == other;
  }

  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      if (CGAL::abs(m_transform.m(i, j) - other.m_transform.m(i, j)) >
          epsilon) {
        return false;
      }
    }
  }

  for (const auto &[key, value] : other.parameters()) {
    PrimitiveParameter primValue1 = parameter(key);
    PrimitiveParameter primValue2 = other.parameter(key);

    bool equal = std::visit(
        [&](auto &&value1) -> bool {
          using T = std::decay_t<decltype(value1)>;

          if (!std::holds_alternative<T>(primValue2)) {
            return false;
          }
          const T &value2 = std::get<T>(primValue2);

          if constexpr (std::is_same_v<T, FT> ||
                        std::is_same_v<T, unsigned int>) {
            return SFCGAL::almostEqual(value1, value2, epsilon);
          } else {
            return value1 == value2; // fallback
          }
        },
        primValue1);

    if (!equal) {
      return false;
    }
  }

  return true;
}

void
Primitive::invalidateCache()
{
  m_polyhedral_surface.reset();
}

auto
Primitive::toString() const -> std::string
{
  using FT = Kernel::FT;

  std::ostringstream stringStream;
  stringStream << "[Primitive type: " << primitiveType() << ", parameters: [";

  for (auto elem = m_parameters.cbegin(); elem != m_parameters.cend();) {
    //  for (const auto &elem : m_parameters) {
    stringStream << elem->first << ": ";
    std::visit(
        [&stringStream](auto &&value) {
          using T = std::decay_t<decltype(value)>;
          if constexpr (std::is_same_v<T, FT> ||
                        std::is_same_v<T, unsigned int>) {
            stringStream << value;
          } else {
            stringStream << "{ unknown alternative }";
          }
        },
        elem->second);

    ++elem;
    if (elem != m_parameters.cend()) {
      stringStream << ", ";
    }
  }
  stringStream << "]]";
  return stringStream.str();
}

auto
operator==(const Primitive &prim1, const Primitive &prim2) -> bool
{
  return prim1.primitiveTypeId() == prim2.primitiveTypeId() &&
         prim1.transformation() == prim2.transformation() &&
         prim1.parameters() == prim2.parameters();
}

} // namespace SFCGAL
