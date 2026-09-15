// Copyright (c) 2025-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_DEREFERENCEITERATOR_H_
#define SFCGAL_DEREFERENCEITERATOR_H_

#include <iterator>
#include <type_traits>

namespace SFCGAL {

/**
 * @brief Iterator adapter that dereferences unique_ptr elements automatically.
 *
 * This class wraps a base iterator (e.g.,
 * std::vector<std::unique_ptr<T>>::iterator) and returns references to the
 * pointed objects instead of unique_ptrs.
 *
 * @tparam BaseIterator the underlying iterator type
 */
template <class BaseIterator>
class DereferenceIterator {
public:
  using value_type =
      typename BaseIterator::value_type::element_type; ///< Dereferenced
  using pointer         = value_type *; ///< Pointer to value type
  using reference       = value_type &; ///< Reference to value type
  using difference_type = typename std::iterator_traits<
      BaseIterator>::difference_type; ///<  Distance between two iterators
  using iterator_concept =
      std::forward_iterator_tag; ///< Forward iterator concept

  /** @brief Default constructor. */
  DereferenceIterator() = default;

  /**
   * @brief Construct from a base iterator.
   * @param it The underlying base iterator to wrap.
   */
  explicit DereferenceIterator(BaseIterator it) : base_(it) {}

  /**
   * @brief Conversion constructor from another compatible DereferenceIterator.
   *
   * Allows implicit conversion from DereferenceIterator<iterator>
   * to DereferenceIterator<const_iterator>.
   *
   * @tparam OtherIterator Another iterator type convertible to BaseIterator.
   * @param other The other DereferenceIterator to copy.
   */
  template <typename OtherIterator>
  DereferenceIterator(const DereferenceIterator<OtherIterator> &other)
    requires(std::is_convertible_v<OtherIterator, BaseIterator>)
      : base_(other.base())
  {
  }

  /**
   * @brief base iterator access
   * @return base iterator
   */
  [[nodiscard]] auto
  base() const -> BaseIterator
  {
    return base_;
  }

  /**
   * @brief Dereference operator returning a reference to the pointed value.
   * @return Reference to the object pointed by the current iterator.
   */
  auto
  operator*() const -> reference
  {
    return *(*base_);
  }

  /**
   * @brief Member access operator returning a pointer to the pointed value.
   * @return Pointer to the object pointed by the current iterator.
   */
  auto
  operator->() const -> pointer
  {
    return base_->get();
  }

  /**
   * @brief Random access operator returning a reference to the nth value.
   * @param idx The offset from the current iterator position.
   * @return Reference to the nth object from the current iterator.
   */
  auto
  operator[](size_t idx) const -> reference
  {
    return *(base_[idx]);
  }

  /**
   * @brief Pre-increment operator.
   * @return Reference to this iterator, advanced by one position.
   */
  auto
  operator++() -> DereferenceIterator &
  {
    ++base_;
    return *this;
  }

  /**
   * @brief Post-increment operator.
   * @return A copy of this iterator before it was advanced.
   */
  auto
  operator++(int) -> DereferenceIterator
  {
    auto tmp = *this;
    ++base_;
    return tmp;
  }

  /**
   * @brief Pre-decrement operator.
   * @return Reference to this iterator, moved back by one position.
   */
  auto
  operator--() -> DereferenceIterator &
  {
    --base_;
    return *this;
  }

  /**
   * @brief Post-decrement operator.
   * @return A copy of this iterator before it was moved back.
   */
  auto
  operator--(int) -> DereferenceIterator
  {
    auto tmp = *this;
    --base_;
    return tmp;
  }

  /**
   * @brief Equality comparison.
   */
  friend auto
  operator==(const DereferenceIterator &, const DereferenceIterator &)
      -> bool = default;

private:
  BaseIterator base_{};
};

/**
 * @brief Helper function to create a DereferenceIterator from a base iterator.
 *
 * @tparam Iterator The underlying iterator type.
 * @param iterator The base iterator to wrap.
 * @return A DereferenceIterator instance wrapping @p t.
 */
template <typename Iterator>
auto
dereference_iterator(Iterator iterator) -> DereferenceIterator<Iterator>
{
  return DereferenceIterator<Iterator>(iterator);
}

} // namespace SFCGAL

#endif // SFCGAL_DEREFERENCEITERATOR_H_
