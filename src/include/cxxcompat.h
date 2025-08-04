#ifndef _CXXCOMPAT_H
#define _CXXCOMPAT_H

#include "config.h"

#if HAVE_PRINT_HEADER
#include <print>
#else
#include <format>
#include <iostream>

namespace std
{
template <typename... Args>
inline void
print (const std::format_string<Args...> fmt, Args &&...args)
{
  std::cout << std::vformat (fmt.get (), std::make_format_args (args...));
}

template <typename... Args>
inline void
println (const std::format_string<Args...> fmt, Args &&...args)
{
  std::cout << std::vformat (fmt.get (), std::make_format_args (args...))
            << '\n';
}
};
#endif

#if !HAVE_STD_EXPECTED
namespace std
{

template <class E> class unexpected
{
public:
  constexpr
  unexpected (const E &e)
      : val_ (e)
  {
  }
  constexpr
  unexpected (E &&e)
      : val_ (std::move (e))
  {
  }

  constexpr const E &
  value () const &
  {
    return val_;
  }
  constexpr E &
  value () &
  {
    return val_;
  }
  constexpr const E &&
  value () const &&
  {
    return std::move (val_);
  }
  constexpr E &&
  value () &&
  {
    return std::move (val_);
  }

private:
  E val_;
};

template <class T, class E> class expected
{
public:
  using value_type = T;
  using error_type = E;
  using unexpected_type = unexpected<E>;

  constexpr
  expected ()
      : has_val_ (true)
  {
    new (&val_) T{};
  }

  constexpr
  expected (const T &val)
      : has_val_ (true)
  {
    new (&val_) T (val);
  }

  constexpr
  expected (T &&val)
      : has_val_ (true)
  {
    new (&val_) T (std::move (val));
  }

  constexpr
  expected (const unexpected<E> &unex)
      : has_val_ (false)
  {
    new (&unex_) unexpected<E> (unex);
  }

  constexpr
  expected (unexpected<E> &&unex)
      : has_val_ (false)
  {
    new (&unex_) unexpected<E> (std::move (unex));
  }

  constexpr
  expected (const expected &other)
      : has_val_ (other.has_val_)
  {
    if (has_val_)
      new (&val_) T (other.val_);
    else
      new (&unex_) unexpected<E> (other.unex_);
  }

  constexpr
  expected (expected &&other)
      : has_val_ (other.has_val_)
  {
    if (has_val_)
      new (&val_) T (std::move (other.val_));
    else
      new (&unex_) unexpected<E> (std::move (other.unex_));
  }

  ~expected ()
  {
    if (has_val_)
      val_.~T ();
    else
      unex_.~unexpected<E> ();
  }

  // Assignment operators
  expected &
  operator= (const expected &other)
  {
    if (this != &other)
      {
        this->~expected ();
        new (this) expected (other);
      }
    return *this;
  }

  expected &
  operator= (expected &&other)
  {
    if (this != &other)
      {
        this->~expected ();
        new (this) expected (std::move (other));
      }
    return *this;
  }

  // Observers
  constexpr bool
  has_value () const noexcept
  {
    return has_val_;
  }
  constexpr explicit
  operator bool () const noexcept
  {
    return has_val_;
  }

  // Value access
  constexpr T &
  value () &
  {
    if (!has_val_)
      throw std::runtime_error ("bad expected access");
    return val_;
  }

  constexpr const T &
  value () const &
  {
    if (!has_val_)
      throw std::runtime_error ("bad expected access");
    return val_;
  }

  constexpr T &&
  value () &&
  {
    if (!has_val_)
      throw std::runtime_error ("bad expected access");
    return std::move (val_);
  }

  constexpr const T &&
  value () const &&
  {
    if (!has_val_)
      throw std::runtime_error ("bad expected access");
    return std::move (val_);
  }

  constexpr T &
  operator* () &
  {
    return val_;
  }
  constexpr const T &
  operator* () const &
  {
    return val_;
  }
  constexpr T &&
  operator* () &&
  {
    return std::move (val_);
  }
  constexpr const T &&
  operator* () const &&
  {
    return std::move (val_);
  }

  constexpr T *
  operator->()
  {
    return &val_;
  }
  constexpr const T *
  operator->() const
  {
    return &val_;
  }

  // Error access
  constexpr E &
  error () &
  {
    return unex_.value ();
  }
  constexpr const E &
  error () const &
  {
    return unex_.value ();
  }
  constexpr E &&
  error () &&
  {
    return std::move (unex_).value ();
  }
  constexpr const E &&
  error () const &&
  {
    return std::move (unex_).value ();
  }

  // Value or alternative
  template <class U>
  constexpr T
  value_or (U &&default_value) const &
  {
    return has_val_ ? val_ : static_cast<T> (std::forward<U> (default_value));
  }

  template <class U>
  constexpr T
  value_or (U &&default_value) &&
  {
    return has_val_ ? std::move (val_)
                    : static_cast<T> (std::forward<U> (default_value));
  }

private:
  bool has_val_;
  union
  {
    T val_;
    unexpected<E> unex_;
  };
};

// Specialization for void
template <class E> class expected<void, E>
{
public:
  using value_type = void;
  using error_type = E;
  using unexpected_type = unexpected<E>;

  constexpr
  expected ()
      : void_ (), has_val_ (true)
  {
  }

  constexpr
  expected (const unexpected<E> &unex)
      : has_val_ (false), unex_ (unex)
  {
  }
  constexpr
  expected (unexpected<E> &&unex)
      : has_val_ (false), unex_ (std::move (unex))
  {
  }

  constexpr ~expected (){};

  constexpr bool
  has_value () const noexcept
  {
    return has_val_;
  }

  constexpr explicit
  operator bool () const noexcept
  {
    return has_val_;
  }

  constexpr void
  value () const
  {
    if (!has_val_)
      throw std::runtime_error ("bad expected access");
  }

  constexpr E &
  error () &
  {
    return unex_.value ();
  }

  constexpr const E &
  error () const &
  {
    return unex_.value ();
  }

  constexpr E &&
  error () &&
  {
    return std::move (unex_).value ();
  }

  constexpr const E &&
  error () const &&
  {
    return std::move (unex_).value ();
  }

private:
  bool has_val_;
  union
  {
    struct
    {
    } void_;
    unexpected<E> unex_;
  };
};
}
#endif

#endif
