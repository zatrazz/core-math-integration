//
// Copyright (c) Adhemerval Zanella. All rights reserved.
// Licensed under the MIT license. See LICENSE file in the project root for
// details.
//

#ifndef _CXXCOMPAT_H
#define _CXXCOMPAT_H

#include "config.h"

// The fallbacks for the case the compiler does not provide required C++
// suppo do not fully conform to the standard, neither provide all
// support (including concepts support).

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
  E val;

public:
  constexpr
  unexpected (const E &e)
      : val (e)
  {
  }
  constexpr
  unexpected (E &&e)
      : val (std::move (e))
  {
  }

  constexpr const E &
  value () const &
  {
    return val;
  }
  constexpr E &
  value () &
  {
    return val;
  }
  constexpr const E &&
  value () const &&
  {
    return std::move (val);
  }
  constexpr E &&
  value () &&
  {
    return std::move (val);
  }
};

template <class T, class E> class expected
{
public:
  using value_type = T;
  using error_type = E;
  using unexpected_type = unexpected<E>;

  constexpr
  expected ()
      : has_val (true)
  {
    new (&val) T{};
  }

  constexpr
  expected (const T &v)
      : has_val (true)
  {
    new (&val) T (v);
  }

  constexpr
  expected (T &&v)
      : has_val (true)
  {
    new (&val) T (std::move (v));
  }

  constexpr
  expected (const unexpected<E> &u)
      : has_val (false)
  {
    new (&unex) unexpected<E> (u);
  }

  constexpr
  expected (unexpected<E> &&un)
      : has_val (false)
  {
    new (&unex) unexpected<E> (std::move (un));
  }

  constexpr
  expected (const expected &other)
      : has_val (other.has_val)
  {
    if (has_val)
      new (&val) T (other.val);
    else
      new (&unex) unexpected<E> (other.unex);
  }

  constexpr
  expected (expected &&other)
      : has_val (other.has_val)
  {
    if (has_val)
      new (&val) T (std::move (other.val));
    else
      new (&unex) unexpected<E> (std::move (other.unex));
  }

  ~expected ()
  {
    if (has_val)
      val.~T ();
    else
      unex.~unexpected<E> ();
  }

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

  constexpr bool
  has_value () const noexcept
  {
    return has_val;
  }
  constexpr explicit
  operator bool () const noexcept
  {
    return has_val;
  }

  constexpr T &
  value () &
  {
    if (!has_val)
      throw std::runtime_error ("bad expected access");
    return val;
  }

  constexpr const T &
  value () const &
  {
    if (!has_val)
      throw std::runtime_error ("bad expected access");
    return val;
  }

  constexpr T &&
  value () &&
  {
    if (!has_val)
      throw std::runtime_error ("bad expected access");
    return std::move (val);
  }

  constexpr const T &&
  value () const &&
  {
    if (!has_val)
      throw std::runtime_error ("bad expected access");
    return std::move (val);
  }

  constexpr T &
  operator* () &
  {
    return val;
  }
  constexpr const T &
  operator* () const &
  {
    return val;
  }
  constexpr T &&
  operator* () &&
  {
    return std::move (val);
  }
  constexpr const T &&
  operator* () const &&
  {
    return std::move (val);
  }

  constexpr T *
  operator->()
  {
    return &val;
  }
  constexpr const T *
  operator->() const
  {
    return &val;
  }

  constexpr E &
  error () &
  {
    return unex.value ();
  }
  constexpr const E &
  error () const &
  {
    return unex.value ();
  }
  constexpr E &&
  error () &&
  {
    return std::move (unex).value ();
  }
  constexpr const E &&
  error () const &&
  {
    return std::move (unex).value ();
  }

  template <class U>
  constexpr T
  value_or (U &&default_value) const &
  {
    return has_val ? val : static_cast<T> (std::forward<U> (default_value));
  }

  template <class U>
  constexpr T
  value_or (U &&default_value) &&
  {
    return has_val ? std::move (val)
                    : static_cast<T> (std::forward<U> (default_value));
  }

private:
  bool has_val;
  union
  {
    T val;
    unexpected<E> unex;
  };
};

// template specialization for void, required by checkulps
template <class E> class expected<void, E>
{
public:
  using value_type = void;
  using error_type = E;
  using unexpected_type = unexpected<E>;

  constexpr
  expected ()
      : void_ (), has_val (true)
  {
  }

  constexpr
  expected (const unexpected<E> &unex)
      : has_val (false), unex (unex)
  {
  }
  constexpr
  expected (unexpected<E> &&unex)
      : has_val (false), unex (std::move (unex))
  {
  }

  constexpr ~expected (){};

  constexpr bool
  has_value () const noexcept
  {
    return has_val;
  }

  constexpr explicit
  operator bool () const noexcept
  {
    return has_val;
  }

  constexpr void
  value () const
  {
    if (!has_val)
      throw std::runtime_error ("bad expected access");
  }

  constexpr E &
  error () &
  {
    return unex.value ();
  }

  constexpr const E &
  error () const &
  {
    return unex.value ();
  }

  constexpr E &&
  error () &&
  {
    return std::move (unex).value ();
  }

  constexpr const E &&
  error () const &&
  {
    return std::move (unex).value ();
  }

private:
  bool has_val;
  union
  {
    struct
    {
    } void_;
    unexpected<E> unex;
  };
};
}
#endif

#endif
