/*------------------------------------------------------------------------------
-- This file is a part of the FFT C++ library
-- Copyright (C) 2021, Plasma Physics Laboratory - CNRS
--
-- This program is free software; you can redistribute it and/or modify
-- it under the terms of the GNU General Public License as published by
-- the Free Software Foundation; either version 2 of the License, or
-- (at your option) any later version.
--
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
-- GNU General Public License for more details.
--
-- You should have received a copy of the GNU General Public License
-- along with this program; if not, write to the Free Software
-- Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
-------------------------------------------------------------------------------*/
/*-- Author : Alexis Jeandet
-- Mail : alexis.jeandet@member.fsf.org
----------------------------------------------------------------------------*/
#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

namespace FFT
{

template <typename T>
inline constexpr T log2(T v)
{
    int res = 0;
    for (auto i = 0UL; i < 64; i++)
    {
        if (T(1UL << i) == v)
        {
            res = i;
            break;
        }
    }
    return res;
}

template <std::size_t size, typename T>
inline constexpr std::array<std::complex<T>, size> generate_twiddle_factors()
{
    std::array<std::complex<T>, size> result;
    for (auto i = 1UL; i <= size; i++)
    {
        auto phi = -2. * M_PI / (1 << i);
        result[i - 1] = std::complex<T> { cos(phi), sin(phi) };
    }
    return result;
}

template <std::size_t size>
inline constexpr unsigned int bit_rev(unsigned int v)
{
    unsigned int result = 0;
    for (auto i = 0UL; i < size / 2; i++)
    {
        const auto shift = (size - 1 - 2 * i);
        const auto lower_mask = (1UL << i);
        const auto upper_mask = (1UL << (size - i - 1));
        result = result | (((v & lower_mask) << shift) | ((v & upper_mask) >> shift));
    }
    if (size / 2 * 2 != size)
    {
        result |= v & (1 << size / 2);
    }
    return result;
}

template <std::size_t size>
inline constexpr std::array<unsigned int, size> generate_reversed_indexes()
{
    std::array<unsigned int, size> result { 0 };
    for (auto i = 0UL; i < size; i++)
    {
        result[i] = bit_rev<log2(size)>(i);
    }
    return result;
}

template <std::size_t size, typename T>
struct FFT
{
    static constexpr auto twiddle_size = log2(size);
    static constexpr auto fft_size = size;

    constexpr FFT() { }

    std::vector<std::complex<T>> operator()(const std::vector<T>& input)
    {
        return fft(input);
    }

    // Cooley–Tukey FFT impl from https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
    static inline std::vector<std::complex<T>> fft(const std::vector<T>& input)
    {
        assert(std::size(input) == size);
        std::vector<std::complex<T>> result(size);
        for (auto i = 0UL; i < size; i++)
        {
            result[reversed_indexes[i]] = input[i];
        }
        for (auto s = 1UL; s <= twiddle_size; s++)
        {
            const auto m = 1UL << s;
            const auto w_m = twiddle_factors[s - 1];
            for (auto k = 0UL; k < size; k += m)
            {
                std::complex<T> w = 1.;
                for (auto j = 0UL; j < m / 2; j++)
                {
                    auto t = w * result[k + j + m / 2];
                    auto u = result[k + j];
                    result[k + j] = u + t;
                    result[k + j + m / 2] = u - t;
                    w = w * w_m;
                }
            }
        }
        return result;
    }

    static inline std::vector<T> mod(const std::vector<std::complex<T>>& input, bool remove_symetry=false)
    {
        assert(std::size(input) == size);
        std::vector<T> result(size/2+1);
        std::transform(std::cbegin(input), remove_symetry?std::cend(input)-(size/2-1):std::cend(input), std::begin(result),
            [](const auto& v) { return std::abs(v) / size; });
        return result;
    }


    static constexpr std::array<std::complex<T>, twiddle_size> twiddle_factors {
        generate_twiddle_factors<twiddle_size, T>()
    };
    static constexpr std::array<unsigned int, size> reversed_indexes {
        generate_reversed_indexes<size>()
    };

private:
};

} // namespace FFT
