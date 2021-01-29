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
#include <FFT/FFT.hpp>

#define CATCH_CONFIG_MAIN
#if __has_include(<catch2/catch.hpp>)
#include <catch2/catch.hpp>
#else
#include <catch.hpp>
#endif

SCENARIO("FFT", "[real to img 8 points]")
{
    for (auto i = 0UL; i < 64; i++)
        REQUIRE(FFT::log2(1UL << i) == i);

    GIVEN("a direct 8 points FFT")
    {
        using fft_t = FFT::FFT<8, double>;
        fft_t fft;
        const std::vector<unsigned int> reversed_indexes { std::cbegin(fft_t::reversed_indexes),
            std::cend(fft_t::reversed_indexes) };
        REQUIRE_THAT(reversed_indexes, Catch::Equals<unsigned int>({ 0, 4, 2, 6, 1, 5, 3, 7 }));
        WHEN("doing a fft of 0")
        {
            auto fft_of_0 = fft({ 0., 0., 0., 0., 0., 0., 0., 0. });
            THEN("result should be all 0+0j")
            {
                // f'{{ {",".join([ ("{" + str(s.real) + "," + str(s.imag) + "}") for s in
                // np.fft.fft([0.]*8) ])}  }}'
                REQUIRE_THAT(fft_of_0,
                    Catch::Equals<decltype(fft_of_0)::value_type>(
                        { { 0.0, 0.0 }, { 0.0, 0.0 }, { 0.0, 0.0 }, { 0.0, 0.0 }, { 0.0, 0.0 },
                            { 0.0, 0.0 }, { 0.0, 0.0 }, { 0.0, 0.0 } }));
            }
        }
        WHEN("doing a fft of a dirac")
        {
            auto fft_of_dirac = fft({ 1., 0., 0., 0., 0., 0., 0., 0. });
            THEN("result should be all 1+0j")
            {
                // f'{{ {",".join([ ("{" + str(s.real) + "," + str(s.imag) + "}") for s in
                // np.fft.fft([1.]+[0.]*7) ])}  }}'
                REQUIRE_THAT(fft_of_dirac,
                    Catch::Equals<decltype(fft_of_dirac)::value_type>(
                        { { 1.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 0.0 },
                            { 1.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 0.0 } }));
            }
            THEN("modulus should be all 1./8.")
            {
                REQUIRE_THAT(fft.mod(fft_of_dirac),
                    Catch::Equals<double>({ 1. / 8., 1. / 8., 1. / 8.,
                        1. / 8., 1. / 8., 1. / 8., 1. / 8., 1. / 8. }));
            }
        }

        WHEN("doing a fft of ones")
        {
            auto fft_of_ones = fft({ 1., 1., 1., 1., 1., 1., 1., 1. });
            THEN("result should be fft'size + 0j only at index 0 IE 0Hz")
            {
                // f'{{ {",".join([ ("{" + str(s.real) + "," + str(s.imag) + "}") for s in
                // np.fft.fft([1.]*8) ])}  }}'
                REQUIRE_THAT(fft_of_ones,
                    Catch::Equals<decltype(fft_of_ones)::value_type>(
                        { { 8.0, 0.0 }, { 0.0, 0.0 }, { 0.0, 0.0 }, { 0.0, 0.0 }, { 0.0, 0.0 },
                            { 0.0, 0.0 }, { 0.0, 0.0 }, { 0.0, 0.0 } }));
            }
        }
    }
}
