// The series for the CIO locator, s + XY/2 (IAU 2006/2000A; Capitaine, Wallace and Chapront
// 2005, A&A 432, 355; IERS Conventions 2010, Table 5.2d), and the complementary terms of the
// equation of the equinoxes (Capitaine and Gontier 1993; IERS Conventions 2010, Table 5.2e).
// Coefficients in arcseconds; arguments are integer combinations of l, l', F, D, Ω, L_Ve, L_E, p_A.
//
// Transcribed mechanically by domains/metrology/tools/transcribe_erfa.py from ERFA's s06.c and eect00.c
// (https://github.com/liberfa/erfa, retrieved 2026-10-28). ERFA is distributed under the following
// terms, which this transcription keeps:
//
//   Copyright (C) 2013-2021, NumFOCUS Foundation. All rights reserved. This library is derived,
//   with permission, from the International Astronomical Union's "Standards of Fundamental
//   Astronomy" library, available from http://www.iausofa.org. Redistribution and use in source
//   and binary forms, with or without modification, are permitted provided that the following
//   conditions are met: (1) redistributions of source code must retain the above copyright
//   notice, this list of conditions and the following disclaimer; (2) redistributions in binary
//   form must reproduce the above copyright notice, this list of conditions and the following
//   disclaimer in the documentation and/or other materials provided with the distribution;
//   (3) neither the name of the Standards Of Fundamental Astronomy Board, the International
//   Astronomical Union nor the names of its contributors may be used to endorse or promote
//   products derived from this software without specific prior written permission. THIS
//   SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
//   IMPLIED WARRANTIES ARE DISCLAIMED.

#pragma once

namespace metrology::frames::series {

struct ArgumentTerm {
  int nfa[8];
  double s, c;  // sine and cosine coefficients, arcseconds
};

// The polynomial part of s + XY/2, arcseconds, in powers 0..5 of t (Julian centuries TT).
inline constexpr double s06_polynomial[6] = {94.00e-6, 3808.65e-6, -122.68e-6, -72574.11e-6, 27.98e-6, 15.62e-6};

inline constexpr ArgumentTerm s06_s0[33] = {
    {{0, 0, 0, 0, 1, 0, 0, 0}, -2640.73e-6, 0.39e-6},
    {{0, 0, 0, 0, 2, 0, 0, 0}, -63.53e-6, 0.02e-6},
    {{0, 0, 2, -2, 3, 0, 0, 0}, -11.75e-6, -0.01e-6},
    {{0, 0, 2, -2, 1, 0, 0, 0}, -11.21e-6, -0.01e-6},
    {{0, 0, 2, -2, 2, 0, 0, 0}, 4.57e-6, 0.00e-6},
    {{0, 0, 2, 0, 3, 0, 0, 0}, -2.02e-6, 0.00e-6},
    {{0, 0, 2, 0, 1, 0, 0, 0}, -1.98e-6, 0.00e-6},
    {{0, 0, 0, 0, 3, 0, 0, 0}, 1.72e-6, 0.00e-6},
    {{0, 1, 0, 0, 1, 0, 0, 0}, 1.41e-6, 0.01e-6},
    {{0, 1, 0, 0, -1, 0, 0, 0}, 1.26e-6, 0.01e-6},
    {{1, 0, 0, 0, -1, 0, 0, 0}, 0.63e-6, 0.00e-6},
    {{1, 0, 0, 0, 1, 0, 0, 0}, 0.63e-6, 0.00e-6},
    {{0, 1, 2, -2, 3, 0, 0, 0}, -0.46e-6, 0.00e-6},
    {{0, 1, 2, -2, 1, 0, 0, 0}, -0.45e-6, 0.00e-6},
    {{0, 0, 4, -4, 4, 0, 0, 0}, -0.36e-6, 0.00e-6},
    {{0, 0, 1, -1, 1, -8, 12, 0}, 0.24e-6, 0.12e-6},
    {{0, 0, 2, 0, 0, 0, 0, 0}, -0.32e-6, 0.00e-6},
    {{0, 0, 2, 0, 2, 0, 0, 0}, -0.28e-6, 0.00e-6},
    {{1, 0, 2, 0, 3, 0, 0, 0}, -0.27e-6, 0.00e-6},
    {{1, 0, 2, 0, 1, 0, 0, 0}, -0.26e-6, 0.00e-6},
    {{0, 0, 2, -2, 0, 0, 0, 0}, 0.21e-6, 0.00e-6},
    {{0, 1, -2, 2, -3, 0, 0, 0}, -0.19e-6, 0.00e-6},
    {{0, 1, -2, 2, -1, 0, 0, 0}, -0.18e-6, 0.00e-6},
    {{0, 0, 0, 0, 0, 8, -13, -1}, 0.10e-6, -0.05e-6},
    {{0, 0, 0, 2, 0, 0, 0, 0}, -0.15e-6, 0.00e-6},
    {{2, 0, -2, 0, -1, 0, 0, 0}, 0.14e-6, 0.00e-6},
    {{0, 1, 2, -2, 2, 0, 0, 0}, 0.14e-6, 0.00e-6},
    {{1, 0, 0, -2, 1, 0, 0, 0}, -0.14e-6, 0.00e-6},
    {{1, 0, 0, -2, -1, 0, 0, 0}, -0.14e-6, 0.00e-6},
    {{0, 0, 4, -2, 4, 0, 0, 0}, -0.13e-6, 0.00e-6},
    {{0, 0, 2, -2, 4, 0, 0, 0}, 0.11e-6, 0.00e-6},
    {{1, 0, -2, 0, -3, 0, 0, 0}, -0.11e-6, 0.00e-6},
    {{1, 0, -2, 0, -1, 0, 0, 0}, -0.11e-6, 0.00e-6},
};

inline constexpr ArgumentTerm s06_s1[3] = {
    {{0, 0, 0, 0, 2, 0, 0, 0}, -0.07e-6, 3.57e-6},
    {{0, 0, 0, 0, 1, 0, 0, 0}, 1.73e-6, -0.03e-6},
    {{0, 0, 2, -2, 3, 0, 0, 0}, 0.00e-6, 0.48e-6},
};

inline constexpr ArgumentTerm s06_s2[25] = {
    {{0, 0, 0, 0, 1, 0, 0, 0}, 743.52e-6, -0.17e-6},
    {{0, 0, 2, -2, 2, 0, 0, 0}, 56.91e-6, 0.06e-6},
    {{0, 0, 2, 0, 2, 0, 0, 0}, 9.84e-6, -0.01e-6},
    {{0, 0, 0, 0, 2, 0, 0, 0}, -8.85e-6, 0.01e-6},
    {{0, 1, 0, 0, 0, 0, 0, 0}, -6.38e-6, -0.05e-6},
    {{1, 0, 0, 0, 0, 0, 0, 0}, -3.07e-6, 0.00e-6},
    {{0, 1, 2, -2, 2, 0, 0, 0}, 2.23e-6, 0.00e-6},
    {{0, 0, 2, 0, 1, 0, 0, 0}, 1.67e-6, 0.00e-6},
    {{1, 0, 2, 0, 2, 0, 0, 0}, 1.30e-6, 0.00e-6},
    {{0, 1, -2, 2, -2, 0, 0, 0}, 0.93e-6, 0.00e-6},
    {{1, 0, 0, -2, 0, 0, 0, 0}, 0.68e-6, 0.00e-6},
    {{0, 0, 2, -2, 1, 0, 0, 0}, -0.55e-6, 0.00e-6},
    {{1, 0, -2, 0, -2, 0, 0, 0}, 0.53e-6, 0.00e-6},
    {{0, 0, 0, 2, 0, 0, 0, 0}, -0.27e-6, 0.00e-6},
    {{1, 0, 0, 0, 1, 0, 0, 0}, -0.27e-6, 0.00e-6},
    {{1, 0, -2, -2, -2, 0, 0, 0}, -0.26e-6, 0.00e-6},
    {{1, 0, 0, 0, -1, 0, 0, 0}, -0.25e-6, 0.00e-6},
    {{1, 0, 2, 0, 1, 0, 0, 0}, 0.22e-6, 0.00e-6},
    {{2, 0, 0, -2, 0, 0, 0, 0}, -0.21e-6, 0.00e-6},
    {{2, 0, -2, 0, -1, 0, 0, 0}, 0.20e-6, 0.00e-6},
    {{0, 0, 2, 2, 2, 0, 0, 0}, 0.17e-6, 0.00e-6},
    {{2, 0, 2, 0, 2, 0, 0, 0}, 0.13e-6, 0.00e-6},
    {{2, 0, 0, 0, 0, 0, 0, 0}, -0.13e-6, 0.00e-6},
    {{1, 0, 2, -2, 2, 0, 0, 0}, -0.12e-6, 0.00e-6},
    {{0, 0, 2, 0, 0, 0, 0, 0}, -0.11e-6, 0.00e-6},
};

inline constexpr ArgumentTerm s06_s3[4] = {
    {{0, 0, 0, 0, 1, 0, 0, 0}, 0.30e-6, -23.42e-6},
    {{0, 0, 2, -2, 2, 0, 0, 0}, -0.03e-6, -1.46e-6},
    {{0, 0, 2, 0, 2, 0, 0, 0}, -0.01e-6, -0.25e-6},
    {{0, 0, 0, 0, 2, 0, 0, 0}, 0.00e-6, 0.23e-6},
};

inline constexpr ArgumentTerm s06_s4[1] = {
    {{0, 0, 0, 0, 1, 0, 0, 0}, -0.26e-6, -0.01e-6},
};

inline constexpr ArgumentTerm eect00_e0[33] = {
    {{0, 0, 0, 0, 1, 0, 0, 0}, 2640.96e-6, -0.39e-6},
    {{0, 0, 0, 0, 2, 0, 0, 0}, 63.52e-6, -0.02e-6},
    {{0, 0, 2, -2, 3, 0, 0, 0}, 11.75e-6, 0.01e-6},
    {{0, 0, 2, -2, 1, 0, 0, 0}, 11.21e-6, 0.01e-6},
    {{0, 0, 2, -2, 2, 0, 0, 0}, -4.55e-6, 0.00e-6},
    {{0, 0, 2, 0, 3, 0, 0, 0}, 2.02e-6, 0.00e-6},
    {{0, 0, 2, 0, 1, 0, 0, 0}, 1.98e-6, 0.00e-6},
    {{0, 0, 0, 0, 3, 0, 0, 0}, -1.72e-6, 0.00e-6},
    {{0, 1, 0, 0, 1, 0, 0, 0}, -1.41e-6, -0.01e-6},
    {{0, 1, 0, 0, -1, 0, 0, 0}, -1.26e-6, -0.01e-6},
    {{1, 0, 0, 0, -1, 0, 0, 0}, -0.63e-6, 0.00e-6},
    {{1, 0, 0, 0, 1, 0, 0, 0}, -0.63e-6, 0.00e-6},
    {{0, 1, 2, -2, 3, 0, 0, 0}, 0.46e-6, 0.00e-6},
    {{0, 1, 2, -2, 1, 0, 0, 0}, 0.45e-6, 0.00e-6},
    {{0, 0, 4, -4, 4, 0, 0, 0}, 0.36e-6, 0.00e-6},
    {{0, 0, 1, -1, 1, -8, 12, 0}, -0.24e-6, -0.12e-6},
    {{0, 0, 2, 0, 0, 0, 0, 0}, 0.32e-6, 0.00e-6},
    {{0, 0, 2, 0, 2, 0, 0, 0}, 0.28e-6, 0.00e-6},
    {{1, 0, 2, 0, 3, 0, 0, 0}, 0.27e-6, 0.00e-6},
    {{1, 0, 2, 0, 1, 0, 0, 0}, 0.26e-6, 0.00e-6},
    {{0, 0, 2, -2, 0, 0, 0, 0}, -0.21e-6, 0.00e-6},
    {{0, 1, -2, 2, -3, 0, 0, 0}, 0.19e-6, 0.00e-6},
    {{0, 1, -2, 2, -1, 0, 0, 0}, 0.18e-6, 0.00e-6},
    {{0, 0, 0, 0, 0, 8, -13, -1}, -0.10e-6, 0.05e-6},
    {{0, 0, 0, 2, 0, 0, 0, 0}, 0.15e-6, 0.00e-6},
    {{2, 0, -2, 0, -1, 0, 0, 0}, -0.14e-6, 0.00e-6},
    {{1, 0, 0, -2, 1, 0, 0, 0}, 0.14e-6, 0.00e-6},
    {{0, 1, 2, -2, 2, 0, 0, 0}, -0.14e-6, 0.00e-6},
    {{1, 0, 0, -2, -1, 0, 0, 0}, 0.14e-6, 0.00e-6},
    {{0, 0, 4, -2, 4, 0, 0, 0}, 0.13e-6, 0.00e-6},
    {{0, 0, 2, -2, 4, 0, 0, 0}, -0.11e-6, 0.00e-6},
    {{1, 0, -2, 0, -3, 0, 0, 0}, 0.11e-6, 0.00e-6},
    {{1, 0, -2, 0, -1, 0, 0, 0}, 0.11e-6, 0.00e-6},
};

inline constexpr ArgumentTerm eect00_e1[1] = {
    {{0, 0, 0, 0, 1, 0, 0, 0}, -0.87e-6, 0.00e-6},
};

}  // namespace metrology::frames::series
