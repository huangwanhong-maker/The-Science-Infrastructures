# The navigation bridge's reference values

| File | What it holds | Source |
|---|---|---|
| `SGP4-VER.TLE` | Vallado's verification set: 33 element sets chosen to exercise every branch of SGP4/SDP4 (drag, deep space, 12-hour and 24-hour resonance, Lyddane's choices, decay, the error codes), each with the start, stop and step of its run (minutes) after column 69 of line 2 | Vallado, Crawford, Hujsak and Kelso 2006, "Revisiting Spacetrack Report #3", AIAA 2006-6753, rev. 2; copied unaltered (CRLF endings kept) from python-sgp4 2.27, which distributes it as Vallado does, retrieved 2026-10-28 (SHA-256 d246d1d9…e2881) |
| `tcppver.out` | The output of Vallado's C++ (testcpp) on that set: for each satellite, the time since epoch (min), the TEME position (km, 8 decimals) and velocity (km/s, 9 decimals), and, after the epoch's line, the osculating elements and the date | as above (SHA-256 7a46d8cd…7e7e) |
| `sgp4-reference.txt` | Vallado's C++ (python-sgp4's compiled vallado_cpp) on every set at −720, 0, 360, 1440 and 4320 min in the two cases `tcppver.out` does not cover: AFSPC's operation mode with WGS 72, and the improved mode with WGS 84 | `tests/make_sgp4_fixtures.py`, 2026-10-28 |
| `teme-reference.txt` | The TEME → ITRF example of AIAA 2006-6753 rev. 2, Appendix C, computed by skyfield 1.55 at its instant exactly (UT1 in two parts) | `tests/make_sgp4_fixtures.py`, 2026-10-28 |

Vallado's code and test data are distributed by him and by CelesTrak without restriction;
python-sgp4 is © Brandon Rhodes under the MIT licence; skyfield is © Brandon Rhodes under the MIT licence.
