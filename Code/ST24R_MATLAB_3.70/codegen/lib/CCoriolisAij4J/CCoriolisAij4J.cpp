//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: CCoriolisAij4J.cpp
//
// MATLAB Coder version            : 5.1
// C/C++ source code generated on  : 14-Nov-2020 16:02:13
//

// Include Files
#include "CCoriolisAij4J.h"
#include "Christoffel.h"
#include "rt_nonfinite.h"
#include <cstring>

// Function Definitions
//
// Arguments    : const double TwMag[28]
//                const double LiMas[28]
//                const double Thetap[4]
//                double Ctdt[16]
// Return Type  : void
//
void CCoriolisAij4J(const double TwMag[28], const double LiMas[28], const double
                    Thetap[4], double Ctdt[16])
{
  double Cosij;
  double d;
  double d1;
  double d2;
  double d3;
  double d4;

  //  "CoriolisAij" CORIOLIS matrix C(t,dt) for an open chain manipulator.
  //  computation based on the Adjoint Transformation Aij.
  //  Use in SE(3).
  //
  //  	Ctdt = CCoriolisAij(TwMag,LiMas,Thetap)
  //
  //  Gives the CORIOLIS MATRIX C(t,dt) for the Lagangian's equations:
  //  M(t)*ddt + C(t,dt)*dt + N(t,dt) = T
  //  of the dynamics of the robot formed by links on an open chain.
  //
  //  INPUTS:
  //  TwMag = [Tw1; Mag1, ..., Twn; Magn] (7xn)
  //  for each rigid body joint (1..n).
  //  Tw1..Twn: The TWIST components for the joint movement.
  //  Mag1..Magn: The MAGNITUDE component for the joint SCREW movement.
  //  LiMas = [CM1; IT1; Mas1, ..., CMn; ITn; Masn] (7xn)
  //  for each rigid body link (1..n).
  //  CM1..CMn: Center of Mass x-y-z Position components to S for each Link.
  //  IT1..ITn: Inertia x-y-z components for each Link refered to its CM.
  //  Mas1..Masn: The Mass for each Link.
  //  Thetap = Vector differential Theta dt (1xn) of the joint VELOCITIES.
  //
  //      |C11...C1n|
  //  C = |         |, With Cij = 1/2 * Sum(1,n)[( dMij/dtk + dMik/dtj - dMkj/dti ) * dtk]   
  //      |Cn1...Cnn|
  //  where:
  //  Cristoffel Symbols are defined by:
  //  dMij/dtk = Sum(l=max(i,j),n)[[Aki*Ei,Ek]'*Alk'*Ml*Alj*Ej + Ei'*Ali'*Ml*Alk*[Akj*Ej,Ek]] 
  //
  //  With Ml being the link inertia nxn LinkInertiaS.
  //  With Ei being the twist xi 6x1.
  //  With Aij being an element 6x6 of the adjoint transformation Aij2adjoint
  //
  //  See also: LinkInertiaS, Aij2adjoint.
  //
  //  Copyright (C) 2001-2020, by Dr. Jose M. Pardos-Gotor.
  //
  //  This file is part of The ST24R "Screw Theory Toolbox for Robotics" MATLAB
  //
  //  ST24R is free software: you can redistribute it and/or modify
  //  it under the terms of the GNU Lesser General Public License as published
  //  by the Free Software Foundation, either version 3 of the License, or
  //  (at your option) any later version.
  //
  //  ST24R is distributed in the hope that it will be useful,
  //  but WITHOUT ANY WARRANTY; without even the implied warranty of
  //  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  //  GNU Lesser General Public License for more details.
  //
  //  You should have received a copy of the GNU Leser General Public License
  //  along with ST24R.  If not, see <http://www.gnu.org/licenses/>.
  //
  //  http://www.
  //
  //  CHANGES:
  //  Revision 1.1  2020/02/11 00:00:01
  //  General cleanup of code: help comments, see also, copyright
  //  references, clarification of functions.
  //
  //  Ctdt = CCoriolisAij(TwMag,LiMas,Thetap)
  //
  std::memset(&Ctdt[0], 0, 16U * sizeof(double));
  d = Thetap[0];
  d1 = Thetap[1];
  d2 = Thetap[2];
  d3 = Thetap[3];
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      int b_i;
      b_i = i + (j << 2);
      Cosij = Christoffel(TwMag, LiMas, static_cast<double>(i) + 1.0,
                          static_cast<double>(j) + 1.0, 1.0);
      Cosij += Christoffel(TwMag, LiMas, static_cast<double>(i) + 1.0, 1.0,
                           static_cast<double>(j) + 1.0);
      Cosij -= Christoffel(TwMag, LiMas, 1.0, static_cast<double>(j) + 1.0,
                           static_cast<double>(i) + 1.0);
      Cosij *= d;
      d4 = Ctdt[b_i] + Cosij;
      Cosij = Christoffel(TwMag, LiMas, static_cast<double>(i) + 1.0,
                          static_cast<double>(j) + 1.0, 2.0);
      Cosij += Christoffel(TwMag, LiMas, static_cast<double>(i) + 1.0, 2.0,
                           static_cast<double>(j) + 1.0);
      Cosij -= Christoffel(TwMag, LiMas, 2.0, static_cast<double>(j) + 1.0,
                           static_cast<double>(i) + 1.0);
      Cosij *= d1;
      d4 += Cosij;
      Cosij = Christoffel(TwMag, LiMas, static_cast<double>(i) + 1.0,
                          static_cast<double>(j) + 1.0, 3.0);
      Cosij += Christoffel(TwMag, LiMas, static_cast<double>(i) + 1.0, 3.0,
                           static_cast<double>(j) + 1.0);
      Cosij -= Christoffel(TwMag, LiMas, 3.0, static_cast<double>(j) + 1.0,
                           static_cast<double>(i) + 1.0);
      Cosij *= d2;
      d4 += Cosij;
      Cosij = Christoffel(TwMag, LiMas, static_cast<double>(i) + 1.0,
                          static_cast<double>(j) + 1.0, 4.0);
      Cosij += Christoffel(TwMag, LiMas, static_cast<double>(i) + 1.0, 4.0,
                           static_cast<double>(j) + 1.0);
      Cosij -= Christoffel(TwMag, LiMas, 4.0, static_cast<double>(j) + 1.0,
                           static_cast<double>(i) + 1.0);
      Cosij *= d3;
      d4 += Cosij;
      Ctdt[b_i] = 0.5 * d4;
    }
  }

  //
}

//
// File trailer for CCoriolisAij4J.cpp
//
// [EOF]
//