%% "F846a_ST24R_IK_UR16e_PG53PK1PG8" UR16e Home POSE.
%
% It solves the INVERSE KINEMATICS for any desired position & orientation
% of the TCP (noap goal) of the Robot, this is tool REFerence.
% by Dr. Pardos-Gotor ST24R "Screw Theory Toolbox for Robotics" MATLAB.
%
% function ThetaOut = F846a_ST24R_IK_UR16e_PG53PK1PG8(u)
%
% The inputs "u" (7x1) are composed by the following vectors.
% "traXYZ" (3x1) desired translations for the TcP (noap - "p" goal).
% "rotXYZ" (3x1) desired rotations for TcP (noap - "noa" goal order (X+Y+Z)
% "Solutions" (1x1) is the value "Theta index" for choosing one out of 8
% possible results (Solutions = 1 to 8) for Robot Joint values
% Solutions = 9 is for sending the robot to the HOME POSITION (ThetaOut=0)
% "ThetaOut" (t1..t6)are the magnitudes solution for the Robot Joints1..6.
%
%
% The SIMSCAPE MULTIBODY SYSTEM used to build the Digital-Twin of the robot
% does not support direxcty Screw Theory and is still based on the DH
% formalism. For this reason, the revolute joints only can rotate on "Z"
% axis. Therefore, to construct the manipulator body we must use the
% Denavit-Hartenberg parameters (in this case the Classical DH - "DHC").
% Afterwards, Screw Theory POE is used to perform forward kinematics.
% For this robot the DHC parameters are:
% Joint-dtra
% J0       0      0     0       0
% J1-J6 = [0.181 t1     0     -90;
%          0     t2     0.478   0;
%          0     t3     0.36    0;
%          0.174 t4     0     -90;
%          0.12  t5     0      90;
%          0.19  t6+180 0       0];
% 
% Example of random trajectory generation for AUTOMATIC input testing.
% F110_CreateToolTra_Trig([0.3 -0.5 0.3 0 0 pi/4], [0.55 0.5 0.7 pi/4 pi/4 pi/2], 0.001, 10)
%
% Copyright (C) 2003-2020, by Dr. Jose M. Pardos-Gotor.
%
% This file is part of The ST24R "Screw Theory Toolbox for Robotics" MATLAB
% 
% ST24R is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% ST24R is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Leser General Public License
% along with ST24R. If not, see <http://www.gnu.org/licenses/>.
%
% http://www.
%
% CHANGES:
% Revision 1.1  2020/02/11 00:00:01
% General cleanup of code: help comments, see also, copyright
% references, clarification of functions.
%
%% F846a_ST24R_IK_UR16e_PG53PK1PG8
%
function ThetaOut = F846a_ST24R_IK_UR16e_PG53PK1PG8(u) %#codegen
%
%
% Mechanical characteristics of the Robot (AT REF HOME POSITION):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Joints TWISTS definition and TcP at home.
pg=[0.838; 0.174; 0.061]; pp=[0.838; 0.364; 0.061];
%
Twist = [0   -0.1810   -0.1810   -0.1810   -0.1740   -0.0610;
         0         0         0         0    0.8380         0;
         0         0    0.4780    0.8380         0    0.8380;
         0         0         0         0         0         0;
         0    1.0000    1.0000    1.0000         0    1.0000;
    1.0000         0         0         0   -1.0000         0];
Hst0 = [-1 0 0 0.838; 0 0 1 0.364; 0 1 0 0.061; 0 0 0 1];
%
% Motion RANGE for the robot joints POSITION rad, (by catalog).
% Thmax = pi/180*[360 360 360 360 360 360];
% Thmin = -pi/180*[360 360 360 360 360 360];
% Maximum SPEED for the robot joints rad/sec, (by catalog).
% Thpmax = pi/180*[120 120 180 180 180 180];
%
%CM1 = [0; 0.01; 0.15]; CM2 = [0.2; 0.174; 0.181]; CM3 = [0.628; 0.05; 0.181];
%CM4 = [0.838; 0.174; 0.144]; CM5 = [0.838; 0.194; 0.061];
%CM6 = [0.838; 0.274; 0.061];
%IT1 = [0.2; 0.2; 0.3]; IT2 = [0.1; 0.1; 0.2]; IT3 = [0.2; 0.1; 0.1];
%IT4 = [0.1; 0.1; 0.1]; IT5 = [0.1; 0.1; 0.1]; IT6 = [0.1; 0.1; 0.1];
%mass = [7.369 10.45 4.321 2.18 2.033 0.907];
%
%
% Calculate Homogeneous transformation for the GOAL "noap"
traXYZ = u(1:3); rotXYZ = u(4:6);
noap = rotX2tform(rotXYZ(1))*rotY2tform(rotXYZ(2))*rotZ2tform(rotXYZ(3));
noap(1,4)= traXYZ(1); noap(2,4)= traXYZ(2); noap(3,4)= traXYZ(3);
%
% number of joints.
n = 6;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IK solution approach PG5+PG3+PK1+PG8 subproblems cosecutively.
% ATTENTION, because this algorithm can give two, four or eight correct
% solutions, as a function of the target and singularities, then
% You must check with the FK which one is right out of the eight.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the IK solutions Theta using the SCREW THEORY.
Theta_STR4 = zeros(8,n);
%
% STEP1: Calculate Theta1.
% With "pg" on the axis of E5, E6. We apply (noap*hs0^-1) to "pg",
% because the screws E5,E6 do not affect "pg" for being on their axes
% and the E2,E3,E4 do not change the plane where "pg" moves, and so do not
% affect the calculation for Theta1 resulting the problem 
% "exp(E1^theta1)*pg = noap*hs0^-1*pg"
% which has one solution for t1 by PARDOS-GOTOR-FIVE.
% Pay attention to the mechanical configuration of the robot, because
% to comply with t12, we must take into consideration the displacement of
% point "pg" (i.e., d4) out of the robot plane of motion (the plane which 
% includes points o, k and r. Besides, for t12 we must add "pi" to the 
% calculation, otherwise it is impossible to reach the point g.
noapHst0ig = noap*(Hst0\[pg; 1]); pk1 = noapHst0ig(1:3);
t1 = PardosGotorFive(Twist(:,1), pg, pk1);
v1 = Twist(1:3,1); w1 = Twist(4:6,1); r1 = cross(w1,v1)/(norm(w1)^2);
v = pk1 - r1; vw1 = w1*w1'*v; vp1 = v-vw1; nvp = norm(vp1);
w = pg - r1; ww1 = w1*w1'*w; up1 = w-ww1; nup = norm(up1);
t11 = t1(1) - asin(pg(2)/nvp) + asin(pg(2)/nup);
t12 = t1(2) + asin(pg(2)/nvp) + asin(pg(2)/nup);
Theta_STR4(1:4,1) = real(t11);
Theta_STR4(5:8,1) = real(t12);
%
% STEP2: Calculate Theta5.
% With "pp" not in the axis of E5 apply E1^-1*noap*gs0^-1 to "pp"
% E6 does not affect "pp" because is in its axis. Then applying E5 to "pp" 
% knowing already Theta1 gives a new point "k2p", but we must consider the
% effect of E2, E3 and E4. These three parallel rotations only make point
% "k2p" move along axis "X" a certain amount. To calculate this magnitude
% we use the Pardos-Gotor-Three (point translation to a given distance to 
% another point). Where the point is "k2p" the distance is the radius of 
% the joint rotation Theta5 (i.e., norm(pp-pg)) to the point "pg"). Solving
% this PG3, we obtain point "k2",  resulting exactly a Canonic problem 
% PADEN-KAHAN-ONE, which has none or one solution for any Th1 known.
% In this case, pay attention to the fact that also -Th1 can be a valid.
%
noapHst0ip = noap*(Hst0\[pp; 1]);
for i = 1:4:5
    pk2ph = (expScrew([Twist(:,1);Theta_STR4(i,1)]))\noapHst0ip;
    pk2p = pk2ph(1:3);
    w7 = [1 0 0]; x7 = [w7 0 0 0]';
    t7 = PardosGotorThree(x7, [pk2p(1:2); 0], [pg(1:2); 0], norm(pp-pg));
    pk2 = pk2p+w7'*t7(2);
    t51 = PadenKahanOne(Twist(:,5), pp, pk2);
    Theta_STR4(i:i+1,5) = real(t51);
    Theta_STR4(i+2:i+3,5) = real(-t51);
end
%
% STEP3: Calculate Theta6.
% Another geometric formulation for the IK to get t6a (alternative).
ox = noap(1,2); oy = noap(2,2); nx = noap(1,1); ny = noap(2,1);
for i = 1:2:7
    s1 = sin(Theta_STR4(i,1)); c1 = cos(Theta_STR4(i,1));
    s5 = sin(Theta_STR4(i,5));
    t61a = atan2((ox*s1-oy*c1)/s5,(ny*c1-nx*s1)/s5);
    Theta_STR4(i:i+1,6) = real(t61a);
end
%
% STEP4: Calculate Theta2, Theta3, Theta4.
% We pass the exponential of t1 to the right-hand side of the 
% kinematics expression, resulting the formula: 
% E2 * E3 * E4 * E5 * E6 * Hst0 = E1^-1 * Hstt => E2 * E3 * E4 * Hp = Hk
% which is the expression for the PARDOS-GOTOR-EIGHT PG8 canonical problem.
% which has none, one or two triple solutions.
for i = 1:2:7
    Hp = (expScrew([Twist(:,6);Theta_STR4(i,6)]))*Hst0;
    Hp = (expScrew([Twist(:,5);Theta_STR4(i,5)]))*Hp;
    Hk = (expScrew([Twist(:,1);Theta_STR4(i,1)]))\noap;
    t234 = PardosGotorEight(Twist(:,2),Twist(:,3),Twist(:,4),Hp,Hk);
    Theta_STR4(i:i+1,2:4) = t234;
end
%
if u(7)==9
    ThetaOut = [0 0 0 0 0 0];
else
    ThetaOut = Theta_STR4(u(7),:);
end
%
end
%