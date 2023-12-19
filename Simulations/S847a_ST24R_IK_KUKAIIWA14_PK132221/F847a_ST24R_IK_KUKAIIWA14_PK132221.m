%% "F857a_ST24R_IK_KUKAIIWA14_PK132221" KUKA IIWA14 Straight POSE.
%
% It solves the INVERSE KINEMATICS for any desired position & orientation
% of the TCP (noap goal) of the Robot, this is tool REFerence.
% by Dr. Pardos-Gotor ST24R "Screw Theory Toolbox for Robotics" MATLAB.
%
% function ThetaOut = F857a_ST24R_IK_KUKAIIWA14_PK132221(u)
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
% J1-J6 = [0.36  t1     0     -90;
%          0     t2     0      90;
%          0.42  t3     0      90;
%          0     t4     0     -90;
%          0.4   t5     0     -90;
%          0     t6     0      90;
%          0.2   t7     0       0];
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
%% F857a_ST24R_IK_KUKAIIWA14_PK132221
%
function ThetaOut = F857a_ST24R_IK_KUKAIIWA14_PK132221(u) %#codegen
%
%
% Mechanical characteristics of the Robot (AT REF HOME POSITION):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Joints TWISTS definition and TcP at home.
pk=[0;0;0.36]; pf=[0;0;1.18]; pp=[0;0;1.38];
% Joints TWISTS definition and TcP at home.
Twist = [0   -0.3600         0    0.7800         0   -1.1800         0;
         0         0         0         0         0         0         0;
         0         0         0         0         0         0         0;
         0         0         0         0         0         0         0;
         0    1.0000         0   -1.0000         0    1.0000         0;
    1.0000         0    1.0000         0    1.0000         0     1.000];
Hst0 = [1 0 0 0; 0 1 0 0; 0 0 1 1.38; 0 0 0 1];
% 
%Motion RANGE for the robot joints POSITION rad, (by catalog).
%Thmax = pi/180*[170 120 170 120 170 120 175];
%Thmin = -pi/180*[170 120 170 120 170 120 175];
% Maximum SPEED for the robot joints rad/sec, (by catalog).
%Thpmax = pi/180*[85 85 100 75 130 135 135];
%
%CM1 = [0; -0.03; 0.2775]; CM2 = [0; 0.042; 0.419]; CM3 = [0; 0.03; 0.6945];
%CM4 = [0; -0.034; 0.847]; CM5 = [0; -0.021; 1];
%CM6 = [0; 0.001; 1.18]; CM7 = [0; 0; 1.28];
%IT1 = [0.1; 0.09; 0.02]; IT2 = [0.018; 0.05; 0.044];
%IT3 = [0.08; 0.075; 0.01];
%IT4 = [0.03; 0.029; 0.01]; IT5 = [0.02; 0.018; 0.005];
%IT6 = [0.005; 0.0036; 0.0047]; IT7 = [0.001; 0.001; 0.001];
%mass = [4 4 3 2.7 1.7 1.8 0.3];
%
% Calculate Homogeneous transformation for the GOAL "noap"
traXYZ = u(1:3); rotXYZ = u(4:6);
noap = rotX2tform(rotXYZ(1))*rotY2tform(rotXYZ(2))*rotZ2tform(rotXYZ(3));
noap(1,4)= traXYZ(1); noap(2,4)= traXYZ(2); noap(3,4)= traXYZ(3);
%
% number of joints.
n = 7;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ABBIRB120 IK solution approach PK3+PK2+PK2+PK1 subproblems cosecutively.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the IK solutions Theta using the SCREW THEORY.
Theta = zeros(16,n);
%
% First, as this is REDUNDANT robot with 7DoF and theoretical infinite
% solutions, we reduce the number of possibilies to 16, and for doing so
% we give some solution INPUTS for t1 and t3 based on some other criteria,
% and then well get 8 solutions for each of these INPUTS for the rest 6DoF.
% For this function we choose t1 and t3 to be the magnitudes for orientate
% these joint towards the TcP, and for doing so we use PADEN-KAHAN-ONE.
t1 = PadenKahanOne(Twist(:,1), [1; 0; 0], noap(1:3,4));
t3 = PadenKahanOne(Twist(:,3), [1; 0; 0], noap(1:3,4));
Theta(1:8,3) = t3; Theta(9:16,1) = t1;
%
% STEP2-1: Calculate Theta4.
% With "pf" on the axis of E5, E6, E7 and "pk" on the axis of E1, E2, E3.
% We apply (noap*hst0^-1) to "pf" and take the norm of the diffence of that
% resulting point and "pk". Doing so we can calculate Theta4 applying the
% Canonic problem PADEN-KAHAN-THREE, because the screws E4,E5,E6 do not 
% affect "pf" and the E1,E2,E3 do not affect the norm of a vector with an
% end on "pk" resulting the problem 
% ||exp(E4^t4)*pf-pk|| = ||noap*hst0^-1*pf-pk|| = ||pk1p-pk|| = de
% which by PADEN-KAHAN-THREE has none, one or two solutions
% t401 & t402.
noapHst0if = noap*(Hst0\[pf; 1]); pk1p = noapHst0if(1:3);
de = norm(pk1p - pk);
t4 = PadenKahanThree(Twist(:,4), pf, pk, de);
Theta(1:4,4) = t4(1); Theta(9:12,4) = t4(1);
Theta(5:8,4) = t4(2); Theta(13:16,4) = t4(2);
%
% STEP2-2: Calculate Theta1 & Theta2, knowing Theta3 (t3in).
% With "pf" on the axis of E5, E6, E7 we apply (noap*hsts0^-1) to "pf" and
% the POE E1..E7 also to "pf" having already known the value for t4 & t3
% resulting exactly a Canonic problem PADEN-KAHAN-TWO, because the screws
% E5,E6,E7 do not affect "pf" and the E4 & E3 are known,resulting
% the problem exp(E1^theta1)*exp(E2^theta2)*pf' = noap*hst0^-1*pf = pk1p
% which by PADEN-KAHAN-TWO has none, one or two DOUBLE solutions
% t401 & t3in => t201-t101 & t202-t102.
% t402 & t3in => t203-t103 & t204-t104. 
for i = 1:4:5
    pf1pt = expScrew([Twist(:,3);Theta(i,3)])*expScrew([Twist(:,4);Theta(i,4)])*[pf; 1];
    pf1p = pf1pt(1:3);
    t1t2 = PadenKahanTwo(Twist(:,1),Twist(:,2),pf1p,pk1p);
    Theta(i,1:2) = t1t2(1,1:2); Theta(i+1,1:2) = t1t2(1,1:2);
    Theta(i+2,1:2) = t1t2(2,1:2); Theta(i+3,1:2) = t1t2(2,1:2);
end
%
% STEP2-3: Calculate Theta2 & Theta3, knowing Theta1 (t1in).
% With "pf" on the axis of E5,E6,E7 we apply E1^-1*(noap*hsts0^-1) to "pf"
% and POE E2..E7 also to "pf" having already known the value for t4 & t1
% resulting exactly a Canonic problem PADEN-KAHAN-TWO, because the screws
% E5,E6,E7 do not affect "pf" and the E4 is known,resulting
% the problem exp(E2^t2)*exp(E3^t3)*pf'' = E1^-1*noap*hst0^-1*pf = pk2p
% which by PADEN-KAHAN-TWO has none, one or two DOUBLE solutions
% t401 & t1in => t301-t205 & t302-t206.
% t402 & t1in => t303-t207 & t304-t208.
%
for i = 9:4:13
    pf2pt = expScrew([Twist(:,4);Theta(i,4)])*[pf; 1];
    pf2p = pf2pt(1:3);
    pk3pt = (expScrew([Twist(:,1);Theta(i,1)]))\noapHst0if;
    pk3p = pk3pt(1:3);
    t2t3 = PadenKahanTwo(Twist(:,2),Twist(:,3),pf2p,pk3p);
    Theta(i,2:3) = t2t3(1,1:2); Theta(i+1,2:3) = t2t3(1,1:2);
    Theta(i+2,2:3) = t2t3(2,1:2); Theta(i+3,2:3) = t2t3(2,1:2);
end
%
% STEP2-4: Calculate Theta5 & Theta6.
% With "pp" on the axis of E7 apply E4^-1*E3^-1**E2^-1*E1^-1*noap*hst0^-1
% to "pp" and also to the POE E5*E6*E7 knowing already t4,t3,t2,t1.
% resulting exactly a Canonic problem PADEN-KAHAN-TWO, because the screws
% E7 does not affect "pp" and the result is the problem
% exp(E5^theta5)*exp(E6^theta6)*pp = pk3p ; with
% pk3p = exp(E4^t4)^-1*exp(E3^t3)^-1*exp(E2^t2)^-1*exp(E1^t1)^-1*noap*hst0^-1*pp 
% which by PADEN-KAHAN-TWO has none, one or two DOUBLE solutions:
% t101-t201-t3in-t401 => t501-t601 & t502-t602
% t102-t202-t3in-t401 => t503-t603 & t504-t604
% t103-t203-t3in-t402 => t505-t605 & t506-t606
% t104-t204-t3in-t402 => t507-t607 & t508-t608
% t1in-t205-t301-t401 => t509-t609 & t510-t610
% t1in-t206-t302-t401 => t511-t611 & t512-t612
% t1in-t207-t303-t402 => t513-t613 & t514-t614
% t1in-t208-t304-t402 => t515-t615 & t516-t616
%
noapHst0ip = noap*(Hst0\[pp; 1]); 
for i = 1:2:15
    pk3pt = (expScrew([Twist(:,1);Theta(i,1)]))\noapHst0ip;
    pk3pt = (expScrew([Twist(:,2);Theta(i,2)]))\pk3pt;
    pk3pt = (expScrew([Twist(:,3);Theta(i,3)]))\pk3pt;
    pk3pt = (expScrew([Twist(:,4);Theta(i,4)]))\pk3pt;
    pk3p = pk3pt(1:3);
    Theta(i:i+1,5:6) = PadenKahanTwo(Twist(:,5),Twist(:,6),pp,pk3p); 
end
%
% STEP2-5: Calculate Theta7.
% With "[1; 0; 0]" not in the axis of E7 apply E6^-1...*E1^-1*noap*gs0^-1 to "po"
% and applying E7 to "po" knowing already Theta6...Theta1,
% resulting exactly a Canonic problem PADEN-KAHAN-ONE, the problem:
% exp(E7^theta7)*po = pk4p ; with
% pk4p = exp(E6^Th6)^-1*...*exp(E1^Th1)^-1*noap*hst0^-1*po 
% which by PADEN-KAHAN-ONE has none or one solution. Then for all
% Th6-Th5-Th4-Th3-Th2-Th1 we get a Th7 = t701...t716:
%
noapHst0io = noap*(Hst0\[1; 0; 0; 1]);
for i = 1:16
    pk4pt = (expScrew([Twist(:,1);Theta(i,1)]))\noapHst0io;
    pk4pt = (expScrew([Twist(:,2);Theta(i,2)]))\pk4pt;
    pk4pt = (expScrew([Twist(:,3);Theta(i,3)]))\pk4pt;
    pk4pt = (expScrew([Twist(:,4);Theta(i,4)]))\pk4pt;
    pk4pt = (expScrew([Twist(:,5);Theta(i,5)]))\pk4pt;
    pk4pt = (expScrew([Twist(:,6);Theta(i,6)]))\pk4pt;
    pk4p = pk4pt(1:3);
    Theta(i,7) = PadenKahanOne(Twist(:,7), [1; 0; 0], pk4p);
end
%
%
if u(7)==17
    ThetaOut = [0 0 0 0 0 0];
else
    ThetaOut = Theta(u(7),:);
end
%
end
%