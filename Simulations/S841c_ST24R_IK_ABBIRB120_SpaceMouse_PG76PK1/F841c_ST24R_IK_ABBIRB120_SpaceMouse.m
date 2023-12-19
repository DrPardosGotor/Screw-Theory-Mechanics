%% "F841c_ST24R_IK_ABBIRB120_SpaceMouse" IRB120 ToolDown POSE.
%
% It solves the INVERSE KINEMATICS for any desired position & orientation
% of the TCP (noap goal) of the Robot, this is tool REFerence.
% by Dr. Pardos-Gotor ST24R "Screw Theory Toolbox for Robotics" MATLAB.
%
% function ThetaOut = F841c_ST24R_IK_ABBIRB120_SpaceMouse(u)
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
% J0       0      0    0      0
% J1-J6 = [0.29  t1    0    -90;
%          0     t2-90 0.27   0;
%          0     t3    0.07 -90;
%          0.302 t4    0     90;
%          0     t5-90 0     90;
%          0.16  t6    0      0];
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
%% F841c_ST24R_IK_ABBIRB120_SpaceMouse
%
function ThetaOut = F841c_ST24R_IK_ABBIRB120_SpaceMouse(u) %#codegen
%
%
% Mechanical characteristics of the Robot (AT REF HOME POSITION):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Joints TWISTS definition and TcP at home.
po=[0;0;0]; pk=[0; 0; 0.290]; pr=[0; 0; 0.560];
pf=[0.302; 0; 0.630]; pp=[0.302; 0; 0.470];
Twist = [0   -0.2900   -0.5600         0   -0.6300         0;
         0         0         0    0.6300         0    0.3020;
         0         0         0         0    0.3020         0;
         0         0         0    1.0000         0         0;
         0    1.0000    1.0000         0    1.0000         0;
    1.0000         0         0         0         0   -1.0000];
Hst0 = [-1 0 0 0.302; 0 1 0 0; 0 0 -1 0.47; 0 0 0 1];
%
% Motion RANGE for the robot joints POSITION rad, (by catalog).
% Themax = pi/180*[165 110 70 160 120 400];
% Themin = -pi/180*[165 110 110 160 120 400];
% Maximum SPEED for the robot joints rad/sec, (by catalog).
% Thpmax = pi/180*[250 250 250 320 320 420]
%
%CM1 = [0; 0; 0.29]; CM2 = [0; 0; 0.425]; CM3 = [0; 0; 0.63];
%CM4 = [0.2; 0; 0.63]; CM5 = [0.302; 0; 0.63]; CM6 = [0.302; 0; 0.53];
%IT1 = [0.1; 0.3; 0.2]; IT2 = [0.3; 0.5; 0.1]; IT3 = [0.1; 0.1; 0.1];
%IT4 = [0.1; 0.3; 0.2]; IT5 = [0.1; 0.1; 0.1]; IT6 = [0.1; 0.1; 0.1];
% mass = [7 6 5 4 2 1];
%
%
% Calculate Homogeneous transformation for the GOAL "noap"
traXYZ = u(1:3); rotXYZ = u(4:6);
noap = rotX2tform(rotXYZ(1))*rotY2tform(rotXYZ(2))*rotZ2tform(rotXYZ(3));
noap(1,4)= traXYZ(1); noap(2,4)= traXYZ(2); noap(3,4)= traXYZ(3);
%
n = 6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ABBIRB120 IK solution approach PG7+PG6+PK1 subproblems cosecutively.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the IK solutions Theta using the SCREW THEORY.
Theta_STR2 = zeros(8,n);
tic % start the ticking for calcule the performance of this algorithm.
%
% STEP1: Calculate Theta3.
% With "pf" on the axis of E4, E5, E6 and "pk" on the axis of E1, E2.
% We apply (noap*gs0^-1) to "pf" and take the norm of the diffence of that
% resulting point and "pk". Doing so we can calculate Theta3 applying the
% Canonic problem PADEN-KAHAN-THREE, because the screws E4,E5,E6 do not affect
% "pf" and the E1,E2 do not affect the norm of a vector with an end on "pk"
% resulting the problem ||exp(E3^theta3)*pf-pk||=||noap*gs0^-1*pf-pk||
% which by PARDOS-GOTOR-SEVEN has none, two or four solutions for t123.
noapHst0if = noap*(Hst0\[pf; 1]); pkp = noapHst0if(1:3);
t123 = PardosGotorSeven(Twist(:,1), Twist(:,2), Twist(:,3), pf, pkp);
Theta_STR2(1,1:3) = t123(1,:);
Theta_STR2(2,1:3) = t123(1,:);
Theta_STR2(3,1:3) = t123(2,:);
Theta_STR2(4,1:3) = t123(2,:);
Theta_STR2(5,1:3) = t123(3,:);
Theta_STR2(6,1:3) = t123(3,:);
Theta_STR2(7,1:3) = t123(4,:);
Theta_STR2(8,1:3) = t123(4,:);
%
% STEP2: Calculate Theta4 & Theta5.
% With "pp" on the axis of E6 apply E3^-1*E2^-1*E1^-1*noap*gs0^-1 to "pp"
% and also the POE E4*E5*E6 to "pp" knowing already Theta3-Theta2-Theta1,
% resulting exactly a Canonic problem PADEN-KAHAN-TWO, because the screws
% E6 does not affect "pp" & Th3-Th2-Th1 known (four solutions), the problem
% exp(E4^theta4)*exp(E5^theta5)*pp = pk2p ; with
% pk2p = exp(E3^Th3)^-1*exp(E2^Th2)^-1*exp(E1^Th1)^-1*noap*gs0^-1*pp 
% which by PADEN-KAHAN-TWO has none, one or two DOUBLE solutions:
%
noapHst0ip = noap*(Hst0\[pp; 1]); 
for i = 1:2:7                     % for the 4 values of t3-t2-t1.
    pk2pt = (expScrew([Twist(:,1);Theta_STR2(i,1)]))\noapHst0ip;
    pk2pt = (expScrew([Twist(:,2);Theta_STR2(i,2)]))\pk2pt;
    pk2pt = (expScrew([Twist(:,3);Theta_STR2(i,3)]))\pk2pt;
    pk2p = pk2pt(1:3);
    t4t5 = PadenKahanTwo(Twist(:,4),Twist(:,5),pp,pk2p);
    Theta_STR2(i:i+1,4:5) = t4t5;
end
%
% STEP3: Calculate Theta6.
% With "po" not in the axis of E6 apply E5^-1...*E1^-1*noap*gs0^-1 to "po"
% and applying E6 to "po" knowing already Theta5...Theta1 (8 solutions),
% resulting exactly a Canonic problem PADEN-KAHAN-ONE, the problem:
% exp(E6^theta6)*po = pk3p ; with
% pk3p = exp(E5^Th5)^-1*...*exp(E1^Th1)^-1*noap*gs0^-1*po 
% which by PADEN-KAHAN-ONE has none or one solution. Then for all
% Th5-Th4-Th3-Th2-Th1 known (eight solutions) we get t61...t68:
%
noapHst0io = noap*(Hst0\[po; 1]);
for i = 1:size(Theta_STR2,1)
    pk2pt = (expScrew([Twist(:,1);Theta_STR2(i,1)]))\noapHst0io;
    pk2pt = (expScrew([Twist(:,2);Theta_STR2(i,2)]))\pk2pt;
    pk2pt = (expScrew([Twist(:,3);Theta_STR2(i,3)]))\pk2pt;
    pk2pt = (expScrew([Twist(:,4);Theta_STR2(i,4)]))\pk2pt;
    pk2pt = (expScrew([Twist(:,5);Theta_STR2(i,5)]))\pk2pt;
    pk3p = pk2pt(1:3);
    Theta_STR2(i,6) = PadenKahanOne(Twist(:,6), po, pk3p);
end
%
if u(7)==9
    ThetaOut = [0 0 0 0 0 0];
else
    ThetaOut = Theta_STR2(u(7),:);
end
%
end
%