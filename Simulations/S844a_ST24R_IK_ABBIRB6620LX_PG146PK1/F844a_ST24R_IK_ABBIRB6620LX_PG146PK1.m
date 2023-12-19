%% Function "F844a_ST24R_IK_ABBIRB6620LX_PG146PK1" IK SIMULATION
% ABB IRB6620LX - Robot ZERO position with horizontal.
%
% It solves the INVERSE KINEMATICS for any desired position & orientation
% of the TCP (noap goal) of the Robot, this is tool REFerence.
% by Dr. Pardos-Gotor ST24R "Screw Theory Toolbox for Robotics" MATLAB.
%
% function ThetaOut = F844a_ST24R_IK_ABBIRB6620LX_PG146PK1(u)
%
% The inputs "u" (7x1) are composed by the following vectors.
% "traXYZ" (3x1) desired translations for the TcP (noap - "p" goal).
% "rotXYZ" (3x1) desired rotations for TcP (noap - "noa" goal order (X+Y+Z)
% "Solutions" (1x1) is the value "Theta index" for choosing one out of 4
% possible results (Solutions = 1 to 4) for Robot Joint values
% Solutions = 5 is for sending the robot to the HOME POSITION (ThetaOut=0)
% "ThetaOut" (t1..t6)are the magnitudes solution for the Robot Joints1..6.
%
% The SIMSCAPE MULTIBODY SYSTEM used to build the Digital-Twin of the robot
% does not support direxcty Screw Theory and is still based on the DH
% formalism. For this reason, the revolute joints only can rotate on "Z"
% axis. Therefore, to construct the manipulator body we must use the
% Denavit-Hartenberg parameters (in this case the Classical DH - "DHC").
% Afterwards, Screw Theory POE is used to perform forward kinematics.
% For this robot the DHC parameters are:
% Joint-dtra
% J0       0     90     2.5     0
% J1-J6 = [0     t1-90  1.468 180;
%          0     t2     0.975   0;
%          0     t3     0.2   -90;
%          0.887 t4+180 0     -90;
%          0     t5-90  0      90;
%          0.357 t6     0       0];
% 
% Example of random trajectory generation for AUTOMATIC input testing.
% F110_CreateToolTra_Sin([2.5 2 1 0 pi/4 0], [3 3 3 pi/4 pi/2 pi/4], 0.001, 10)
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
%% F844a_ST24R_IK_ABBIRB6620LX_PG146PK1
%
function ThetaOut = F844a_ST24R_IK_ABBIRB6620LX_PG146PK1(u) % #codegen
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mechanical characteristics of the Robot (AT REF HOME POSITION):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Joints TWISTS definition and TcP at home.
po=[0;0;0]; pf=[2.643;1.613;0]; pp=[3.000;1.613;0];
%
Twist = [0   -2.5000   -2.5000         0   -1.6130         0;
         0    1.4680    2.4430         0    2.6430         0;
    1.0000         0         0   -2.6430         0   -1.6130;
         0         0         0         0         0    1.0000;
         0         0         0   -1.0000         0         0;
         0   -1.0000   -1.0000         0   -1.0000         0];
Hst0 = [0 0 1 3; -1 0 0 1.613; 0 -1 0 0; 0 0 0 1];
%
%Motion RANGE for the robot joints POSITION rad, (by catalog).
%Thmax = [33 pi/180*[125 70 300 130 300]]; Th1max limited to 2.5m
%Thmin = [1.8 -pi/180*[125 180 300 130 300]]; Th1min limited to 0m.
%Maximum SPEED for the robot joints rad/sec, (by catalog).
%Thpmax = [3.3 pi/180*[90 90 150 120 190]];
%
%CM1 = [1.2; 2.5; 0]; CM2 = [2; 2.5; 0.3]; CM3 = [2.5; 2.5; -0.15];
%CM4 = [2.5; 2; 0]; CM5 = [2.65; 1.6; 0]; CM6 = [2.8; 1.6; 0];
%IT1 = [15; 10; 15]; IT2 = [5; 15; 15]; IT3 = [5; 5; 5];
%IT4 = [15; 5; 15]; IT5 = [3; 5; 5]; IT6 = [0.1; 0.1; 0.1];
%mass = [125 75 50 35 20 10];
%
% Calculate Homogeneous transformation for the GOAL "noap"
traXYZ = u(1:3); rotXYZ = u(4:6);
noap = rotX2tform(rotXYZ(1))*rotY2tform(rotXYZ(2))*rotZ2tform(rotXYZ(3));
noap(1,4)= traXYZ(1); noap(2,4)= traXYZ(2); noap(3,4)= traXYZ(3);
%
n = 6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IK solution approach PG1+PG4+PG6+PK1 subproblems cosecutively.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the IK solutions Theta using the SCREW THEORY.
Theta_STR4 = zeros(4,n);
%
% With "pf" on the axis of E4, E5, E6. We apply (noap*hs0^-1) to "pf"
% Doing so we get Theta1 applying the Canonic problem PARDOS-ONE,
% because the screws E4,E5,E6 do not affect "pf" for being on their axes
% and the E2,E3 do not change the plane where "pf" moves (perpendicular to
% the axis of those screws, and so do not affect the calculation for Theta1
% resulting the problem "exp(E1^theta1)*pf = noap*hs0^-1*pf" by PARDOS-ONE
% which has one solution for t1.
noapHst0if = noap*(Hst0\[pf; 1]); pk1 = noapHst0if(1:3);
t1 = PardosGotorOne(Twist(:,1), pf, pk1);
% prepare Theta for next calculation
Theta_STR4(1:4,1) = t1;
%
% STEP2: Calculate Theta2 & Theta3.
% With "pf" on the axis of E4, E5, E6 we apply (noap*hs0^-1) to "pf" and
% the POE E1..E6 also to "pf" having already known the value for Theta1
% resulting exactly a Canonic problem PARDOS-FOUR, because the screws
% E4,E5,E6 do not affect "pf" and the E1 is known,resulting the problem
% exp(E2^theta2)*exp(E3^theta3)*pf = exp(E1^Th1)^-1*noap*gs0^-1*pf = pk1p
% which by PARDOS-FOUR has none, one or two DOUBLE solutions.
% t21-t31 & t22-t32 for each value of t11
%
E1inoapHst0if = (expScrew([Twist(:,1);t1]))\noapHst0if;
pk4 = E1inoapHst0if(1:3);
t2t3 = PardosGotorFour(Twist(:,2),Twist(:,3),pf,pk4);
Theta_STR4(1,2:3) = t2t3(1,:);
Theta_STR4(2,2:3) = t2t3(1,:);
Theta_STR4(3,2:3) = t2t3(2,:);
Theta_STR4(4,2:3) = t2t3(2,:);
%
% STEP3: Calculate Theta4 & Theta5.
% With "pp" on the axis of E6 apply E3^-1*E2^-1*E1^-1*noap*gs0^-1 to "pp"
% and also the POE E4*E5*E6 to "pp" knowing already Theta3-Theta2-Theta1,
% resulting exactly a Canonic problem PADEN-KAHAN-TWO, because the screws
% E6 does not affect "pp" & Th3-Th2-Th1 known (four solutions), the problem
% exp(E4^theta4)*exp(E5^theta5)*pp = pk2p ; with
% pk2p = exp(E3^Th3)^-1*exp(E2^Th2)^-1*exp(E1^Th1)^-1*noap*gs0^-1*pp 
% which by PARDOS-GOTOR-SIX has none, one or two DOUBLE solutions:
%
noapHst0ip = noap*(Hst0\[pp; 1]); 
for i = 1:2:3                     % for the 2 values of t3-t2-t1.
    pk2pt = (expScrew([Twist(:,1);Theta_STR4(i,1)]))\noapHst0ip;
    pk2pt = (expScrew([Twist(:,2);Theta_STR4(i,2)]))\pk2pt;
    pk2pt = (expScrew([Twist(:,3);Theta_STR4(i,3)]))\pk2pt;
    pk2p = pk2pt(1:3);
    t4t5 = PardosGotorSix(Twist(:,4),Twist(:,5),pp,pk2p);
    Theta_STR4(i:i+1,4:5) = t4t5; 
end
%
% STEP4: Calculate Theta6.
% With "po" not in the axis of E6 apply E5^-1...*E1^-1*noap*gs0^-1 to "po"
% and applying E6 to "po" knowing already Theta5...Theta1 (8 solutions),
% resulting exactly a Canonic problem PADEN-KAHAN-ONE, the problem:
% exp(E6^theta6)*po = pk3p ; with
% pk3p = exp(E5^Th5)^-1*...*exp(E1^Th1)^-1*noap*gs0^-1*po 
% which by PADEN-KAHAN-ONE has none or one solution. Then for all
% Th5-Th4-Th3-Th2-Th1 known (eight solutions) we get t61...t68:
%
noapHst0io = noap*(Hst0\[po; 1]);
for i = 1:size(Theta_STR4,1)
    pk2pt = (expScrew([Twist(:,1);Theta_STR4(i,1)]))\noapHst0io;
    pk2pt = (expScrew([Twist(:,2);Theta_STR4(i,2)]))\pk2pt;
    pk2pt = (expScrew([Twist(:,3);Theta_STR4(i,3)]))\pk2pt;
    pk2pt = (expScrew([Twist(:,4);Theta_STR4(i,4)]))\pk2pt;
    pk2pt = (expScrew([Twist(:,5);Theta_STR4(i,5)]))\pk2pt;
    pk3p = pk2pt(1:3);
    Theta_STR4(i,6) = PadenKahanOne(Twist(:,6), po, pk3p);
end
%
if u(7)==5
    ThetaOut = [0 0 0 0 0 0];
else
    ThetaOut = Theta_STR4(u(7),:);
end
%
end
%