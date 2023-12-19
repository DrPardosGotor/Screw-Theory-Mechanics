%% Function "F845a_ST24R_IK_ABBIR910SC_PG14PK1" IK SIMULATION
% ABB IRB910SC - Robot ZERO position straight on "X".
%
% It solves the INVERSE KINEMATICS for any desired position & orientation
% of the TCP (noap goal) of the Robot, this is tool REFerence.
% by Dr. Pardos-Gotor ST24R "Screw Theory Toolbox for Robotics" MATLAB.
%
% function ThetaOut = F845a_ST24R_IK_ABBIR910SC_PG14PK1(u)
%
% The inputs "u" (7x1) are composed by the following vectors.
% "traXYZ" (3x1) desired translations for the TcP (noap - "p" goal).
% "rotXYZ" (3x1) desired rotations for TcP (noap - "noa" goal order (X+Y+Z)
% "Solutions" (1x1) is the value "Theta index" for choosing one out of 2
% possible results (Solutions = 1 to 2) for Robot Joint values
% Solutions = 3 is for sending the robot to the HOME POSITION (ThetaOut=0)
% "ThetaOut" (t1..t4)are the magnitudes solution for the Robot Joints1..4.
%
% The SIMSCAPE MULTIBODY SYSTEM used to build the Digital-Twin of the robot
% does not support direxcty Screw Theory and is still based on the DH
% formalism. For this reason, the revolute joints only can rotate on "Z"
% axis. Therefore, to construct the manipulator body we must use the
% Denavit-Hartenberg parameters (in this case the Classical DH - "DHC").
% Afterwards, Screw Theory POE is used to perform forward kinematics.
% For this robot the DHC parameters are:
% Joint-dtra
% J0       0      0    0    -90
% J1-J6 = [0.258 t1    0.4    0;
%          0     t2    0.25   0;
%          t3   180    0    180;
%          0.133 t4    0      0];
% 
% Example of random trajectory generation for AUTOMATIC input testing.
% F110_CreateToolTra_Trig([0.3 0 -0.3 pi/2 0 0], [0.6 0.12 0.4 pi/2 0 pi/4], 0.001, 10)
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
%% F845a_ST24R_IK_ABBIR910SC_PG14PK1
%
function ThetaOut = F845a_ST24R_IK_ABBIR910SC_PG14PK1(u) % #codegen
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mechanical characteristics of the Robot (AT REF HOME POSITION):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Joints TWISTS definition and TcP at home.
po=[0;0;0]; pf=[0.65;0;0];
%
Twist = [0         0         0         0;
         0         0    1.0000         0;
         0    0.4000         0   -0.6500;
         0         0         0         0;
    1.0000    1.0000         0   -1.0000;
         0         0         0         0];
Hst0 = [-1 0 0 0.65; 0 0 -1 0.125; 0 -1 0 0; 0 0 0 1];
%
% Motion RANGE for the robot joints POSITION rad, (by catalog).
% Thmax = [pi/180*140 pi/180*150 0,18 pi/180*400];
% Thmin = [-pi/180*140 -pi/180*150 0 -pi/180*400];
% Maximum SPEED for the robot joints m/s and rad/sec, (by catalog).
% Thpmax = [7.58 7.58 1.02 pi/180*2400];
%
%CM1 = [0.2; 0.2; 0]; CM2 = [0.5; 0.258; 0];
%CM3 = [0.65; 0.258; 0]; CM4 = [0.65; 0.208; 0];
%IT1 = [0.1; 0.3; 0.2]; IT2 = [0.1; 0.5; 0.3];
%IT3 = [0.1; 0.1; 0.1]; IT4 = [0.1; 0.1; 0.1];
%mass = [7 5 1 0.5];
%
% Calculate Homogeneous transformation for the GOAL "noap"
traXYZ = u(1:3); rotXYZ = u(4:6);
noap = rotX2tform(rotXYZ(1))*rotY2tform(rotXYZ(2))*rotZ2tform(rotXYZ(3));
noap(1,4)= traXYZ(1); noap(2,4)= traXYZ(2); noap(3,4)= traXYZ(3);
%
% number of joints.
n = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IK solution approach PG1+PG4+PK1 subproblems cosecutively.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the IK solutions Theta using the SCREW THEORY.
Theta = zeros(2,n);
%
% STEP2: Calculate Theta3.
% With "pf" on the axis of E4. We apply (noap*hs0^-1) to "pf"
% doing so we get Theta3 applying the Canonic problem PARDOS-ONE,
% because the screws E4 do not affect "pf" for being on its axis
% and the E1,E2 do not change the plane where "pf" moves (perpendicular to
% the axis of those screws, and so do not affect the calculation for Theta3
% resulting the problem "exp(E3^theta3)*pf = noap*hst0^-1*pf" by PARDOS-ONE
% which has one solution for t31.
pkp = noap*(Hst0\[pf; 1]);
Theta(1,3) = PardosGotorOne(Twist(:,3), pf, pkp(1:3));
% prepare Theta for next calculation
Theta(2,3) = Theta(1,3);
%
% STEP2: Calculate Theta1 & Theta2.
% With "pf" on the axis of E4, we apply (noap*hs0^-1) to "pf" and
% the POE E1..E4 also to "pf" having already known the value for Theta3
% resulting exactly a Canonic problem PARDOS-FOUR, because the screw
% E4 do not affect "pf" and the E3 is known,resulting the problem
% exp(E1^theta1)*exp(E2^theta2)*exp(E3^theta3)*pf = 
% exp(E1^theta1)*exp(E2^theta2)*pfp = noap*hst0^-1*pf = pkp
% which by PARDOS-FOUR has none, one or two DOUBLE solutions.
% t11-t21 & t12-t22 for each value of t31
%
pfp = expScrew([Twist(:,3);Theta(1,3)])*[pf; 1];
Theta(1:2,1:2) = PardosGotorFour(Twist(:,1),Twist(:,2),pfp(1:3),pkp(1:3));
%
% STEP2: Calculate Theta4.
% With "po" not in the axis of E4 apply E3^-1...*E1^-1*noap*hst0^-1 to "po"
% and applying E4 to "po" knowing already Theta3...Theta1 solutions,
% resulting exactly a Canonic problem PADEN-KAHAN-ONE, the problem:
% exp(E4^theta4)*po = pk3p ; with
% pk3p = exp(E3^Th3)^-1*...*exp(E1^Th1)^-1*noap*hst0^-1*po 
% which by PADEN-KAHAN-ONE has none or one solution. Then for all
% Th3-Th2-Th1 known (two solutions) we get t4:
noapHst0io = noap*(Hst0\[po; 1]);
for i = 1:2
    pk3p = (expScrew([Twist(:,1);Theta(i,1)]))\noapHst0io;
    pk3p = (expScrew([Twist(:,2);Theta(i,2)]))\pk3p;
    pk3p = (expScrew([Twist(:,3);Theta(i,3)]))\pk3p;
    Theta(i,4) = PadenKahanOne(Twist(:,4), po, pk3p(1:3));
end
%
if u(7)==3
    ThetaOut = [0 0 0 0];
else
    ThetaOut = Theta(u(7),:);
end
%
end
%