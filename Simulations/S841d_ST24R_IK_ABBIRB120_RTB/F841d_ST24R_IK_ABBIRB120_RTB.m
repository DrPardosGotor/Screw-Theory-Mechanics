%% "F080_STR_Simulation_421_ABBIRB120_ST24RandRTB_IK" IRB120 ToolDown POSE.
%
% Solves rhe INVERSE KINEMATICS for any desired position & orientation
% of the TCP (noap goal) of the Robot.
% by P. Corke RTB "Robotics Toolbox v10.2 MATLAB.
% by Dr. Pardos-Gotor ST24R "Screw Theory Toolbox for Robotics" MATLAB.
%
% function ThetaOut = F080_STR_Simulation_421_ABBIRB120_ST24RandRTB_IK(u)
%
% The inputs "u" (7x1) are composed by the following vectors.
% "traXYZ" (3x1) desired translations for the TcP (noap - "p" goal).
% "rotXYZ" (3x1) desired rotations for TcP (noap - "noa" goal order (X+Y+Z)
% "Solutions" (1x1) is the value "Theta index" for choosing one out of 8
% possible results (Solutions = 1 to 8) for Robot Joint values
% Solutions = 9 is for sending the robot to the HOME POSITION (ThetaOut=0)
% Solutions = 12 is for sending the robot to the FREEZE POSITION.
% "ThetaOut" (t1..t6)are the magnitudes solution for the Robot Joints1..6.
%
% PUMA type Robot mechanical characteristics at Home (ZERO REFERENCE POSE):
% "TcpInit" (4x4) Initial Pose for TcP (noap tform) at robot zero pose.
% "Dimensi" (3x4) [po pk pr pf] dimensions of the robot structure
% defined by the following points at robot home (zero reference) pose:
% po = Origen for he STATIONARY system of reference.
% pk = point in the crossing of the DOF Th1(rot) & Th2(rot).
% pr = point in the axis of Th3(rot).
% pf = point in the crossing of the DOF Th4(rot), Th5(rot), Th6(rot).
% The considered HOME ZERO REFERENCE POSE is with the ToolUp and the robot
% arm on the X-Y plane:
%                                          Tcp
%                                           ^
%           ->W4[1 0 0] W5[0 0 1]/ W6[0 1 0]|
% W3[0 0 1]/
%
% W2[0 0 1]/
%           ^
%  W1[0 1 0]|
%           ^
%          Y|
%            ->X
%         Z/
%
% Copyright (C) 2003-2018, by Dr. Jose M. Pardos-Gotor.
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
% Revision 1.1  2018/02/11 00:00:01
% General cleanup of code: help comments, see also, copyright
% references, clarification of functions.
%
%% F080_STR_Simulation_421_ABBIRB120_ST24RandRTB_IK
%
function ThetaOut = F080_STR_Simulation_421_ABBIRB120_ST24RandRTB_IK(u) %#codegen
%
persistent Theta ThetaNUM RobotMOV;
if isempty(Theta)
   Theta = zeros(10,6);
end
if isempty(ThetaNUM)
   ThetaNUM = 5;
end
if isempty(RobotMOV)
   RobotMOV = 1;
end
%
% Now we read the DashBoard for checking the solution choosen
% If "Solutions" = u(7) = 1..8 then choosing Theta1..Theta8 Solutions
% If "Solutions" = u(7) = 9 then robot to HOME position and Theta(9) = 0
% If "Solutions" = u(7) = 12 then robot FREEZE to last position Theta(10).
if u(7)<=9
    ThetaNUM = u(7);
    RobotMOV = 1;
end
%
if u(7)==12
    Theta(10,:) = Theta(ThetaNUM,:);
    ThetaNUM = 10;
    RobotMOV = 0;
end
%
% Mechanical characteristics of the IRB120 Robot:
Twist6R = zeros(6,6);
po=[0;0;0]; pk=[0;0.290;0]; pr=[0;0.560;0];
pf=[0.302;0.630;0]; pp=[0.302;0.790;0];
AxisX = [1; 0; 0]; AxisY = [0; 1; 0]; AxisZ = [0; 0; 1]; 
Axis = [AxisY AxisZ AxisZ AxisX AxisZ AxisY];
Point = [pk pk pr pf pf pp];
Joint = ['R'; 'R'; 'R'; 'R'; 'R'; 'R'];
for i = 1:6
    Twist6R(:,i) = Twist(Joint(i,:), Axis(:,i), Point(:,i));
end
%
Hst0 = [rotx(-pi/2)*rotz(pi/2) pp; 0 0 0 1];
%
% Calculate Homogeneous transformation for the TcpGOAL "noap"
traXYZ = u(1:3); rotXYZ = u(4:6);
noap = [rotx(rotXYZ(1))*roty(rotXYZ(2))*rotz(rotXYZ(3)) traXYZ'; 0 0 0 1];
%
% Calculate the IK solutions Theta using the SCREW THEORY techniques
% and basically the PADEN-KAHAN-PARDOS Canonic Subproblems.
%
% STEP1: Calculate Theta3.
% With "pf" on the axis of E4, E5, E6 and "pk" on the axis of E1, E2.
% We apply (noap*gs0^-1) to "pf" and take the norm of the diffence of that
% resulting point and "pk". Doing so we can calculate Theta3 applying the
% Canonic problem PADEN-KAHAN-THREE, because the screws E4,E5,E6 do not affect
% "pf" and the E1,E2 do not affect the norm of a vector with an end on "pk"
% resulting the problem ||exp(E3^theta3)*pf-pk||=||noap*gs0^-1*pf-pk||
% which by PADEN-KAHAN-THREE has none, one or two solutions for t31 t32.
noapHst0if = noap*(Hst0\[pf; 1]); pkp = noapHst0if(1:3);
de = norm(pkp - pk);
t3 = PadenKahanThree(Twist6R(:,3), pf, pk, de);
Theta(1,3) = t3(1); Theta(5,3) = t3(2); % put results into Theta
Theta(3,:) = Theta(1,:); Theta(7,:) = Theta(5,:); % prepare Theta for next
%
% STEP2: Calculate Theta1 & Theta2.
% With "pf" on the axis of E4, E5, E6 we apply (noap*gs0^-1) to "pf" and
% the POE E1..E6 also to "pf" having already known the value for Theta3
% resulting exactly a Canonic problem PADEN-KAHAN-TWO, because the screws
% E4,E5,E6 do not affect "pf" and the E3 is known (two values),resulting
% the problem exp(E1^theta1)*exp(E2^theta2)*pf' = noap*gs0^-1*pf
% which by PADEN-KAHAN-TWO has none, one or two DOUBLE solutions
% t11-t21 & t12-t22 for each value of t3, but we have two, then consider
% for t31 we get t11-t21 & t12-t22 & for t32 we get t13-t23 & t14-t24.
t1t2 = zeros(4,2);
for i = 1:2                   % for the TWO values of t3.
    j = i+fix(i/2);
    pfpt = trexp(Twist6R(:,3)',t3(i))*[pf; 1];
    pfp = pfpt(1:3);
    t1t2(j:j+1,1:2) = PadenKahanTwo(Twist6R(:,1),Twist6R(:,2),pfp,pkp);
end
for i = 1:2:7
    j = i-fix(i/2);
    Theta(i,1:2) = t1t2(j,1:2); % put t1t2 values into Theta
    Theta(i+1,:) = Theta(i,:);  % prepare Theta for next step
end
%
% STEP3: Calculate Theta4 & Theta5.
% With "pp" on the axis of E6 apply E3^-1*E2^-1*E1^-1*noap*gs0^-1 to "pp"
% and also the POE E4*E5*E6 to "pp" knowing already Theta3-Theta2-Theta1,
% resulting exactly a Canonic problem PADEN-KAHAN-TWO, because the screws
% E6 does not affect "pp" & Th3-Th2-Th1 known (four solutions), the problem
% exp(E4^theta4)*exp(E5^theta5)*pp = pk2p ; with
% pk2p = exp(E3^Th3)^-1*exp(E2^Th2)^-1*exp(E1^Th1)^-1*noap*gs0^-1*pp 
% which by PADEN-KAHAN-TWO has none, one or two DOUBLE solutions:
% t31,t21,t11 to t41-t51 & t42-t52 ; t31,t22,t12 to t43-t53 & t44-t54
% t32,t23,t13 to t45-t55 & t46-t56 ; t32,t24,t14 to t47-t57 & t48-t58
%
t4t5 = zeros(8,2);
noapHst0ip = noap*(Hst0\[pp; 1]); 
for i = 1:2:7                     % for the 4 values of t3-t2-t1.
    pk2pt = (trexp(Twist6R(:,1)',Theta(i,1)))\noapHst0ip;
    pk2pt = (trexp(Twist6R(:,2)',Theta(i,2)))\pk2pt;
    pk2pt = (trexp(Twist6R(:,3)',Theta(i,3)))\pk2pt;
    pk2p = pk2pt(1:3);
    t4t5(i:i+1,1:2) = PadenKahanTwo(Twist6R(:,4),Twist6R(:,5),pp,pk2p); 
end
    Theta(1:8,4:5) = t4t5;          % put t4t5 values into Theta
% STEP4: Calculate Theta6.
% With "po" not in the axis of E6 apply E5^-1...*E1^-1*noap*gs0^-1 to "po"
% and applying E6 to "po" knowing already Theta5...Theta1 (8 solutions),
% resulting exactly a Canonic problem PADEN-KAHAN-ONE, the problem:
% exp(E6^theta6)*po = pk3p ; with
% pk3p = exp(E5^Th5)^-1*...*exp(E1^Th1)^-1*noap*gs0^-1*po 
% which by PADEN-KAHAN-ONE has none or one solution. Then for all
% Th5-Th4-Th3-Th2-Th1 known (eight solutions) we get t61...t68:
noapHst0io = noap*(Hst0\[po; 1]);
for i = 1:8
    pk2pt = (trexp(Twist6R(:,1)',Theta(i,1)))\noapHst0io;
    pk2pt = (trexp(Twist6R(:,2)',Theta(i,2)))\pk2pt;
    pk2pt = (trexp(Twist6R(:,3)',Theta(i,3)))\pk2pt;
    pk2pt = (trexp(Twist6R(:,4)',Theta(i,4)))\pk2pt;
    pk2pt = (trexp(Twist6R(:,5)',Theta(i,5)))\pk2pt;
    pk3p = pk2pt(1:3);
    Theta(i,6) = PadenKahanOne(Twist6R(:,6), po, pk3p);
end
%
ThetaOut = Theta(ThetaNUM,:);
%
end
%