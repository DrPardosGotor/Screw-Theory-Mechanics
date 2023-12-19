%% Function "F837_ST24R_FK_KUKAIIWA14"
% KUKA IIWA 14 R820 - Robot HOME Straight Up with complete dynamics.
% Function solves FORWARD KINEMATICS in REAL TIME to robot joint.
%
% The aim is simply to apply motion to the Joints positions
% and check that the FK formulation by POE of Screw Theory coincides
% with the actual motion and pose for the robot tool
% by Dr. Pardos-Gotor ST24R "Screw Theory Toolbox for Robotics" MATLAB.
%
% function ToolPose = F847_ST24R_FK_KUKAIIWA14a(Theta)
%
% INPUTS "Theta" (7x1) are composed by the following vectors.
% Theta (1:7) Joint positions in DEGREES Value
% means to calculate the FORWARD KINEMATICS.
%
% OUTPUTS "ToolPose"(1:6):
% ToolPose(1:3)POSITION X-Y-Z position (m) for the Tool in Spatial frame.
% ToolPose(4:6)ORIENTATION by Euler Rotation XYZ position (rad) for Tool.
%
% Mechanical characteristics of the Robot (AT REF HOME POSITION):
% The S Spatial system has the "Z" axis up (i.e., -g direction).
% po = Origen for he STATIONARY system of reference.
% pk = point in the crossing of the DOF Th1(rot) & Th2(rot) & Th3(rot).
% pr = point in the axis of Th4(rot) Th5(rot).
% pf = point in the crossing of the DOF Th5(rot), Th6(rot), Th7(rot).
% pp = TcP Tool Center Point
% hst0 = Tool (TcP) POSE configuration (rot+tra) at reference position. 
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
%% F837_ST24R_FK_KUKAIIWA14
%
function ToolPose = F837_ST24R_FK_KUKAIIWA14(Theta) % #codegen
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mechanical characteristics of the Robot (AT REF HOME POSITION):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
TwMag = [Twist; Theta];
%
% from the Forward Kinematics with the current joint position values.
HstThe = ForwardKinematicsPOE(TwMag)*Hst0;
TargetREF = [HstThe(1:3,4)' rotm2eul(HstThe(1:3,1:3), 'XYZ')];
%
ToolPose = TargetREF;
%
end
%