%% Function F811_ST24R_FK_ABBIRB120" - FORWARD KINEMATICS SIMULATION
% ABB IRB120 - Robot ZERO position with Tool-Down.
% Function solves FORWARD KINEMATICS in REAL TIME to robot joint with two
% ways: MANUAL (commanding the joints) & AUTOMATIC (input trajectory).
%
% The aim is simply to apply motion to the Joints positions
% and check that the FK formulation by POE of Screw Theory coincides
% with the actual motion and pose for the robot tool
% by Dr. Pardos-Gotor ST24R "Screw Theory Toolbox for Robotics" MATLAB.
%
% function ToolPose = F841_ST24R_FK_ABBIRB120(Theta)
%
% INPUTS "Theta" (6x1) are composed by the following vectors.
% Theta (1:6) Joint positions in DEGREES Value
% means to calculate the FORWARD KINEMATICS.
%
% OUTPUTS "ToolPose"(1:6):
% ToolPose(1:3)POSITION X-Y-Z position (m) for the Tool in Spatial frame.
% ToolPose(4:6)ORIENTATION by Euler Rotation XYZ position (rad) for Tool. 
%
% The SIMSCAPE MULTIBODY SYSTEM used to build the Digital-Twin of the robot
% does not support directy Screw Theory and is still based on the DH
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
%% F831_ST24R_FK_ABBIRB120
%
function ToolPose = F831_ST24R_FK_ABBIRB120(Theta) % #codegen
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mechanical characteristics of the Robot (AT REF HOME POSITION):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Joints TWISTS definition and TcP at home.
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