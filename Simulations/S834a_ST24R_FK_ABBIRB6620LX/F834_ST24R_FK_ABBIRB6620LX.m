%% Function "F834_ST24R_FK_ABBIRB6620LX" - FORWARD KINEMATICS SIMULATION
% ABB IRB6620LX - Robot ZERO position with Tool on "X".
% Function solves FORWARD KINEMATICS in REAL TIME to robot joint with two
% ways: MANUAL (commanding the joints) & AUTOMATIC (input trajectory).
%
% The aim is simply to apply motion to the Joints positions
% and check that the FK formulation by POE of Screw Theory coincides
% with the actual motion and pose for the robot tool
% by Dr. Pardos-Gotor ST24R "Screw Theory Toolbox for Robotics" MATLAB.
%
% function ToolPose = F844_ST24R_FK_ABBIRB6620LX(Theta)
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
%% F834_ST24R_FK_ABBIRB6620LX
%
function ToolPose = F834_ST24R_FK_ABBIRB6620LX(Theta) % #codegen
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mechanical characteristics of the Robot (AT REF HOME POSITION):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Joints TWISTS definition and TcP at home.  
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