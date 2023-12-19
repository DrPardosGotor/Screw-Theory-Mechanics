%% Function "F833_ST24R_FK_ABBIRB1600" - FORWARD KINEMATICS SIMULATION
% ABB IRB1600 - Robot ZERO position with Tool horizontal on "X".
% Function solves FORWARD KINEMATICS in REAL TIME to robot joint with two
% ways: MANUAL (commanding the joints) & AUTOMATIC (input trajectory).
%
% The aim is simply to apply motion to the Joints positions
% and check that the FK formulation by POE of Screw Theory coincides
% with the actual motion and pose for the robot tool
% by Dr. Pardos-Gotor ST24R "Screw Theory Toolbox for Robotics" MATLAB.
%
% function ToolPose = F843_ST24R_FK_ABBIRB1600(Theta)
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
% J0       0      0      0      0
% J1-J6 = [0.4865 t1     0.15 -90;
%          0      t2-90  0.475  0;
%          0      t3     0    -90;
%          0.6    t4     0     90;
%          0      t5+180 0     90;
%          0.15   t6     0      0];
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
%% F833_ST24R_FK_ABBIRB1600
%
function ToolPose = F833_ST24R_FK_ABBIRB1600(Theta) % #codegen
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mechanical characteristics of the Robot (AT REF HOME POSITION):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Joints TWISTS definition and TcP at home.
Twist = [0   -0.4865   -0.9615         0   -0.9615         0;
         0         0         0    0.9615         0    0.9615;
         0    0.1500    0.1500         0    0.7500         0;
         0         0         0    1.0000         0    1.0000;
         0    1.0000    1.0000         0    1.0000         0;
    1.0000         0         0         0         0         0];
Hst0 = [0 0 1 0.9; 0 1 0 0; -1 0 0 0.9615; 0 0 0 1];
%
% Motion RANGE for the robot joints POSITION rad, (by catalog).
% Thmax = pi/180*[180 150 65 190 115 400];
% Thmin = -pi/180*[180 90 245 190 115 400];
% Maximum SPEED for the robot joints rad/sec, (by catalog).
% Thpmax = pi/180*[180 180 185 385 400 460];
%
%CM1 = [0.15; 0; 0.5]; CM2 = [0.15; -0.15; 0.75]; CM3 = [0.2; 0; 0.96];
%CM4 = [0.55; 0; 0.96]; CM5 = [0.75; 0; 0.96]; CM6 = [0.8; 0; 0.96];
%IT1 = [0.2; 0.2; 0.3]; IT2 = [0.1; 0.1; 0.2]; IT3 = [0.2; 0.1; 0.1];
%IT4 = [0.1; 0.1; 0.1]; IT5 = [0.1; 0.1; 0.1]; IT6 = [0.1; 0.3; 0.3];
%mass = [75 35 25 20 10 5];
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