%% "F854a_ST24R_DK_ABBIRB6620LX" Gantry Robot Differential Kinematics
% ABB IRB 6620LX - Robot with complete dynamics.
%
% Function solves DIFFERENTIAL KINEMATICS in REAL TIME for robot joints.
% It applies the EULER EXPLICIT method for integration.
%
% It solves the DIFFERENTIAL KINEMATICS for desired position & orientation
% velocities of the TCP (noap goal) of the Robot.
% Actually, the aim is to calculate de Joints position without solving
% the Inverse Kinematics, but using the DIFFERENTIAL INVERSE KINEMATICS
% from the GEOMETRIC JACOBIAN (inverse) and the Tool Velocities.
%
% function ThetaOut = F854a_ST24R_DK_ABBIRB6620LX(u)
%
% The inputs "u" (7x1) are composed by the following vectors.
% TargetVAL - Target Value composed by "traXYZ" and "rotXYZ" as
% "traXYZ" (3x1) translations for the TcP (noap - "p" goal) in S frame.
% "rotXYZ" (3x1) rotations for Tool (noap - "noa" (X+Y+Z)) in S frame.
% "Solutions" (1x1) is the value for choosing either track the target (1)
% using screw differential kinematics or sending the robot Home (2).
%
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
% F110_CreateToolTra_Trig([1.5 1.5 0.5 0 pi/4 0], [3 3 2.5 pi/4 pi/2 pi/4], 0.001, 10)
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
%% F854a_ST24R_DK_ABBIRB6620LX
%
function ThetaOut = F854a_ST24R_DK_ABBIRB6620LX(u) %#codegen
%
persistent ThetaVAL;
if isempty(ThetaVAL)
    ThetaVAL = zeros(1,6);
end
%
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
Themax = [3 pi/180*[125 70 300 130 300]];
% Th1max is 33 but limited to 3m
Themin = [0 -pi/180*[125 180 300 130 300]]
% Th1min 1.8 but extended down to 0m by the Spatial frame definition.
%Maximum SPEED for the robot joints rad/sec, (by catalog).
Thepmax = [3.3 pi/180*[90 90 150 120 190]];
Thepmin = -[3.3 pi/180*[90 90 150 120 190]];
%
%CM1 = [1.2; 2.5; 0]; CM2 = [2; 2.5; 0.3]; CM3 = [2.5; 2.5; -0.15];
%CM4 = [2.5; 2; 0]; CM5 = [2.65; 1.6; 0]; CM6 = [2.8; 1.6; 0];
%IT1 = [15; 10; 15]; IT2 = [5; 15; 15]; IT3 = [5; 5; 5];
%IT4 = [15; 5; 15]; IT5 = [3; 5; 5]; IT6 = [0.1; 0.1; 0.1];
%mass = [125 75 50 35 20 10];
%
% stamp is the integration step size (seconds).
stamp = 0.001;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROBOT DIFFERENTIAL KINEMATICS by Screw Theory.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Robot CONFIGURATION Matrix for the possible Joint Magnitudes Solutions.
% Target VALUE is the current cartesian trajectory point.
% it is expressed as a Tool pose [trvX trvY trvZ rotX rotY rotZ] in Euler
% coordinates with translation X-Y-Z and orientation with scheme X-Y-Z:
% from the Forward Kinematics with the current joint position values.
TwMag = [Twist; ThetaVAL];
HstThe = ForwardKinematicsPOE(TwMag)*Hst0;
TargetVAL = [HstThe(1:3,4)' rotm2eul(HstThe(1:3,1:3), 'XYZ')];
%
% Target REFERENCE is the desired next cartesian trajectory point INPUT
% it is expressed as a Tool pose [trvX trvY trvZ rotX rotY rotZ] in Euler
% coordinates with translation X-Y-Z and orientation with scheme X-Y-Z:
TargetREF = u(1:6);
%
% We now solve for the THETAP VELOCITIES values.
% VtS is the classical velocity for Tool Pose in spatial frame (S) rad/s.
% consider the differentiation step size.
VtS = minusposeEul(TargetVAL, TargetREF) / stamp;
%
% GEOMETRIC JACOBIAN JstS and SPATIAL TWIST VELOCITY "VstS" at current pose
JstS = GeoJacobianS(TwMag);
VstS = [VtS(1:3)-axis2skew(VtS(4:6))*TargetVAL(1:3)'; VtS(4:6)];
%
% Using Moore-Penrose generalized inverse for getting Theta velocities.
% Thetap = JstS'*((JstS*JstS')\VstS);
% Next formulation is slower, but works too for Non-Squarre matrices.
% Thetap = (pinv(JstS)*VstS)'; % it is giving worse results.
% Next formulation is faster, but only works for Square matrices 
Thetap = (JstS\VstS)';
%
% The Theta VELOCITIES values are limited by the joints spped limits.
Thetap = jointmag2limits(Thetap, Thepmax, Thepmin);
%
% from Inverse DK we get the incremental joint coordinates
% and then integrating with EULER Explicit Method the Theta VALUE 
% with the new joint positions vector (1x7) OUTPUT.
% consider the integration step size.
ThetaVALThe = ThetaVAL + (Thetap * stamp);
%
% The Theta POSITION values are limited by the joints position limits.
ThetaVALnew = jointmag2limits(ThetaVALThe, Themax, Themin);
%
%
%
if u(7)==2
    ThetaVAL = [0 0 0 0 0 0];
else
    ThetaVAL = ThetaVALnew;
end
ThetaOut = ThetaVAL;
%
end
%