%% Function "F861a_ST24R_ID_ABBIRB120_CLA" - Classical Closed-Form Lagrange
% ABB IRB 120 - Robot HOME TOOL down.
% Function solves INVERSE DYNAMICS with CONTROL, which
% takes into account a feedback for POSITION and VELOCITY.
%
% It is necessary to calculate de JOINTS magnitudes (position),
% to implement the INVERSE DYNAMICS algorithm.
% To do so, we use the DIFFERENTIAL KINEMATICS algorithm with the
% GEOMETRIC JACOBIAN (inverse) and the Tool Velocities, to define the
% joint trajectory target (i.e., position, velocity and acceleration).
% Then the classical screw theory closed-form solution for the Lagrange
% formulation is applied.
% by Dr. Pardos-Gotor ST24R "Screw Theory Toolbox for Robotics" MATLAB.
%
% function TorOut = F861a_ST24R_ID_ABBIRB120_CLA(u)
%
% INPUTS: "u" (1x20) = are composed by the following vectors.
% "traXYZ" (1..3) = translations for the TcP (noap - "p" goal) in S frame.
% "rotXYZ" (4..6) = rotations for Tool (noap - "noa" (X+Y+Z)) in S frame.
% "Solutions" (7) =  sending the robot Home (0) or track the target (1)
% "Safety" (8) = Stop de robot (0) or On (1) to track the target 
% "JointVAL"(9..14) = TheVAL Feedback of actual joints POSITIONS rad.
% "JointpVAL"(15..20) = ThepVAL Feedback actual joints VELOCITIES rad/sec
%
% OUTPUTS:
% "TorOut" (t1..t6) are the TORQUES solution for the Robot Joints1..6.
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
% Example of random trajectory generation for AUTOMATIC input testing.
% F110_CreateToolTra_Trig([0.4 -0.2 0.4 0 pi/2 0], [0.5 0.2 0.5 pi/4 3*pi/4 pi/4], 0.001, 10)
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
%% F861a_ST24R_ID_ABBIRB120_CLA
%
function TorOut = F861a_ST24R_ID_ABBIRB120_CLA(u) % #codegen
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mechanical characteristics of the Robot (AT REF HOME POSITION):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
persistent ThetaVAL ThetapVAL ThetappVAL
if isempty(ThetaVAL)
    ThetaVAL = zeros(1,6);
    ThetapVAL = zeros(1,6);
end
%
% The S Spatial system has the "Z" axis up (i.e. -g direction).
PoAcc = [0 0 -9.80665]';
%
% stamp is the integration step size (seconds).
stamp = 0.01;
%
% Robot DOF
DoF = 6;
%
% Joints TWISTS definition and TcP at home.
po=[0;0;0]; pk=[0; 0; 0.29]; pr=[0; 0; 0.56]; ps=[0; 0; 0.63];
pf=[0.302; 0; 0.63]; pu=[0.302; 0; 0.558]; pp=[0.302; 0; 0.47];
AxisX = [1 0 0]'; AxisY = [0 1 0]'; AxisZ = [0 0 1]'; 
Point = [po pk pr ps pf pu];
Joint = ['rot'; 'rot'; 'rot'; 'rot'; 'rot'; 'rot'];
Axis = [AxisZ AxisY AxisY AxisX AxisY -AxisZ];
Twist = zeros(6,DoF);
for i = 1:DoF
    Twist(:,i) = joint2twist(Axis(:,i), Point(:,i), Joint(i,:));
end
Hst0 = trvP2tform(pp)*rotY2tform(pi);
%
% LiMas matrix stands for "Link Mass" and is a 7x7 matrix.
CM = [0 0 0.29; 0 0 0.425; 0 0 0.63; 0.2 0 0.63; 0.302 0 0.63; 0.302 0 0.53]';
IT = [0.1 0.3 0.2; 0.3 0.5 0.1; 0.1 0.1 0.1; 0.1 0.3 0.2; 0.1 0.1 0.1; 0.1 0.1 0.1]';
mass = [7 6 5 4 2 1];
LiMas = [CM; IT; mass];
%
%
% Motion RANGE for the robot joints POSITION rad, (by catalog).
Themax = pi/180*[165 110 70 160 120 400];
Themin = -pi/180*[165 110 110 160 120 400];
% Maximum SPEED for the robot joints rad/sec, (by catalog).
Thepmax = pi/180*[250 250 250 320 320 420];
Thepmin = -pi/180*[250 250 250 320 320 420];
% Maximum Magnitue for the robot joints TORQUE Nm.
Tormax = [195 195 195 195 48 48];
Tormin = -[195 195 195 195 48 48];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading INPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Target REFERENCE is the desired next cartesian trajectory point INPUT
% it is expressed as a Tool pose [trvX trvY trvZ rotX rotY rotZ] in Euler
% coordinates with translation X-Y-Z and orientation with scheme X-Y-Z:
TargetREF = u(1:6);
%
Home = u(7);
Safety = u(8);
%
ThetaFbc = u(9:14);
ThetapFbc = u(15:20);
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INVERSE DIFFERENTIAL KINEMATICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Target VALUE is the current cartesian trajectory point.
% it is expressed as a Tool pose [trvX trvY trvZ rotX rotY rotZ] in Euler
% coordinates with translation X-Y-Z and orientation with scheme X-Y-Z:.
TwMag = [Twist; ThetaVAL];
HstThe = ForwardKinematicsPOE(TwMag)*Hst0;
TargetVAL = [HstThe(1:3,4)' rotm2eul(HstThe(1:3,1:3), 'XYZ')];
%
% We now solve for the THETAP VELOCITIES values.
% VtS is the velocity for the Tool Pose in spatial frame (S) rad/s.
% consider the differentiartion step size.
VtS = minusposeEul(TargetVAL, TargetREF) / stamp;
%
% GEOMETRIC JACOBIAN JstS and SPATIAL TWIST VELOCITY "VstS" at current pose
JstS = GeoJacobianS(TwMag);
VstS = [VtS(1:3)-axis2skew(VtS(4:6))*TargetVAL(1:3)'; VtS(4:6)];
%
% JOINT VELOCITIES:
% Using Moore-Penrose generalized inverse for getting Theta velocities.
% Thetap = JstS'*((JstS*JstS')\VstS);
% Next formulation is slower, but works too for Non-Squarre matrices.
% ThetapNew = (pinv(JstS)*VstS)';
% Next formulation is faster, but only works for Square matrices 
ThetapNew = (JstS\VstS)';
% The Theta VELOCITIES values are limited by the joints spped limits.
ThetapNew = jointmag2limits(ThetapNew, Thepmax, Thepmin);
%
% The Theta ACCELERATION is calculated with reference to Theta VELOCITIES.
ThetappNew = (ThetapNew - ThetapVAL) / 2;
%ThetappNew = (ThetapNew - ThetapVAL) / stamp;
%
% Now we solve for the NEW THETA POSITION values.
% DIFFERENTIAL KINEMATICS is applied to get the incremental joint positions
% and then integrating with EULER Explicit Method the Theta VALUE 
% with the new joint positions vector & consider the integration step size.
ThetaNew = ThetaVAL + (ThetapNew * stamp);
% The Theta POSITION values are limited by the joints position limits.
ThetaNew = jointmag2limits(ThetaNew, Themax, Themin);
%
%
if (Home == 0)
    ThetaNew = zeros(1,DoF);
    ThetapNew = zeros(1,DoF);
    ThetappNew = zeros(1,DoF);
end
if (Safety == 0)
    ThetaNew = ThetaFbc;
    ThetapNew = zeros(1,DoF);
    ThetappNew = zeros(1,Dof);
end
ThetaVAL = ThetaNew;
ThetapVAL = ThetapNew;
ThetappVAL = ThetappNew;
%   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INVERSE DYNAMICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
TwMagFbc = [Twist; ThetaFbc];
%
% complete algebraic solution with ST24R.
% M(t) Inertia matrix by the use of Jsl LINK Jacobian.
%ST24RJsl = MInertiaJsl(TwMagFbc,LiMas);
MtST24RJsl = MInertiaJsl6J_mex(TwMagFbc,LiMas);
% C(t,dt) Coriolis matrix by the use of Aij Adjoint transformation.
%CtdtST24RAij = CCoriolisAij(TwMagFbc,LiMas,ThetapFbc);
CtdtST24RAij = CCoriolisAij6J_mex(TwMagFbc,LiMas,ThetapFbc);
% N(t) NEW Potential Matrix by the use of the Spatial Jacobian.
%NtST24RWre = NPotentialWre(TwMagFbc,LiMas,PoAcc);
NtST24RWre = NPotentialWre6J_mex(TwMagFbc,LiMas,PoAcc);
% Inverse Dynamics solution for the joint TORQUES T Direct (no control).
% FEED-FORWARD contribution to joint TORQUES T.
TorFF = MtST24RJsl*ThetappVAL' + CtdtST24RAij*ThetapFbc' + NtST24RWre;
%
% CONTROL: 
% FEED-BACK contribution to joint TORQUES T.
% POSITION ERROR ThErr and position error velocity.
ThErr = ThetaFbc - ThetaVAL;
KpWeight = [0.4, 0.9, 0.7, 0.4, 0.3, 0.2];
Kp = 950 * KpWeight;
%
% VELOCITY ERROR ThpErr
ThpErr = ThetapFbc - ThetapVAL;
KvWeight = [0.4, 0.9, 0.7, 0.4, 0.3, 0.2];
Kv = 125 * KvWeight;
%
% Torque FeedBack
TorFB = MtST24RJsl*(Kp.*ThErr + Kv.*ThpErr)';
%
% The total TORQUE is the sum of INVERSE DYNAMICS and CONTROL calculations.
TorTot = (TorFF - TorFB)';
%
% The Theta POSITION values are limited by the joints position limits.
TorOut = jointmag2limits(TorTot, Tormax, Tormin);
%
end
%