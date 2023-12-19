%% Function "F865a_ST24R_ID_ABBIR910SC_CLA"
% ABB IRB 910SC - Robot HOME Straight ahead.
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
%
% function TorOut = F865a_ST24R_ID_ABBIR910SC_CLA(u)
%
% INPUTS: "u" (1x20) = are composed by the following vectors.
% "traXYZ" (1..3) = translations for the TcP (noap - "p" goal) in S frame.
% "rotXYZ" (4..6) = rotations for Tool (noap - "noa" (X+Y+Z)) in S frame.
% "Solutions" (7) =  sending the robot Home (0) or track the target (1)
% "Safety" (8) = Stop de robot (0) or On (1) to track the target 
% "JointVAL"(9..12) = TheVAL Feedback of actual joints POSITIONS.
% "JointpVAL"(13..16) = ThepVAL Feedback actual joints VELOCITIES.
%
% OUTPUTS:
% "TorOut" (t1..t4) TORQUES (t1,t2,t4) & FORCE (t3)solution for Joints1..4.
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
%% F865a_ST24R_ID_ABBIRB910SC_CLA
%
function TorOut = F865a_ST24R_ID_ABBIRB910SC_CLA(u) % #codegen
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mechanical characteristics of the Robot (AT REF HOME POSITION):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
persistent ThetaVAL ThetapVAL ThetappVAL
if isempty(ThetaVAL)
    ThetaVAL = zeros(1,4);
    ThetapVAL = zeros(1,4);
end
%
% The S Spatial system has the "Z" axis up (i.e. -g direction).
PoAcc = [0 -9.80665 0]';
%
% stamp is the integration step size (seconds).
stamp = 0.01;
%
% Robot DOF
n = 4;
%
% Joints TWISTS definition and TcP at home.
po=[0;0;0]; pr=[0.4;0;0]; pf=[0.65;0;0]; pp=[0.65;0.125;0]; 
%AxisX = [1 0 0]'; AxisZ = [0 0 1]'; 
AxisY = [0 1 0]';
Point = [po pr pf pp];
Joint = ['rot'; 'rot'; 'tra'; 'rot'];
Axis = [AxisY AxisY AxisY -AxisY];
Twist = zeros(6,n);
for i = 1:4
    Twist(:,i) = joint2twist(Axis(:,i), Point(:,i), Joint(i,:));
end
Hst0 = trvP2tform(pp)*rotX2tform(pi/2)*rotZ2tform(pi);
%
% LiMas matrix stands for "Link Mass" and is a 7x7 matrix.
CM1 = [0.2; 0.2; 0]; CM2 = [0.5; 0.258; 0];
CM3 = [0.65; 0.258; 0]; CM4 = [0.65; 0.208; 0];
IT1 = [0.1; 0.3; 0.2]; IT2 = [0.1; 0.5; 0.3];
IT3 = [0.1; 0.1; 0.1]; IT4 = [0.1; 0.1; 0.1];
mass = [7 5 1 0.5];
LiMas = [CM1 CM2 CM3 CM4; IT1 IT2 IT3 IT4; mass];
%
%Motion RANGE for the robot joints POSITION rad, (by catalog).
Themax = [pi/180*140 pi/180*150 0 pi/180*400];
Themin = [-pi/180*140 -pi/180*150 -0.125 -pi/180*400];
% Maximum SPEED for the robot joints m/s and rad/sec, (by catalog).
Thepmax = [7.58 7.58 1.02 pi/180*2400];
Thepmin = -[7.58 7.58 1.02 pi/180*2400];
% Maximum Magnitue for the robot joints TORQUE Nm.
Tormax = [250 250 25 25];
Tormin = -[250 250 25 25];
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
ThetaFbc = u(9:12);
ThetapFbc = u(13:16);
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
ThetapNew = (pinv(JstS)*VstS)';
% Next formulation is faster, but only works for Square matrices 
%ThetapNew = (JstS\VstS)';
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
    ThetaNew = zeros(1,n);
    ThetapNew = zeros(1,n);
    ThetappNew = zeros(1,n);
end
if (Safety == 0)
    ThetaNew = ThetaFbc;
    ThetapNew = zeros(1,n);
    ThetappNew = zeros(1,n);
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
%MtST24RJsl = MInertiaJsl(TwMagFbc,LiMas);
MtST24RJsl = MInertiaJsl4J_mex(TwMagFbc,LiMas);
% C(t,dt) Coriolis matrix by the use of Aij Adjoint transformation.
%CtdtST24RAij = CCoriolisAij(TwMagFbc,LiMas,ThetapFbc);
CtdtST24RAij = CCoriolisAij4J_mex(TwMagFbc,LiMas,ThetapFbc);
% N(t) NEW Potential Matrix by the use of the Spatial Jacobian.
%NtST24RWre = NPotentialWre(TwMagFbc,LiMas,PoAcc);
NtST24RWre = NPotentialWre4J_mex(TwMagFbc,LiMas,PoAcc);
% Inverse Dynamics solution for the joint TORQUES T Direct (no control).
% FEED-FORWARD contribution to joint TORQUES T.
TorFF = MtST24RJsl*ThetappVAL' + CtdtST24RAij*ThetapFbc' + NtST24RWre;
%
% CONTROL: 
% FEED-BACK contribution to joint TORQUES T.
% POSITION ERROR ThErr and position error velocity.
ThErr = ThetaFbc - ThetaVAL;
KpWeight = [0.4, 0.4, 0.3, 0.1];
Kp = 950 * KpWeight;
%Kp = 950 * KpWeight;
%
% VELOCITY ERROR ThpErr
ThpErr = ThetapFbc - ThetapVAL;
KvWeight = [0.4, 0.4, 0.3, 0.1];
Kv = 125 * KvWeight;
%Kv = 125 * KvWeight;
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