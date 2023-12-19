%% Function "F867a_ST24R_ID_KUKAIIWA14_CLA"
% KUKA IIWA 14 R820 - Robot HOME Straight Up.
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
% function TdynOut = F867a_ST24R_ID_KUKAIIWA14_CLA(Trajectory)
%
% INPUTS: "u" (1x22) = are composed by the following vectors.
% "traXYZ" (1..3) = translations for the TcP (noap - "p" goal) in S frame.
% "rotXYZ" (4..6) = rotations for Tool (noap - "noa" (X+Y+Z)) in S frame.
% "Solutions" (7) =  sending the robot Home (0) or track the target (1)
% "Safety" (8) = Stop de robot (0) or On (1) to track the target 
% "JointVAL"(9..15) = TheVAL Feedback of actual joints POSITIONS rad.
% "JointpVAL"(16..22) = ThepVAL Feedback actual joints VELOCITIES rad/sec
%
% OUTPUTS:
% "TdynOut" (t1..t7) are the TORQUES solution for the Robot Joints1..7.
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
% J1-J7 = [0.36  t1     0     -90;
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
%% F867a_ST24R_ID_KUKAIIWA14_CLA
%
function TorOut = F867a_ST24R_ID_KUKAIIWA14_CLA(u) % #codegen
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mechanical characteristics of the Robot (AT REF HOME POSITION):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
persistent ThetaVAL ThetapVAL ThetappVAL
if isempty(ThetaVAL)
    ThetaVAL = zeros(1,7);
    ThetapVAL = zeros(1,7);
end
%
% The S Spatial system has the "Z" axis up (i.e. -g direction).
PoAcc = [0 0 -9.80665]';
%
% stamp is the integration step size (seconds).
stamp = 0.01;
%
% Robot DOF
DoF = 7;
%
% Joints TWISTS definition and TcP at home.
po=[0;0;0]; pk=[0;0;0.36]; pr=[0;0;0.78];
pf=[0;0;1.18]; pp=[0;0;1.38];
%AxisX = [1 0 0]';
AxisY = [0 1 0]'; AxisZ = [0 0 1]'; 
Point = [po pk po pr po pf po];
Joint = ['rot'; 'rot'; 'rot'; 'rot'; 'rot'; 'rot'; 'rot'];
Axis = [AxisZ AxisY AxisZ -AxisY AxisZ AxisY AxisZ];
Twist = zeros(6,7);
for i = 1:7
    Twist(:,i) = joint2twist(Axis(:,i), Point(:,i), Joint(i,:));
end
Hst0 = trvP2tform(pp);
%
% LiMas matrix stands for "Link Mass" and is a 7x7 matrix.
CM1 = [0; -0.03; 0.2775]; CM2 = [0; 0.042; 0.419]; CM3 = [0; 0.03; 0.6945];
CM4 = [0; -0.034; 0.847]; CM5 = [0; -0.021; 1];
CM6 = [0; 0.001; 1.18]; CM7 = [0; 0; 1.28];
IT1 = [0.1; 0.09; 0.02]; IT2 = [0.018; 0.05; 0.044];
IT3 = [0.08; 0.075; 0.01];
IT4 = [0.03; 0.029; 0.01]; IT5 = [0.02; 0.018; 0.005];
IT6 = [0.005; 0.0036; 0.0047]; IT7 = [0.001; 0.001; 0.001];
mass = [4 4 3 2.7 1.7 1.8 0.3];
LiCMIT = [CM1 CM2 CM3 CM4 CM5 CM6 CM7; IT1 IT2 IT3 IT4 IT5 IT6 IT7];
LiMas = [LiCMIT; mass];
% Motion RANGE for the robot joints POSITION rad, (by catalog).
Themax = pi/180*[170 120 170 120 170 120 175];
Themin = -pi/180*[170 120 170 120 170 120 175];
% Maximum SPEED for the robot joints rad/sec, (by catalog).
Thepmax = pi/180*[85 85 100 75 130 135 135];
Thepmin = -pi/180*[85 85 100 75 130 135 135];
% Maximum Magnitue for the robot joints TORQUE Nm.
Tormax = [320 320 176 176 110 40 40];
Tormin = -[320 320 176 176 110 40 40];
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
ThetaFbc = u(9:15);
ThetapFbc = u(16:22);
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
%Thetap = (JstS\VstS)';
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
    ThetappNew = zeros(1,DoF);
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
% MtST24RJsl = MInertiaJsl(TwMagFbc,LiMas);
MtST24RJsl = MInertiaJsl7J_mex(TwMagFbc,LiMas);
% C(t,dt) Coriolis matrix by the use of Aij Adjoint transformation.
%CtdtST24RAij = CCoriolisAij(TwMagFbc,LiMas,ThetapFbc);
CtdtST24RAij = CCoriolisAij7J_mex(TwMagFbc,LiMas,ThetapFbc);
%CtdtST24RAij = zeros(7,7);
% N(t) NEW Potential Matrix by the use of the Spatial Jacobian.
%NtST24RWre = NPotentialWre(TwMagFbc,LiMas,PoAcc);
NtST24RWre = NPotentialWre7J_mex(TwMagFbc,LiMas,PoAcc);
% Inverse Dynamics solution for the joint TORQUES T Direct (no control).
% FEED-FORWARD contribution to joint TORQUES T.
TorFF = MtST24RJsl*ThetappVAL' + CtdtST24RAij*ThetapFbc' + NtST24RWre;
%
% CONTROL: 
% FEED-BACK contribution to joint TORQUES T.
% POSITION ERROR ThErr and position error velocity.
ThErr = ThetaFbc - ThetaVAL;
KpWeight = [0.3, 0.8, 0.6, 0.6, 0.3, 0.2, 0.1];
Kp = 950 * KpWeight; 
%Kp = 950 * KpWeight;
%
% VELOCITY ERROR ThpErr
ThpErr = ThetapFbc - ThetapVAL;
KvWeight = [0.3, 0.8, 0.6, 0.6, 0.3, 0.2, 0.1];
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