%% Function "F866b_ST24R_ID_UR16e_SVA" - Spatial Vector Algebra.
% UNIVERSAL ROBOT UR16e - Robot HOME Extended.
% Function solves INVERSE DYNAMICS with CONTROL, which
% takes into account a feedback for POSITION and VELOCITY.
%
% It is necessary to calculate de JOINTS magnitudes (position),
% to implement the INVERSE DYNAMICS algorithm.
% To do so, we use the DIFFERENTIAL KINEMATICS algorithm with the
% GEOMETRIC JACOBIAN (inverse) and the Tool Velocities, to define the
% joint trajectory target (i.e., position, velocity and acceleration).
% Then the SVA (Spatial Vector Algebra) is used implementing the
% RNEA Recursive Newton-Euler Algorithm with use of POE for kinematics.
% by Dr. Pardos-Gotor ST24R "Screw Theory Toolbox for Robotics" MATLAB.
%
% function TorOut = F866b_ST24R_ID_UR16e_SVA(u)
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
% J0       0      0     0       0
% J1-J6 = [0.181 t1     0     -90;
%          0     t2     0.478   0;
%          0     t3     0.36    0;
%          0.174 t4     0     -90;
%          0.12  t5     0      90;
%          0.19  t6+180 0       0];
%
% Example of random trajectory generation for AUTOMATIC input testing.
% F110_CreateToolTra_Trig([0.3 -0.4 0.3 0 0 pi/2], [0.55 0.4 0.7 pi/4 pi/4 5*pi/4], 0.001, 10)
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
%% F866b_ST24R_ID_UR16e_SVA
%
function TorOut = F866b_ST24R_ID_UR16e_SVA(u) % #codegen
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
%po=[0;0;0];
pg=[0; 0; 0.1]; pk=[0; 0; 0.181]; pr=[0.478; 0; 0.181];
ps=[0.838; 0; 0.181]; pf=[0.838; 0.174; 0.181];
pu=[0.838; 0.174; 0.061]; pp=[0.838; 0.364; 0.061];
%AxisX = [1 0 0]';
AxisY = [0 1 0]'; AxisZ = [0 0 1]'; 
Point = [pg pk pr ps pf pu];
Joint = ['rot'; 'rot'; 'rot'; 'rot'; 'rot'; 'rot'];
Axis = [AxisZ AxisY AxisY AxisY -AxisZ AxisY];
Twist = zeros(6,DoF);
for i = 1:DoF
    Twist(:,i) = joint2twist(Axis(:,i), Point(:,i), Joint(i,:));
end
Hst0 = trvP2tform(pp)*rotX2tform(-pi/2)*rotZ2tform(pi);
%
% DYNAMIC Parameters of the Robot at REF HOME POSITION - Only aproximation
CM13 = [0 0.01 0.15; 0.2 0.174 0.181; 0.628 0.05 0.181]';
CM46 = [0.838 0.174 0.144; 0.838 0.194 0.061; 0.838 0.274 0.061]';
CM = [CM13 CM46];
IT = [0.2 0.2 0.3; 0.1 0.1 0.2; 0.2 0.1 0.1; 0.1 0.1 0.1; 0.1 0.1 0.1; 0.1 0.1 0.1]';
mass = [7.369 10.45 4.321 2.18 2.033 0.907];
%LiMas = [CM; IT; mass];
%
% Motion RANGE for the robot joints POSITION rad, (by catalog).
Themax = pi/180*[360 360 360 360 360 360];
Themin = -pi/180*[360 360 360 360 360 360];
% Maximum SPEED for the robot joints rad/sec, (by catalog).
Thepmax = pi/180*[120 120 180 180 180 180];
Thepmin = -pi/180*[120 120 180 180 180 180];
% Maximum Magnitue for the robot joints TORQUE Nm.
%Tormax = [195 195 195 195 48 48];
%Tormin = -[195 195 195 195 48 48];
Tormax = [950 2100 950 200 200 200];
Tormin = -[950 2100 950 200 200 200];
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
%TwMagFbc = [Twist; ThetaVAL];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ID with the Spatial Vector Algebra.
% by Featherstone Spatial Toolbox 2015.
% RNEA-POE - Recursive Newton-Euler Algorithm with Screw Theory POE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Xts is the transformation from the Tool to the Spatial Systems
% using the Spatial Vector Featherstone nomenclature
% the magnitudes are always measured in the predecesor link
% In fact, Featherstone nomenclature moves the succesor link in order to
% make it coincide wigth the predecessor. In such a way, the X(6x6)
% Featherstone matrix, multiplies magnitudes from the predecesor to give
% as a result the transformed magnitude in the successor frame. 
%
% Matrices of Inertia (I1...I3) defined in LINK Frame with S orientation.
% the postion on the CM is defined in relation to each LINK Frame.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This block of variables is useful for the preparation of the
% first recursive OUTWARDS PASS
%
% Gravitational field is modelled as a fictitious acceleration of the base.
% To do this, we replace the initial value with a0 = - ag = - PoAcc
% With this trick the values for ai and fi (net force body) are no longer
% the true values, as they are offset by gravity acceleration and force.
% Nonetheless, the result of the ID is correct, as we are interested only
% in the Joint Torques.
%
% Gravity definition: PoAcc = [0 0 -9.81]';
ai = [0;0;0; -PoAcc];
%
% Motion Subspace for the Joints.
S = [Axis; zeros(3,DoF)];
%
% Initial values for the recursive algorithm.
PoE = eye(4); % Product of Exponentials.
Hs0 = eye(4); % Homogeneous transformation for the Link Frame.
Pre = [0; 0; 0]; % Origin of Base Frame.
% Xst stores the Spatial Vector transformation to the base (zero), plus the
% Links Frames and besides the X for the Tool.
Xst = zeros(6,6,DoF+2);
Xst(:,:,1) = eye(6); % the initial zero transformation for the base.
% Link velocity, initial value.
vli = zeros(6,1);
% Body or Link Forces, inital value.
fli = zeros(6,DoF);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recursive Algorithm OUTWARDS PASS from Base (Spatial Frame) to Tool
% to calculate Velocities and Accelerations for each Link.
% Besides, the Link forces required are aso calculated.
%
for i = 1:DoF
    % The spatial vector transformation X is obtained with the POE instead
    % of the DH parameters. All Link Coordinate frames have the same
    % orientation of the Spatial frame and are positioned on the points of
    % the Twist for each Joint
    PoE = PoE * ForwardKinematicsPOE(TwMagFbc(:,i));
    Hs0 = Hs0 * trvP2tform(Point(:,i) - Pre);
    Pre = Point(:,i); % Position to the previous Link Frame
    Hsi = PoE * Hs0;
    % Once we have the homogeneous transformation for all Link frames, we
    % store them in form of Spatial Vector Plücker transformation X.
    Xst(:,:,i+1) = tform2xpluc(Hsi);
    % Xst(:,:,i+1) \ Xst(:,:,i) is the X transform spatial MOTION vector
    % between Xi to Xi-1, which is calculated from the FK expressions as
    % Xi_i-1 = inv(X0i) * X0_i-1 = Xi0 * X0_i-1.
    Xi_im1 =  Xst(:,:,i+1) \ Xst(:,:,i);
    % Joint transmitted velocity is obtained with the Joint Motion Subspace
    % & the Joint velocites in the Joint space Thp.
    vji = S(:,i) * ThetapVAL(i);
    % Link velocity is obtained with the Joint transmited velocity and the 
    % link velocity of the predecessor link, trasformed to succersor link
    % with spatial vector MOTION transformation X.
    vli = vji + Xi_im1 * vli;
    % Link acceleration with the Featherstone formulation.
    ai =  Xi_im1 * ai + S(:,i) * ThetappVAL(i) + spavec2crm(vli) * vji;
    % fli are the necessary forces applied to the links are calculated
    Ii = LinkInertiaB(CM(:,i)-Point(:,i),diag(IT(:,i)),mass(i));
    fli(:,i) = Ii * ai + spavec2crf(vli) * Ii * vli; 
end
% This last step is not necessary to solve the ID for Joint torques (Force)
% but is useful to complete the FK analysis of the robot
% besides is used for the next implementation of the inwards pass.
% Position of the Tool to the previous Link Frame
Hs0 = Hs0 * trvP2tform(pp - Pre) * rotY2tform(pi);
Hsi = PoE * Hs0;
Xst(:,:,i+2) = tform2xpluc(Hsi);
%
% It is possible to test the FK with both the direct POE and SVA evolution
% from base to tool, thoughout all the Link frames.
%Hst_POE = ForwardKinematicsPOE(TwMag) * Hst0
%Hst_SVA = Hsi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recursive Algorithm INWARDS PASS from the Tool to Base (Spatial Frame)
% to calculate Joint Torques (& Forces) required are calculated.
%
ID_SVA = zeros (DoF,1); % 
fji = zeros(6,1); % the Joint in the Tool has zero Force value
%
for i = DoF:-1:1
    % fji is the value of the Joint Force. Input i+1 & output i.
    % fli is the value of the Link Force
    % (Xst(:,:,i+2) \ Xst(:,:,i+1))' is the value of X transform vector
    % between Xi+1 to Xi, which is calculated from the FK expressions as
    % Xi+1_i = inv(X0_i+1) * X0_i = Xi+1_0 * X0i.
    % Attention because the Plücker transformation for forcers works
    % in a complementary way to the transformation for motion & velocities.
    % Joint force is obtained with the Link force and the Joint force of 
    % the succesor link, transformed to predecessor link with spatial
    % vector FORCE transformation X.
    % In fact, you can manage transformations without any distinction of 
    % MOTION and FORCE spatial vector, if you follow the logic of the
    % Plücker definition of coordinate transformations.
    % Then, with the TRANSPOSE spatial vector transformation X for motion,
    % it is possible to pass forces from the origin to the end of the X.
    %fji = fli(:,i) + inv(Xst(:,:,i+1) \ Xst(:,:,i+2))' * fji;
    fji = fli(:,i) + (Xst(:,:,i+2) \ Xst(:,:,i+1))' * fji;
    % The Joint Motion Subspace moves the Joint forces to the Joint-Space.
    ID_SVA(i) = S(:,i)' * fji; 
end
%
% Inverse Dynamics solution for the joint TORQUES T Direct (no control).
% FEED-FORWARD contribution to joint TORQUES T.
TorFF = ID_SVA;
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
%TorFB = MtST24RJsl*(Kp.*ThErr + Kv.*ThpErr)';
TorFB = (Kp.*ThErr + Kv.*ThpErr)';
%
% The total TORQUE is the sum of INVERSE DYNAMICS and CONTROL calculations.
TorTot = (TorFF - TorFB)';
%
% The Theta POSITION values are limited by the joints position limits.
TorOut = jointmag2limits(TorTot, Tormax, Tormin);
%
end
%