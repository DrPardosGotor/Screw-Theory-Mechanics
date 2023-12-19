%% Screw Theory in Robotics - FORWARD DYNAMICS
% ABB IRB120 Home position Elbow & Tool down
% & Gravity acting in direction -Z (gz).
%
% FIRST: The goal is to prove the FORWARD DYNAMICS Classical
% Lagrange to obtain joint Accelerations, knowing joint Torques (or forces)
% ddt= inv(M(t)) * (T  - C(t,dt) * dt + N(t)) 
% by Dr. Pardos-Gotor ST24R "Screw Theory Toolbox for Robotics" MATLAB.
%
% SECOND with the Spatial Vector Algebra.
% and the Articulated Body Inertia Algorithm by Featherstone
% but with the screw theory POE for the management of the robot kinematics
%
% Copyright (C) 2003-2021, by Dr. Jose M. Pardos-Gotor.
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
% You should have received a copy of the GNU Lesser General Public License
% along with ST24R.  If not, see <http://www.gnu.org/licenses/>.
%
% http://www.
%
% CHANGES:
% Revision 1.1  2021/02/11 00:00:01
% General cleanup of code: help comments, see also, copyright
% references, clarification of functions.
%
%%
%
clear;
clc;
%
% Potential Action Vector - Gravity definition (i.e., -g direction).
PoAcc = [0 0 -9.81]';
%
% Degress of Freedon of the Robot
DoF = 6;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mechanical characteristics of the Robot (AT REF HOME POSITION):
% kinematics defined with the screw theory POE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
po=[0;0;0]; pg=[0; 0; 0.19]; pk=[0; 0; 0.29]; pr=[0; 0; 0.56]; ps=[0; 0; 0.63];
pf=[0.302; 0; 0.63]; pu=[0.302; 0; 0.558]; pp=[0.302; 0; 0.47];
AxisX = [1 0 0]'; AxisY = [0 1 0]'; AxisZ = [0 0 1]'; 
Point = [pg pk pr ps pf pu];
Joint = ['rot'; 'rot'; 'rot'; 'rot'; 'rot'; 'rot'];
Axis = [AxisZ AxisY AxisY AxisX AxisY -AxisZ];
Twist = zeros(6,DoF);
for i = 1:DoF
    Twist(:,i) = joint2twist(Axis(:,i), Point(:,i), Joint(i,:));
end
Hst0 = trvP2tform(pp)*rotY2tform(pi);
%
% Motion RANGE for the robot joints POSITION rad, (by catalog).
Thmax = pi/180*[165 110 70 160 120 400];
Thmin = -pi/180*[165 110 110 160 120 400];
% Maximum SPEED for the robot joints rad/sec, (by catalog).
Thpmax = pi/180*[250 250 250 320 320 420];
Thpmin = pi/180*[250 250 250 320 320 420];
%
% DYNAMIC Parameters of the Robot at REF HOME POSITION - Only aproximation
CM = [0 0 0.29; 0 0 0.425; 0 0 0.63; 0.2 0 0.63; 0.302 0 0.63; 0.302 0 0.53]';
IT = [0.1 0.3 0.2; 0.3 0.5 0.1; 0.1 0.1 0.1; 0.1 0.3 0.2; 0.1 0.1 0.1; 0.1 0.1 0.1]';
mass = [7 6 5 4 2 1];
LiMas = [CM; IT; mass];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TARGET RANDOM definition of POSITION, VELOCITY and TORQUES for the
% Joints of the Robot
% It is only one random target point and the differentiability of the
% position and velocity trajectory is given for granted. Here we are
% concerned with the Dynamic solution for a single trajectory point.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
Th = zeros(1,DoF); Thp = zeros(1,DoF)'; Tor = zeros(1,DoF)';
for i = 1:DoF
    Th(i) = rand*Thmax(i)-rand*Thmin(i); % for testing various Theta POS
    Thp(i) = rand*Thpmax(i)-rand*Thpmin(i); % for testing various Theta VEL
    Tor(i) = rand*Thpmax(i)-rand*Thpmin(i); % for testing various Theta ACC
end
%
% Twist and Magnitude for the Joint position (Th).
TwMag = [Twist; Th];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLASSICAL SCREW THEORY FORWARD DYNAMICS solution
% Third with the WRENCH algorithm.
% This is with the new gravity wrench matrix for N(t).
tic;
%
Tor
%
% M(t) Inertia matrix by the use of Jsl LINK TOOL Jacobian.
MtST24RJsl = MInertiaJsl(TwMag,LiMas);
%
% C(t,dt) Coriolis matrix by the use of Aij Adjoint transformation.
CtdtST24RAij = CCoriolisAij(TwMag,LiMas,Thp);
%
% N(t) Potential by the use of the new GRAVITY WRENCH Matrix.
NtST24RWre = NPotentialWre(TwMag,LiMas,PoAcc);
%
% Forward Dynamics to get Joint accelerations Thpp from joint torques Tor.
Thpp = MtST24RJsl \ (Tor - CtdtST24RAij*Thp - NtST24RWre)
%
toc
%
% Check the result solving the Classical ID
%TdynST24RWre = MtST24RJsl*Thpp + CtdtST24RAij*Thp + NtST24RWre
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FD with the Spatial Vector Algebra.
% with the Articulated Body Inertia Algorithm by Featherstone but 
% with Screw Theory POE to define the kinematics transformations X(6x6)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matrices of Inertia (I1...I3) defined in LINK Frame with S orientation.
% the postion on the CM is defined in relation to each LINK Frame.
%
% Define the spatial vector motion Subspace “Si” for each joint,
S = [Axis; zeros(3,DoF)];
%
% Xst stores the Spatial Vector TRANSFORMATION to the LINKS Frames
Xst = zeros(6,6,DoF);
%
% vli stores the Spatial Vector VELOCICITIES for the LINKS
vli = zeros(6,DoF);
%
% Ali stores the Spatial Vector ACCELERATIONS for the LINKS
ali = zeros(6,DoF);
%
% Ilii Inertia Matrix for the LINKS
Ili = zeros(6,6,DoF);
%
% fli stores the Spatial Vector FORCES for the LINKS
fli = zeros(6,DoF);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This block of variables is useful for the preparation of the
% first recursive OUTWARDS PASS
%
% Initial values for the BASE or Link Zero.
% Define the spatial coordinate system “S” with a homogeneous matrix H0(0),
% according to the free selection for its orientation matrix “R0” and 
% position “q0”.
R0 = eye(3);
q0 = [0 0 0]';
H0Ze = [R0 q0; 0 0 0 1];
% Initial value Product of Exponentials.
POE0 = eye(4);
% the initial transformation, no rotation no translation.
Xli0 = eye(6);
% Velocity for Link Zero is null.
vli0 = zeros(6,1);
% Gravitational field is modelled as a fictitious acceleration of the base.
% To do this, we replace the initial value with a0 = - ag = - PoAcc
% With this trick the values for ai and fi (net force body) are no longer
% the true values, as they are offset by gravity acceleration and force.
% Nonetheless, the result of the ID is correct, as we are interested only
% in the Joint Torques.
ali0 = [0;0;0; -PoAcc];
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recursive Algorithm OUTWARDS PASS from Base (Spatial Frame) to Tool
% to calculate Velocities and Accelerations for each Link.
% Besides, the Link forces required are aso calculated.
%
tic;
%
qim1 = q0;
Him10 = H0Ze;
POEim1 = POE0;
Xim1 = Xli0;
vlim1 = vli0;
alim1 = ali0;
%
for i = 1:DoF
    % The spatial vector transformation X is obtained with the POE instead
    % of the DH parameters. All Link Coordinate frames have the same
    % orientation of the Spatial frame and are positioned on the points of
    % the Twist for each Joint
    qi = Point(:,i);
    Hi0 = trvP2tform(qi - qim1) * Him10;
    POEi = POEim1 * ForwardKinematicsPOE(TwMag(:,i));
    HiThe = POEi * Hi0;
    % Once we have the homogeneous transformation for all Link frames, we
    % store them in form of Spatial Vector Plücker transformation X.
    Xst(:,:,i) = tform2xpluc(HiThe);
    % Preparing the value for the next loop pass
    qim1 = qi;
    Him10 = Hi0;
    POEim1 = POEi;
    %
    % Xi_im1 is the X transform spatial vector MOTION vector
    % between thee link_i to link_i-1, using the absolute definitions
    % of Xi & Xi-1 in respect to the Spatial system.
    Xi_im1 =  Xst(:,:,i) \ Xim1;
    % Preparing the value for the next loop pass
    Xim1 = Xst(:,:,i);
    %
    % Joint transmitted velocity is obtained with the Joint Motion Subspace
    % & the Joint velocites in the Joint space Thp.
    vji = S(:,i) * Thp(i);
    %
    % vli Spatial vector Link Velocity.
    vlii = Xi_im1 * vlim1 + vji;
    vli(:,i) = vlii;
    % Preparing the value for the next loop pass
    vlim1 = vlii;
    %
    % ali Spatial vector Link Acceleration.
    alii =  Xi_im1 * alim1 + S(:,i) * Thpp(i) + spavec2crm(vlii) * vji;
    ali(:,i) = alii;
    % Preparing the value for the next loop pass
    alim1 = alii;
    %
    % fli Spatial vector Link Forces.
    Ili(:,:,i) = LinkInertiaB(CM(:,i)-qi,diag(IT(:,i)),mass(i));
    fli(:,i) = Ili(:,:,i) * alii + spavec2crf(vlii) * Ili(:,:,i) * vlii;
    %
end
% This last step is not necessary to solve the ID for Joint torques (Force)
% but is useful to complete the FK analysis of the robot
% besides is used for the next implementation of the inwards pass.
% Position of the Tool to the previous Link Frame
Hsp0 = trvP2tform(pp - qim1) * Him10 * rotY2tform(pi);
HiThe = POEim1 * Hsp0;
Xip1 = tform2xpluc(HiThe);
%
% It is possible to test the FK with both the direct POE and SVA evolution
% from base to tool, thoughout all the Link frames.
%Hst_POE = ForwardKinematicsPOE(TwMag) * Hsp0
%Hst_SVA = HiThe
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recursive Algorithm INWARDS PASS from the Tool to Base (Spatial Frame)
% to calculate Joint Torques (& Forces) required are calculated.
%
% Initial values for the Tool
% fjip1 Force for the successor joint to the last link.
fjip1 = zeros(6,1);
%
ID_SVA = zeros (DoF,1);
%
for i = DoF:-1:1
    % fji is the value of the Joint Force. Input i+1 & output i.
    % fli is the value of the Link Force
    % (Xst(:,:,i+1) \ Xst(:,:,i))' is the value of X transform vector
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
    % fji = fli(:,:,i) + inv(Xst(:,:,i) \ Xip1)' * fji;
    fji = fli(:,i) + (Xip1 \ Xst(:,:,i))' * fjip1;
    % Preparing the value for the next loop pass
    fjip1 = fji;
    Xip1 = Xst(:,:,i);
    % The Joint Motion Subspace moves the Joint forces to the Joint-Space.
    ID_SVA(i) = S(:,i)' * fji; 
end
ID_SVA
%
toc
%