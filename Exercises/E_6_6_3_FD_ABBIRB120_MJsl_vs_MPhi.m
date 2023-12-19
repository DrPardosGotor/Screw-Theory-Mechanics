%% Screw Theory - INERTIA MATRIX "M" for DYNAMICS - Classical STR vs. SVA.
% "M" Geometric Link Jacobian Inertia Matrix - classical STR.
% "M" Physical Interpretation Inertia Matrix - SVA.
% Both are "M" Joint-Space Inertia Matrix.
%
% ABB IRB120 Home position Elbow & Tool down
% & Gravity acting in direction -Z (gz).
%
% The goal of this exercise is to calculate the INERTIA MATRIX "M" with two
% different approaches.
%
% First with the classical Screw Theory for Robotics.
% by Dr. Pardos-Gotor ST24R "Screw Theory Toolbox for Robotics" MATLAB.
% M(t)*ddt + C(t,dt)*dt + N(t,dt) = T
%
% Second with the Spatial Vector Algebra.
% by the COMPOSITE-RIGID-BODY Algorithm by Featherstone
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
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
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
%% E721a_ST24R_IM_ABBIRB120_CLAvSVA
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
% RANDOM definition of POSITION, VELOCITY and ACCELERANTIONS for the
% Joints of the Robot
% TARGET defined by JOINT Position, Velocity & Acceleration
% It is only one random target point and the differentiability of the
% position and velocity trajectory is given for granted. Here we are
% concerned with the Dynamic solution for a single trajectory target.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
Th = zeros(1,DoF); Thp = zeros(1,DoF); Thpp = zeros(1,DoF);
for i = 1:DoF
    Th(i) = rand*Thmax(i)-rand*Thmin(i); % for testing various Theta POS
    Thp(i) = rand*Thpmax(i)-rand*Thpmin(i); % for testing various Theta VEL
    Thpp(i) = rand*Thpmax(i)-rand*Thpmin(i); % for testing various Theta ACC
end
%
% To apply the Physical interpretation of the M Inertia Matrix.
Thp = zeros(1,6);
Thpp = ones(1,6);
%
% Twist and Magnitude for the Joint position (Th).
TwMag = [Twist; Th];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "M" with the classical Screw Theory for Robotics, Closed-Solution.
% by Dr. Pardos-Gotor ST24R "Screw Theory Toolbox for Robotics" MATLAB.
% M(t)*ddt + C(t,dt)*dt + N(t,dt) = T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
% M(t) Inertia matrix by the use of Jsl LINK TOOL Jacobian.
MtST24RJsl = MInertiaJsl(TwMag,LiMas)
%
toc
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "M" with the Spatial Vector Algebra.
% by the "M" Physical Interpretation Inertia Matrix by Featherstone
% but with the screw theory POE for the management of the robot kinematics
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
% ali stores the Spatial Vector ACCELERATIONS for the LINKS
% For this algorithm THIS IS THE KEY IDEA, fixing the acceleration 
% of the links to its UNITARY VALUE
ali = S;
%
% Ili Inertia Matrix for the LINKS
Ili = zeros(6,6,DoF);
%
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
    % fli Spatial vector Link Forces.
    Ili(:,:,i) = LinkInertiaB(CM(:,i)-qi,diag(IT(:,i)),mass(i));
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
% to calculate the "M" Composite-Rigid-Body Inertia Matrix - SVA.
%
MtSVA = zeros (DoF,DoF); % Inertia Matrix initialized to zero.
Iicp1 = zeros (6,6);
%
for i = DoF:-1:1
    %
    % (Xst(:,:,i+1) \ Xst(:,:,i)' is the value of X transform vector
    % between Xi+1 to Xi, which is calculated from the FK expressions as
    % Xi+1_i = inv(X0_i+1) * X0_i = Xi+1_0 * X0i.
    % In so doing, you can manage transformations without any distinction
    % of MOTION and FORCE spatial vector, if you follow the logic of the
    % Plücker definition of coordinate transformations.
    % Then, with the TRANSPOSE spatial vector transformation X for motion,
    Xip1_i = Xip1 \ Xst(:,:,i);
    Iic = Ili(:,:,i) + (Xip1_i)' * Iicp1 * Xip1_i;
    Fi = Iic * ali(:,i);
    MtSVA(i,i) = S(:,i)' * Fi;
    for j = i:-1:1
        if j > 1
            Xjm1 = Xst(:,:,j-1);
        else
            Xjm1 = Xli0;
        end            
        Fi = (Xst(:,:,j) \ Xjm1)'* Fi;
        if j > 1 
            MtSVA(i,j-1) = Fi' * S(:,j-1);
            MtSVA(j-1,i) = MtSVA(i,j-1)';
        end
    end
    %
    % Preparing the value for the next loop pass
    Xip1 = Xst(:,:,i);
    Iicp1 = Iic;
end
MtSVA
%
toc
%