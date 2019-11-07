%% t_fixBrainstormHeadshapeAlignment

% This tutorial is under construction, but started with the intention to
% show how to fix misalignment of the digitized headshape and the spheres
% made by Brainstorm (extracted from a T1 scan).

% If you have the issue that the headshape xyz points are reflected, you
% can run this script first and then keep the new R matrix and T vector in
% the workspace memory. Then put a stop at in_fopen_kit.m. line 188. This is
% just after Brainstorm has created R and T.

% Next, run brainstorm and import MEG data. When Matlab stops at the break
% point in_fopen_kit.m, I use Matlab?s evalin function to copy the R and T
% (i.e. R = evalin('base','R') variables from the base workspace to the
% in_fopen_kit?s workspace. Make sure they are in the correct orientation
% and I also checked whether one of the R matrices computed by my code is
% the same as the R matrix computed by brainstorm.



%  Example with: 10_SSMEG_08_12_2014_wl_subj004
% subjFolder = '/Volumes/server-1/Projects/MEG/SSMEG/10_SSMEG_08_12_2014_wl_subj004/raw/';
% markerSQDFile   = fullfile(subjFolder, 'R0774_Marker1_8.12.14.sqd');
% HSTxtFile       = fullfile(subjFolder, 'R0774_8.12.14_HS.txt');
% PointsTxtFile   = fullfile(subjFolder, 'R0774_8.12.14_Points.txt');

% subjFolder = '/Volumes/server/Projects/MEG/SSMEG/fullOnly/03_SSMEG_R1374_03.14.2018/raw/';
% markerSQDFile   = fullfile(subjFolder, 'R1374_Marker1_03.14.2018.sqd');
% HSTxtFile       = fullfile(subjFolder, 'R1374_03.14.18_HS.txt');
% PointsTxtFile   = fullfile(subjFolder, 'R1374_03.14.18_Points.txt');

subjFolder = '/Volumes/server/Projects/MEG/SSMEG/fullOnly/04_SSMEG_R1021_08.30.2018/raw/';
markerSQDFile   = fullfile(subjFolder, 'R1021_SSMEG2_run1-8.sqd');
HSTxtFile       = fullfile(subjFolder, 'R1021_8.30.18_HS.txt');
PointsTxtFile   = fullfile(subjFolder, 'R1021_8.30.18_Points.txt');

% subjFolder      = '/Volumes/server/Projects/MEG/Retinotopy/Data/MEG/wl_subj040/wl_subj040_20170406/Raw/';
% markerSQDFile   = fullfile(subjFolder, 'R1151_Marker1_04.06.17.sqd');
% HSTxtFile       = fullfile(subjFolder, 'R1151_4.6.2017_HS.txt');
% PointsTxtFile   = fullfile(subjFolder, 'R1151_4.6.2017_Points.txt');

%% Get points from the three files
headShapeLabels = {'Nasion','Left Tragus', 'Right Tragus', 'Left PA', 'Right PA', ...
    'Central frontal', 'Left frontal', 'Right frontal'};

labels = {'LPA', 'RPA', 'CPF', 'LPF', 'RPF'}; % Obtained by Brainstorm after loading in the SQD Marker1 file

% Open Points textfile
fid = fopen(PointsTxtFile, 'r');

txtCell = textscan(fid, '%f%f%f', 'Delimiter', '\n', 'CommentStyle', '%');
fclose(fid);
points_xyz = [txtCell{1}, txtCell{2}, txtCell{3}]' ./ 1000;
points_xyz = points_xyz';
points_xyz = points_xyz(4:end,:);

% Open Head Sweeps txt file
fid = fopen(HSTxtFile, 'r');
txtCell = textscan(fid, '%f%f%f', 'Delimiter', '\n', 'CommentStyle', '%');
fclose(fid);
hs_xyz = [txtCell{1}, txtCell{2}, txtCell{3}]' ./ 1000;

step_size = round(size(hs_xyz,2) / 1000);
hs_xyz = hs_xyz(:,1:step_size:end);

% Open Marker SQD file
header.coreg = getYkgwHdrCoregist(markerSQDFile);
header.sensors = getYkgwHdrChannel(markerSQDFile);
meg_xyz = cat(1, header.coreg.hpi.meg_pos)';

% Get XYZ position of MEG Sensors:
meg_sensors = [];
for ii = 1:157
    meg_sensors(ii,:) = [header.sensors.channel(ii).data.x,header.sensors.channel(ii).data.y, header.sensors.channel(ii).data.z];
end

%% Original rotation + Visualize

% Compute the transformation matrix for data set:
[R_wl, T_wl] = rot3dfit(points_xyz, meg_xyz');
Tmat_wl = [R_wl', T_wl'; 0 0 0 1];

% Visualize original rotation/translation
% 1. Apply new rotation and transformation
reorient_hs_xyz = R_wl * hs_xyz + T_wl'.*ones(3,size(hs_xyz,2));
reorient_points_xyz = R_wl * points_xyz' + T_wl'.*ones(3, size(points_xyz',2));
reorient_meg_sensors = meg_sensors';

% 2. Visualize
figure(1); clf; hold all; title('Reoriented headshape and points using extra point 2')

% Plot the MEG markers
plot3(reorient_points_xyz(1,:) ,reorient_points_xyz(2,:), reorient_points_xyz(3,:),...
    'R.','MarkerSize',24)

% Plot the head shape:
plot3(reorient_hs_xyz(1,:),reorient_hs_xyz(2,:), reorient_hs_xyz(3,:), 'k.');

% Plot the sensor's of the meg
plot3(reorient_meg_sensors(1,:),reorient_meg_sensors(2,:),reorient_meg_sensors(3,:), 'gd', 'MarkerSize',20);

xlabel('x')
ylabel('y')
zlabel('z')

%% Find plane and add extra point

% Use affine fit to find the normal of the plane through the head points (p and n are 1x3):
[n, v, p] = affine_fit(points_xyz);

% Normal is normalized, so rescale:
n = n'./100;

% Transpose points (from 3x5 to 5x3)
% points_xyz = points_xyz';

% Plot the head points:
figure;
plot3(points_xyz(:,1) ,points_xyz(:,2), points_xyz(:,3),...
    'R.','MarkerSize',24)
hold on

% Plot the normal of the surface through the plane, in both directions:
quiver3(p(1),p(2),p(3),n(1),n(2),n(3),'b')
quiver3(p(1),p(2),p(3),-n(1),-n(2),-n(3),'r')

% Plot the head shape:
plot3(hs_xyz(1,:),hs_xyz(2,:), hs_xyz(3,:), 'k.');

% Compute the mean position of the head shape (93382x3) --> (1x3):
m_hs = mean(hs_xyz,2);

% Plot mean position of the head shape:
plot3(m_hs(1), m_hs(2), m_hs(3), 'bd','MarkerFaceColor','b','MarkerSize',18)

% Determine on which side of the plane the head shape is by finding the
% normal that closest the the head shape's mean:
dist_01 = sum(((p(1:3) + -n(1:3)') - m_hs).^2);
dist_02 = sum(((p(1:3) + n(1:3)') - m_hs).^2);

% Pick the best:
if dist_01 < dist_02
    n = -n;
end

% Add this point to the head points:
hndl = p + n; %(1x3)
wl_subj_points_hndl = [points_xyz;  hndl];


% Do the same thing for the MEG marker points:
[n2, v2, p2] = affine_fit(meg_xyz');
n2 = n2./100;

% Add points along the normal orthogonal to the plane in both directions:
hndl_01 = p2 + -n2';
hndl_02 = p2 + n2';

wl_subj_markers_hndl_01 = [meg_xyz' ; hndl_01];
wl_subj_markers_hndl_02 = [meg_xyz' ; hndl_02];

% Recompute the transformation:
% This one is almost the same as the initial transformation:
[R_01, T_01] = rot3dfit(wl_subj_points_hndl, wl_subj_markers_hndl_01);

% This one is not, and when one replaces it with the transformation matrix
% Brainstorm computes, the alignment has a much better result:
[R_02, T_02] = rot3dfit(wl_subj_points_hndl, wl_subj_markers_hndl_02);





%% Apply rotations + visualize

% These are the old R and T matrices
R = R_01;
T = T_01;

% 1. Apply new rotation and transformation
reorient_hs_xyz = R' * hs_xyz  + T'.*ones(3,size(hs_xyz,2));
reorient_points_xyz = R' * points_xyz' + T'.*ones(3, size(points_xyz',2));
reorient_meg_sensors = meg_sensors';


% 2. Visualize
figure(1); clf; hold all;

% Plot the MEG markers
subplot(211); hold on; title('Reoriented headshape and points using extra point 1')
plot3(reorient_points_xyz(1,:) ,reorient_points_xyz(2,:), reorient_points_xyz(3,:),...
    'R.','MarkerSize',24)

% Plot the head shape:
plot3(reorient_hs_xyz(1,:),reorient_hs_xyz(2,:), reorient_hs_xyz(3,:), 'k.');

% Plot the sensor's of the meg
plot3(reorient_meg_sensors(1,:),reorient_meg_sensors(2,:),reorient_meg_sensors(3,:), 'gd', 'MarkerSize',20);

xlabel('x')
ylabel('y')
zlabel('z')


%% For comparison: the new R and T matrices 
R = R_02;
T = T_02;

% 1. Apply new rotation and transformation
reorient_hs_xyz = R * hs_xyz + T'.*ones(3,size(hs_xyz,2));
reorient_points_xyz = R * points_xyz' + T'.*ones(3, size(points_xyz',2));
reorient_meg_sensors = meg_sensors';

% 2. Visualize
subplot(212); hold all; title('Reoriented headshape and points using extra point 2')

% Plot the MEG markers
plot3(reorient_points_xyz(1,:) ,reorient_points_xyz(2,:), reorient_points_xyz(3,:),...
    'R.','MarkerSize',24)

% Plot the head shape:
plot3(reorient_hs_xyz(1,:),reorient_hs_xyz(2,:), reorient_hs_xyz(3,:), 'k.');

% Plot the sensor's of the meg
plot3(reorient_meg_sensors(1,:),reorient_meg_sensors(2,:),reorient_meg_sensors(3,:), 'gd', 'MarkerSize',20);

xlabel('x')
ylabel('y')
zlabel('z')

return


%% Some figures for debugging

names = {'LPA','RPA','CPF','LPF','RPF'};

left_idx = [1 4];
right_idx = [2 5];
central_idx = 3;


figure;
plot3(points_xyz(:,1),...
    points_xyz(:,2),...
    points_xyz(:,3),...
    'k.','MarkerSize',24);

text(points_xyz(:,1),...
    points_xyz(:,2),...
    points_xyz(:,3),...
    names)


%% Compare Markers handel 1 and 2
figure;
plot3(wl_subj_004_markers_hndl_02(left_idx,1), ...
    wl_subj_004_markers_hndl_02(left_idx,2), ...
    wl_subj_004_markers_hndl_02(left_idx,3),...
    'k.','MarkerSize', 24)

hold on;
plot3(wl_subj_004_markers_hndl_02(right_idx,1), ...
    wl_subj_004_markers_hndl_02(right_idx,2), ...
    wl_subj_004_markers_hndl_02(right_idx,3),...
    'ks','MarkerSize', 12,...
    'MarkerFaceColor','k')

plot3(wl_subj_004_markers_hndl_02(central_idx,1), ...
    wl_subj_004_markers_hndl_02(central_idx,2), ...
    wl_subj_004_markers_hndl_02(central_idx,3),...
    'kd','MarkerSize', 12,...
    'MarkerFaceColor','k')


plot3(wl_subj_004_markers_hndl_01(left_idx,1), ...
    wl_subj_004_markers_hndl_01(left_idx,2), ...
    wl_subj_004_markers_hndl_01(left_idx,3),...
    'r.','MarkerSize', 12,...
    'MarkerFaceColor','r')

plot3(wl_subj_004_markers_hndl_01(right_idx,1), ...
    wl_subj_004_markers_hndl_01(right_idx,2), ...
    wl_subj_004_markers_hndl_01(right_idx,3),...
    'r.','MarkerSize', 12,...
    'MarkerFaceColor','r')

plot3(wl_subj_004_markers_hndl_01(central_idx,1), ...
    wl_subj_004_markers_hndl_01(central_idx,2), ...
    wl_subj_004_markers_hndl_01(central_idx,3),...
    'r.','MarkerSize', 12,...
    'MarkerFaceColor','r')


figure; hold on;
plot3(hs_xyz(1,:), ...
    hs_xyz(2,:), ...
    hs_xyz(3,:), 'k.','MarkerSize', 12)




[R, T; 0 0 0 1]






function [n,V,p] = affine_fit(X)
%Computes the plane that fits best (lest square of the normal distance
%to the plane) a set of sample points.
%INPUTS:
%
%X: a N by 3 matrix where each line is a sample point
%
%OUTPUTS:
%
%n : a unit (column) vector normal to the plane
%V : a 3 by 2 matrix. The columns of V form an orthonormal basis of the
%plane
%p : a point belonging to the plane
%
%NB: this code actually works in any dimension (2,3,4,...)
%Author: Adrien Leygue
%Date: August 30 2013

%the mean of the samples belongs to the plane
p = mean(X,1);

%The samples are reduced:
R = bsxfun(@minus,X,p);
%Computation of the principal directions if the samples cloud
[V,D] = eig(R'*R);
%Extract the output from the eigenvectors
n = V(:,1);
V = V(:,2:end);
end


function [R,T,Yf,Err] = rot3dfit(X,Y)
%ROT3DFIT Determine least-square rigid rotation and translation.
% [R,T,Yf] = ROT3DFIT(X,Y) permforms a least-square fit for the
% linear form
%
% Y = X*R + T
%
% where R is a 3 x 3 orthogonal rotation matrix, T is a 1 x 3
% translation vector, and X and Y are 3D points sets defined as
% N x 3 matrices. Yf is the best-fit matrix.
%
% See also SVD, NORM.
%
% rot3dfit: Frank Evans, NHLBI/NIH, 30 November 2001
%

% ROT3DFIT uses the method described by K. S. Arun, T. S. Huang, and
% S. D. Blostein, "Least-Squares Fitting of Two 3-D Point Sets",
% IEEE Transactions on Pattern Analysis and Machine Intelligence,
% PAMI-9(5): 698 - 700, 1987.
%
% A better theoretical development is found in B. K. P. Horn,
% H. M. Hilden, and S. Negahdaripour, "Closed-form solution of
% absolute orientation using orthonormal matrices", Journal of the
% Optical Society of America A, 5(7): 1127 - 1135, 1988.
%
% Special cases, e.g. colinear and coplanar points, are not implemented.

% error(nargchk(2,2,nargin));
if size(X,2) ~= 3, error('X must be N x 3'); end;
if size(Y,2) ~= 3, error('Y must be N x 3'); end;
if size(X,1) ~= size(Y,1), error('X and Y must be the samesize'); end;

% mean correct

Xm = mean(X,1); X1 = X - ones(size(X,1),1)*Xm;
Ym = mean(Y,1); Y1 = Y - ones(size(Y,1),1)*Ym;

% calculate best rotation using algorithm 12.4.1 from
% G. H. Golub and C. F. van Loan, "Matrix Computations"
% 2nd Edition, Baltimore: Johns Hopkins, 1989, p. 582.

XtY = (X1')*Y1;
[U,S,V] = svd(XtY);
R = U*(V');

% solve for the translation vector

T = Ym - Xm*R;

% calculate fit points

Yf = X*R + ones(size(X,1),1)*T;

% calculate the error

dY = Y - Yf;
Err = norm(dY,'fro'); % must use Frobenius norm
end