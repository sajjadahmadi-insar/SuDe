clear;close all;clc
%% Script to produce simulated interferograms using okada formulations

% Sajjad Ahmadi, University of Tehran, October 2023
% ahmadi.sajjad@ut.ac.ir

% Requires the following additional file:
% disloc.m

%% Definition of geometric parameters of fault plane

Quake.Slip = 3;                                  % fault Slip in m
Quake.Length = 46;                               % fault Length in Km
Quake.Top = 12;                                  % depth to Top of fault in km
Quake.Bottom = 20;                               % depth to Bottom of fault in km

%% Define Source Parameters from focal mechanism

Quake.Strike = 335;                              % Strike in degrees
Quake.Dip = 16;                                  % Dip in degrees
Quake.Rake = 120;                                % Rake in degrees

%% Set Wavelength, and Heading and Incidence Angles for Satellite

Radar_Wavelength = 0.056;                        % Wavelength in m
Heading = -11;                                   % Heading (azimuth) of satellite measured clockwise from North, in degrees
Incidence = 36;                                  % Incidence angle of satellite in degrees

%% Set the regular grid parameters for visualization

x = -100; % in km
y = 100; % in km

row = 500; % number of pixels in X direction (Shoud be even)
col = 500; % number of pixels in Y direction (Shoud be even)


%%                                                                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                   DO NOT EDIT ANY TEXT BELOW HERE                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                                                       %%

%% Define Elastic Lame parameters

lambda = 2.3e10;                                         %units = pascals
mu = 2.3e10; 
v = lambda / (2*(lambda + mu));                          %calculate poisson's ration

%% okada model

Quake.Depth = (Quake.Bottom+Quake.Top)/2;                % Depth of earthquake in Km
Quake.Width = (Quake.Bottom-Quake.Top)/sind(Quake.Dip);  % fault Width in Km
xx=linspace(x,y,row);yy=linspace(x,y,col);
[E,N] = meshgrid(xx,yy);
[uE,uN,uZ] = disloc(E,N,Quake.Depth,Quake.Strike,Quake.Dip,Quake.Length,Quake.Width,Quake.Rake,Quake.Slip,v);

%% fault plane

DEG2RAD = pi/180;
alpha = pi/2 - Quake.Strike*DEG2RAD;
x_fault = Quake.Length/2*cos(alpha)*[-1,1,1,-1] + sin(alpha)*cosd(Quake.Dip)*Quake.Width/2*[-1,-1,1,1];
y_fault = Quake.Length/2*sin(alpha)*[-1,1,1,-1] + cos(alpha)*cosd(Quake.Dip)*Quake.Width/2*[1,1,-1,-1];

%% fault trace

xt = -Quake.Top/tand(Quake.Dip)-(Quake.Bottom-Quake.Top)/(2*tand(Quake.Dip));
yt = -Quake.Length/2;
yt1 = -yt;
rotM = [cosd(Quake.Strike) sind(Quake.Strike)
       -sind(Quake.Strike) cosd(Quake.Strike)];
xy1 = rotM*[xt;yt];
xy2 = rotM*[xt;yt1];
x11 = xy1(1);
y11 = xy1(2);
x22 = xy2(1);
y22 = xy2(2);

%% LOS_vector

LOS_vector = [sind(Incidence)*sind(Heading) -sind(Incidence)*cosd(Heading) cosd(Incidence)];
unwrap = LOS_vector(1)*uN+LOS_vector(2)*uE+LOS_vector(3)*uZ;
wrap = wrapToPi(-4*pi/Radar_Wavelength*unwrap);

%% Wrapped LOS

figure(1)
imagesc(E(:)',N(:)',(wrap')')
hold on
patch(x_fault,y_fault,0.3*[1,1,0],'EdgeColor','k','FaceColor','none','LineWidth',2)
plot(x11,y11,'ko');
plot([x11,x22],[y11,y22],'k--','LineWidth',2)

axis xy
axis image
title('Wrapped LOS')
xlabel('Easting (km)')
ylabel('Northing (km)')
h=colorbar('vert');
ylabel(h,'radians')
colormap(flipud(jet))
hold off

%% Unwrapped LOS

figure(2)
imagesc(E(:)',N(:)',unwrap)
hold on
patch(x_fault,y_fault,0.3*[1,1,1],'EdgeColor','k','FaceColor','none','LineWidth',1)
plot(x11,y11,'ko');
plot([x11,x22],[y11,y22],'k--','LineWidth',1)

axis xy
axis image
title('Unwrapped LOS')
xlabel('Easting (km)')
ylabel('Northing (km)')
h=colorbar('vert');
ylabel(h,'metres')
colormap(jet)
hold off

%% 3D Deformation

figure(3)
imagesc(E(:)',N(:)',uZ)
hold on
[E,N] = meshgrid(xx,yy);
stepc=ceil(col/40);
stepr=ceil(row/40);
quiver(E(1:stepc:end,1:stepr:end),N(1:stepc:end,1:stepr:end), ...
       uE(1:stepc:end,1:stepr:end),uN(1:stepc:end,1:stepr:end),'w')
patch(x_fault,y_fault,0.3*[1,1,0],'EdgeColor','k','FaceColor','none','LineWidth',1)
plot(x11,y11,'ko');
plot([x11,x22],[y11,y22],'k--','LineWidth',1)

axis xy
axis image
title('3D Deformation')
xlabel('Easting (km)')
ylabel('Northing (km)')
h=colorbar('vert');
ylabel(h,'metres')
colormap(jet)
hold off

%% N-S Deformation 

figure(4)
imagesc(E(:)',N(:)',uN)
hold on
[E,N] = meshgrid(xx,yy);
patch(x_fault,y_fault,0.3*[1,1,0],'EdgeColor','k','FaceColor','none','LineWidth',1)
plot(x11,y11,'ko');
plot([x11,x22],[y11,y22],'k--','LineWidth',1)

axis xy
axis image
title('N-S Deformation')
xlabel('Easting (km)')
ylabel('Northing (km)')
h=colorbar('vert');
ylabel(h,'metres')
colormap(flipud(jet))
hold off

%% E-W Deformation 

figure(5)
imagesc(E(:)',N(:)',uE)
hold on
[E,N] = meshgrid(xx,yy);
patch(x_fault,y_fault,0.3*[1,1,0],'EdgeColor','k','FaceColor','none','LineWidth',1)
plot(x11,y11,'ko');
plot([x11,x22],[y11,y22],'k--','LineWidth',1)

axis xy
axis image
title('E-W Deformation')
xlabel('Easting (km)')
ylabel('Northing (km)')
h=colorbar('vert');
ylabel(h,'metres')
colormap(flipud(jet))
hold off

%% Displacement Profiles of 3D Deformation along Northing=0

figure (6)
plot(E(col/2+1,:)',uZ(col/2+1,:),'b')
xlabel('Easting (km)')
ylabel('Displacement (m)')
hold on
plot(E(col/2+1,:)',uE(col/2+1,:),'r')
plot(E(col/2+1,:)',uN(col/2+1,:),'g')
legend ('Vert','E-W','N-S');
title('Displacement Profiles of 3D Deformation')
hold off

%% All plots

figure (7)
subplot(2,3,1)
imagesc(E(:)',N(:)',(wrap')')
hold on
patch(x_fault,y_fault,0.3*[1,1,0],'EdgeColor','k','FaceColor','none','LineWidth',2)
plot(x11,y11,'ko');
plot([x11,x22],[y11,y22],'k--','LineWidth',2)

axis xy
axis image
title('Wrapped LOS')
xlabel('Easting (km)')
ylabel('Northing (km)')
h=colorbar('vert');
ylabel(h,'radians')
colormap(flipud(jet))
hold off

subplot(2,3,2)
imagesc(E(:)',N(:)',unwrap)
hold on
patch(x_fault,y_fault,0.3*[1,1,1],'EdgeColor','k','FaceColor','none','LineWidth',1)
plot(x11,y11,'ko');
plot([x11,x22],[y11,y22],'k--','LineWidth',1)

axis xy
axis image
title('Unwrapped LOS')
xlabel('Easting (km)')
ylabel('Northing (km)')
h=colorbar('vert');
ylabel(h,'metres')
colormap(jet)
hold off

subplot(2,3,3)
plot(E(col/2+1,:)',uZ(col/2+1,:),'b')
xlabel('Easting (km)')
ylabel('Displacement (m)')
hold on
plot(E(col/2+1,:)',uE(col/2+1,:),'r')
plot(E(col/2+1,:)',uN(col/2+1,:),'g')
legend ('Vert','E-W','N-S');
title('Displacement Profiles of 3D Deformation')
hold off

subplot(2,3,4)
imagesc(E(:)',N(:)',uZ)
hold on
[E,N] = meshgrid(xx,yy);
quiver(E(1:stepc:end,1:stepr:end),N(1:stepc:end,1:stepr:end), ...
       uE(1:stepc:end,1:stepr:end),uN(1:stepc:end,1:stepr:end),'w')
patch(x_fault,y_fault,0.3*[1,1,0],'EdgeColor','k','FaceColor','none','LineWidth',1)
plot(x11,y11,'ko');
plot([x11,x22],[y11,y22],'k--','LineWidth',1)

axis xy
axis image
title('3D Deformation')
xlabel('Easting (km)')
ylabel('Northing (km)')
h=colorbar('vert');
ylabel(h,'metres')
colormap(jet)
hold off

subplot(2,3,5)
imagesc(E(:)',N(:)',uN)
hold on
[E,N] = meshgrid(xx,yy);
patch(x_fault,y_fault,0.3*[1,1,0],'EdgeColor','k','FaceColor','none','LineWidth',1)
plot(x11,y11,'ko');
plot([x11,x22],[y11,y22],'k--','LineWidth',1)

axis xy
axis image
title('N-S Deformation')
xlabel('Easting (km)')
ylabel('Northing (km)')
h=colorbar('vert');
ylabel(h,'metres')
colormap(flipud(jet))
hold off

subplot(2,3,6)
imagesc(E(:)',N(:)',uE)
hold on
[E,N] = meshgrid(xx,yy);
patch(x_fault,y_fault,0.3*[1,1,0],'EdgeColor','k','FaceColor','none','LineWidth',1)
plot(x11,y11,'ko');
plot([x11,x22],[y11,y22],'k--','LineWidth',1)

axis xy
axis image
title('E-W Deformation')
xlabel('Easting (km)')
ylabel('Northing (km)')
h=colorbar('vert');
ylabel(h,'metres')
colormap(flipud(jet))
hold off
