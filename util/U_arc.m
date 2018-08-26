function [xx,yy] = U_arc(x,y,rx,ry,t,N,th)
%Number of points, arc start*end angles (in radians)
if~exist('N','var');N = 100;end        %default: 100 points
if~exist('th','var');th = [0 2*pi];end  %default: full ellipse

%distribute N points between arc start & end angles
th = linspace(th(1),th(2),N);

%calculate x and y points
xx = x + rx*cos(th)*cos(t) - ry*sin(th)*sin(t);
yy = y + rx*cos(th)*sin(t) + ry*sin(th)*cos(t);

