% MIT License
%
% Copyright (c) 2017 Roman Szewczyk
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
% 
%
% DESCRIPTION:
% Demonstration script to show the results of calculations
%
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% R. Szewczyk "Magnetostatic modelling of thin layers using the method of moments 
%              and its implementation in Octave/Matlab", Springer, 2018. 
%              DOI 10.1007/978-3-319-77985-0, ISBN 978-3-319-77984-3
%
%

load('res_calc.mat');

mi0=4.*pi.*1e-7;

[xx,yy]=meshgrid(1:N,1:N);

xx=xx.*dL.*1000;
yy=yy.*dL.*1000;


figure (1,'Position',[100,100,720,550]);
colormap('copper');
pcolor(xx,yy,mi);
set(gca, "fontsize", 14);
xlabel('x (mm)',"fontweight", "bold");
ylabel('y(mm)',"fontweight", "bold");
h=colorbar;
set(h,'fontsize',14)


figure (2,'Position',[100,100,720,550]);
contourf(xx,yy,Bx(:,:,5));
set(gca, "fontsize", 14);
xlabel('x (mm)',"fontweight", "bold");
ylabel('y(mm)',"fontweight", "bold");
h=colorbar;
set(h,'fontsize',14);

figure (3,'Position',[100,100,720,550]);
contourf(xx,yy,By(:,:,5));
set(gca, "fontsize", 14);
xlabel('x (mm)',"fontweight", "bold");
ylabel('y(mm)',"fontweight", "bold");
h=colorbar;
set(h,'fontsize',14);


figure (4,'Position',[100,100,720,550]);
contourf(xx,yy,Bx(:,:,1));
set(gca, "fontsize", 14);
xlabel('x (mm)',"fontweight", "bold");
ylabel('y(mm)',"fontweight", "bold");
h=colorbar;
set(h,'fontsize',14);

figure (5,'Position',[100,100,720,550]);
contourf(xx,yy,By(:,:,1));
set(gca, "fontsize", 14);
xlabel('x (mm)',"fontweight", "bold");
ylabel('y(mm)',"fontweight", "bold");
h=colorbar;
set(h,'fontsize',14);


