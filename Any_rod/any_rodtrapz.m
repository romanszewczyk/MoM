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
% Demonstration script for magnetostatic modelling of rod with any crossection rod.
% Script uses the method of moments and trapezoid integration method
%
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% R. Szewczyk "Magnetostatic modelling of thin layers using the method of moments 
%              and its implementation in Octave/Matlab".
%


t_ver=ver;
if strcmp(t_ver(1).Name,'Octave')==1
   page_screen_output(0);
   page_output_immediately(1);
end

fprintf('--- Any rod calculation --- \n\n');

L=0.25;               % (m) length of the rod

N=30;                 % number of cells

dL = L./N;

mi=50;                % relative magnetic permeability
mi0=4.*pi.*1e-7;      % magnetic constant

Hext0=1000;

% Cross-section
y_range=1e-2;         % thickness: 50 um
y_step=1e-4;

z_range=1e-2;
z_step=1e-4;

nI=1000;

fprintf('Calculate baricenter of (y,z) cross-section\n');

[y,z]=meshgrid(-1.*y_range:y_step:y_range,-1.*z_range:z_step:z_range);

k=isinshape(y,z);


f0=@(y,z) isinshape(y,z);
m=dbltrapz(f0,-1.*y_range,y_range,10.*N,-1.*z_range,z_range,10.*N);

f=@(y,z) isinshape(y,z).*y;

y0=dbltrapz(f,-1.*y_range,y_range,10.*N,-1.*z_range,z_range,10.*N)./m;

f=@(y,z) isinshape(y,z).*z;

z0=dbltrapz(f,-1.*y_range,y_range,10.*N,-1.*z_range,z_range,10.*N)./m;


fprintf('Baricenter: y0=%f, z0=%f \nDone.\n\n',y0,z0);

fprintf('Calculate A matrix\n');

Hext = ones(N,1).*Hext0.*(mi-1);
A = zeros(N,N);

for k=1:N
  fprintf('%i of %i:',k,N);
  tic;
  for i=1:N

    fc0=@(ry,rz) isinshape(ry,rz).*(2.*(dL.*(i-k+0.5)).^2-(ry-y0).^2-(rz-z0).^2)./(((dL.*(i-k+0.5)).^2+(ry-y0).^2+(rz-z0).^2).^(5/2));
    fc1=@(ry,rz) isinshape(ry,rz).*(2.*(dL.*(i-1-k+0.5)).^2-(ry-y0).^2-(rz-z0).^2)./(((dL.*(i-1-k+0.5)).^2+(ry-y0).^2+(rz-z0).^2).^(5/2));
 
    A(i,k)=(mi-1).*dL./(4.*pi).*(sign(i-k+0.5).*dbltrapz(fc0,-1.*y_range,y_range,nI,-1.*z_range,z_range,nI)-sign(i-1-k+0.5).*dbltrapz(fc1,-1.*y_range,y_range,nI,-1.*z_range,z_range,nI));
    if k==i
       A(i,k)=A(i,k)+1;   
    end
    fprintf('.');
  end
  toc
end


fprintf('Done.\n');

fprintf('Solve linear equations\n');

Mt=A \ Hext;

fprintf('Done.\n');

Bt=Hext0.*mi0+Mt.*mi0;

x_g=0:L./300:L;
x_t=(1:N).*dL-0.5.*dL;
B1=ones(size(x_g)).*Hext0.*mi0.*mi;
B2=ones(size(x_g)).*Hext0.*mi0;

Bi=interp1(x_t,Bt,x_g,'nearest','extrap');


if strcmp(t_ver(1).Name,'Octave')==1
  % version for Octave
  plot(x_t,Bt,'ok','linewidth',1,'markersize', 7, x_g,Bi,'r','linewidth',2,x_g,B1,'k','linewidth',2,x_g,B2,'k','linewidth',2);
  set(gca, 'fontsize', 12);
  xlabel('x (m)','fontweight', 'bold');
  ylabel('B (T)','fontweight', 'bold');
  grid;
else
  % version for MATLAB
  plot(x_t,Bt,'ok', x_g,Bi,'r',x_g,B1,'k',x_g,B2,'k');
  xlabel('x (m)');
  ylabel('B (T)');
  grid;
end

 

 
