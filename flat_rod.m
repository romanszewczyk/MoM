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
% Demonstration script for magnetostatic modelling of flat rod.
% Script uses the method of moments.
%
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% R. Szewczyk "Magnetostatic modelling of thin layers using the method of moments 
%              and its implementation in Octave/Matlab".
%

t=ver;
if strcmp(t(1).Name,'Octave')==1
   page_screen_output(0);
   page_output_immediately(1);
end % if Octave, than configure screen for Octave


L=0.25;             % (m) the length of thin rod

N=30;               % number of cells in the rod

dL = L./N;

mi=50;              % relative magnetic permeability of material
mi0=4.*pi.*1e-7;    % magnetic constant

Hext0=1000;         % (A/m) external magnetiing field in the direction of rod's axis

% Flat rod (ribbon)
g=5.*10e-5;     % thickness 50 um
w=1e-2;         % width 1 cm



Hext = ones(N,1).*Hext0.*(mi-1);
A = zeros(N,N);

for k=1:N
  for i=1:N

    fc0=@(ry) (2.*(dL.*(i-k+0.5)).^2-ry.^2)./((dL.*(i-k+0.5)).^2+ry.^2).^(5/2);
    fc1=@(ry) (2.*(dL.*(i-1-k+0.5)).^2-ry.^2)./((dL.*(i-1-k+0.5)).^2+ry.^2).^(5/2);
 
    A(i,k)=(mi-1).*g.*dL./(2.*pi).*(sign(i-k+0.5).*quad(fc0,0,w/2)-sign(i-1-k+0.5).*quad(fc1,0,w/2));
                  % numerical integration
    if k==i
       A(i,k)=A(i,k)+1;   
    end
  end
end

Mt=A \ Hext;

Bt=Hext0.*mi0+Mt.*mi0;                % (T) Flux density in the rod


x_g=0:L./300:L;
x_t=(1:N).*dL-0.5.*dL;
B1=ones(size(x_g)).*Hext0.*mi0.*mi;
B2=ones(size(x_g)).*Hext0.*mi0;       % interpolation of results

Bi=interp1(x_t,Bt,x_g,'nearest','extrap');

plot(x_t,Bt,'ok','linewidth',1,'markersize', 7, x_g,Bi,'r','linewidth',2,x_g,B1,'k','linewidth',2,x_g,B2,'k','linewidth',2);
set(gca, 'fontsize', 12);
xlabel('x (m)','fontweight', 'bold');
ylabel('B(T)','fontweight', 'bold');
grid;

