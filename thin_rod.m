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
% Demonstration script for magnetostatic modelling of thin rod.
% Script uses the method of moments.
%
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% R. Szewczyk "Magnetostatic modelling of thin layers using the method of moments 
%              and its implementation in Octave/Matlab".
%


L=0.2;              % (m) the length of thin rod

s=pi.*(0.002).^2;   % (m^2) rod's crossection

N=50;               % number of cells in the rod

dL = L./N;

mi=30;              % relative magnetic permeability of material
mi0=4.*pi.*1e-7;    % magnetic constant

Hext0=1;


a = zeros(N,N);                   % matrix for a coefficents
Hext = ones(N,1).*Hext0.*(mi-1);  % matrix for external field Hext

for k=1:N       % loop of equations
  for i=1:N     % loop of elements
    
    b1=(mi-1).*s./(2.*pi.*dL.^2);
    b2=i-k+0.5;
    b3=(i-1)-k+0.5;               % calculate parameters to simplify
    
    a(i,k)=b1.*(sign(b2)./(abs(b2)).^3-sign(b3)./(abs(b3)).^3);
                                  % calculate a coefficient
    if k==i
       a(i,k)=a(i,k)+1;           % add 1 for Mk in its equation
    endif
    
  endfor        % end of the elements's loop
endfor          % end of the equation's loop

M=a \ Hext;     % solve set of linear equations

B=Hext0.*mi0+M.*mi0;  % calculate flux density B

x_t=(1:N).*dL-0.5.*dL;

x_g=0:L./300:L;
B1=ones(size(x_g)).*Hext0.*mi0.*mi;
B2=ones(size(x_g)).*Hext0.*mi0;

Bi=interp1(x_t,B,x_g,'nearest','extrap');

plot(x_t,B,'or','linewidth',2,x_g,Bi,'r','linewidth',2,x_g,B1,'k','linewidth',2,x_g,B2,'k','linewidth',2);
set(gca, 'fontsize', 12);
xlabel('x (m)','fontweight', 'bold');
ylabel('B(T)','fontweight', 'bold');
grid;

