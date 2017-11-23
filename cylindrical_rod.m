% MIT License
%
% Copyright (c) 2017 Roman Szewczyk
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the 'Software'), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
% 
%
% DESCRIPTION:
% Demonstration script for magnetostatic modelling of cylindrical rod.
% Script uses the method of moments.
%
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% R. Szewczyk 'Magnetostatic modelling of thin layers using the method of moments 
%              and its implementation in Octave/Matlab'.
%

L=0.25;             % (m) the length of cylindrical rod

R=0.0022;           % (m) rod's radius

s=pi.*R.^2;         % (m^2) rod's crossection

N=25;               % number of cells in the rod

dL = L./N;

mi=30;             % relative magnetic permeability of material
mi0=4.*pi.*1e-7;   % magnetic constant

Hext0=100;         % (A/m) external magnetiing field in the direction of rod's axis



a = zeros(N,N);                   % matrix for a coefficents
Hext = ones(N,1).*Hext0.*(mi-1);  % matrix for external field Hext

g=@(iw,kw,Rw,dLw) Rw.^2.*sign(iw-kw+0.5)./(Rw.^2+(dLw.*(iw-kw+0.5)).^2).^(3/2);
                % anonymous function to simplify calculations
                
for k=1:N       % loop of equations
  for i=1:N     % loop of elements
        
    a(i,k)=(mi-1).*dL./2.*(g(i,k,R,dL)-g(i-1,k,R,dL));
                            % calculate a coefficient
     if k==i
       a(i,k)=a(i,k)+1;     % add 1 for Mk in its equation
    end
  endfor        % end of the elements's loop
endfor          % end of the equation's loop

M=a \ Hext;     % solve set of linear equations

B=Hext0.*mi0+M.*mi0;  % calculate flux density B


Bantder=B;    % (T) flux density B calculated with antiderivative


a = zeros(N,N);                   % matrix for a coefficents
Hext = ones(N,1).*Hext0.*(mi-1);  % matrix for external field Hext

for k=1:N       % loop of equations
  for i=1:N     % loop of elements

    gc=@(ry) ry.*(2.*(dL.*(i-k+0.5)).^2-ry.^2)./((dL.*(i-k+0.5)).^2+ry.^2).^(5/2);
    gc1=@(ry) ry.*(2.*(dL.*(i-1-k+0.5)).^2-ry.^2)./((dL.*(i-1-k+0.5)).^2+ry.^2).^(5/2);
                % two ananymous functions for numerical integration
                
    a(i,k)=(mi-1).*dL./2.*(sign(i-k+0.5).*quad(gc,0,R)-sign(i-1-k+0.5).*quad(gc1,0,R));
                % calculate a coefficient
    if k==i
       a(i,k)=a(i,k)+1;   % add 1 for Mk in its equation
    end
  endfor        % end of the elements's loop
endfor          % end of the equation's loop

M=a \ Hext;     % solve set of linear equations

B=Hext0.*mi0+M.*mi0;  % calculate flux density B

Bnum=B;         % (T) flux density B calculated with numerical integration

x_g=0:L./300:L;
x_t=(1:N).*dL-0.5.*dL;
B1=ones(size(x_g)).*Hext0.*mi0.*mi;
B2=ones(size(x_g)).*Hext0.*mi0;

Bi=interp1(x_t,Bnum,x_g,'cubic','extrap');

plot(x_t,Bantder,'+k','linewidth',1,'markersize',10,x_t,Bnum,'xk','linewidth',1,'markersize', 7, x_g,Bi,'g','linewidth',2,x_g,B1,'k','linewidth',2,x_g,B2,'k','linewidth',2);
set(gca, 'fontsize', 12);
xlabel('x (m)','fontweight', 'bold');
ylabel('B(T)','fontweight', 'bold');
grid;

