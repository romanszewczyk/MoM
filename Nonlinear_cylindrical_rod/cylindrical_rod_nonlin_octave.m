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
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% R. Szewczyk 'Magnetostatic modelling of thin layers using the method of moments 
%              and its implementation in Octave/Matlab', Springer, 2018. 
%              DOI 10.1007/978-3-319-77985-0, ISBN 978-3-319-77984-3
%

page_screen_output(0);
page_output_immediately(1);

L=0.25;               % (m) rod length

R=0.0015;     % rdius

s=pi.*R.^2;   % (m^2) rod crosssection

N=25;               % number of points

dL = L./N;

mi=1100;             % materials relative permeability 
mi0=4.*pi.*1e-7;    % magnetic constant

Hext0=5000;


% --- Calculation of linear model from antiderivative dependence ----


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

Blin_antder=Hext0.*mi0+M.*mi0;  % calculate flux density B


% --- Calculation of nonlinear model from numerical integration ----

a = zeros(N,N);       % matrix for a coefficents without adding ones
b = eye (N);          % diagonal matrix


for k=1:N       % loop of equations
  for i=1:N     % loop of elements

    gc=@(ry) ry.*(2.*(dL.*(i-k+0.5)).^2-ry.^2)./((dL.*(i-k+0.5)).^2+ry.^2).^(5/2);
    gc1=@(ry) ry.*(2.*(dL.*(i-1-k+0.5)).^2-ry.^2)./((dL.*(i-1-k+0.5)).^2+ry.^2).^(5/2);
                % two ananymous functions for numerical integration
                
    a(i,k)=dL./2.*(sign(i-k+0.5).*quad(gc,0,R)-sign(i-1-k+0.5).*quad(gc1,0,R));
                % calculate a coefficient
    %if k==i
    %   a(i,k)=a(i,k)+1;   % add 1 for Mk in its equation
    %end
  endfor        % end of the elements's loop
endfor          % end of the equation's loop


%M=a \ Hext;     % solve set of linear equations
%B=Hext0.*mi0+M.*mi0;  % calculate flux density B
%Bnum=B;


% - JA model - 

load('JAresR.mat');

mi_M=@(Mnl) interp1(M_JA,mi_JA,Mnl,'nearest','extrap');

Ftarget=@(M_) ((a.*(mi_M(M_)-1)+b)*M_-(mi_M(M_)-1).*Hext0);

o=optimset('MaxIter',10000,'MaxFunEvals',10000,'TolX',1e-10,'TolFun',1e-10);
[M_JA_wyn, fval, exitflag]=fsolve(Ftarget,M,o);

%[M_JA_wyn, fval, exitflag]=fsolve(Ftarget,M)

BJA=Hext0.*mi0+M_JA_wyn.*mi0;  % calculate flux density B

% --- Draw results ---
figure (1);

plot(H_JA,M_JA.*mi0,'-ob',H_JA,H_JA.*mi.*mi0,'-xr');
set(gca, 'fontsize', 14);
title('Jiles-Atherton B(H) dependence');
xlabel('H (A/m)','fontweight', 'bold');
ylabel('B (T)','fontweight', 'bold');
grid;

figure (2);
plot(H_JA,mi_JA,'-o',H_JA,(H_JA.*0+1).*mi,'-xr');
set(gca, 'fontsize', 14);
title('Jiles-Atherton \mu(H) dependence');
xlabel('\mu','fontweight', 'bold');
ylabel('M (A/m)','fontweight', 'bold');
grid;

figure (3);
plot(M_JA,mi_JA,'-o',H_JA(1:200).*mi,(H_JA(1:200).*0+1).*mi,'-xr');
set(gca, 'fontsize', 14);
title('Jiles-Atherton \mu(M) dependence');
xlabel('M (A/m)','fontweight', 'bold');
ylabel('\mu','fontweight', 'bold');
grid;


figure (4);

x_g=0:L./300:L;   % high density length scale 
x_t=(1:N).*dL-0.5.*dL; % results length scale

B_mimax=ones(size(x_g)).*Hext0.*mi0.*mi;  % maximal from mi 
B_mi0=ones(size(x_g)).*Hext0.*mi0;        % calculated for mi0

plot(x_g,B_mimax,'k','linewidth',2,x_g,B_mi0,'k','linewidth',2, ... 
     x_t,BJA,'ob','linewidth',2,x_g,interp1(x_t,BJA,x_g,'cubic','extrap'),'-b','linewidth',2,
     x_t,Blin_antder,'xr','linewidth',1,'markersize',10,x_g,interp1(x_t,Blin_antder,x_g,'cubic','extrap'),'-r','linewidth',1,'markersize', 7);


set(gca, 'fontsize', 12);
xlabel('x(m)','fontweight', 'bold');
ylabel('B(T)','fontweight', 'bold');
grid;

