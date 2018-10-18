function [Bx, By]= MoM2Dvect (mi, dL, g, Hxext, Hyext)
%
% MIT License
%
% Copyright (c) 2018 Roman Szewczyk
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
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% DESCRIPTION:
% Function for modelling the magnetic field distribution in the racetrack-shaped
% core of fluxgate sensors.
%
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION:
% R. Szewczyk "Shape optimization of the race-track core for fluxgate sensor"
% (under review process)
%
% For detailed explanation of the Method of moments please refer to:
% R. Szewczyk "Magnetostatic modelling of thin layers using the method of moments 
%              and its implementation in Octave/Matlab", Springer, 2018. 
%              DOI 10.1007/978-3-319-77985-0, ISBN 978-3-319-77984-3
%
%

% dL = 1e-2./(N-6);  % sizeof the cell
% g = 3e-5;          % layer thickness
% Hxext = 1000;      % external magnetizing field - x axis
% Hyext = 0;         % external magnetizing field - y axis
% mi=ones(N,N);      % relative permeability matrix

ShowTime;
fprintf(' Initialization... ');

if (size(mi,1)~=size(mi,2))
   ShowEmpty;
   fprintf(' mi has to be the square matrix \n\n');
   Bx=[];
   By=[];
   return;
end

mi0=4.*pi.*1e-7;    % magnetic constant

N=size(mi,1);       % number of cells on each dimmension

ShowEmpty;
fprintf(' mi matrix: %i x %i, dL= %f (m), g= %3.1f (um)',N,N,dL,g.*1e6);

ShowEmpty;
fprintf(' Hxext= %4.2f, Hyext= %4.2f ',Hxext,Hyext);

ShowEmpty;
fprintf(' Done.');

% -----------------------------------------------------------------
% --- Matrix of linear equations and external magnetizing field ---
% -----------------------------------------------------------------
ShowTime;
fprintf(' Generation of matrixes for linear equations and external magnetizing field. ');

aM = zeros(2.*N.^2,2.*N.^2);

aH = zeros(2.*N.^2,1);

ShowEmpty;
fprintf(' Number of points to generate: %i ...',size(aM,1).*size(aM,2));

c = g.*dL.^2./(4.*pi);


CTablesReady=0;

if exist('Cf.mat')==2
  load('Cf.mat');
  ShowEmpty;
  fprintf(' C matrixes loaded'); 
  if (N==Nf) && (abs(dLf-dL)<1e-9)
      CTablesReady=1;
      ShowEmpty;
      fprintf(' C matrixes ok.'); 
  else
  ShowEmpty;
  fprintf(' C matrixes cleared'); 
  end
end

if ~CTablesReady

   ShowEmpty;
   fprintf(' Calculate C matrixes...'); 
   x=linspace(-1.*N,N,2.*N+1);
   y=linspace(-1.*N,N,2.*N+1);

   [xx,yy]=meshgrid(x,y);

   TCxx=arrayfun(@(a,b) Cxx(a,b,dL),yy,xx);
   TCxy=arrayfun(@(a,b) Cxy(a,b,dL),yy,xx);
   TCyx=arrayfun(@(a,b) Cyx(a,b,dL),yy,xx);
   TCyy=arrayfun(@(a,b) Cyy(a,b,dL),yy,xx);

   Nf=N;
   dLf=dL;
   save -v6 Cf.mat TCxx TCyy TCxy TCyx Nf dLf  
   
end

ShowTime;
fprintf(' Done.'); 

dp=0;
tic;

ShowTime;
fprintf(' Calculate mi matrix....'); 

x=1:2.*N.^2;
y=1:2.*N.^2;
[xx,yy]=meshgrid(x,y);

mixx=mod(floor((yy-1)./2),N)+1;
miyy=floor((floor((yy-1)./2))./N)+1;

TAmi=mi(mixx+(miyy-1).*N);  % table of mi - acces in vectorial form

ShowTime;
fprintf(' Done.'); 

ShowTime;
fprintf(' Calculate aM matrix...'); 

maskCxx=xx.*0;
maskCxx(1:2:end,1:2:end)=1;

maskCxy=xx.*0;
maskCxy(1:2:end,2:2:end)=1;

maskCyx=xx.*0;
maskCyx(2:2:end,1:2:end)=1;

maskCyy=xx.*0;
maskCyy(2:2:end,2:2:end)=1;

aMCxx=aM;
aMCxy=aM;
aMCyx=aM;
aMCyy=aM;

Lx=xx;
Ly=yy;

aMCxx = TCxx(1+N+mod(floor((Lx-1)./2),N)-mod(floor((Ly-1)./2),N)+ ...
        (2.*N+1).*(((1+N+floor((floor((Lx-1)./2))./N)-floor((floor((Ly-1)./2))./N)))-1))- ...
        TCxx(N+mod(floor((Lx-1)./2),N)-mod(floor((Ly-1)./2),N)+ ... 
        (2.*N+1).*(((1+N+floor((floor((Lx-1)./2))./N)-floor((floor((Ly-1)./2))./N)))-1));

aMCxy = TCxy(1+N+mod(floor((Lx-1)./2),N)-mod(floor((Ly-1)./2),N)+ ... 
        (2.*N+1).*(((1+N+floor((floor((Lx-1)./2))./N)-floor((floor((Ly-1)./2))./N)))-1))- ...
        TCxy(1+N+mod(floor((Lx-1)./2),N)-mod(floor((Ly-1)./2),N)+ ...
        (2.*N+1).*(((N+floor((floor((Lx-1)./2))./N)-floor((floor((Ly-1)./2))./N)))-1));
           
aMCyx = TCyx(1+N+mod(floor((Lx-1)./2),N)-mod(floor((Ly-1)./2),N)+ ...
        (2.*N+1).*(((1+N+floor((floor((Lx-1)./2))./N)-floor((floor((Ly-1)./2))./N)))-1))- ...
        TCyx(N+mod(floor((Lx-1)./2),N)-mod(floor((Ly-1)./2),N)+ ...
        (2.*N+1).*(((1+N+floor((floor((Lx-1)./2))./N)-floor((floor((Ly-1)./2))./N)))-1));
                       
aMCyy = TCyy(1+N+mod(floor((Lx-1)./2),N)-mod(floor((Ly-1)./2),N)+ ...
        (2.*N+1).*(((1+N+floor((floor((Lx-1)./2))./N)-floor((floor((Ly-1)./2))./N)))-1))- ...
        TCyy(1+N+mod(floor((Lx-1)./2),N)-mod(floor((Ly-1)./2),N)+ ...
        (2.*N+1).*(((N+floor((floor((Lx-1)./2))./N)-floor((floor((Ly-1)./2))./N)))-1));


aM=maskCxx.*aMCxx+maskCxy.*aMCxy+maskCyx.*aMCyx+maskCyy.*aMCyy;

aM= TAmi.*c.*aM+eye(2.*N.^2);

ShowTime;
fprintf(' Done.'); 

% ------------------------------------
% --- Prepare H vector via mi mask ---
% ------------------------------------

ShowTime;
fprintf(' Calculate aH matrix...'); 

aH(1:2:end)=Hxext;
aH(2:2:end)=Hyext;

mix=mod(floor((y'-1)./2),N)+1;
miy=floor((floor((y'-1)./2))./N)+1;

aH=aH.*mi(mix+(miy-1).*N);  % table of mi - acces in vectorial form;


ShowTime;
fprintf(' Done.');

% -----------------------------
% --- Solve linear eqations ---
% -----------------------------
ShowTime;
fprintf(' Solving the linear equations... ');

M = aM \ aH;

ShowTime;
fprintf(' Done.');
% -----------------------
% --- Present results ---
% -----------------------

Mx=M(1:2:end);
My=M(2:2:end);

Mx=reshape(Mx,N,N);
My=reshape(My,N,N);

Bx=Mx.*mi0+Hxext.*mi0;
By=My.*mi0+Hyext.*mi0;   % TRANSPOZYCJA MACIERZY ???

ShowTime;
fprintf(' *** Calculations done sucessfully. *** \n\n');

