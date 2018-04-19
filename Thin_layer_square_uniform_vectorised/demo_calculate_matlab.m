%
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
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% R. Szewczyk "Magnetostatic modelling of thin layers using the method of moments 
%              and its implementation in Octave/Matlab", Springer, 2018. 
%              DOI 10.1007/978-3-319-77985-0, ISBN 978-3-319-77984-3
%

clear all;
diary on;

mi0=4.*pi.*1e-7;

clc;
fprintf('\n\n\n');
 
fprintf('\n *************************');
 
fprintf('\n *** MoM 2D vectorized ***');
 
fprintf('\n *************************\n');
 
fprintf('\n Initialization... ');

N = 25;        % number of cells on each dimmension

dL = 1e-2./(N-6);         % sizeof the cell

g = 1e-6;         % layer thickness

fprintf('\nN= %i , mi= %i, g= %f ',N,1000,g);

Hxext = 1000;      % external magnetizing field - x axis
Hyext = 0;        % external magnetizing field - y axis

mi=ones(N,N);   % relative permeability matrix

fprintf('\nok.');
% --------------------------------
% --- Ring on the permeability ---
% --------------------------------
 
fprintf('\n Generation of permeability matrix. Number of points = %i ...',N.^2);

% Shape of the layer

for ix=1:N
  for iy=1:N
    if (ix>2) && (ix<N-2) && (iy>2) && (iy<N-2)
      mi(iy,ix)=1000;
    end
      
    if (ix==10) && (iy==10)
      mi(iy,ix)=1;
    end 
  end
end


 
fprintf('\n Done.');

% -----------------------------------------------------------------
% --- Matrix of linear equations and external magnetizing field ---
% -----------------------------------------------------------------
 
fprintf('\n Generation of matrixes for linear equations and external magnetizing field. ');

aM = zeros(2.*N.^2,2.*N.^2);

aH = zeros(2.*N.^2,1);

 
fprintf('\n Number of points to generate: %i ...',size(aM,1).*size(aM,2));

c = g.*dL.^2./(4.*pi);


CTablesReady=0;

if exist('Cf.mat')==2
  load('Cf.mat');
   
  fprintf('\n C matrixes loaded'); 
  if (N==Nf) && (abs(dLf-dL)<1e-9)
      CTablesReady=1;
       
      fprintf('\n C matrixes ok.'); 
  else
   
  fprintf('\n C matrixes cleared'); 
  end
end

if ~CTablesReady

    
   fprintf('\n Calculate C matrixes...'); 
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

 
fprintf('\n Done.'); 



dp=0;
tic;

 
fprintf('\n Calculate mi matrix....'); 

x=1:2.*N.^2;
y=1:2.*N.^2;
[xx,yy]=meshgrid(x,y);

mixx=mod(floor((yy-1)./2),N)+1;
miyy=floor((floor((yy-1)./2))./N)+1;

TAmi=mi(mixx+(miyy-1).*N);  % table of mi - acces in vectorial form

 
fprintf('\n Done.'); 

 
fprintf('\n Calculate aM matrix...'); 

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

 
fprintf('\n Done.'); 

% ------------------------------------
% --- Prepare H vector via mi mask ---
% ------------------------------------

 
fprintf('\n Calculate aH matrix...'); 

aH(1:2:end)=Hxext;
aH(2:2:end)=Hyext;

mix=mod(floor((y'-1)./2),N)+1;
miy=floor((floor((y'-1)./2))./N)+1;

aH=aH.*mi(mix+(miy-1).*N);  % table of mi - acces in vectorial form;


 
fprintf('\n Done.');

% -----------------------------
% --- Solve linear eqations ---
% -----------------------------
 
fprintf('\n Solving the linear equations... ');

M = aM \ aH;

 
fprintf('\n Done. \n');
% -----------------------
% --- Present results ---
% -----------------------

Mx=M(1:2:end);
My=M(2:2:end);

Mx=reshape(Mx,N,N);
My=reshape(My,N,N);

Bx=Mx'.*mi0+Hxext.*mi0;
By=My'.*mi0+Hyext.*mi0;   % transposition of matrixes

save -v6 Results.mat Bx By Mx My mi N dL g Hxext Hyext

 

fprintf('\n *** All done. *** \n\n\n');

diary off;

