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
% Demonstration script for magnetostatic modelling of thin layer with 2D uniform mesh.
% Script calculates 2D flux density distribution for different thickneses of the layer
% Script uses the method of moments.
%
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% R. Szewczyk "Magnetostatic modelling of thin layers using the method of moments 
%              and its implementation in Octave/Matlab", Springer, 2018. 
%              DOI 10.1007/978-3-319-77985-0, ISBN 978-3-319-77984-3
%
%

clear all;

t_ver=ver;
if strcmp(t_ver(1).Name,'Octave')==1
   page_screen_output(0);
   page_output_immediately(1);
end

diary on;

mi0=4.*pi.*1e-7;

fprintf('\n\n Initialization... ');

Tk=[1e-6 5e-6 1e-5 3e-5 5e-5 1e-4 2e-4 5e-4 1e-3 3e-3]; % the list of thickneses for calculations

N = 50;                         % number of cells on each dimmension

Bx=zeros(N,N,numel(Tk));
By=zeros(N,N,numel(Tk));
Tmi=zeros(N,N,numel(Tk));

for iTk=1:numel(Tk)             % the main loop of the script

dL = 1e-2./(N-6);               % sizeof the cell

g = Tk(iTk);                    % set layer thickness

fprintf('N= %i , g= %f ',N,g);

Hxext = 1000;                   % external magnetizing field - x axis
Hyext = 0;                      % external magnetizing field - y axis


mi=ones(N,N);                   % relative permeability matrix

fprintf('ok. \n');

% ----------------------------------
% --- Prepare 2D mi distribution ---
% ----------------------------------
fprintf('\n ');
 
fprintf(' Generation of permeability matrix. Number of points = %i ...',N.^2);


center_x=(N./2+0.5).*dL;
center_y=(N./2+0.5).*dL;

r_o=N./2.*dL.*0.92;
r_i=N./2.*dL.*0.72;

for ix=1:N
  for iy=1:N
    if (ix>2) && (ix<N-2) && (iy>2) && (iy<N-2)
      mi(iy,ix)=1000;
      end
      
    if (ix>28) && (ix<30) && (iy>28) && (iy<30)
      mi(iy,ix)=1;
      end 
    end
end

fprintf('\n ');
 
fprintf(' ok.\n');

% -----------------------------------------------------------------
% --- Matrix of linear equations and external magnetizing field ---
% -----------------------------------------------------------------
fprintf('\n ');
 
fprintf(' Generation of matrixes for linear equations and external magnetizing field. ');

aM = zeros(2.*N.^2,2.*N.^2);

aH = zeros(2.*N.^2,1);

fprintf('\n                     Number of points to generate: %i ...',size(aM,1).*size(aM,2));

c = g.*dL.^2./(4.*pi);

if exist('Cf.mat')==2

   load('Cf.mat');
   fprintf('\n C matrixes loaded'); 
   
   if N~=Nf
      TCxx=zeros(2.*N,2.*N).*NaN;
      TCxy=zeros(2.*N,2.*N).*NaN;
      TCyx=zeros(2.*N,2.*N).*NaN;
      TCyy=zeros(2.*N,2.*N).*NaN;         % matrixes for speed up 
      fprintf('\n C matrixes cleared. N<>Nf ');
      end
      
   if abs(dLf-dL)>1e-9
      TCxx=zeros(2.*N,2.*N).*NaN;
      TCxy=zeros(2.*N,2.*N).*NaN;
      TCyx=zeros(2.*N,2.*N).*NaN;
      TCyy=zeros(2.*N,2.*N).*NaN;        % matrixes for speed up 
      fprintf('\n C matrixes cleared. dL<>dLf');
      end
      
else

   TCxx=zeros(2.*N,2.*N).*NaN;
   TCxy=zeros(2.*N,2.*N).*NaN;
   TCyx=zeros(2.*N,2.*N).*NaN;
   TCyy=zeros(2.*N,2.*N).*NaN;          % matrixes for speed up 
   fprintf('\n C matrixes empty');
end

dp=0;
tic;

for Ly=1:2.*N.^2

    for Lx=1:2.*N.^2
       
       [ix,iy,kx,ky]= L2ik(Lx,Ly,N);
     
   
       if mod(Ly,2)==1 		  % => Mx in equations (equation 1) 	=> cx?
         if mod(Lx,2)==1 		% => Mx in calculation of aM(iy,ix) => cxx
            
            if isnan(TCxx(1+N+ix-kx,1+N+iy-ky))
               Cxx_=Cxx(ix,iy,kx,ky,dL);
               TCxx(1+N+ix-kx,1+N+iy-ky)=Cxx_;
               else
               Cxx_=TCxx(1+N+ix-kx,1+N+iy-ky);
               end

            if isnan(TCxx(1+N+ix-1-kx,1+N+iy-ky))
               Cxx1_=Cxx(ix-1,iy,kx,ky,dL);
               TCxx(1+N+ix-1-kx,1+N+iy-ky)=Cxx1_;
               else
               Cxx1_=TCxx(1+N+ix-1-kx,1+N+iy-ky);
               end   
     
            aM(Ly,Lx)=c.*(mi(kx,ky)-1).*(Cxx_- Cxx1_);

         end
         
         if mod(Lx,2)==0 		% => My in calculation of aM(iy,ix) => cxy
            
            if isnan(TCxy(1+N+ix-kx,1+N+iy-ky))
               Cxy_=Cxy(ix,iy,kx,ky,dL);
               TCxy(1+N+ix-kx,1+N+iy-ky)=Cxy_;
               else
               Cxy_=TCxy(1+N+ix-kx,1+N+iy-ky);
               end

            if isnan(TCxy(1+N+ix-kx,1+N+iy-1-ky))
               Cxy1_=Cxy(ix,iy-1,kx,ky,dL);
               TCxy(1+N+ix-kx,1+N+iy-1-ky)=Cxy1_;
               else
               Cxy1_=TCxy(1+N+ix-kx,1+N+iy-1-ky);
               end  
            
            aM(Ly,Lx)=c.*(mi(kx,ky)-1).*(Cxy_- Cxy1_);
            
            end
       end
      
      
       if mod(Ly,2)==0  		  % => My in equations (equation 2)   => cy?
          if mod(Lx,2)==1 		% => Mx in calculation of aM(iy,ix) => cyx 

            if isnan(TCyx(1+N+ix-kx,1+N+iy-ky))
               Cyx_=Cyx(ix,iy,kx,ky,dL);
               TCyx(1+N+ix-kx,1+N+iy-ky)=Cyx_;
               else
               Cyx_=TCyx(1+N+ix-kx,1+N+iy-ky);
               end

            if isnan(TCyx(1+N+ix-1-kx,1+N+iy-ky))
               Cyx1_=Cyx(ix-1,iy,kx,ky,dL);
               TCyx(1+N+ix-1-kx,1+N+iy-ky)=Cyx1_;
               else
               Cyx1_=TCyx(1+N+ix-1-kx,1+N+iy-ky);
               end
            
            aM(Ly,Lx)=c.*(mi(kx,ky)-1).*(Cyx_- Cyx1_);
          
          end
            
          if mod(Lx,2)==0 		% => My in calculation of aM(iy,ix)     => cyy

            if isnan(TCyy(1+N+ix-kx,1+N+iy-ky))
               Cyy_=Cyy(ix,iy,kx,ky,dL);
               TCyy(1+N+ix-kx,1+N+iy-ky)=Cyy_;
               else
               Cyy_=TCyy(1+N+ix-kx,1+N+iy-ky);
               end

            if isnan(TCyy(1+N+ix-kx,1+N+iy-1-ky))
               Cyy1_=Cyy(ix,iy-1,kx,ky,dL);
               TCyy(1+N+ix-kx,1+N+iy-1-ky)=Cyy1_;
               else
               Cyy1_=TCyy(1+N+ix-kx,1+N+iy-1-ky);
               end
               
            aM(Ly,Lx)=c.*(mi(kx,ky)-1).*(Cyy_- Cyy1_);

            end 
       end      
      
    end

aM(Ly,Ly)= aM(Ly,Ly)+1;

if mod(Ly,2)==1  
   aH(Ly,1)=Hxext.*(mi(kx,ky)-1);
   else aH(Ly,1)=Hyext.*(mi(kx,ky)-1);		% prepare aH matrix
   end 

if floor(Ly./(2.*N.^2).*100)>dp
   dp=floor(Ly./(2.*N.^2).*100);
   czas=toc;
   fprintf('\n done: %i proc., working: %3.2f h., to end: %3.2f h.',dp, czas./3600, czas./(3600.*dp)*(100-dp));
   end
   
end


fprintf('\n ');
 
fprintf(' ok.\n');

% -----------------------------
% --- Solve linear eqations ---
% -----------------------------

fprintf('\n ');
 
fprintf(' Solving the linear equations... ');

M = aM \ aH;

fprintf('\n ');
 
fprintf(' ok. \n');

% -----------------------
% --- Present results ---
% -----------------------

Mx=zeros(N,N);
My=zeros(N,N);

for Ly=1:2.*N.^2

[ix,iy,kx,ky]= L2ik(1,Ly,N);

if mod(Ly,2)==1
   Mx(ky,kx)=M(Ly,1);
   else
   My(ky,kx)=M(Ly,1);
   
   end

end 

Bx(:,:,iTk)=Mx.*mi0+Hxext.*mi0;
By(:,:,iTk)=My.*mi0+Hyext.*mi0;
Tmi(:,:,iTk)=mi;

Nf=N;
dLf=dL;

save -v6 Cf.mat TCxx TCyy TCxy TCyx Nf dLf

save -v6 res_calc.mat Bx By Tmi iTk Tk mi Mx My N dL g Hxext Hyext

fprintf('\n ');
 
fprintf(' *** Cycle done. *** \n\n\n');

end

fprintf('\n ');
 
fprintf(' *** All done. *** \n\n\n');

diary off;

