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
% Demonstration script for magnetostatic modelling of cylindrical rod.
% Script uses the method of moments.
%
% AUTHOR: Roman Szewczyk, rszewczyk@onet.pl
%
% RELATED PUBLICATION(S):
% R. Szewczyk "Magnetostatic modelling of thin layers using the method of moments 
%              and its implementation in Octave/Matlab", Springer, 2018. 
%              DOI 10.1007/978-3-319-77985-0, ISBN 978-3-319-77984-3
%

clear all
clc

[err, msg] = unlink ('diary');
clear err msg
diary on

page_screen_output(0);
page_output_immediately(1);

PresentPlots=true;

tic

fprintf('*** Magnetostatic calculations for thin layer described by the 2D adaptive mesh ***');

mesh_dir='romb_02/';

fprintf('\n\nLoading mesh from directory: %s...', mesh_dir);

[nodes, triangles, borders_external]=LoadMesh(mesh_dir);

fprintf('done.\nnodes: %i, triangles: %i, external border lines: %i \n\n',size(nodes,1),size(triangles,1),size(borders_external,1));

clear mesh_dir

% --- Formats ---

% nodes = [nx ny] 
% where
% nx - x coordinate of the node
% ny - y coordinate of the node

% triangles = [p1 p2 p3 bx by s mi Mx My]
% where
% p1, p2, p3 - vertices of the triangle
% bx, by - x,y coordinates of the barycentre
% s - area of the triangle
% mi - relative permeability
% Mx, My - value of the magnetisation in the trangle-shaped cell

% borders_external = [p1 p2]
% where
% p1, p2 - points of the border

% borders = [p1 p2 t1 t2]
% where
% p1, p2 - points of the border
% t1, t2 - neighbour triangles connected by the border


fprintf('Calculations for triangles...');

% --- calculation for triangles ---

triangles_=zeros(size(triangles,1),9);
triangles_(:,1:3)=triangles;

% calculate baricentres
triangles_(:,4)=(nodes(triangles_(:,1),1)+nodes(triangles_(:,2),1)+nodes(triangles_(:,3),1))./3;
triangles_(:,5)=(nodes(triangles_(:,1),2)+nodes(triangles_(:,2),2)+nodes(triangles_(:,3),2))./3;

% calculate fields

triangles_(:,6)=TriangleField(nodes(triangles_(:,1),1), nodes(triangles_(:,1),2), nodes(triangles_(:,2,1)), nodes(triangles_(:,2),2), ...
                  nodes(triangles_(:,3),1), nodes(triangles_(:,3),2));

triangles=triangles_;
clear triangles_
  
fprintf('done.\n'); 
 
% --- calculation for borders ---
fprintf('Create table of borders...');

borders=zeros(size(triangles,1).*3,4);
b_count=1;

for i=1:size(triangles,1)

% find borders

p1=triangles(i,1);
p2=triangles(i,2);

pos1=sum(ismember(borders(1:b_count,1:2),[p1 p2],'rows'));
pos2=sum(ismember(borders(1:b_count,1:2),[p2 p1],'rows'));

if (pos1==0) && (pos2==0)
    borders(b_count,1:2)=[p1 p2];
    b_count=b_count+1;
end

p1=triangles(i,1);
p2=triangles(i,3);

pos1=sum(ismember(borders(1:b_count,1:2),[p1 p2],'rows'));
pos2=sum(ismember(borders(1:b_count,1:2),[p2 p1],'rows'));

if (pos1==0) && (pos2==0)
    borders(b_count,1:2)=[p1 p2];
    b_count=b_count+1;
end

p1=triangles(i,2);
p2=triangles(i,3);

pos1=sum(ismember(borders(1:b_count,1:2),[p1 p2],'rows'));
pos2=sum(ismember(borders(1:b_count,1:2),[p2 p1],'rows'));

if (pos1==0) && (pos2==0)
    borders(b_count,1:2)=[p1 p2];
    b_count=b_count+1;
end

end

borders=borders(1:b_count-1,:);

fprintf(' Found: %i borders. Done.\n',b_count-1); 

clear b_count p1 p2 pos1 pos2 i

% find neighbour domains

for i=1:size(borders,1)

    x1=nodes(borders(i,1),1);
    y1=nodes(borders(i,1),2);
    
    x2=nodes(borders(i,2),1);
    y2=nodes(borders(i,2),2);

    px=x1+(x2-x1)./2;
    py=y1+(y2-y1)./2; 

    wx=0.01.*(x2-x1);
    wy=0.01.*(y2-y1);
    
    wx_=wx.*cos(pi./2)-wy.*sin(pi./2);
    wy_=wx.*sin(pi./2)+wy.*cos(pi./2);
    
    borders(i,3)=FindTriangle(nodes, triangles, px+wx_, py+wy_);
 
    wx_=wx.*cos(-1.*pi./2)-wy.*sin(-1.*pi./2);
    wy_=wx.*sin(-1.*pi./2)+wy.*cos(-1.*pi./2);
    
    borders(i,4)=FindTriangle(nodes, triangles, px+wx_, py+wy_); 
    
end

clear x1 y1 x2 y2 px py wx wy wx_ wy_ i


% --- Draw the results of import of the mesh ---

if PresentPlots

  figure(1);
  hold on;

  plot(nodes(:,1),nodes(:,2),'ok');  % plot nodes

  for i=1:size(triangles,1)
      plot([nodes(triangles(i,1),1) nodes(triangles(i,2),1) nodes(triangles(i,3),1) nodes(triangles(i,1),1)], ...
        [nodes(triangles(i,1),2) nodes(triangles(i,2),2) nodes(triangles(i,3),2) nodes(triangles(i,1),2)],'k');
  end                               % draw triangles

  for i=1:size(borders,1)

    if ((borders(i,3)==0) || (borders(i,4)==0))
        plot([nodes(borders(i,1),1) nodes(borders(i,2),1)], ...
             [nodes(borders(i,1),2) nodes(borders(i,2),2)],'r-','LineWidth',4);
    end   
  end                               % draw external borders from M=0


  for i=1:size(borders_external,1)

  plot([nodes(borders_external(i,1),1) nodes(borders_external(i,2),1)], ...
       [nodes(borders_external(i,1),2) nodes(borders_external(i,2),2)],'g-','LineWidth',1);
  
  end                             % draw external borders from mesh


  for i=1:size(triangles,1)
    % plot(triangles(i,4),triangles(i,5),'*k');
  end                               % draw baricentres

  set(gca, "fontsize", 14);
  
  xlabel('x',"fontweight", "bold");
  ylabel('y',"fontweight", "bold");
  hold off;

end


% --- Prepare environent for calculations ---

x_min=min(nodes(:,1))-abs(min(nodes(:,1))).*0.065;
x_max=max(nodes(:,1))+abs(max(nodes(:,1))).*0.12;
y_min=min(nodes(:,2))-abs(min(nodes(:,2))).*0.065;
y_max=max(nodes(:,2))+abs(max(nodes(:,2))).*0.12;

Nx=20;        % number of cells in x direction
dL = (x_max-x_min)./Nx;         % sizeof the cell
Ny = round((y_max-y_min)./dL)+1;

mi0 = 4.*pi.*1e-7;
mi_mat = 1000;              % magnetic relative permeability of the layer

Hextx = 0;
Hexty = 100;


fprintf('\n\n Initialization... ');


fprintf('Nx = %i , Nx = %i , dL = %2.2e ',Nx,Ny,dL);

g = 7e-3;         % layer thickness

Hxext = 0;      % external magnetizing field - x axis
Hyext = 100;        % external magnetizing field - y axis

mi=ones(Nx,Ny);   % relative permeability matrix

[xx,yy]=meshgrid(1:Ny,1:Nx);

xx=x_min+(xx-1).*dL;
yy=y_min+(yy-1).*dL; 

for i=1:size(xx,2)
  for j=1:size(xx,1)
  if FindTriangle(nodes,triangles,xx(j,i),yy(j,i))>0
    mi(j,i)=mi_mat;
    else
    mi(j,i)=1;
    end
end
end

fprintf('ok. \n');

figure(2);
pcolor(xx,yy,mi);
colorbar;
title('mi');



fprintf('\n ');
fprintf(strftime ("%Y-%m-%d %H:%M:%S", localtime (time ()))); 
fprintf(' ok.\n');

% -----------------------------------------------------------------
% --- Matrix of linear equations and external magnetizing field ---
% -----------------------------------------------------------------
fprintf('\n ');
fprintf(strftime ("%Y-%m-%d %H:%M:%S", localtime (time ()))); 
fprintf(' Generation of matrixes for linear equations and external magnetizing field. ');

aM = zeros(2.*Nx.*Ny,2.*Nx.*Ny);

aH = zeros(2.*Nx.*Ny,1);

fprintf('\n                     Number of points to generate: %i ...',size(aM,1).*size(aM,2));

c = g.*dL.^2./(4.*pi);

if exist('Cf.mat')==2

   load('Cf.mat');
   fprintf('\n C matrixes loaded'); 
   
   if (Nx~=Nfx)
      TCxx=zeros(2.*Nx,2.*Ny).*NaN;
      TCxy=zeros(2.*Nx,2.*Ny).*NaN;
      TCyx=zeros(2.*Nx,2.*Ny).*NaN;
      TCyy=zeros(2.*Nx,2.*Ny).*NaN;        % matrixes for speed up 
      fprintf('\n C matrixes cleared. N<>Nf ');
      end

   if (Ny~=Nfy)
      TCxx=zeros(2.*Nx,2.*Ny).*NaN;
      TCxy=zeros(2.*Nx,2.*Ny).*NaN;
      TCyx=zeros(2.*Nx,2.*Ny).*NaN;
      TCyy=zeros(2.*Nx,2.*Ny).*NaN;        % matrixes for speed up 
      fprintf('\n C matrixes cleared. N<>Nf ');
      end
      
   if abs(dLf-dL)>1e-9
      TCxx=zeros(2.*Nx,2.*Ny).*NaN;
      TCxy=zeros(2.*Nx,2.*Ny).*NaN;
      TCyx=zeros(2.*Nx,2.*Ny).*NaN;
      TCyy=zeros(2.*Nx,2.*Ny).*NaN;        % matrixes for speed up 
      fprintf('\n C matrixes cleared. dL<>dLf');
      end
      
else

   TCxx=zeros(2.*Nx,2.*Ny).*NaN;
   TCxy=zeros(2.*Nx,2.*Ny).*NaN;
   TCyx=zeros(2.*Nx,2.*Ny).*NaN;
   TCyy=zeros(2.*Nx,2.*Ny).*NaN;        % matrixes for speed up 
   fprintf('\n C matrixes empty');
end

dp=0;
tic;

for Ly=1:2.*Nx.*Ny

    for Lx=1:2.*Nx.*Ny
       
       [ix,iy,kx,ky]= L2ik(Lx,Ly,Nx);
     
   
       if mod(Ly,2)==1 		  % => Mx w typie równania (równanie 1) 	=> cx?
         if mod(Lx,2)==1 		% => Mx w obliczeniu aM(iy,ix)         	=> cxx
            
            if isnan(TCxx(1+Nx+ix-kx,1+Ny+iy-ky))
               Cxx_=Cxx(ix,iy,kx,ky,dL);
               TCxx(1+Nx+ix-kx,1+Ny+iy-ky)=Cxx_;
               else
               Cxx_=TCxx(1+Nx+ix-kx,1+Ny+iy-ky);
               end

            if isnan(TCxx(1+Nx+ix-1-kx,1+Ny+iy-ky))
               Cxx1_=Cxx(ix-1,iy,kx,ky,dL);
               TCxx(1+Nx+ix-1-kx,1+Ny+iy-ky)=Cxx1_;
               else
               Cxx1_=TCxx(1+Nx+ix-1-kx,1+Ny+iy-ky);
               end   
     
            aM(Ly,Lx)=c.*(mi(kx,ky)-1).*(Cxx_- Cxx1_);

         end
         
         if mod(Lx,2)==0 		% => My w obliczeniu aM(iy,ix)         	=> cxy
            
            if isnan(TCxy(1+Nx+ix-kx,1+Ny+iy-ky))
               Cxy_=Cxy(ix,iy,kx,ky,dL);
               TCxy(1+Nx+ix-kx,1+Ny+iy-ky)=Cxy_;
               else
               Cxy_=TCxy(1+Nx+ix-kx,1+Ny+iy-ky);
               end

            if isnan(TCxy(1+Nx+ix-kx,1+Ny+iy-1-ky))
               Cxy1_=Cxy(ix,iy-1,kx,ky,dL);
               TCxy(1+Nx+ix-kx,1+Ny+iy-1-ky)=Cxy1_;
               else
               Cxy1_=TCxy(1+Nx+ix-kx,1+Ny+iy-1-ky);
               end  
            
            aM(Ly,Lx)=c.*(mi(kx,ky)-1).*(Cxy_- Cxy1_);
            
            end
       end
      
      
       if mod(Ly,2)==0  		% => My w typie równania (równanie 2)	    => cy?
          if mod(Lx,2)==1 		% => Mx w obliczeniu aM(iy,ix)   		    => cyx 

            if isnan(TCyx(1+Nx+ix-kx,1+Ny+iy-ky))
               Cyx_=Cyx(ix,iy,kx,ky,dL);
               TCyx(1+Nx+ix-kx,1+Ny+iy-ky)=Cyx_;
               else
               Cyx_=TCyx(1+Nx+ix-kx,1+Ny+iy-ky);
               end

            if isnan(TCyx(1+Nx+ix-1-kx,1+Ny+iy-ky))
               Cyx1_=Cyx(ix-1,iy,kx,ky,dL);
               TCyx(1+Nx+ix-1-kx,1+Ny+iy-ky)=Cyx1_;
               else
               Cyx1_=TCyx(1+Nx+ix-1-kx,1+Ny+iy-ky);
               end
            
            aM(Ly,Lx)=c.*(mi(kx,ky)-1).*(Cyx_- Cyx1_);
          
          end
            
          if mod(Lx,2)==0 		% => My w obliczeniu aM(iy,ix)          => cyy

            if isnan(TCyy(1+Nx+ix-kx,1+Ny+iy-ky))
               Cyy_=Cyy(ix,iy,kx,ky,dL);
               TCyy(1+Nx+ix-kx,1+Ny+iy-ky)=Cyy_;
               else
               Cyy_=TCyy(1+Nx+ix-kx,1+Ny+iy-ky);
               end

            if isnan(TCyy(1+Nx+ix-kx,1+Ny+iy-1-ky))
               Cyy1_=Cyy(ix,iy-1,kx,ky,dL);
               TCyy(1+Nx+ix-kx,1+Ny+iy-1-ky)=Cyy1_;
               else
               Cyy1_=TCyy(1+Nx+ix-kx,1+Ny+iy-1-ky);
               end
               
            aM(Ly,Lx)=c.*(mi(kx,ky)-1).*(Cyy_- Cyy1_);

            end 
       end      
      
    end

aM(Ly,Ly)= aM(Ly,Ly)+1;

if mod(Ly,2)==1  
   aH(Ly,1)=Hxext.*(mi(kx,ky)-1);
   else aH(Ly,1)=Hyext.*(mi(kx,ky)-1);		% od razu zrób macierz aH
   end 

if floor(Ly./(2.*Nx.*Ny).*100)>dp
   dp=floor(Ly./(2.*Nx.*Ny).*100);
   czas=toc;
   % fprintf('\n done: %i proc., working: %3.2f h., to end: %3.2f h.',dp, czas./3600, czas./(3600.*dp)*(100-dp));

   end

   
end


fprintf('\n ');
fprintf(strftime ("%Y-%m-%d %H:%M:%S", localtime (time ()))); 
fprintf(' ok.\n');
% -----------------------------
% --- Solve linear eqations ---
% -----------------------------
fprintf('\n ');
fprintf(strftime ("%Y-%m-%d %H:%M:%S", localtime (time ()))); 
fprintf(' Solving the linear equations... ');

M = aM \ aH;

fprintf('\n ');
fprintf(strftime ("%Y-%m-%d %H:%M:%S", localtime (time ()))); 
fprintf(' ok. \n');
% -----------------------
% --- Present results ---
% -----------------------

Mx=zeros(Nx,Ny);
My=zeros(Nx,Ny);

for Lx=1:2.*Nx.*Ny

[ix,iy,kx,ky]= L2ik(1,Lx,Nx);

if mod(Lx,2)==1
   Mx(kx,ky)=M(Lx,1);
   else
   My(kx,ky)=M(Lx,1);
   
   end

end 

Bx=Mx.*mi0+Hxext.*mi0;
By=My.*mi0+Hyext.*mi0;

% save -v6 dane2.mat mi Mx My Nx Ny dL g Hxext Hyext

Nfx=Nx;
Nfy=Ny;
dLf=dL;

save -v6 Cf.mat TCxx TCyy TCxy TCyx Nfx Nfy dLf

% save -v6 daneB.mat Bx By mi 

figure (3);
pcolor(xx,yy,Bx);
colorbar;
title('Bx');

figure (4);
pcolor(xx,yy,By);
colorbar;
title('By');

fprintf('\n ');
fprintf(strftime ("%Y-%m-%d %H:%M:%S", localtime (time ()))); 
fprintf(' *** All done. *** \n\n\n');

% save -binary results.mat *

diary off;

