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
% Script for calculations in modelling the magnetic field distribution 
% in the racetrack-shaped core of fluxgate sensors.
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

page_screen_output(0);

page_output_immediately(1);


a=15; % dlugosc prostego odcinka: 2a
% b=12;  % szerokosc sciezki: b = 1:12
% c=10;  % odleglosc miedzy prostymi odcinkami: 2c

mi=1000;
dL=1e-3;

g=2e-5;

Hyext=100;
Hxext=0;

b_max=12;
c_max=10;

TBx=zeros(b_max.*c_max,80,80);
TBy=zeros(b_max.*c_max,80,80);
Tmi=zeros(b_max.*c_max,80,80);
Tco=zeros(b_max.*c_max,2);

n=1;

for b=1:b_max
   for c=1:c_max

fprintf('\n\n ### Calculation b=%i of %i, c=%i of %i \n\n',b,b_max,c,c_max);

      mi_res=MiRacetrack (a,b,c,mi);

      [Bx, By]= MoM2Dvect (mi_res, dL, g, Hxext, Hyext);

      TBx(n,:,:)=Bx;
      TBy(n,:,:)=By;
      Tmi(n,:,:)=mi_res;
      Tco(n,1)=b;
      Tco(n,2)=c;
      
      n=n+1;
      
    end
end


save -v7 Results.mat TBx TBy Tco Tmi a n

