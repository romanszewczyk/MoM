function mi_res=MiRacetrack (a,b,c,mi)
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

mi_res=ones(80,80);


if 2.*(a+b+c)>76
  fprintf('\n\n Racetrack is to long \n\n');
  return
end

mi_res(40-(c+b-1):40-c,40-a+1:40+a)=mi;
mi_res(40+c+1:40+(c+b),40-a+1:40+a)=mi;

for i=1:b
  for fi=0:180
    mi_res(round(40.5+(c+i-0.49).*cos(fi./360.*2.*pi)), ... 
       round(41+a+(c+i-0.49).*sin(fi./360.*2.*pi)))=mi;    
  end     
end

for i=1:b
  for fi=180:-1:0
    mi_res(round(40.5+(c+i-0.49).*cos(fi./360.*2.*pi)), ... 
       round(41-a-(c+i-0.49).*sin(fi./360.*2.*pi)))=mi;
  end     
end

end


