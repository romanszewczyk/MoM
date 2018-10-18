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
% Script for visualization of the results of modelling the sensitivity 
% of racetrack-shaped core fluxgate sensor.
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

clear all

load ('Results.mat');

RozklB=zeros(max(Tco(:,1)),max(Tco(:,2)));

for i=1:n-1

    By=squeeze(TBy(i,:,:));
    SigBy=By(:,31:51);
    SigSum=sum(sum(SigBy));
    RozklB(Tco(i,1),Tco(i,2))=SigSum;

end   

RozklB=RozklB./max(max(RozklB));

surf(RozklB);
set(gca, 'fontsize', 20)
xlabel ('c (mm)');
ylabel ('b (mm)');
zlabel ('sensitivity (a.u.)');

