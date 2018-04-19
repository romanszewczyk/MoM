function res=IsInTriangle(x1, y1, x2, y2, x3, y3, x, y)
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

d=(x1.*(y2-y3)+y1.*(x3-x2)+x2.*y3-y2.*x3);

t1=(x.*(y3-y1)+y.*(x1-x3)-x1.*y3+y1.*x3)./d;
t2=(x.*(y2-y1)+y.*(x1-x2)-x1.*y2+y1.*x2)./(-1.*d);

% if (0<=t1)&&(t1<=1)&&(0<=t2)&&(t2<=1)&&(t1+t2<=1)
%    res=1;
%    else
%    res=0;
%    end
% end
 
 res=(0<=t1)&(t1<=1)&(0<=t2)&(t2<=1)&(t1+t2<=1);
    
% http://totologic.blogspot.fr/2014/01/accurate-point-in-triangle-test.html
