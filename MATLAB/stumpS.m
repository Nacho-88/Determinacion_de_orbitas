% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
function s = stumpS(z)
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
%
% This function evaluates the Stumpff function S(z) according
% to Equation 3.49.
%
% z- input argument
% s- value of S(z)
%
% User M-functions required: none
%-----------------------------------------------------------
if z > 0
    s = (sqrt(z)- sin(sqrt(z)))/(sqrt(z))ˆ3;
elseif z < 0
    s = (sinh(sqrt(-z))- sqrt(-z))/(sqrt(-z))ˆ3;
else
    s = 1/6;
end
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜

