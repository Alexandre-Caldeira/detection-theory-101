%Marc Reese
%shows the graphical convolution of two discrete signals, a and b.
%
% computed as: a(n)**b(n) = sum(a(m)b(n-m))
% 
% "m" can be thought of as a "shift time"
%
% usage: [conv,Frames] = graphConvolve(a,b,pausespacing)
%
% -----INPUT------
% a and b are vectors, can be of different length
%
% pausespacing is the pause betweeen each step in the convolution in seconds. Use it
% to control the speed of the animation. Default is 0.1 seconds if not
% included.
%
% ------OUTPUT------
%
% conv:     the result of the convolution operation
% Frames:   captured frames for playback with: movie(Frames);
%
% ------NOTES------
%
% - a and b must be zero padded to work properly, this is done automatically
%   by graphConvolve.
% - it is much easier and faster for the program to run if b is shorter
%   than a.
%
function [conv,M] = graphConvolve(a,b,pausespacing)
if nargin < 3
    pausespacing = 0.1;
end
%if b is a column vector, make it a row vector
sizeb = size(b);
if sizeb(1) ~=1
    if sizeb(2) ~=1
        error('make sure you have vector inputs')
    else
        b = b';
    end
end
%Length of the time shifting function
mlen = length(a)+length(b);
%zero pad both signals to mlen
zeropadda = zeros(1,(mlen-length(a)));
zeropaddb = zeros(1,(mlen-length(b)));
b = [b zeropaddb]; a = [a zeropadda];
%x-axes for a vector
xa = 0:length(a)-1;
%x-axes for b vector
xb = 0:length(b)-1;
%plot the two functions on the same graph; supblot top
h = figure;
red = [0.8 0 0];    %specify less intense red color for plots.
subplot(211);
stem(xa,a); hold on; stem(xb,b,'Color',red);
xlim([-length(b) mlen]);
M(1) = getframe;
pause(1);         %pause for a second to let the user see the image
%Flip the b function, with the goal b[k] = b[-k]
b = fliplr(b);
%shift the b axes accordingly so that b[k] = b[-k]
xb = -fliplr(xb);
%plot these changes
hold off;
stem(xa,a); hold on; stem(xb,b,'Color',red);
xlim([-length(b) mlen]);
M(2) = getframe;
pause(1);    
for m = 1:mlen-1
    conv(m) = sum((a(1:m).*b((end-m+1):end)));
    subplot(212); hold off; stem(conv,'k'); xlim([-length(b) mlen]);
    xb = xb+1;
    pause(pausespacing);
    subplot(211); hold off; stem(xa,a); hold on;stem(xb,b,'Color',red);xlim([-length(b) mlen]);
    M(m+2) = getframe;
end