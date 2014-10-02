function h=labeltext(string, p)
if nargin<2, p=0; end

aL=axis;
% Set a little above the current plane
h=text(p, 1-p, ...
     string, 'units', 'normalized', ...
     'verticalalignment', 'top', ...
     'BackgroundColor', 'w', 'EdgeColor', 'k');

