% calculate rastrigin
function r = rastrigin(coord)
   r = 20 + (coord(:,1)./10).^2 + (coord(:,2)./10).^2  - 10.*(cos((2*pi*coord(:,1))/10) + cos((2*pi*coord(:,2))/10));
end