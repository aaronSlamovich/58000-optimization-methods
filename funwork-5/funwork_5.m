% --------problem 1---------
% plotting rastrigin's function
axis = linspace(-11,11,100);
x1 = (axis./10).^2 - 10 * cos((2*pi*axis)/10);
x2 = (axis./10).^2 - 10 * cos((2*pi*axis)/10);
M = 20 + x1 + x2';
surfc(axis,axis,M)
% calculate rastrigin
function r = rastrigin(coord)
   r = 20 + (coord(1)./10).^2 + (coord(2)./10).^2  - 10.*(cos((2*pi*coord(1))/10) + cos((2*pi*coord(2))/10));
end
% convert phenotype to genotype
order = 16; % order for each coordinate
function g = genotype(order, coord)
   coord = coord - 8.5;
   coord = coord ./ 2;
   coord = coord * (2.^order - 1);
   coord = round(coord);
   coord = dec2bin(coord);
   g = cat(2, coord(1,:), coord(2,:)); % will truncate to make as short as possible
end
g = genotype(16, [9,10])
% convert genotype to phenotype
function p = phenotype(order, genotype)
   p(1) = bin2dec(genotype(1:16));
   p(2) = bin2dec(genotype(17:32));
   p = p./(2.^order - 1);
   p = p * 2;
   p = p + 8.5;
end
p = phenotype(order, g)

