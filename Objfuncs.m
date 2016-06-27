function Y = Objfuncs(x,y,state);

Y = (state(1,1)*(1000/sqrt(4*pi).*exp(-((x-5).^2/2)-((y-9).^2/2))))+(state(1,2)*(400/sqrt(pi/2).*exp(-((x-12).^2/2)-((y-12).^2/2))))+(state(1,4)*(1000/sqrt(pi).*exp(-((x-12).^2/2)-((y-3).^2/2))))+(state(1,3)*(1000/sqrt(2*pi).*exp(-((x-9).^2/2)-((y-5).^2/2))));
end