function [best]=firefly_simple(instr)
% n=number of fireflies
% MaxGeneration=number of pseudo time steps
if nargin<1,   instr=[40 1000];     end
n=instr(1);  MaxGeneration=instr(2);
% Show info
help firefly_simple.m
rand('state',0);  % Reset the random generator
% ------ Four peak functions ---------------------
str1='exp(-(x-4)^2-(y-4)^2)+exp(-(x+4)^2-(y-4)^2)';
str2='+2*exp(-x^2-(y+4)^2)+2*exp(-x^2-y^2)';
str3 = '(1000/sqrt(2*pi).*exp(-((x-2).^2/2)-((y-2).^2/2))+(200/sqrt(pi/2).*exp(-((x+1).^2/2)-((y-1).^2/2))+(1000/sqrt(pi).*exp(-((x-4).^2/2)-((y+3).^2/2)))))+(1000/sqrt(2*pi).*exp(-((x+9).^2/2)-((y+5).^2/2)))';
funstr=strcat(str3);
% Converting to an inline function
f=vectorize(inline(funstr));

% range=[xmin xmax ymin ymax];  a
range=[-15 15 -10 10];

% ------------------------------------------------
alpha=0.6;      % Randomness 0--1 (highly random)
gamma=2.0;      % Absorption coefficient
delta=1;      % Randomness reduction (similar to 
                % an annealing schedule)
% ------------------------------------------------
% Grid values are used for display only
Ngrid=100;
dx=(range(2)-range(1))/Ngrid;
dy=(range(4)-range(3))/Ngrid;
[x,y]=meshgrid(range(1):dx:range(2),...
               range(3):dy:range(4));
z=f(x,y);
% Display the shape of the objective function
figure(1);    surfc(x,y,z);

% ------------------------------------------------
% generating the initial locations of n fireflies
[xn,yn,Lightn]=init_ffa(n,range);
% Display the paths of fireflies in a figure with
% contours of the function to be optimized
writerObj = VideoWriter('Firefly.avi'); % Name it.
writerObj.FrameRate = 20; % How many frames per second.
open(writerObj);

 figure(2);
% Iterations or pseudo time marching
for i=1:MaxGeneration,     %%%%% start iterations
% Show the contours of the function
 contour(x,y,z,15); hold on;
% Evaluate new solutions
zn=f(xn,yn);

% Ranking the fireflies by their light intensity
[Lightn,Index]=sort(zn);
xn=xn(Index); yn=yn(Index);
xo=xn;   yo=yn;    Lighto=Lightn;
% Trace the paths of all roaming  fireflies
plot(xn,yn,'.','markersize',10,'markerfacecolor','g');
% Move all fireflies to the better locations
[xn,yn]=ffa_move(xn,yn,Lightn,xo,yo,Lighto,alpha,gamma,range);
drawnow;
% Use "hold on" to show the paths of fireflies
    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
    writeVideo(writerObj, frame);
    hold off;
% Reduce randomness as iterations proceed
if i<100
    delta = 1.0;
else
    delta = 0.97;
end
    alpha=newalpha(alpha,delta);
    
end   %%%%% end of iterations
close(writerObj);
best(:,1)=xo'; best(:,2)=yo'; best(:,3)=Lighto';

% ----- All subfunctions are listed here ---------
% The initial locations of n fireflies
function [xn,yn,Lightn]=init_ffa(n,range)
xrange=range(2)-range(1);
yrange=range(4)-range(3);
% xn=rand(1,n)*xrange+range(1);
% yn=rand(1,n)*yrange+range(3);

xn = range(1)*ones(1,n/4);
yn = range(3)*ones(1,n/4);
xn((n/4)+1:n/2) = range(2)*ones(1,n/4);
yn((n/4)+1:n/2) = range(4)*ones(1,n/4);
xn((n/2)+1:3*n/4) = range(2)*ones(1,n/4);
yn((n/2)+1:3*n/4) = range(3)*ones(1,n/4);
xn((3*n/4)+1:n) = range(1)*ones(1,n/4);
yn((3*n/4)+1:n) = range(4)*ones(1,n/4);
Lightn=zeros(size(yn));

% Move all fireflies toward brighter ones
function [xn,yn]=ffa_move(xn,yn,Lightn,xo,yo,...
    Lighto,alpha,gamma,range)
ni=size(yn,2); nj=size(yo,2);
str3 = '(1000/sqrt(2*pi).*exp(-((x-2).^2/2)-((y-2).^2/2))+(200/sqrt(pi/2).*exp(-((x+1).^2/2)-((y-1).^2/2))+(1000/sqrt(pi).*exp(-((x-4).^2/2)-((y+3).^2/2)))))';
funstr=strcat(str3);
% Converting to an inline function
f=vectorize(inline(funstr));
for i=1:ni,
% The attractiveness parameter beta=exp(-gamma*r)
    for j=1:nj,
r=sqrt((xn(i)-xo(j))^2+(yn(i)-yo(j))^2);
if Lightn(i)<Lighto(j), % Brighter and more attractive
beta0=1;     beta=beta0*exp(-gamma*r.^2);
xn(i)=xn(i).*(1-beta)+xo(j).*beta+alpha.*(rand-0.5);
yn(i)=yn(i).*(1-beta)+yo(j).*beta+alpha.*(rand-0.5);
% elseif f(xn(i),yn(i))<=0.5
% xn(i) = xn(i)+(rand);
% yn(i) = yn(i)+(rand);
end
    end % end for j
end % end for i
[xn,yn]=findrange(xn,yn,range);

% Reduce the randomness during iterations
function alpha=newalpha(alpha,delta)
alpha=alpha*delta;

% Make sure the fireflies are within the range
function [xn,yn]=findrange(xn,yn,range)
for i=1:length(yn),
   if xn(i)<=range(1), xn(i)=xn(i)+rand(1); end
   if xn(i)>=range(2), xn(i)=xn(i)-rand(1); end
   if yn(i)<=range(3), yn(i)=range(3); end
   if yn(i)>=range(4), yn(i)=range(4); end
end
