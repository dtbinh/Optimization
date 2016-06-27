clc; clear all;
N_grp = 5;
N_grp_popu = 5;
N_popu = N_grp*N_grp_popu;
alpha = 0.92;
betal = 0.1;
betag = 0.1;
for i = 1:N_grp
    
    point(:,:,i) = zeros(N_grp_popu,2);
    theta = (((i-1)*90/N_grp)+(90/N_grp)*lhsdesign(N_grp_popu,1))*pi/180;
    vel_mag = 0.1;
    in_vel(:,:,i) = vel_mag*rand(N_grp_popu,2).*[cos(theta) sin(theta)];
end

Target_locs = [5 9; 12 12; 9 5; 12 3];
state = [1 1 1 1];
% To display the function
range = [0 15 0 15];
Ngrid=100;
dx=(range(2)-range(1))/Ngrid;
dy=(range(4)-range(3))/Ngrid;
[xx,yy]=meshgrid(range(1):dx:range(2),...
               range(3):dy:range(4));
zz=Objfuncs(xx,yy,state);
% Display the shape of the objective function
figure(1);    surfc(xx,yy,zz);

threshold = 5;
figure(2);
exp_popu = zeros(N_grp_popu,1,N_grp);
loc_best = zeros(N_grp_popu,2,N_grp);
p_gbest = [];
count = 0;
% Making a movie
% writerObj = VideoWriter('Group-PSO1.avi'); % Name it.
% writerObj.FrameRate = 20; % How many frames per second.
% open(writerObj);
for iter = 1:1000
    count1 =1;
    for i = 1:N_grp
        F(:,:,i) = Objfuncs(point(:,1,i),point(:,2,i),state);
        var =[];
        var = find(F(:,:,i)>=threshold);
        var1 = find(F(:,:,i)<threshold);
        Func_value(:,iter,i) = F(:,:,i);
        Func_value(var1,iter,i)= 0;
        if length(var)==0
            point(:,:,i) = point(:,:,i)+in_vel(:,:,i);
            for jjj = 1:N_grp_popu
                if point(jjj,1,i)<range(1)
                    point(jjj,1,i) = range(1)
                    in_vel(jjj,1,i) = rand/4;
                end
                if point(jjj,2,i)<range(3)
                    point(jjj,2,i) = range(3)
                    in_vel(jjj,2,i) = rand/4;
                end
                if point(jjj,1,i)>range(2)
                    point(jjj,1,i) = range(2)
                    in_vel(jjj,1,i) = -rand/4;
                end
                if point(jjj,2,i)>range(4)
                    point(jjj,2,i) = range(4)
                    in_vel(jjj,2,i) = -rand/4;
                end
            end
        else
            count = count+count1
            max_F(count,:,i) = max(F(:,:,i));
            exp_popu(var,1,i) = 1;
            
            if count>1
                if max_F(end,:,i)>=max_F(end-1,:,i)
                    [ee gbest] = max(F(:,:,i));
                    p_gbest(:,:,i) = point(gbest,:,i);
                else
                    p_gbest(:,:,i) = p_gbest(:,:,i);
                end
            else
                [ee gbest] = max(F(:,:,i));
                p_gbest(:,:,i) = point(gbest,:,i);
            end

            
            if length(find(exp_popu(:,:,i)==1))>0
                var1 = find(exp_popu(:,:,i)==1);
                for ii = 1:length(find(exp_popu(:,:,i)==1))
                    if max(F(var1(ii),:,i))> Objfuncs(loc_best(var1(ii),1,i),loc_best(var1(ii),1,i),state)
                        loc_best(var1(ii),:,i) =  point(var1(ii),:,i);
                    end
                    in_vel(var1(ii),:,i) = alpha*in_vel(var1(ii),:,i)+betag*rand*(p_gbest(:,:,i)-point(var1(ii),:,i))/2+betal*rand*(loc_best(var1(ii),:,i)-point(var1(ii),:,i))/2;
                    
                end
            end
            if length(find(exp_popu(:,:,i)==0))>0
                var2 = find(exp_popu(:,:,i)==0);
                for jj = 1:length(find(exp_popu(:,:,i)==0))
                    in_vel(var2(jj),:,i) = alpha*in_vel(var2(jj),:,i)+betag*rand*(p_gbest(:,:,i)-point(var2(jj),:,i));
                end
            end
            point(:,:,i) = point(:,:,i) + in_vel(:,:,i);
            count1 = 0;
            if iter>5
                m = mean(Func_value(:,iter-10:iter,i),2);
                std_dev = std(m);
                if std_dev <0.5
                    cnt = 1;
                    for iii = 1:4
                        if p_gbest(:,1,i)<Target_locs(iii,1)+0.3 && p_gbest(:,1,i)>Target_locs(iii,1)-0.3 && p_gbest(:,2,i)<Target_locs(iii,2)+0.3 && p_gbest(:,2,i)>Target_locs(iii,2)-0.3
                            fprintf('Target Found \n');
                            state(iii) = 0;
                            break
                        end
                    end
                    p_gbest = [];
                    loc_best = zeros(N_grp_popu,2,N_grp);
                    in_vel(:,:,i) = (1/4)*(rand(5,2)*2-1);
                    max_F = [];
                    count = 0;
                    count1 = 1;
                end
            end
                
        end
        for jjj = 1:N_grp_popu
            if point(jjj,1,i)<range(1)
                point(jjj,1,i) = range(1)
                in_vel(jjj,1,i) = rand/4;
            end
            if point(jjj,2,i)<range(3)
                point(jjj,2,i) = range(3)
                in_vel(jjj,2,i) = rand/4;
            end
            if point(jjj,1,i)>range(2)
                point(jjj,1,i) = range(2)
                in_vel(jjj,1,i) = -rand/4;
            end
            if point(jjj,2,i)>range(4)
                point(jjj,2,i) = range(4)
                in_vel(jjj,2,i) = -rand/4;
            end
        end
        xn((i-1)*N_grp_popu+1:i*N_grp_popu,1) = point(:,1,i);
        
        yn((i-1)*N_grp_popu+1:i*N_grp_popu,2) = point(:,2,i);
    
    end
    ZZ = Objfuncs(xx,yy,state);
    contour(xx,yy,ZZ,15); hold on;
    plot(xn,yn,'.','markersize',10,'markerfacecolor','g');
    hold off;
    if length(find(ZZ==0))>0
        fprintf('All the targets are found\n')
        break
    end
%     frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
%     writeVideo(writerObj, frame);
    pause(0.001)
end

% close(writerObj);