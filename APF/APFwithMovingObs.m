clc
clear
close all

%% global path
x = 0:0.01:20;
y = 0.7.*x.*(x>=0&x<10)+(0.3.*x+4).*(x>=10&x<=20);
plot(x,y)
axis([0 20 0 10])
set(gca,'XTick',0:1:20);
set(gca,'YTick',0:1:10);
title('Paths in a 2D grid map')
grid on
hold on

% local goal point 
for i = 1:21
    x_lg(i) = i-1;
    if x_lg(i)>=0 && x_lg(i)<10
        y_lg(i) = 0.7*x_lg(i);
    else
        y_lg(i) = 0.3*x_lg(i)+4;
    end
    plot(x_lg(i),y_lg(i),'*');
    hold on
end
lg = [x_lg;y_lg;1:i]';

%% obstacles
O = [0  4;   2 1.5;4 3;6 5;5 5.5;8 5.5;9 7;10 7;11 6.8;12 8;13 7;15 7;17 9;19 9];
% O = [100 100];
v = [0.041 0;0 0;0 0;0 0;0 0;  0   0;0 0; 0 0; 0   0; 0 0; 0 0; 0 0; 0 0; 0 0];
% v = [0 0;0 0;0 0;0 0;0 0;  0   0;0 0; 0 0; 0   0; 0 0; 0 0; 0 0; 0 0; 0 0];
plot(O(:,1),O(:,2),'o');
hold on

%% initialize
sp = [0 0];                                          % position of start ponit
vr = [0 0];
n = size(lg,1);                                      % number of local points
gp = lg(n,1:2);
% k = 5;                                               % gravitational gain coeffecient
eta = 1000;
eta_v = 500;
rou0 = 1;
len_step = 0.05;
no = size(O,1);                                      % number of obstacles
Num_iter = 400;
A = 1000;
A0 = 5;
B = 100;
B0 = 5;
pic_num = 1;

%% main cyclic
Xr = sp;
i = 0;

while sqrt((Xr(1)-lg(n,1)).^2+(Xr(2)-lg(n,2)).^2) > 0.1  % tolerance
    i = i + 1;
    Path(i,:) = Xr;                                  % record the path point
    
    if O(1,1) < 12
        O(1,1:2) = [v(1,1)*i,O(1,2)];
        plot(O(1,1),O(1,2),'o');
%         text(O(1,1)+0.01,O(1,2)+0.01,num2str(i)) ;
        hold on
    end
    
    for j = 1:n
        delta(j,1:2) = lg(j,1:2)-Xr(1:2);
        delta(j,3) = lg(j,3);
        dist(j,1) = norm(delta(j,1:2));
        dist(j,2) = lg(j,3);
    end
    
    dist_ = sortrows(dist,1);
    [m2,o2] = max(dist_(1:2,2));
    [m1,o1] = min(dist_(1:2,2));
    dist_lg = dist_(o2,1);                               % distance between current & local goal point
    id2 = dist_(o2,2);
    id1 = dist_(o1,2);
    x1 = lg(find(lg(:,3)==id1),1:2);
    x2 = lg(find(lg(:,3)==id2),1:2);
    unitVector_lg = [delta(find(delta(:,3)==id2),1)/dist_lg,delta(find(delta(:,3)==id2),2)/dist_lg];
    d1 = abs(det([x2-x1;Xr-x1]))/norm(x2-x1);            % d = abs(det([Q2-Q1;P-Q1]))/norm(Q2-Q1);
%     p = 24;
    p = A*d1 + A0; %adaptive attracting factor 
    
    delta_gg = gp-Xr;                  % >0
    dist_gg = norm(delta_gg);                    % distance between current point & global goal point
    unitVector_gg(1:2)=[delta_gg(1)/dist_gg,delta_gg(2)/dist_gg];
%     g = 0;
    g = B/(dist_lg+B0);
    
    F_att_gg = [2*g*dist_gg*unitVector_gg(1), 2*g*dist_gg*unitVector_gg(2)];
    F_att_lg = [2*p*dist_lg*unitVector_lg(1), 2*p*dist_lg*unitVector_lg(2)];
    F_att = F_att_gg + F_att_lg;
    
    for k = 1:no
        del(k,1:2) = Xr(1:2)-O(k,1:2); % vector pointing from obstacle position to robot position
        rou(k) = norm(del(k,:));
        if rou(k) > rou0
            F_rep(k,:) = [0 0];
        else
            uv(k,1:2) = [del(k,1)/rou(k),del(k,2)/rou(k)];  % umit vector pointing from obstacle position to robot position
            F_rep1_norm = eta*(1/rou(k)-1/rou0)*(p/(p+g)*dist_lg^2+g/(p+g)*rou(k)^2)/rou(k)^2;
            F_rep1 = [F_rep1_norm*uv(k,1), F_rep1_norm*uv(k,2)];                       % point from obstacles to robot
            F_rep2 = eta*(1/rou(k) - 1/rou0)^2*(p/(p+g)*(x2-Xr)+g/(p+g)*(gp-Xr));         % point from robot to goal
            
            vri = vr-v(k,1:2);    % vector pointing from obstacle to robot
            sigma = acos(dot(vri,uv(k,1:2))/(norm(vri)*norm(uv(k,1:2)))); % rad
            if sigma>-pi/2 && sigma<pi/2 && norm(v(k,1:2))>0
                F_rep3_norm = norm(eta_v.*vri/rou(k));
                F_rep3 = [F_rep3_norm*uv(k,2), F_rep3_norm*(-uv(k,1))];
                F_rep4 = eta_v*uv(k,1:2);
                F_rep(k,:) = F_rep1+F_rep2+F_rep3+F_rep4;
            else
                F_rep(k,:) = F_rep1+F_rep2;
            end
        end
    end
    
    F_rep = [sum(F_rep(:,1)),sum(F_rep(:,2))]; 
    F_sum = [F_rep(1,1)+F_att(1,1),F_rep(1,2)+F_att(1,2)];
    UnitVec_Fsum(i,:) = 1/norm(F_sum) * F_sum;
    
    vr = len_step*UnitVec_Fsum(i,:);
    Xr(1,1:2) = Xr(1,1:2)+vr;
    figure(1);
    set(gcf,'outerposition',get(0,'screensize'));
    plot(Path(i,1),Path(i,2),'.');
%     pause(0.1);
    drawnow;
    F=getframe(gcf);
    I=frame2im(F);
    [I,map]=rgb2ind(I,256);
    if pic_num == 1
        imwrite(I,map,'test.gif','gif', 'Loopcount',inf,'DelayTime',0.03);
    else
        imwrite(I,map,'test.gif','gif','WriteMode','append','DelayTime',0.03);
    end
    pic_num = pic_num + 1;
%     text(Path(i,1)+0.01,Path(i,2)+0.01,num2str(i)) ;
%     if i>Num_iter
%         break
%     end
end
