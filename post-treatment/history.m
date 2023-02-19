


clear
hold off
clf
hold on
history_dat;

for j = 0:1:round(nsteps/hist_step)-1;
wrt_inc = j * nf;
for i = 1:1:nf;
%worm3(X(i+wrt_inc,:),Y(i+wrt_inc,:),Z(i+wrt_inc,:),0.2,12);
%plot(X(i+wrt_inc,:),Z(i+wrt_inc,:));
%plot(X(5,:),Z(5,:));
%plot(X(i+wrt_inc,:)+j/(round(nsteps/hist_step)-1)*2,Y(i+wrt_inc,:)) % vortex ring
%plot(Y(i+wrt_inc,:)+j/(round(nsteps/hist_step)-1)*2,X(i+wrt_inc,:))  % pair of vortices
plot3(X(i+wrt_inc,:),Y(i+wrt_inc,:),Z(i+wrt_inc,:));

hold on
end;
axis equal
box on

plot(X(1+wrt_inc,:),(sqrt((Y(1+wrt_inc,:)-Y(2+wrt_inc,:)).^2+(Z(1+wrt_inc,:)-Z(2+wrt_inc,:)).^2)))


end;
hold off

%v = [-2.0,8,-2.0,2.0,-2.0,2.0];
%axis(v);
% plot(Y(:,2)-1)

###################################
clear
close all
hold off
clf
hold on
history_dat;

for j = 0:1:round(nsteps/hist_step)-1;
wrt_inc = j * nf;
tot(j+1)=min((sqrt((Y(1+wrt_inc,:)-Y(2+wrt_inc,:)).^2+(Z(1+wrt_inc,:)-Z(2+wrt_inc,:)).^2)));
end;

ss=[0:1:round(nsteps/hist_step)-1];
plot(ss,tot)
hold on 
plot(ss,ss.^0*4*0.0112*0.540)
