


clear
hold off
clf
hold on
fc_dat;

for i = 1:1:nf;
worm3(Ux(i,:),Uy(i,:),Uz(i,:),0.2,12);
hold on
end;
v = [-2.0,2.0,-2.0,2.0,-2.0,2.0];
axis(v);


%v = [6,10,-2.0,2.0,-2.0,2.0];
%axis(v)

