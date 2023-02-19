


clear
hold off
clf
hold on
ic_dat;

for i = 1:1:nf;
worm3(Ux(i,:),Uy(i,:),Uz(i,:),0.2,12);
%line(Ux(i,:),Uy(i,:),Uz(i,:));
hold on
end;
v = [-2.0,2.0,-2.0,2.0,-2.0,2.0];
axis(v);


%          s_step =  2*pi/(np-1);
%	  radius = 0.5;
	    
%	  for i=1:1:np;
          
%	  Ux(1,i) = sin(2 * (i-1) * s_step - pi/4);
%	  Ux(1,i) = radius * Ux(1,i);
%	  Uy(1,i) = 2. * cos((i-1) * s_step);
%	  Uy(1,i) = radius * Uy(1,i); 
%	  Uz(1,i) = 1.5 * sin((i-1) * s_step);
%          Uz(1,i) = radius * Uz(1,i); 
%          end;


