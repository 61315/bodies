% Test driver for odeBoxBox

rng(1);
density = 1;
s1 = [1 1 1];%[3 2 2];
s2 = [1 1 1];%[3 2 1];
E1 = eye(4); se3.randE();
E2 = [eye(3), [0 0.1 0.99]'; 0 0 0 1]%se3.randE();
% E1 = [rand(3),[0 0 0]'; 0 0 0 1];
% E2 = [rand(3),[0 0 -0.45]'; 0 0 0 1];
clf;
hold on;
axis equal;
view(3);
grid on;
se3.drawCuboid(E1,s1);
se3.drawCuboid(E2,s2);

collisions = odeBoxBox_mex(E1,s1,E2,s2)
%collisions = AffineBoxBox(E1,s1,E2,s2);

for i = 1 : collisions.count
	xw = collisions.pos(:,i);
	nw = collisions.nor;
	xnw = [xw,xw+nw];
	plot3(xw(1,:),xw(2,:),xw(3,:),'rs')
	plot3(xnw(1,:),xnw(2,:),xnw(3,:),'r-')
end
