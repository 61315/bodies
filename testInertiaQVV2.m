clear all; clc; close all;

%% Test with maximal coordinates

density = 1;
l = 5;
w = 1;
sides = [l w w];

E = eye(4);
%E(1:3,1:3) = se3.aaToMat([1 1 1],pi/4);
%E(1:3,4) = [0.5*l 0 2*l]';

% phi = [3 -4 5 -3 2 -5]';
phi = [1 2 3 0.1 0.2 0.3]'; % [angular velocity; linear velocity]

E0 = E;
phi0 = phi;

clf
hold on;
axis equal;
axis(5*[-1 1 -1 1 -1 1]);
se3.drawCuboid(E,sides);
view(3);
grid on;
ax = gca;
ax.Clipping = 'off';

% 6x1, 1:3 is rotational (principal), 4:6 is just mass repeated (or translational)
m = se3.inertiaCuboid(sides,density) ;
% `M` here is a 6x6 diagonal
M = diag(m);

tEnd = 3;
h = 1e-3;
drawHz = 30;
nsteps = ceil(tEnd/h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test with affine
E = E0;
phi = phi0;

R = E(1:3,1:3); % linear transform
p = E(1:3,4); % translation
r = reshape(R',9,1); % linear transform (the paper use this row major order)

% Edot means, you can take the time derivative of a transform matrix `E`.
% i can also put this way. the time derivative of transform matrix `E` is, `E` experiencing generalized velocity `phi`.
% similar to the time derivative of a vector is that vector experiencing an angular velocity omega.
% meaning, you can apply a generalized velocity [w; v] directly to the SE(3) `E`using that brac thing.
Edot = E*se3.brac(phi);
% Compute time derivative of E using body-fixed twist
% Geometric Insight:
% 1. E*se3.brac(phi) represents the transformation of body-fixed velocity to world frame
% 2. This right multiplication (E*[phi]) is equivalent to left multiplication ([Ad_E(phi)]*E)
%    where Ad_E is the Adjoint transformation
% 3. It preserves the group structure of SE(3) and is useful for on-manifold integration
% 4. This formulation directly relates to the Lie group exponential map: E(t+dt) ≈ E(t)*exp([phi]*dt)

% E = |R t|  Ad_E = |R [t]×R|
%     |0 1|         |0     R|

% E is SE(3), phi is se(3) ?? Edot and phi are isomorphic???



rdot = reshape(Edot(1:3,1:3)',9,1);
pdot = Edot(1:3,4);

R2A = 0.5*[
	-1  1  1
	1 -1  1
	1  1 -1
	];

M;
M(1:3,1:3); % this is rotational inertia (principal axised)
diag(M(1:3,1:3)); % first three diagonal entries as 3x1 vector
repmat(diag(M(1:3,1:3)),3,1); % above repeated for each column vector (axis) in R

I = R2A*diag(M(1:3,1:3)); % linear transformational?? inertia.
W = repmat(1./I,3,1);

iters = 1;
t = 0;
for k = 1 : nsteps
	% Unconstrained step
	p0 = p; % set previous translation
	p = p + h*pdot; % just add like velocity for translation
	r0 = r; % set previous linear transform
	r = r + h*rdot; % just add like velocity for linear transform. it IS linear. cool. you just add two.
	% Apply affine constraint
	for i = 1 : iters
		for j = 1 : 6
			[C,dC] = conAffine(r,j); % get the current violation and the direction of the maximum change
			WdC = W.*dC';
			numerator = -C;
			denominator = dC*WdC;
			dlambda = numerator/denominator;
			r = r + dlambda*WdC;
		end
	end
	% Compute velocities
	pdot = (p - p0)/h;
	rdot = (r - r0)/h;
	t = t + h;
	if drawHz > 0 && (floor(t*drawHz) > floor((t-h)*drawHz))
		cla;
		hold on;
		R = reshape(r,3,3)';
		E(1:3,1:3) = R;
		E(1:3,4) = p;
		se3.drawCuboid(E,sides);
		Rdot = reshape(rdot,3,3)';
		%phi = se3.unbrac(R'*Rdot);
		phibrac = R'*Rdot;
		phi = se3.unbrac(0.5*(phibrac - phibrac'));
		w = R*phi(1:3);
		I = R*M(1:3,1:3)*R';
		l = I*w;
		ppw = [p,p+w];
		ppl = [p,p+l];
		plot3(ppw(1,:),ppw(2,:),ppw(3,:),'r');
		plot3(ppl(1,:),ppl(2,:),ppl(3,:),'g');
		title(sprintf('%f',t));
		drawnow;
	end
end


%% Test with affine without matrix exponential
E = E0;
phi = phi0;

R = E(1:3,1:3);
p = E(1:3,4);
r = reshape(R',9,1);

% let's get pdot and rdot without Edot
% Rdot is R*[w] and pdot is R*v
omega = phi(1:3); % and this omega is body frame omega. if you want the world omega, you need to do R*w
v = phi(4:6); % is `v` local or global? it's local!! (at least in this file)
Rdot = R*se3.brac(omega);
rdot = reshape(Rdot',9,1);

% `p` is in world frame and also is `pdot`

% if v is body frame velocity, you multiply R and make it world frame velocity
pdot = R*v;

% if v is in world frame, 
% pdot = v + cross(omega, p)

R2A = 0.5*[
	-1  1  1
	1 -1  1
	1  1 -1
	];

I = R2A*diag(M(1:3,1:3));
W = repmat(1./I,3,1);

iter = 1;
t = 0;
for k = 1 : nsteps
	p0 = p;
	p = p + h*pdot;
	r0 = r;
	r = r + h*rdot;
	
	for i = 1:iters
		for j = 1:6
			[C, dC] = conAffine(r,j);
			WdC = W.*dC';
			numerator = -C;
			denominator = dC*WdC;
			dlambda = numerator/denominator;
			r = r + dlambda*WdC;
		end
	end

	pdot = (p - p0)/h;
	rdot = (r - r0)/h;
	t = t + h;
	if drawHz > 0 && (floor(t*drawHz) > floor((t-h)*drawHz))
		cla;
		hold on;
		R = reshape(r,3,3)';
		E(1:3,1:3) = R;
		E(1:3,4) = p;
		se3.drawCuboid(E,sides);
		Rdot = reshape(rdot,3,3)';
		
		% i still don't get this phibrac part
		% this is to get the phi back from (corrected q vector)
		phibrac = R'*Rdot; % why is this R has to be transposed? > it's the inverse, dumbass!! R' == inv(R)

		% i think i know what this thing is!!!
		% it's making a skew-symmetric matrix by substract a matrix by it's transpose
		% and 0.5 is there bc, the matrix appeared twice!! it all makes sense now
		phi = se3.unbrac(0.5*(phibrac - phibrac'));

		w = R*phi(1:3); % is this omega? `world omega` shoud i say? > this means the omega in phi is local. > then is v local too?
		I = R*M(1:3,1:3)*R'; % this should be the world inertia
		l = I*w; % and this is the (world) angular momentum?
		ppw = [p, p+w]; % 3x2
		ppl = [p, p+l];
		plot3(ppw(1,:),ppw(2,:),ppw(3,:),'r');
		plot3(ppl(1,:),ppl(2,:),ppl(3,:),'g');
		title(sprintf('%f',t));
		drawnow;
		
	end
end

%%
function [C,dC] = conAffine(r,i)
A11 = r(1);
A12 = r(2);
A13 = r(3);
A21 = r(4);
A22 = r(5);
A23 = r(6);
A31 = r(7);
A32 = r(8);
A33 = r(9);
% i know R'R is symmetric, but is it really ok to use only 6 components?
% C = 6x1, []
% |2*A11 0 0
% |0 2*A12 0
% |0 0 2*A13
switch i
	case 1
		C = A11^2 + A21^2 + A31^2 - 1;
		dC = [2*A11, 0, 0, 2*A21, 0, 0, 2*A31, 0, 0];
	case 2
		C = A12^2 + A22^2 + A32^2 - 1;
		dC = [0, 2*A12, 0, 0, 2*A22, 0, 0, 2*A32, 0];
	case 3
		C = A13^2 + A23^2 + A33^2 - 1;
		dC = [0, 0, 2*A13, 0, 0, 2*A23, 0, 0, 2*A33];
	case 4
		C = A11*A12 + A21*A22 + A31*A32;
		dC = [A12, A11, 0, A22, A21, 0, A32, A31, 0];
	case 5
		C = A12*A13 + A22*A23 + A32*A33;
		dC = [0, A13, A12, 0, A23, A22, 0, A33, A32];
	case 6
		C = A11*A13 + A21*A23 + A31*A33;
		dC = [A13, 0, A11, A23, 0, A21, A33, 0, A31];
end
end
