classdef BodyCuboid < Body
	%%
	properties
		sides       % Side lengths
		collision2D % Apply 2D collision along the X=0, Y=0, or Z=0 plane
	end
	
	%%
	methods
		%%
		function this = BodyCuboid(density,sides)
			this = this@Body(density);
			this.sides = sides;
			this.collision2D = '';
		end
		
		%%
		function computeInertia_(this)
			% I is the 6x1 diagonal rigid inertia, assuming that the frame origin is at the
			% center of mass and the axes are oriented along the principal axes. We store
			% the rotation on top of translations, so that I(1:3) is the rotational inertia
			% and I(4:6) is the translational inertia.
			I = se3.inertiaCuboid(this.sides, this.density);
			this.mass = I(4);
			% This matrix transforms from rigid to affine.
			R2A = [
				-1  1  1
				1 -1  1
				1  1 -1
				]/2;
			% Transform translational inertia
			this.M(1:3) = I(4:6);
			% Transform rotational inertia, repeating 3 times
			% each axis gets its own mass, hence the 3.
			this.M(4:6) = R2A*I(1:3);
			this.M(7:9) = this.M(4:6);
			this.M(10:12) = this.M(4:6);
			this.M = this.M';
			
			% computed once, which is great
			% assert(false)
		end
		
		%%
		function collisions = collideGround_(this,groundE,planar,collisions)
			% assert(false)
			% caller: this@Body.collideGround
			%
			% Scene.solvePositions > ConstraintGroundContact.solvePosition > ConstraintGroundContact.collide > Body.collideGround > BodyCuboid.collideGround_
			
			% Collision detection
			xg = groundE(1:3,4); % ground origin
			ng = groundE(1:3,3); % ground normal
			S = eye(4);
			S(1:3,1:3) = diag(0.5*this.sides); % 0.5 for half extent
			% switch this.collision2D
			% 	case 'X'
			% 		xl = S*[
			% 			 0 -1 -1  1
			% 			 0 -1  1  1
			% 			 0  1 -1  1
			% 			 0  1  1  1
			% 			]';
			% 	case 'Y'
			% 		xl = S*[
			% 			-1  0 -1  1
			% 			-1  0  1  1
			% 			 1  0 -1  1
			% 			 1  0  1  1
			% 			]';
			% 	case 'Z'
			% 		xl = S*[
			% 			-1 -1  0  1
			% 			-1  1  0  1
			% 			 1 -1  0  1
			% 			 1  1  0  1
			% 			]';
			% 	otherwise
			xl = S*[
				-1 -1 -1  1
				-1 -1  1  1
				-1  1 -1  1
				-1  1  1  1
				 1 -1 -1  1
				 1 -1  1  1
				 1  1 -1  1
				 1  1  1  1
				]';
			% end

			% interpretation: if any of the vertices for each side is below the ground plane, there is a collision??
			
			% what is `xl`??? (xl:xw = x_local:x_world) is it the case??
			
			%xl = xl(:,1); % for debugging JointPlanar
			%xl = xl(:,end); % for debugging
			%xl = xl(:,[1,5]); % for debugging
			
			xw = this.E_wi*xl; % it is the case!
			for i = 1 : size(xl,2)
				xli = xl(1:3,i); % i-th x_local
				xwi = xw(1:3,i); % i-th x_world
				% penetration depth
				d = ng'*(xwi - xg);
				if d > 0
					% No collision
					continue
				end
				
				% Contact point world velocity
				collision.bodies{1} = this;
				collision.xl{1} = xli; % local pos
				collision.xw = xwi; % world pos
				collision.nw = ng; % world nor from ground
				%collision.Tw = this.getTangentBasis(collision,planar); % world tangent based on velocity
				collision.d = d; % penetration depth
				%collision.l = []; % initial Lagrange multiplier
				%collision.idxB = []; % index into the global
				collisions{end+1} = collision; %#ok<AGROW>
			end
		end
		
		%%
		function collisions = collideBody_(this,that,planar,collisions)
			if ~isa(that,'BodyCuboid')
				% Only Cuboid-Cuboid is supported
				return;
			end
			c = odeBoxBox_mex(this.E_wi,this.sides,that.E_wi,that.sides);
			for i = 1 : c.count
				xw = c.pos(:,i);
				% if this.collision2D == 'Z'
				% 	% Assume that the 2D scene is on the z=0 plane
				% 	if xw(3) < 0
				% 		% Remove collisions with negative z
				% 		continue;
				% 	end
				% 	% Make the positive collisions happen at z=0
				% 	xw(3) = 0;
				% end
				nw = -c.nor; % sign depends on the convention of the collision detector
				d = -c.depth(i);
				collision.bodies{1} = this;
				collision.bodies{2} = that;
				collision.xl{1} = this.E_iw(1:3,:)*[xw;1];
				collision.xl{2} = that.E_iw(1:3,:)*[xw;1];
				collision.xw = xw;
				collision.nw = nw;
				%collision.Tw = this.getTangentBasis(collision,planar); % world tangent based on velocity
				collision.d = d;
				collisions{end+1} = collision; %#ok<AGROW>
			end
		end
		
		%%
		function [F,V] = draw_(this)
			[F,V] = se3.patchCuboid(this.E_wi,this.sides);
			patch('Faces',F,'Vertices',V,'FaceColor',this.color);
		end
		
		%%
		function s = getAxisSize(this)
			s = min(this.sides);
		end
	end
end