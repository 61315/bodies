
classdef (Abstract) Body < handle
	%Body An affine body connected through a joint to a parent
	
	%%
	properties
		name     % Optional name
		color    % Color for rendering
		density  % Mass/volume
		damping  % Viscous damping
		collide  % Whether collision is enabled
		I_i      % Inertia at body center
		I_inv    % Inverese of Inertia at body center
		M        % Affine Inertia
		mass     % Mass
		q        % states
		q_prev   % Previous states
		qdot     % Derivative of states
		phi      % Rigid velocity in local space
		E_wi     % Where the body is wrt world, world transform
		E_wi_dot
		E_iw     % The inverse of E_wi, world inverse transform
		next     % Next body in traversal order
		
		isFixed
	end
	
	%%
	methods
		%%
		function this = Body(density)
			this.name = ['body'];
			this.color = [0.8 0.8 0.3];
			this.density = density;
			this.damping = 0;
			this.collide = false;
			this.isFixed = false;
			% Note: whenever we use reshape(), don't forget that we want q to be
			% row-wise, so we need the transpose.
			this.E_wi = eye(4);
			this.q = [this.E_wi(1:3,4); reshape(this.E_wi(1:3,1:3)',9,1)];
			this.phi = [0 0 0 0 0 0]'; % rigid velocity in local space
			this.E_wi_dot = this.E_wi*se3.brac(this.phi);
			this.qdot = zeros(12,1);
			this.qdot(1:3) = this.E_wi_dot(1:3,4);
			this.qdot(4:12) = reshape(this.E_wi_dot(1:3,1:3)',9,1);
		end
		%%
		function updateE(this)
			this.E_wi(1:3,4) = this.q(1:3);
			this.E_wi(1:3,1:3) = reshape(this.q(4:12),3,3)'; % row-by-row
			this.E_iw = se3.inv(this.E_wi);
		end
		%%
		function updateQdot(this)
			if ~this.isFixed
				% Sets the transform of this body wrt parent joint
				this.E_wi_dot = this.E_wi*se3.brac(this.phi);
				% this.WorldTransformTimeDerivative = this.WorldTransform * ???
				
				this.qdot = zeros(12,1);
				this.qdot(1:3) = this.E_wi_dot(1:3,4);
				this.qdot(4:12) = reshape(this.E_wi_dot(1:3,1:3)',9,1);
			else
				this.qdot = zeros(12,1);
			end
		end
		%%
		function setBodyTransform(this, x, A)
			% Sets the transform of this body wrt parent joint
			this.E_wi = eye(4);
			this.E_wi(1:3,1:3) = A;
			this.E_wi(1:3,4) = x;
			this.q = [this.E_wi(1:3,4); reshape(this.E_wi(1:3,1:3)',9,1)];
		end
		
		%%
		function computeInertia(this)
			% Computes inertia at body and joint
			this.computeInertia_();
		end
		%%
		function [Fs,Vs] = draw(this,Fs,Vs)
			if nargin == 1
				Fs = {};
				Vs = {};
			end
			%s = this.getAxisSize();
			%if s > 0
			%	se3.drawAxis(this.E_wi,s);
			%end
			%text(this.E_wi(1,4),this.E_wi(2,4),this.E_wi(3,4),this.name);
			[Fs{end+1},Vs{end+1}] = this.draw_();
			% Go to the next body
			if ~isempty(this.next)
				[Fs,Vs] = this.next.draw(Fs,Vs);
			end
		end
		
		%%
		function s = getAxisSize(this) %#ok<MANU>
			s = 1;
		end
		%%
		function collisions = collideGround(this,groundE,planar,collisions)
			% Runs collision detection with the ground
			if this.collide
				collisions = this.collideGround_(groundE,planar,collisions);
			end
			% Go to the next body
			if ~isempty(this.next)
				%collisions = this.next.collideGround(groundE,planar,collisions);
			end
		end
		
		%%
		function collisions = collideBody(this,that,planar,collisions)
			if this.collide && that.collide
				collisions = this.collideBody_(that,planar,collisions);
			end
		end
		%%
		function updateWithoutConstraints(this, force_ext, h)
			if ~this.isFixed
				this.q_prev = this.q;
				this.q = this.q + h * this.qdot + h^2*(force_ext./this.M);
				this.updateE();
			else
				this.q_prev = this.q;
				this.updateE();
			end
		end
		%%
		
		function updateAfterSolve(this, h)
			if ~this.isFixed
				this.qdot = (this.q - this.q_prev)/h;
			else
				this.qdot = zeros(12,1);
			end
		end
		
	end
end