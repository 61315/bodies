classdef Scene < handle
	%Scene Test scenes for redmax
	
	properties
		name % scene name
		bodies % list of bodies
		equalities % list of equality constraints
		collisions % list of collisions
		limits % list of joint limits
		constraints %list of constraints
		ground % ground transform, with Z up
		tEnd % end time
		t % current time
		h % time step
		k % current step
		mu % coefficient of friction
		nsteps % number of steps to take
		sub_steps% number of substeps per step
		grav % gravity
	end
	
	methods
		function this = Scene()
			% Default values
			this.name = '';
			this.bodies = {};
			this.constraints = {};
			this.ground.E = zeros(4);
			this.ground.size = 10;
			this.tEnd = 1;
			this.h = 3e-2;
			this.t = 0;
			this.k = 0;
			this.mu = 0;
			this.nsteps = 0;
			this.grav = [0 0 -980]';
			this.sub_steps = 30;
		end
		
		%%
		function init(this)
			if usejava('jvm')
				colormap('default'); % Restore colormap mode
			end
			
			% Initialize bodies
			nbodies = length(this.bodies);
			for i = 1 : nbodies
				this.bodies{i}.computeInertia();
				if i < nbodies
					this.bodies{i}.next = this.bodies{i+1}; %#ok<*SAGROW>
				end
			end
			
			% Other initial values
			this.nsteps = ceil(this.tEnd/this.h);
		end
		
		%%
		function solvePositions(this)
			for i = 1 : length(this.constraints)
				this.constraints{i}.solvePositions();
			end
		end
	end
end

