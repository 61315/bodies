classdef ConstraintGroundContact < Constraint
    %CONSTRAINTGROUNDCONTACT Constraint for ground contact
    
    properties
        collisions
        E_g
        planar % 2D planar scene
    end

    %%
    methods
        %%
        function this = ConstraintGroundContact(body_, E_g_)
            this = this@Constraint(body_);
            this.E_g = E_g_;
            this.planar = false;
        end

        %%
        function collide(this)
			% Clear collisions
			this.collisions = {};
			
			% Ground collisions (only if ground.E isn't the zero matrix)
			if this.E_g(4,4) == 1
				this.collisions = this.bodies{1}.collideGround(this.E_g, this.planar, this.collisions);
            end
        end

        %%
        function solvePositions(this)
            this.collide();
            for i = 1 : length(this.collisions)
				c = this.collisions{i};
		        this.C = c.nw'*c.xw;
		        this.dC = c.nw'* se3.getJacobian(c.xl{1});
                solvePositions@Constraint(this);
            end
        end
    end
end

