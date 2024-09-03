classdef (Abstract) Constraint < handle
    %CONSTRAINTPOSITION Super class for Constraints
    
    properties
        bodies
        C
        dC
        dC2
    end
    
    %%
    methods
        %%
        function this = Constraint(bodies_)
            this.bodies = bodies_;
        end
        
        %%
        function solvePositions(this)
            nbodies = length(this.bodies);
            if nbodies == 1
                body = this.bodies{1};
                lambda = -this.C/(this.dC*(this.dC'./body.M));
                
                if ~body.isFixed
                    body.q = body.q + lambda*(this.dC'./body.M);
                    body.updateE();
                end
            elseif nbodies == 2
                body1 = this.bodies{1};
                body2 = this.bodies{2};
                
                %lambda = -this.C/(this.dC*(this.dC'./body1.M + this.dC'./body2.M));
                lambda = - this.C/(this.dC * (this.dC'./body1.M) + this.dC2 * (this.dC2'./body2.M));
                if ~body1.isFixed
                    body1.q = body1.q + lambda*(this.dC'./body1.M);
                    body1.updateE();
                end
                
                %lambda = -this.C/(this.dC2*(this.dC2'./body1.M + this.dC2'./body2.M));
                if ~body2.isFixed
                    body2.q = body2.q - lambda*(this.dC2'./body2.M);
                    body2.updateE();
                end
            end
        end
    end
end

