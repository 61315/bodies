classdef ConstraintVolume < Constraint
    %CONSTRAINTGROUNDCONTACT Constraint for ground contact
    
    properties
    end
    
    methods
        %%
        function this = ConstraintVolume(body_)
            this = this@Constraint(body_);
        end
        
        %%
        function solvePositions(this)
            body = this.bodies{1};
            A = reshape(body.q(4:12),3,3)';
            
            % phi = det(A) - I
            this.C = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) - A(1,2)*A(2,1)*A(3,3) ...
                    + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1) - 1;
            
            % interpretation: translation part of A(3) does not affect the volumn change
            this.dC = [
                                            0,                             0,                             0, ...
                A(2,2)*A(3,3) - A(2,3)*A(3,2), A(2,3)*A(3,1) - A(2,1)*A(3,3), A(2,1)*A(3,2) - A(2,2)*A(3,1), ...
                A(1,3)*A(3,2) - A(1,2)*A(3,3), A(1,1)*A(3,3) - A(1,3)*A(3,1), A(1,2)*A(3,1) - A(1,1)*A(3,2), ...
                A(1,2)*A(2,3) - A(1,3)*A(2,2), A(1,3)*A(2,1) - A(1,1)*A(2,3), A(1,1)*A(2,2) - A(1,2)*A(2,1)
                ];
            solvePositions@Constraint(this);
            this.bodies{1}.updateE();
        end
    end
end

