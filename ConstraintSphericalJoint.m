classdef ConstraintSphericalJoint < Constraint
    %CONSTRAINTGROUNDCONTACT Constraint for spherical joint
    
    % NOTE: this whole joint implementation looks like it's from an articulated
    % body implementation, which is confusing (to me). let's make it more me
    % friendly and drop the Lie stuff for now if possible.
    
    % joint will always have two bodies
    
    properties
        b0
        b1
        
        r0 % relative attachment point of joint in body0 coordinate frame (local)
        r1 % relative attachment point of joint in body1 coordinate frame (local)
        q0 % relative attachment orientation in body0 coordinate frame (local)
        q1 % relative attachment orientation in body1 coordinate frame (local)
    end
    
    methods    
        function this = ConstraintSphericalJoint(b0, b1, r0, r1)
            this = this@Constraint({b0, b1});

            this.r0 = r0;
            this.r1 = r1;
        end
        
        function solvePositions(this)
            body1 = this.bodies{1};
            R0 = body1.E_wi(1:3,1:3);
            t0 = body1.E_wi(1:3,4);
            r0_world = t0 + R0 * this.r0;
            
            body2 = this.bodies{2};
            R1 = body2.E_wi(1:3,1:3);
            t1 = body2.E_wi(1:3,4);
            r1_world = t1 + R1 * this.r1;
            
            delta_x = r0_world - r1_world; % note the direction
            
            this.C = norm(delta_x); % current violation
            delta_x = delta_x ./ this.C;
            
            this.dC = delta_x'* se3.getJacobian(this.r0);
            this.dC2 = delta_x'* se3.getJacobian(this.r1);
            
            % Jacobian, 3x12
            % |     |x x x 0 0 0 0 0 0|
            % |  I  |0 0 0 x x x 0 0 0|
            % |     |0 0 0 0 0 0 x x x|
            % translation does no contribution (I), work is done by q(4:12)
            % local attachment point repeated 3 times for each column in the
            % linear transform matrix
            
            if this.C > 1e-12
                solvePositions@Constraint(this);
            end
        end
    end
end