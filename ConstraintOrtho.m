classdef ConstraintOrtho < Constraint
    %CONSTRAINTGROUNDCONTACT Constraint for ground contact
    
    properties
    end
    
    methods
        %%
        function this = ConstraintOrtho(body_)
            this = this@Constraint(body_);
        end

        %%
        function solvePositions(this)
            body = this.bodies{1};

            % `q` is just the components of A(3) in 12x1 vector.
            A = reshape(body.q(4:12),3,3)';

            % `C` here is a 6x1 vector.
            C(1,1) = A(1,1)^2 + A(2,1)^2 + A(3,1)^2 - 1; % norm of the first column being 1
            C(2,1) = A(1,2)^2 + A(2,2)^2 + A(3,2)^2 - 1; % norm of the second column being 1
            C(3,1) = A(1,3)^2 + A(2,3)^2 + A(3,3)^2 - 1; % norm of the third column being 1
            C(4,1) = A(1,1)*A(1,2) + A(2,1)*A(2,2) + A(3,1)*A(3,2); % dot product between the first and the second column
            C(5,1) = A(1,2)*A(1,3) + A(2,2)*A(2,3) + A(3,2)*A(3,3); % dot product between the first and the third column
            C(6,1) = A(1,1)*A(1,3) + A(2,1)*A(2,3) + A(3,1)*A(3,3); % dot product between the second and the second column

            % `dC` here is (partical C)/(partial A), 6x3x3 tensor??
            % no, it's (partial C)/(partial q), 6x12 matrix.
            % i think `q` holds translation in the first 3 elements
            % interpretation: translation does not affect the orthogonal potential
            dC = [
                0, 0, 0, 2*A(1,1),        0,        0, 2*A(2,1),        0,        0, 2*A(3,1),        0,        0
                0, 0, 0,        0, 2*A(1,2),        0,        0, 2*A(2,2),        0,        0, 2*A(3,2),        0
                0, 0, 0,        0,        0, 2*A(1,3),        0,        0, 2*A(2,3),        0,        0, 2*A(3,3)
                0, 0, 0,   A(1,2),   A(1,1),        0,   A(2,2),   A(2,1),        0,   A(3,2),   A(3,1),        0
                0, 0, 0,        0,   A(1,3),   A(1,2),        0,   A(2,3),   A(2,2),        0,   A(3,3),   A(3,2)
                0, 0, 0,   A(1,3),        0,   A(1,1),   A(2,3),        0,   A(2,1),   A(3,3),        0,   A(3,1)
            ];
            for i = 1 : length(C)
                this.C = C(i);
                this.dC = dC(i,:);
                solvePositions@Constraint(this);
            end
            this.bodies{1}.updateE();
        end
    end
end

