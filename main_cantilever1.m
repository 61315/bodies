clear; close all;
rng(2);

scene = Scene();
density = 1.0;
sides = [2 1 0.5];
scene.ground.E = eye(4);
scene.tEnd = 5;

% bodies
scene.bodies{1} = BodyCuboid(density,[1 1 1]);
scene.bodies{1}.setBodyTransform([-2 0 5]', eye(3));
scene.bodies{1}.isFixed = true;

scene.bodies{2} = BodyCuboid(density,[1 1 1]);
scene.bodies{2}.setBodyTransform([-1 0 5]', eye(3));

scene.bodies{3} = BodyCuboid(density,[1 1 1]);
scene.bodies{3}.setBodyTransform([0 0 5]', eye(3));

scene.bodies{4} = BodyCuboid(density,[1 1 1]);
scene.bodies{4}.setBodyTransform([1 0 5]', eye(3));

scene.bodies{5} = BodyCuboid(density,[1 1 1]);
scene.bodies{5}.setBodyTransform([2 0 5]', eye(3));

scene.bodies{6} = BodyCuboid(density,[1 1 1]);
scene.bodies{6}.setBodyTransform([3 0 5]', eye(3));

% connect bodies
connectCuboids(scene, scene.bodies{1}, scene.bodies{2}, 'x+', 'x-');
connectCuboids(scene, scene.bodies{2}, scene.bodies{3}, 'x+', 'x-');
connectCuboids(scene, scene.bodies{3}, scene.bodies{4}, 'x+', 'x-');
connectCuboids(scene, scene.bodies{4}, scene.bodies{5}, 'x+', 'x-');
connectCuboids(scene, scene.bodies{5}, scene.bodies{6}, 'x+', 'x-');

% other constraints
for i = 1:length(scene.bodies)
    % scene.constraints{end+1} = ConstraintVolume({scene.bodies{i}});
    scene.constraints{end+1} = ConstraintOrtho({scene.bodies{i}});
    scene.constraints{end+1} = ConstraintGroundContact({scene.bodies{i}}, scene.ground.E);
end

scene.sub_steps = 5;
scene.init();
grav = [0 0 -9.8]';
% grav = [0 0 0]';

f = zeros(12,1);
f(1:3) = scene.bodies{1}.M(1:3).*grav;

figure_handle = figure;

nsteps = scene.nsteps;
for i = 0 : nsteps - 1
    if ~ishandle(figure_handle)
        break;
    end
    h = scene.h / scene.sub_steps;
    for j = 0 : scene.sub_steps -1
        
        for k = 1 : length(scene.bodies)
            scene.bodies{k}.updateWithoutConstraints(f, h);
        end
        scene.solvePositions();
        for k = 1 : length(scene.bodies)
            scene.bodies{k}.updateAfterSolve(h);
        end
        % draw(scene);
    end
    % assert(false)
    draw(scene);
    scene.t = scene.t + scene.h;
end

%%
function draw(scene)
if scene.t == 0
    % clf;
    hold on;
    axis equal;
    axis(3*[-1 1 -1 1 0 2]);
    grid on;
    view(3);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    ax = gca;
    ax.Clipping = 'off';
    zoom('on');
    pan('on');
    rotate3d('on');
end

cla;
for i = 1 : length(scene.bodies)
    E = scene.bodies{i}.E_wi;
    sides = scene.bodies{i}.sides;
    se3.drawAxis(E);
    se3.drawCuboid(E,sides);
end
title(sprintf('%f',scene.t));
drawnow
end

%%
function connectCuboids(scene, body1, body2, face1, face2)
% face1, face2 are strings: 'x+', 'x-', 'y+', 'y-', 'z+', 'z-'

% define the vertices of the faces
vertices = containers.Map();
vertices('x+') = [ 0.5  0.5  0.5;  0.5  0.5 -0.5;  0.5 -0.5  0.5;  0.5 -0.5 -0.5];
vertices('x-') = [-0.5  0.5  0.5; -0.5  0.5 -0.5; -0.5 -0.5  0.5; -0.5 -0.5 -0.5];
vertices('y+') = [ 0.5  0.5  0.5; -0.5  0.5  0.5;  0.5  0.5 -0.5; -0.5  0.5 -0.5];
vertices('y-') = [ 0.5 -0.5  0.5; -0.5 -0.5  0.5;  0.5 -0.5 -0.5; -0.5 -0.5 -0.5];
vertices('z+') = [ 0.5  0.5  0.5; -0.5  0.5  0.5;  0.5 -0.5  0.5; -0.5 -0.5  0.5];
vertices('z-') = [ 0.5  0.5 -0.5; -0.5  0.5 -0.5;  0.5 -0.5 -0.5; -0.5 -0.5 -0.5];

v1 = vertices(face1);
v2 = vertices(face2);

% scale the vertices according to the body sizes
v1 = v1 .* repmat(body1.sides, 4, 1);
v2 = v2 .* repmat(body2.sides, 4, 1);

for i = 1:4
    scene.constraints{end+1} = ConstraintSphericalJoint(body1, body2, v1(i,:)', v2(i,:)');
end
end