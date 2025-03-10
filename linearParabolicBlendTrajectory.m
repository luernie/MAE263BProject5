function [q,qd,qdd, qddd] = linearParabolicBlendTrajectory(t0, tf, q0, qf, tb, qdd, N)% linearParabolicBlendTrajectory:
    % linearParabolicBlend generates joint-space positions from q0 to qf
    %   with a trapezoidal velocity profile (accel-constant-accel).
    %
    % INPUTS:
    %       q0  - initial joint angles/Px initial
    %       qf  - final joint angles/Px final
    %       t0  - initial time
    %       tf  - final time
    %       tb  - blend time (calculate qd)
    %       qdd - desired acceleration (calculate blend time)
    
    %       N   - number of trajectory points
    %
    % Outputs:
    %   pos - [1 x N] positions
    %   vel - [1 x N] velocities
    %   acc - [1 x N] accelerations
    %   jerk   - [1 x N] jerk
    
    T = tf - t0;
    % Time vector
    t = linspace(t0, tf, N);
    
    if(tb == 0 && qdd ~= 0) % calculate tb
        tb = 0.5*T-sqrt(qdd^2*T^2-4*abs(qdd)*abs(qf-q0))/abs(2*qdd);
    elseif(tb ~= 0 && qdd == 0) % calculate qdd
        qdd = abs((qf-q0)/(tb*(T-tb)));
    end
    
    % if qdd < 4 * (qf - q0) / T^2 return  !!!
    
    
    if q0 < qf
        ab0 = [1, t0, t0^2; 0, 1, 2*t0;0, 0, 2] \ [q0 0 qdd]'; % Coefficient for first blend
        abf = [1, tf, tf^2; 0, 1, 2*tf;0, 0, 2] \ [qf 0 -qdd]'; % Coefficient for second blend
    else
        ab0 = [1, t0, t0^2; 0, 1, 2*t0;0, 0, 2] \ [q0 0 -qdd]'; % Coefficient for first blend
        abf = [1, tf, tf^2; 0, 1, 2*tf;0, 0, 2] \ [qf 0 qdd]'; % Coefficient for second blend
    end
    qb1 = ab0(1) + ab0(2)*(t0+tb)+ab0(3)*(t0+tb)^2;
    qb2 = abf(1)+abf(2)*(tf-tb)+abf(3)*(tf-tb)^2;
    
    a = [1, t0+tb;1, tf-tb]\[qb1;qb2]; % Coefficient for linear region
    
    % first parabolic region
    t11 = t((t0<=t) & (t<=t0+tb));
    q = ab0(1)+ab0(2)*t11+ab0(3)*t11.^2;
    qd  = ab0(2)+2*ab0(3)*t11;
    qdd = 2*ab0(3)*ones(size(t11));

    % Linear region
    t22 = t((t0+tb<t) & (t<tf-tb)); % linear region
    q = [q, a(1)+a(2)*t22];
    qd = [qd, a(2).*ones(size(t22))];
    qdd = [qdd, zeros(size(t22))];

    % second parabolic region
    t33 = t((tf-tb<=t) & (t<=tf));
    q   = [q, abf(1)+abf(2)*t33+abf(3)*t33.^2];
    qd  = [qd, abf(2)+2*abf(3)*t33];
    qdd = [qdd, 2*abf(3)*ones(size(t33))];
    qddd = zeros(1, N);
    end
    