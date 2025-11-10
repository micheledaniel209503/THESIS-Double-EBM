% gives points that define an arc, given two points and the sagitta
function [xarc,yarc] = arc_two_pts(A,B,s,NP, bulgeSign)
    % A,B: 2 points
    % s: sagitta (freccia) must be orthogonal to the chord (alpha is just
    % that)
    % NP: number of points to discretize arc
    if nargin < 4, NP = 200; end
    v = B - A;                 % chord vector
    L = hypot(v(1), v(2));     % chord length
    assert(L > 0, 'A and B are coincident.');
    c = 0.5*L; % half chord
    M = 0.5*(A + B); % midpoint on the chord

    % tangent and normal
    t = v / L; % tangent is actually the direction of the chord
    n = [-t(2), t(1)]; % normal is the tangent rotated 90 degrees!!!
    if bulgeSign * n(2) < 0, n = -n; end   % switch the direction of the bulge ("try and error" basically)

    % radius + center of the circumference associated with the arc
    R = (c^2 + s^2) / (2*s);
    d = R - s;
    C = M + d * n;

    % peak bulge point = midpoint on the chord + the sagitta, directed as
    % the normal (see notes)
    Pbulge = M + s * n;

    % BUILD THE ARC, discretized
    % angles of A, B, Pbulge from the CENTER point
    thA = atan2(A(2)-C(2), A(1)-C(1));
    thB = atan2(B(2)-C(2), B(1)-C(1));
    thP = atan2(Pbulge(2)-C(2), Pbulge(1)-C(1));

    % wrap [-pi,pi)
    wrap = @(x) mod(x+pi, 2*pi) - pi;

    % build the two paths A->B (there are 2 possible paths to reach B from A and viceversa), choose the one that goes through thP
    th1a = thA; th1b = thB; if wrap(th1b - th1a) < 0, th1b = th1b + 2*pi; end
    th2a = thA; th2b = thB; if wrap(th2b - th2a) > 0, th2a = th2a + 2*pi; end

    % check that thP is in between [tha, thb], otherwise no good, choose
    % the other path
    in_interval = @(tha,thb,thx) (wrap(thx-tha) >= 0) && (wrap(thx-thb) <= 0);

    if in_interval(th1a, th1b, thP)
        th = linspace(th1a, th1b, NP);
    else
        th = linspace(th2a, th2b, NP);
    end

    % calculate the arc points based on the center, angles and radius
    xarc = C(1) + R*cos(th);
    yarc = C(2) + R*sin(th);
end
