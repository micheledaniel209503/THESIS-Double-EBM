clc;
clear all;
%% 
function double_EBM_interactive()
% PARAMETERS
ri = 8e-3;               % [m]
ro = 15e-3;              % [m]
h0 = 1.5e-3;             % [m]
params.ri = ri; params.ro = ro; params.h0 = h0;

% slider arrays (domain)
lh = 201; lr = 201;
h_vec  = linspace(1e-5, h0-1e-5, lh)';
rc_vec = linspace(ro - 6e-3, ro*0.99, lr)';

% initial values
h  = h_vec(round(lh*0.4));
rc = rc_vec(round(lr*0.6));

% FIGURE
f = figure('Name','DOUBLE EBM cross section - interactive');
ax = axes('Parent',f); hold(ax,'on'); grid(ax,'on'); box(ax,'on'); axis(ax,'equal');
xlabel(ax,'x [mm]'); ylabel(ax,'y [mm]');
xlim(ax, 1e3*[-ro ro]);
ylim(ax, 1e3*[0 2.5*h0]);

% graphical objects
% bottom EBM
hl = gobjects(9,1);   % 9 bottom ebm segments
% top EBM
hTopMid = plot(ax, NaN, NaN, 'b', 'LineWidth', 2.0); % mid line
harc = gobjects(4,1); % 4 top ebm arcs

% reference at h0
y0line = yline(ax, h0*1e3, 'k--', 'DisplayName','h0');

% title
ttl = title(ax,'');

% SLIDERS AND LABELS
% slider h
uicontrol('Style','text','Units','normalized','Position',[0.12 0.02 0.15 0.04],...
          'String','h [mm]');
sH = uicontrol('Style','slider','Units','normalized','Position',[0.12 0.06 0.30 0.04],...
               'Min',h_vec(1),'Max',h_vec(end),'Value',h,...
               'SliderStep',[1/(lh-1) min(5/(lh-1),1)],'Callback',@onChange);
labH = uicontrol('Style','text','Units','normalized','Position',[0.44 0.06 0.10 0.04],...
                 'String',sprintf('%.3f',h*1e3));

% slider rc
uicontrol('Style','text','Units','normalized','Position',[0.58 0.02 0.15 0.04],...
          'String','r_c [mm]');
sRC = uicontrol('Style','slider','Units','normalized','Position',[0.58 0.06 0.30 0.04],...
                'Min',rc_vec(1),'Max',rc_vec(end),'Value',rc,...
                'SliderStep',[1/(lr-1) min(5/(lr-1),1)],'Callback',@onChange);
labRC = uicontrol('Style','text','Units','normalized','Position',[0.90 0.06 0.08 0.04],...
                  'String',sprintf('%.3f',rc*1e3));

% update during slide
addlistener(sH, 'Value', 'PostSet', @(~,~) onChange());
addlistener(sRC,'Value', 'PostSet', @(~,~) onChange());

% First draw
redraw();

% callbacks
    function onChange(~,~)
        h  = sH.Value;
        rc = sRC.Value;
        labH.String  = sprintf('%.3f',h*1e3);
        labRC.String = sprintf('%.3f',rc*1e3);
        redraw();
    end

    function redraw()
        % geometry commons
        h_bot = max(h0 - h, 0);
        h_top = h0 + h;
        h_tot = h_bot + h_top; % = 2*h0

        % BOTTOM EBM points
        P0  = [0,   h_bot];
        P1D = [ri,  h_bot];
        P2D = [rc,  h_bot/2];
        P3D = [ro,  h_bot/2];
        P4D = [ri,  0];

        P1L = [-ri, h_bot];
        P2L = [-rc, h_bot/2];
        P3L = [-ro, h_bot/2];
        P4L = [-ri, 0];

        % TOP EBM points
        P1Dt = [ri, h_tot];
        P3Dt = [ro, h_bot + h_top/2];
        P3Lt = [-ro, h_bot + h_top/2];
        P1Lt = [-ri, h_tot];

        % compute alpha
        alpha = max(0, alpha_sol(h, rc, params));

        % DRAW BOTTOM
        segs = {
            P0,P1D; P1D,P2D; P2D,P3D; P2D,P4D; P4D,P4L; P4L,P2L; P2L,P3L; P2L,P1L; P1L,P0
        };
        for k=1:9
            [xk, yk] = segXY(segs{k,1}, segs{k,2});
            if ~isgraphics(hl(k))
                hl(k) = plot(ax, xk, yk, 'b', 'LineWidth', 2.0);
            else
                set(hl(k), 'XData', xk, 'YData', yk);
            end
        end

        % MID LINE
        [xk, yk] = segXY(P1Lt, P1Dt);
        set(hTopMid, 'XData', xk, 'YData', yk);

        % 4 ARCS
        [x1,y1] = arc_two_pts(P3Dt, P1Dt, alpha, 200, -1);
        [x2,y2] = arc_two_pts(P3Dt, P1D,  alpha, 200, +1);
        [x3,y3] = arc_two_pts(P1Lt, P3Lt, alpha, 200, -1);
        [x4,y4] = arc_two_pts(P3Lt, P1L,  alpha, 200, +1);

        harc(1) = updateArc(ax, harc(1), x1, y1, 'arcTop1');
        harc(2) = updateArc(ax, harc(2), x2, y2, 'arcTop2');
        harc(3) = updateArc(ax, harc(3), x3, y3, 'arcTop3');
        harc(4) = updateArc(ax, harc(4), x4, y4, 'arcTop4');

        % TITLE
        ttl.String = sprintf('h = %.3f mm,  r_c = %.3f mm,  \\alpha = %.3f mm', ...
                             h*1e3, rc*1e3, alpha*1e3);
    end


    function [xmm, ymm] = segXY(A,B)
        xmm = [A(1) B(1)]*1e3;
        ymm = [A(2) B(2)]*1e3;
    end

    

end
% helper for updating the arc
function h = updateArc(ax, h, x, y, tagStr)
    X = x*1e3; Y = y*1e3;
    if ~isgraphics(h)
        h = plot(ax, X, Y, 'b', 'LineWidth', 2.0, 'Tag', tagStr);
    else
        set(h, 'XData', X, 'YData', Y);
    end
end

double_EBM_interactive
