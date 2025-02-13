figure
ax = axes()
a = integrator("method", "middle")



function mygraph(f, binf, bsup, from, to)
    
    % a = integrator("method", "left");
    % b = integrator("method", "right");
    % c = integrator("method", "middle");
    % d = integrator();
    ee = integrator("method", "gauss2");
    ff = integrator("method", "gauss3");

    figure;
    ax = axes();
    hold on;

    for i = from:to
        % plot(ax, [i i+1], [a.integrate(f, binf, bsup, i) a.integrate(f, binf, bsup, 1+i)] , 'r'); % left
        % plot(ax, [i i+1], [b.integrate(f, binf, bsup, i) b.integrate(f, binf, bsup, 1+i)] , 'g'); % right
        % plot(ax, [i i+1], [c.integrate(f, binf, bsup, i) c.integrate(f, binf, bsup, 1+i)] , 'b'); % middle
        % plot(ax, [i i+1], [d.integrate(f, binf, bsup, i) d.integrate(f, binf, bsup, 1+i)] , 'y'); % trapezes
        plot(ax, [i i+1], [ee.integrate(f, binf, bsup, i) ee.integrate(f, binf, bsup, 1+i)] , 'm'); % gauss 2
        plot(ax, [i i+1], [ff.integrate(f, binf, bsup, i) ff.integrate(f, binf, bsup, 1+i)] , 'c'); % gauss 3
        drawnow;
    endfor

endfunction

mygraph(@(x) x.^2, 0, 10, 100, 500)