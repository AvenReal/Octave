
function res = mygraph(f, binf, bsup, from, to)
    a = integrator("method", "left");
    b = integrator("method", "right");
    c = integrator("method", "middle");
    d = integrator();

    figure;
    ax = axes();
    hold on;

    for i = from:to
        plot(ax, i, a.integrate(f, binf, bsup, i), 'r');
        plot(ax, i, b.integrate(f, binf, bsup, i), 'g');
        plot(ax, i, c.integrate(f, binf, bsup, i), 'b');
        plot(ax, i, d.integrate(f, binf, bsup, i));
    endfor

    return;
endfunction
