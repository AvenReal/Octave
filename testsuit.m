figure
ax = axes()
a = integrator("method", "middle")



function static res = mygraph(f, binf, bsup, from, to)



    f = @(x) 1./x;
    binf = -10;
    bsup = 10;
    from = 1;
    to = 500;
    
    a = integrator("method", "left");
    b = integrator("method", "right");
    c = integrator("method", "middle");
    d = integrator();

    figure;
    ax = axes();
    hold on;

    for i = from:to
        plot(ax, [i i+1], [a.integrate(f, binf, bsup, i) a.integrate(f, binf, bsup, 1+i)] , 'r');
        plot(ax, [i i+1], [b.integrate(f, binf, bsup, i) b.integrate(f, binf, bsup, 1+i)] , 'g');
        plot(ax, [i i+1], [c.integrate(f, binf, bsup, i) c.integrate(f, binf, bsup, 1+i)] , 'b');
        plot(ax, [i i+1], [d.integrate(f, binf, bsup, i) d.integrate(f, binf, bsup, 1+i)] );
    endfor

    res = 0
endfunction