
function result = test ()

    f = @(x) x^2;

    a = integrator("method", "left", "dx");

    result = a.integrate(f, 0, 10, 100);
endfunction
