

a = integrator("method", "left");
b = integrator("method", "right");
c = integrator("method", "middle");
d = integrator();

figure;
ax = axes();
hold on;

for i = 0:100
    plot(ax, i, a.integrate(@(x) sin(x.^2)*e^x + log(sqrt(x)), 0, 2*pi, i), 'r');
    plot(ax, i, b.integrate(@(x) sin(x.^2)*e^x + log(sqrt(x)), 0, 2*pi, i), 'g');
    plot(ax, i, c.integrate(@(x) sin(x.^2)*e^x + log(sqrt(x)), 0, 2*pi, i), 'b');
    plot(ax, i, d.integrate(@(x) sin(x.^2)*e^x + log(sqrt(x)), 0, 2*pi, i));
endfor





printf(a.get("method"))
% hold off;
% 
% b = integrator(a);
% 
% b.set("method", "middle");
% 
% hold on;
% a.integrate(@(x) x.^3 -7, -10, 20, 25, ax)