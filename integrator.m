classdef integrator < handle
  properties (Access = protected)
    % A character string. Can be 'left', 'right', 'middle', 'trapezes', 'gauss2'
    % or 'gauss3'
    method;
    xk; % the interpolation points, scaled to the interval [0,1]
    wk; % The weight of each interpolation point, to approximate the integral on [0,1]
    % The integral on a small interval [ai,bi] of width di=bi-ai is approximated by
    % di * ( wk(1)*f(ai + xk(1)*di) + wk(2) * f(ai+xk(2)*di + ... )

    dx; % When computing primitive(f,from,x), the interval [0,x] is split in n
        % subintervals, n such that (x-from)/n ~= dx.
        % More precisely, n = round( (x-from)/dx ) (or n=1 if the latter formula
        % evaluates to 0).
    draw;
  endproperties


  methods (Access = public)
    % Constructor. The arguments can be:
    % - An integrator. In this case, varargin has length 1, it is varargin{1}
    % - Options about the properties method and dx. You call it with pairs :
    %     integrator("method",value,"dx",value)
    %   Default value of method is "trapezes"
    %   Default value of dx is 0.1
    function obj = integrator (varargin)
      obj.method = "trapezes";
      obj.dx = 0.1;
      if nargin > 0
        if isa(varargin{1}, 'integrator')
          obj.method = varargin{1}.method;
          obj.xk = varargin{1}.xk;
          obj.wk = varargin{1}.wk;
          obj.dx = varargin{1}.dx;
        else
          for i = 1: nargin
          switch varargin{i}
            case "method"
              obj.method = varargin{i+1};
            case "dx"
              obj.dx = varargin{i+1};
          endswitch
          endfor
        endif
      endif
      obj.update();
    endfunction

    % Return the property value. Argument prop can be "method" or "dx".
    % The function is already coded, you don't need to change it.
    function retval = get (this, prop)
      switch (prop)
        case "method"
          retval = this.method;
        case "dx"
          retval = this.dx;
      endswitch
    endfunction

    % Set the required properties to the required values. Arguments are a list
    % of pairs property,value
    % The properties can be "method" or "dx". You can call the method with, for
    % example (assuming that itg is an integrator):
    %   itg.set("method","left","dx",0.1)
    %   or itg.set("method","left")
    %   or itg.set("dx",0.1)
    % If varargin is empty, nothing is done.
    function this = set (this, varargin)
      try
        for i = 1: nargin
          switch varargin{i}
            case "method"
              this.method = varargin{i+1};
            case "dx"
              this.dx = varargin{i+1};
          endswitch
        endfor
        this.update();
      end_try_catch
      this.update();
    endfunction

    % This function displays the integrator properties.
    % No need to change it.
    function disp(this)
      printf("method: %s\n", this.method);
      printf("dx = %.2e\n", this.dx);
    endfunction

    function this = update(this)
      switch this.method
      
        case "trapezes"
          this.wk = [0.5 0.5];
          this.xk = [0 1];
          this.draw = @(f, ai, bi, hax) plot(hax, [ai ai bi bi ai], [0 f(ai) f(bi) 0 0], 'r', "linewidth", 1);
          
        case "left"
          this.wk = [1];
          this.xk = [0];
          this.draw = @(f, ai, bi, hax) plot(hax, [ai ai bi bi ai], [0 f(ai) f(ai) 0 0], 'r', "linewidth", 1);
          
        case "right"
          this.wk = [1];
          this.xk = [1];
          this.draw = @(f, ai, bi, hax) plot(hax, [ai ai bi bi ai], [0 f(bi) f(bi) 0 0], 'r', "linewidth", 1);
        
        case "middle"
          this.wk = [1];
          this.xk = [1/2];
          this.draw = @(f, ai, bi, hax) plot(hax, [ai ai bi bi ai], [0 f((bi+ai) / 2) f((bi+ai) / 2) 0 0], 'r', "linewidth", 1);
        
        case "gauss2"
          this.wk = [0.5 0.5];
          this.xk = [(3-sqrt(3))/6 (3+sqrt(3))/6];
        case "gauss3"
          this.wk = [5/18 4/9 5/18]
          this.xk = [(5-sqrt(15))/10 0.5 (5+sqrt(15))/10];
      endswitch
    endfunction



    % This function compute the integral of f between a and b, using the
    % quadrature formula specified by the property method.
    % Arguments are:
    % - f : a function handle, for example @sin or @(x) x.^2
    % - a and b : two scalar numerical values, the bounds of integrations
    % - n : a scalar numerical value (integer value). When doing the integration,
    %   the interval [a,b] is split in n subintervals [ai,bi] of equal width
    %   bi-ai = (b-a)/n
    % - hax is an optional argument. If it is set, it is a graph handle. Use it
    %   to draw an explanaory graphic of the integration method (graph of the
    %   function to integrate, graph of the appromating polynomial on each interval
    %   [ai,bi], etcleft
    % The function returns the integral approximation value (a scalar numerical).

    % A character string. Can be 'left', 'right', 'middle', 'trapezes', 'gauss2'
    % or 'gauss3'
    function I = integrate(this,f,a,b,n,hax)
      range = linspace(a, b, n);
      

      if(nargin == 6)
        drawing = true;
        hold off;
        plot(hax, linspace(a, b, 1000), f(linspace(a, b, 1000)), 'b', "linewidth", 2);
        hold on;
      else
        drawing = false;
      endif

      result = 0;
      for i = 1: length(range)-1
        
        ai = range(i);
        bi = range(i+1);
        if(drawing)
          this.draw(f, ai, bi, hax);
        endif

        sum = 0;
        for k=1:length(this.xk)
          sum += this.wk(k)*f(ai + this.xk(k)*(bi - ai));
        endfor

        result += (bi-ai)*sum;
        drawnow()
        % result += calcul(f, ai, bi);
        
      endfor

      
      
      % hold off;
      I = result;
    endfunction

    % Computes the values of a primitive of f: the integral of f between from,
    % and x.
    % Arguments are:
    % - f : a function handle, for example @sin or @(x) x.^2
    % - from : a scalar numerical value (the bound below in the integration)
    % - x : the values at which the primitive function is evaluated. It can be
    %   a scalar value or a matrix.
    % The returned value y is a matrix of same size as x
    %
    % Note:
    % The function needs to execute a for loop, if argument x has more than 1
    % element. For each x(i,j), compute y(i,j) using the function integrate.
    % The number of subintervals of [from,x(i,j)] (argument n in the integrate
    % function) depends on the length x(i,j)-from and on this.dx. Take
    %   n = round( (x(i,j)-from) / this.dx ) (and set it to 1 if the latter
    % formula results in 0).
    function y = primitive (this, f, from, x)
      y = zeros(size(x)); % Erase this line and write your code
    endfunction


    % This function studies the integration error, and returns a model of the
    % form |Error| = C / n^alpha.
    % when the interval [a,b] is split in n subintervals [ai,bi]
    %
    % Arguments
    % - f: a function handle
    % - a and b: two scalar numerical values, the integration bounds
    % - Ith: theoretical value of the integral. A scalar numerical value.
    % - ns: the values of n that we test. A vector of numerical (with integrer
    %   values). At each value n of ns, call the function integrate and compare
    %   its output with Ith.
    % - hax: an optional argument. If it is given it is a axe handle where you
    %   can do a graphical representation of the error pattern.
    %
    % Output is a 2-element vector [C,alpha] such that the error model
    % is C / n^alpha.
    function retval = integration_error (this, f, a, b, Ith, ns, hax)
      retval = [0,0]; % Erase this line and write your code
    endfunction
  endmethods
endclassdef





