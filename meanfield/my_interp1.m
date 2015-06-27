function y1 = my_interp1(x,y,x1,dlogx,n_range)
    logx1 = log10(x1);
    if logx1 <= n_range(1)
      ind = 1;
    elseif logx1 < n_range(2)
      ind = ceil((logx1-n_range(1))/dlogx)+1;
    else
        ind = length(x)-2; % if x1 is outside the range, 
                           %it jsut continues the linear extrapolation.
    end
    dy = y(ind+1)-y(ind);
    dx = x(ind+1)-x(ind);
    
    y1 = (x1-x(ind))*(dy/dx)+y(ind);
end