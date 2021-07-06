function inv = betainvIQM(x, a, b)
% Quantile function of the Beta distribution
%
% USAGE:
% ======
% inv = betainvIQM(x, a, b)
%
% For each component of 'x', compute the quantile (the inverse of
% the CDF) at 'x' of the Beta distribution with parameters 'a'
% and 'b'.

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>

  if (nargin ~= 3)
    error('Incorrect number of input arguments'); 
  end

  if (~isscalar (a) || ~isscalar(b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
      error ('x, a and b must be of common size or scalars');
    end
  end
  
  sz = size (x);
  inv = zeros (sz);

  k = find ((x < 0) | (x > 1) | ~(a > 0) | ~(b > 0) | isnan (x));
  if (any (k))
    inv (k) = NaN;
  end

  k = find ((x == 1) & (a > 0) & (b > 0));
  if (any (k))
    inv (k) = 1;
  end

  k = find ((x > 0) & (x < 1) & (a > 0) & (b > 0));
  if (any (k))
    if (~isscalar(a) || ~isscalar(b))
      a = a (k);
      b = b (k);
      y = a ./ (a + b);
    else
      y = a / (a + b) * ones (size (k));
    end
    x = x (k);
    l = find (y < eps);
    if (any (l))
      y(l) = sqrt (eps) * ones (length (l), 1);
    end
    l = find (y > 1 - eps);
    if (any (l))
      y(l) = 1 - sqrt (eps) * ones (length (l), 1);
    end

    y_old = y;
    for i = 1 : 10000
      h     = (betacdfIQM(y_old, a, b) - x) ./ betapdfIQM(y_old, a, b);
      y_new = y_old - h;
      ind   = find (y_new <= eps);
      if (any (ind))
        y_new (ind) = y_old (ind) / 10;
      end
      ind = find (y_new >= 1 - eps);
      if (any (ind))
        y_new (ind) = 1 - (1 - y_old (ind)) / 10;
      end
      h = y_old - y_new;
      if (max (abs (h)) < sqrt (eps))
        break;
      end
      y_old = y_new;
    end

    inv (k) = y_new;
  end

end
