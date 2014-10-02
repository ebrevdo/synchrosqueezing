% function x = atrend(t, xp, p)
function x = atrend(t, xp, p)
  x = xp + polyval(p, t);
end