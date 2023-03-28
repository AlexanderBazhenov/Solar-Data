function z = wedge(x, y)
  if (not (isa (x, "kinterval")))
    x = kinterval (x);
  endif
    if (not (isa (y, "kinterval")))
    y = kinterval (y);
  endif
  zinf = max(inf(x), inf(y));
  zsup = min(sup(x), sup(y));

  z = kinterval(zinf, zsup);
endfunction