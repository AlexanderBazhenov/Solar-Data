function z = vee(x, y)
  if (not (isa (x, "kinterval")))
    x = kinterval (x);
  endif
    if (not (isa (y, "kinterval")))
    y = kinterval (y);
  endif
  zinf = min(inf(x), inf(y));
  zsup = max(sup(x), sup(y));

  z = kinterval(zinf, zsup);
endfunction