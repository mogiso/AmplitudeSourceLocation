program make_stationlist
  use nrtype, only : fp



  call getarg(1, outfile)
  call getarg(2, min_epdist_t)
  call getarg(3, max_epdist_t)



