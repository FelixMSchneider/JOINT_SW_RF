c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: CROSS                                                 c
c                                                                     c
c      COPYRIGHT (C)  1997 R. B. Herrmann                             c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      3507 Laclede Avenue                                            c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
        subroutine cross( ix,  iy, c)
        integer*4 ix, iy
        character c*(*)
        call dfcros(ix, iy, c)
        return
        end
