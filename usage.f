      subroutine usage
      write(*,*)'Usage:'
      write(*,*)'      cross_rm -T[max time shift] -W[windows length]
     & -R[reference time] -N[number of control points] -L[fl] -H[fh]
     & file0 file1 ...filen '
      write(*,*)'      -T  maximum time shift to search, default 20.'
      write(*,*)'      -W  time length of the window, default 8.'
      write(*,*)'      -R  reference time e.g. t1..t9 or real number, de
     &fault t1.'
      write(*,*)'      -L  low frequency, default 0.8.'
      write(*,*)'      -H  high frequency, default 2.'
      write(*,*)'      -N  number of control points, default 5'
      write(*,*)'      file0 reference sacfile'
      write(*,*)'      file1..filen other sacfile'
      end subroutine
