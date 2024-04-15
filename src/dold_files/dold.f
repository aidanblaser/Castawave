
      program ww
c==============================================================================
c
c            "ww1.01a"
c
c  A series of routines for following the unsteady motion of a two-dimensional
c  wave-surface, based on a Cauchy-integral boundary formulation of the fluid
c  equations and the use of a high-order approach in space and time.
c
c            This copy is dated:   11th October, 1993
c
c  Two NAG routines "c06faf" and "c06fbf" are used for some purposes in this
c  version, to carry out fast Fourier transforms.  They can, of course, be
c  replaced by any standard numerical routines performing the same tasks.
c==============================================================================
c
c    written by:                   J.W. Dold
c
c             address:             School of Mathematics
c                                  University of Bristol
c                                  Bristol   BS8 1TW
c                                  United Kingdom
c
c          Electronic-mail:        J.W.Dold@bristol.ac.uk
c
c             telephone:         (+44)-(+272)-30-30-30
c              telefax:          (+44)-(+272)-30-33-16
c               telex:             445938 BSUNIV G
c
c    +-------------------------------------------------------------+
c    |                                                             |
c    |          Copyright:         ALL RIGHTS RESERVED             |
c    |                                                             |
c    +-------------------------------------------------------------+
c
c      No commercial exploitation of this program or any part of it or any
c   of its derivatives is allowed except with the written consent of the
c   author.
c
c      Use of the program for purely non-remunerative academic research
c   purposes is welcomed, but please keep the author informed of all use,
c   modifications and any relevant developments.
c
c      Any use whatsoever of the programme will be taken as an indication of
c   acceptance of these conditions.
c
c      There is no guarantee of any kind associated with the use of this
c   software.  However, the programmes making up the package are undergoing
c   a process of development and limited advice concerning problems may be
c   obtainable by contacting the author directly.  Any suggestions and reports
c   of problems will be taken seriously and will help towards producing an
c   effective and widely usable final product.
c
c      The program at the moment appears to be bug-free, but please
c   report any problems identified to the author.
c
c==============================================================================
c
c    Introduction:
c    -------------
c
c  This program follows the unsteady development of a spatially periodic
c  surface wave profile on a two-dimensional inviscid, incompressible and
c  irrotational fluid.
c
c------------------------------------------------------------------------------
c
c  [length], [time] and [mass] scales are taken to be normalised such that:
c     spatial period                  =       wl   x  [length]
c    |gravitational acceleration|     =      gty   x  [length]/[time]/[time]
c     density                         =       1    x  [mass]/[length]/[length]
c  with pressure expressed in units of  [mass]/[length]/[time]/[time].
c
c  The user can provide values for  wl  and  gty.
c
c  The fluid is taken to be bounded below by a flat impermeable horizontal
c  bottom at a specified depth.  The default depth is infinity (bottomless).
c
c  The reference frame is normalised to be the zero mean-current reference
c  frame (i.e. moving with the mean fluid motion at depth).  Because the bottom
c  is flat, there is no generality lost in this and the user is free to add an
c  arbitrary current to the results.
c
c------------------------------------------------------------------------------
c
c  Initial data are required for the horizontal and vertical positions of a
c  finite number of points on the surface and the velocity potential at these
c  points.  Additional data may optionally be provided as outlined later.
c
c  The points MUST be `smoothly' distributed so that a reasonably accurate
c  numerical estimation of derivatives with respect to the "point label"
c  parameter is possible.
c
c  Spatial periodicity of  wl  or  2.py  (the default value) is expected for
c  both the surface elevation and the velocity potential.  The latter amounts
c  to defining the zero mean-current reference frame.
c
c  If no initial data are supplied then the program requests details for
c  the calculation of steadily-propagating waves on given depths.  A least-
c  squares approach is used for this.  [THIS OPTION IS NOW ANULLED---USE "STW"]
c
c------------------------------------------------------------------------------
c
c  Pressures at the surface must be specified to their second Lagrangian time
c  derivative.  The default assumes a constant surface pressure of zero.
c  Any other, for instance non-constant, pressure behaviour can be specified
c  by suitable amendment of the program following entry point "prsure" of
c  subroutine "miscel".
c
c------------------------------------------------------------------------------
c
c  Time step sizes are determined so as to maintain a specified degree of
c  accuracy per period or, still providing this accuracy, until it is possible
c  to step precisely to times at which the surface profile solution is desired.
c  Unless defaults are over-ridden (see later) time-steps are restricted
c  further by a threshold value beyond which strong numerical instabilities
c  are encountered.
c
c  The program normally provides significantly better than the specified
c  degree of accuracy per period (with errors possibly accumulating over
c  many periods), provided the spatial discretisation is "adequate."
c
c  SEE:   Dold, J.W. (1992) Journal of Computational Physics, 103, 90-115.
c
c------------------------------------------------------------------------------
c
c  Output data (to a file called "ww.p," for post-processing) are printed at
c  times which are incremented by a specified value, starting at the initial
c  time, up to a specified final time, at which data are also printed out.
c
c  With "ww.p" renamed to "ww.in," to form the input data, the program will
c  carry out its post-processing offering a menu of types of analyses.  The
c  results of the chosen analysis appear in the file "ww.out" as well as in
c  unformatted data files appropriate to the type of analysis. (See below)
c
c  When the program stops, output is made to a file (called "ww.res")
c  of the same kind and layout as the initial input file.  This allows for
c  subsequent resumption of the processing with possible amendment.
c
c  For resuming processing, one has only to rename "ww.res" to "ww.in",
c  after amendment if desired.  Unless altered, the program finds that the
c  requested final time in this data is not after the initial time and calls
c  for a re-entry or modification of all program control parameters via the
c  terminal.
c
c==============================================================================
c
c  Starting data  (in unformatted records from a file called "ww.in")
c  as follows  (single lines of record between lines, "======";
c       records are double-precision unless stated as "integer"):-
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c  Physical specification of initial motion:
c
c  t = initial time      (normally zero, but optional)
c  n = number of points  (integer)
c  h = depth  (zero for infinite depth, positive referring to distance between
c              bottom and mean level, negative referring to distance between
c              the bottom and the initial minimum level)
c  gty = gravitational acceleration (normally unity, or possibly zero)
c  wl  = spatial period   (defaults to 2.py if zero, |wl|.2.py if negative)
c  uv  = uniform vorticity (horizontal shear) of the fluid
c
c======
c  x1 x2 x3 ... xn       (horizontal position of the surface points)
c======
c  y1 y2 y3 ... yn       (vertical      "     "   "     "      "   )
c======
c  f1 f2 f3 ... fn       (velocity potential at each point)
c======
c
c  Program Control Parameters:
c
c  e = approximate precision per period          (defaults to -0.001 if zero)
c             with   e < 0  =>   relative to initial disturbance amplitude
c                &   e > 0  =>   absolute precision
c  sm = control of smoothing;  0 => no smoothing             (see note below)
c  cs = use routines for conserving mean-level and total energy  (note below)
c  bd = backward differencing & instability threshold action (see note below)
c  pts = time interval between solution printouts            (see note below)
c  tl = time limit for solution  (if tl<=t the program requests a revision of
c                                  the control input parameters from terminal)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c {The following data are optional and, if available, are output to "ww.res".
c  If included in input they provide information for use in initialising the
c  calculated degree of surface "strain" and for providing first estimates
c  of the normal potential gradient solutions for iteration:
c
c======
c  str1 str2 ... strn  (normalisation of surface strain, based on initial data)
c======
c  p1n1 p1n2 ... p1nn  (normal gradient of 'first' potential function)
c======
c  p2n1 p2n2 ... p2nn  (  "       "     "  'second'    "        "    )
c======
c  p3n1 p3n2 ... p3nn  (  "       "     "  'third'     "        "    )      }
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c  Note:  An estimate of the "roughness" of the calculations is made at each
c         time-step by determining the norm of the degree of smoothing which
c         would arise if the THIRD time derivatives of position and potential
c         were smoothed using an 11-point smoothing formula
c
c      :  If  sm  is non-zero then, after each time-step, an |sm|-point
c         smoothing formula is applied to the profile and velocity potential
c         with  |sm|  being 8 or 14 or an odd number between  5  and  15.  If
c         the sign of  sm  is negative then smoothing is applied only if the
c         effect of "roughness" exceeds "|precision|^(1/4)", and then only
c         where the roughness is most significant.
c
c      :  Smoothing is stronger (more radical) the lower the value of  |sm|,
c         and for larger values of |sm| (eg. 15) smoothing is fairly effective
c         at selectively removing two-point or `sawtooth' surface modes (should
c         this be considered necessary) with little loss of other information
c         for a given number of points.  The formulae with |sm| = 8 or 14
c         also remove any 3-point periodic modes as well as the `sawtooth'
c         (or 2-point periodic) modes.
c
c  SEE:   Dold, J.W. (1992) Journal of Computational Physics, 103, 90-115.
c
c----------
c
c  Note:  If  cs < 0,  y  is adjusted to set the initial mean level to zero.
c         At each subsequent time-step:
c         cs  = 0, or  |cs| > 3, has no effect;
c        |cs| = 1, mean-level is restored;
c        |cs| = 2, total energy is restored (very nearly) to its initial value;
c        |cs| = 3, both of the latter.
c
c----------
c
c  Note on bd:     With  0 < |bd| < 6  (default is  |bd| = 5),  |bd|  is the
c                     order of backward differencing to be used;
c       bd >= 0    causes an estimated numerical `strong instability' limit to
c                     be applied as the maximum allowable size of time-step;
c       bd < 0     prevents the application of this restriction;
c       bd = 0     defaults to  bd = 5  (recommended option).
c
c  SEE:   Dold, J.W. (1992) Journal of Computational Physics, 103, 90-115.
c
c----------
c
c  Note:   The number of points  `n'  must not exceed the value of  `NN'  as
c             specified in the `parameter' statements.  If necessary, modify
c             all such values of  `NN'  as well as  `N7'  (which should equal,
c             or exceed, the value  NN + 7) in every `parameter' statement.
c
c----------
c
c  Note: |pts| > 0  demands some time-steps to land on time multiples of |pts|,
c                     causing values of  x, y, f and  pn  to be stored in
c                     "ww.p" for post-processing at these times as well as
c                     the initial and final times;
c        pts < 0    also causes such printout to occur at times of max / min
c                     kinetic energy levels (as well as normal printout times),
c                     useful with studying "standing" waves;
c       otherwise,  printout between the initial and final times is inhibited.
c
c------------------------------------------------------------------------------
c
c  POST PROCESSING:
c
c  Using the file  "ww.p"  as input data for the program (ie running the
c  program with  "ww.p"  renamed to  "ww.in") a menu is displayed for
c  alternative forms of output.  Output data (appearing in  "ww.out") is
c  translated into the form selected.
c
c  The following options are available in this version of the code:
c
c   : A column-wise list of properties at surface points
c   : Fourier modes of  y(x), or of both  y(x)  and  v(x)
c   : Complex amplitude of velocity potential at evenly spaced points moving
c        at the linear group velocity, 1/2.  The carrier wave is expected to
c        have wavelength exactly  2.py  for this to be useful.
c   : Values of  y(x)  at evenly-spaced points (in  x)  moving with the linear
c        group velocity.
c   : Profile and velocity potential values, at surface points
c   : Profile, potential and 1st time derivatives of surface particles
c   : Profile, potential, 1st and 2nd derivatives
c   : Profile, potential, 1st, 2nd and 3rd derivatives
c
c==============================================================================
c
      implicit double precision (a-h,o-z)
      parameter (NN=4100, N7=NN+7)
      common /profil/ x(-6:N7),y(-6:N7),f(-6:N7),
     *                gty,h,t,ym,enk,enp,ent,en0,dt,dtl,
     *                n,ndtn,uv, eps, BW,lbd,mbd
      common /auxlry/ ft(-6:N7),sf(-6:N7),sft(-6:N7),sft2(-6:N7)
      common /scalng/ wl,rx,tx,fx,xr,xt,xf
      common /drvtvs/   u(-6:N7),  v(-6:N7), fd(-6:N7),
     *                 ud(-6:N7), vd(-6:N7),fd2(-6:N7),
     *                ud2(-6:N7),vd2(-6:N7),fd3(-6:N7),
     *                udb(NN,5), vdb(NN,5), fdb(NN,5),
     *                ubd(NN,5), vbd(NN,5), fbd(NN,5),
     *      ut(-6:N7), vt(-6:N7),ut2(-6:N7),vt2(-6:N7),
     *      ux(-6:N7), vx(-6:N7),utx(-6:N7),vtx(-6:N7),
     *     uxx(-6:N7),vxx(-6:N7)
      common /pressr/ p0(-6:N7),pd1(-6:N7),pd2(-6:N7)
      common /potntl/ p(-6:N7),ps(-6:N7),pn(-6:N7),
     *                pno(-6:N7,3),pnl(-6:N7,3),k(3),lbp
      common /matrix/ a(NN,NN),b(NN,NN),
     *                zwa(-6:N7),zwr(-6:N7),zwi(-6:N7)
      common /contrl/ e,erp,lsm,lcs,istop,lsym
      common /rufnes/ err(-6:N7),erm
      common /print/  pt,pts,tl,mpt,lpt,ig
      common /propty/ slp(-6:N7),curv(-6:N7),angvel(-6:N7),
     *                strn(-6:N7),rstr(-6:N7),astr(-6:N7),
     *                Prn(-6:N7),strn0(-6:N7),erl(-6:N7)
      common /utilty/ s(-6:N7),d(-6:N7)
      common /lbldrv/ x1(-6:N7),y1(-6:N7),r1(-6:N7),
     *                x2(-6:N7),y2(-6:N7),
     *                u1(-6:N7),v1(-6:N7),
     *                ut1(-6:N7),vt1(-6:N7),
     *                ux1(-6:N7),vx1(-6:N7)
      logical uflag
      common /under/ uflag
c------------------------------------------------------------------------------
      uflag=.false.
      call input
      if (n.gt.0) then
        call begin
    1     call output    
          call timstp
          call motion        
        go to 1
      else if (n.lt.0) then
        call postps
        call inpost
    2     call output
          call nextps
          call motion
        go to 2
      else
        write(*,*) 'Input data required in "ww.in"'
      end if
      end
c==============================================================================
      subroutine motion
      implicit double precision (a-h,o-z)
      integer ii
      parameter (NN=4100, N7=NN+7, py2=.62831853071795864769d1)
      common /profil/ x(-6:N7),y(-6:N7),f(-6:N7),
     *                gty,h,t,ym,enk,enp,ent,en0,dt,dtl,
     *                n,ndtn,uv,lbd,mbd
      common /auxlry/ ft(-6:N7),sf(-6:N7),sft(-6:N7),sft2(-6:N7)
      common /scalng/ wl,rx,tx,fx,xr,xt,xf
      common /drvtvs/   u(-6:N7),  v(-6:N7), fd(-6:N7),
     *                 ud(-6:N7), vd(-6:N7),fd2(-6:N7),
     *                ud2(-6:N7),vd2(-6:N7),fd3(-6:N7),
     *                udb(NN,5), vdb(NN,5), fdb(NN,5),
     *                ubd(NN,5), vbd(NN,5), fbd(NN,5),
     *      ut(-6:N7), vt(-6:N7),ut2(-6:N7),vt2(-6:N7),
     *      ux(-6:N7), vx(-6:N7),utx(-6:N7),vtx(-6:N7),
     *     uxx(-6:N7),vxx(-6:N7)
      common /pressr/ p0(-6:N7),pd1(-6:N7),pd2(-6:N7)
      common /potntl/ p(-6:N7),ps(-6:N7),pn(-6:N7),
     *                pno(-6:N7,3),pnl(-6:N7,3),k(3),lbp
      common /matrix/ a(NN,NN),b(NN,NN),
     *                zwa(-6:N7),zwr(-6:N7),zwi(-6:N7)
      common /contrl/ e,erp,lsm,lcs,istop,lsym
      common /rufnes/ err(-6:N7),erm
      common /print/  pt,pts,tl,mpt,lpt,ig
      common /propty/ slp(-6:N7),curv(-6:N7),angvel(-6:N7),
     *                strn(-6:N7),rstr(-6:N7),astr(-6:N7),
     *                Prn(-6:N7),strn0(-6:N7),erl(-6:N7)
      common /utilty/ s(-6:N7),d(-6:N7)
      common /lbldrv/ x1(-6:N7),y1(-6:N7),r1(-6:N7),
     *                x2(-6:N7),y2(-6:N7),
     *                u1(-6:N7),v1(-6:N7),
     *                ut1(-6:N7),vt1(-6:N7),
     *                ux1(-6:N7),vx1(-6:N7)
      CHARACTER(LEN=20) :: fname
      logical uflag
      common /under/ uflag
      character*3 nan
      save
c==============================================================================
c  calculating time-derivatives at current time, for profile  <x,y,f>
c------------------------------------------------------------------------------
c  calculation or re-setting of mean-level
      call meany
c  set up matrices for performing Cauchy integral solution in transformed plane
      call kernel
c  surface pressure distributions
      call prsure
c  Lagrangian time derivatives of x, y and f
      call times(ltm,e)
c  11-point roughness on the third derivatives calculated
      if (ltm.ge.0.or.ltm.eq.-3)
     *  call rough(11,ud2,vd2,fd3,n,0.d0,err,erm)
c  calculation of energy
      if (ltm.ge.0) call energy
C      print*, t
      return
c============================================================================
      entry output
C      ii=0
      if (isc.eq.10.or.mpt.eq.1.or.istop.gt.0) then
        if (ltm.lt.0) then
          call ddi(x,n,py2,x1)
          call ddi(y,n,0.d0,y1)
        end if
        call slope
        call hilow(slp,n,slph,slpl)
        call hilow(y,n,yh,yl)
C        write(*,*) xtf*t
C        write(*,90) xtf*t,ndtn,log10(xr*erm+1.002d-10),xr*yl,
C     *    1.d2*(ym-ymo)/yho,xr*yh,slpl,slph,lth*ent/en0,1.d2*enk/ent
c        write(*,90) xtf*t,ndtn,log10(xr*erm+1.002d-10),xr*yl,xr*ym,
c     *              xr*yh,slpl,slph,lth*ent/en0,100.d0*enk/ent
        rewind(11)
        write(11,'(f4.1)') xr*erm
        rewind(11)
C       print*, erm
        read(11,'(a4)') nan
C        print*, nan
C       if (nan.eq.'NaN'
        if (abs(erm).gt.5.d8) stop 'data is no longer intelligible'
        if (nan.eq.'NaN') stop 'data is no longer intelligible'
        isc = 0 
C       ii=ii+1
C      print*,  ii
C       print*, y(1:96)
C       WRITE(fname,'(A,I4.4,A)') 'y.txt'
C       OPEN(11,file=fname,form='formatted')
C       WRITE(11,*) (y(i),i=1,n)
C        open(unit=22,file='y.txt',
C     &  status = 'unknown')
C       write(22,*)(y)
       
C       WRITE(fname,'(A,I4.4,A)') 'x.txt'
C       OPEN(21,file=fname,form='formatted')
C       WRITE(21,*) (x(i),i=1,n)
C      if (mod(ii,5).eq.0) then

C      WRITE(fname,'(A,F0.2,F0.3,A)') 'y2048_s.txt'
C      OPEN(81,file=fname,form='formatted')
C      WRITE(81,*) (y(i),i=1,n)
C      WRITE(fname,'(A,F0.2,F0.3,A)') 'x2048_s.txt'
C      OPEN(71,file=fname,form='formatted')
C      WRITE(71,*) (x(i),i=1,n)
C      WRITE(fname,'(A,F0.2,F0.3,A)') 'phi2048_s.txt'
C      OPEN(61,file=fname,form='formatted')
C      WRITE(61,*) (f(i),i=1,n)        
C 
C      tout=t 
C      WRITE(fname,'(A,F0.2,F0.3,A)') 't1024.txt'
C      OPEN(52,file=fname,form='formatted')
C      WRITE(52,*) t

C      print*, n
C      if (mod(ii,10).eq.0) then
       BW=BW
       WRITE(fname,'(A,F0.4,F0.3,F0.3,A)') 'y.txt'
       OPEN(81,file=fname,form='formatted')
       WRITE(81,*) (y(i),i=1,n)
       WRITE(fname,'(A,F0.4,F0.3,F0.3,A)') 'x.txt'
       OPEN(71,file=fname,form='formatted')
       WRITE(71,*) (x(i),i=1,n)
       WRITE(fname,'(A,F0.4,F0.3,F0.3,A)') 'phi.txt'
       OPEN(61,file=fname,form='formatted')
       WRITE(61,*) (f(i),i=1,n)        
C      print*, ii   
C      print*, t
       WRITE(fname,'(A,F0.4,F0.3,F0.3,A)') 't.txt'
       OPEN(59,file=fname,form='formatted')
       WRITE(59,*) (t,i=1,n)

       WRITE(fname,'(A,F0.4,F0.3,F0.3,A)') 'u.txt'
       OPEN(215,file=fname,form='formatted')
       WRITE(215,*) (u(i),i=1,n)
C      
       WRITE(fname,'(A,F0.4,F0.3,F0.3,A)') 'enk.txt'
       OPEN(34,file=fname,form='formatted')
       WRITE(34,*) enk
C      CLOSE(UNIT=11)
       WRITE(fname,'(A,F0.4,F0.3,F0.3,A)') 'enp.txt'
       OPEN(24,file=fname,form='formatted')
       WRITE(24,*) enp
C      CLOSE(UNIT=11) 
        
C      end if
C      print*, t
C      CLOSE(UNIT=11)    
C      WRITE(fname,'(A,I4.4,A)') 'y_',ii,'.txt'
C      OPEN(81,file=fname,form='formatted')
C      WRITE(81,*) (y(i),i=1,n)
CC      CLOSE(UNIT=11) 
C      WRITE(fname,'(A,I4.4,A)') 'x_',ii,'.txt'
C      OPEN(71,file=fname,form='formatted')
C      WRITE(71,*) (x(i),i=1,n)
CC      CLOSE(UNIT=11)  
C      WRITE(fname,'(A,I4.4,A)') 'phi_',ii,'.txt'
C      OPEN(61,file=fname,form='formatted')
C      WRITE(61,*) (f(i),i=1,n)
CC      CLOSE(UNIT=11)         
CCC       print*, t    
C      WRITE(fname,'(A,I4.4,A)') 't_',ii,'.txt'
C      OPEN(51,file=fname,form='formatted')
C      WRITE(51,*) t
CC      CLOSE(UNIT=11)    
C       WRITE(fname,'(A,I4.4,A)') 'enk_',ii,'.txt'
C       OPEN(11,file=fname,form='formatted')
C       WRITE(11,*) enk
C       CLOSE(UNIT=11)
C        WRITE(fname,'(A,I4.4,A)') 'enp_',ii,'.txt'
C       OPEN(11,file=fname,form='formatted')
C       WRITE(11,*) enp
C       CLOSE(UNIT=11)      
C       WRITE(fname,'(A,I4.4,A)') 'u_',ii,'.txt'
C       OPEN(11,file=fname,form='formatted')
C       WRITE(11,*) (u(i),i=1,n)
C       CLOSE(UNIT=11)  
C       WRITE(fname,'(A,I4.4,A)') 'v_',ii,'.txt'
C       OPEN(11,file=fname,form='formatted')
C       WRITE(11,*) (v(i),i=1,n)
C       CLOSE(UNIT=11)        
        
      end if
      isc = isc+1
c  control of w
c  printout and termination
c----------------------------------------------------------------------------
c  mpt = 1 or istop > 0  trigger relevant printout (see "pts" in introduction)
c  istop > 0 forces a termination
c------------------------------------------------------------------------------
      if (mpt.eq.1.or.istop.gt.0) then
        if (lpt.ne.1.and.lpt.ne.2) then
c  post-processing output:
       if(uflag) then
         if(lpt.eq.3) call wrudat
       else
          write(8,91)
          write(8,92) xtf*t,ndtn,log10(xr*erm+1.002d-10),xr*ym,
     *            ent/en0,100.*enk/ent,xr*yl,xr*yh,slpl,slph
          if (abs(lpt).gt.3) call xfixed
          if (lpt.eq.3) call column
          endif
          if (lpt.le.0.and.lpt.ge.-3) then
            call prtxyf
            if (lpt.le.-1) call prLagr
            write(8,97)
          end if
          if (ltm.lt.0) return
        else
c  normal time-stepping operation:
          i = abs(ndtn)
          write(7) t,n,i,real(erm),real(enk),real(enp)
          write(7) (real(x(i)),i=1,n)
          write(7) (real(y(i)),i=1,n)
          write(7) (real(f(i)),i=1,n)
C          print*, t
          do 1 i=1,3
    1     write(7) (real(pno(j,i)),j=1,n)
        end if
        if (t.gt.tl-e.or.istop.gt.0) then
          open(9,file='ww.res',form='unformatted')
          rewind(9)
          call resump
        else
          mpt = 0
          return
        end if
      else
        return
      end if
      close(7)
      if(uflag) then
       close(11)
      else
       close(8)
      endif
      close(9)
      stop '<information for resumption is stored in ww.res>'
c==============================================================================
      entry timstp
        call bakxyf
c  time-step size determined from results of  "sizedt"
        if (lbd.ge.2) then
c  calculated from  3rd  and  4th  time derivatives
c  with quadratic or higher orders of backward differencing:
          do 3 i=1,n
            vtx(i) = ubd(i,1)
            uxx(i) = vbd(i,1)
    3       vxx(i) = fbd(i,1)
          call sizedt(4,vtx,uxx,vtx,dt1)
          call sizedt(3,ud2,vd2,fd3,dt2)
        else
c  or  2nd  and  3rd  derivatives with linear or zero differencing:
          call sizedt(3,ud2,vd2,fd3,dt1)
          call sizedt(2,ud,vd,fd2,dt2)
          if (lbd.lt.1) then
c  However, for zero backward differencing,  2nd  and  3rd derivatives
c      are used but the resulting step-size is halved:
            dt1 = .5d0*dt1
            dt2 = .5d0*dt2
c  [ normally this is necessary for the first time-step only ]
          end if
        end if
        if (dt1.gt.dt2)  then
c  one should also usually expect that  dt1 > dt2,
c  in which case take the geometric mean time-step:
          dtn = sqrt(dt1*dt2)
          ndtn = abs(ndtn)
        else
c  otherwise, something is likely to be wrong, so make a reduced time-step:
          dtn = .5d0*dt1
c  flag for this:
          ndtn = -abs(ndtn)
        end if
c  restricting the time-step size to below a `strong' instability
c  threshold for infinitesimal disturbances on a flat surface:
c  This is (theoretically) proportional to  min{sqrt(|dR/di|/Pn)},
        dts = 0.d0
        do 5 i=1,n
c  determine  max{ Pn/|dR/di| }   :   Pn = normal pressure gradient,
          tp = abs(x1(i)*(vd(i)+gty)-y1(i)*(ud(i)+uv*v(i)))
     *           /(x1(i)*x1(i)+y1(i)*y1(i))
    5     if (dts.lt.tp) dts = tp
        tp = sqrt(dts)
c  conservatively estimated threshold proportionality constants,
c  [ depending on the order of backward differencing ]:
        if (lbd.eq.5) then
          dts = .6d0/tp
        else if (lbd.eq.2) then
          dts = 1.1d0/tp
        else if (lbd.eq.1) then
          dts = 1.2d0/tp
        else if (lbd.eq.0) then
          dts = 1.d0/tp
        else if (lbd.eq.3) then
          dts = 1.d0/tp
        else if (lbd.eq.4) then
          dts = .9d0/tp
        else
          dts = .6d0/tp
        end if
c  restriction below the strong instability threshold:
        if (dtn.gt.dts) then
          if (mbd.gt.0) dtn = dts
          lth = -1
        else
          lth = 1
        end if
c  a more suitable value of dtn may be defined by "nicedt"
        call nicedt(dtn)
        call nextpn(dtn)
c  storing time-step values:
        dtl = dt
        dt = dtn
c  time-stepping to next set of values of x, y and f
        call Taylor
        if (lsm.gt.0) call smxyf
        if (lsm.lt.0.and.(e+erm)**4.gt.e) call smxyf
      return
c==============================================================================
      entry input
c----------------------------------------
c  section to read in initial data
c------------------------------------------------------------------------------
c  starting flags:
      n = 0
      ltm = -999
      ig = 0
      lth = 1
        open(11,form='formatted',status='scratch')
        open(7,file='ww.in',form='unformatted',
     &   status='old',err=20)
        rewind(7)
        read(7,end=20,err=20) t,n,h,gty,wl,uv, eps, BW
C        print*,n
      if (n.lt.0) read(7,end=20,err=20) erp,sm,cs,pts,tl
      if (n.le.0) return
        read(7,end=20,err=20) (x(i),i=1,n)
        read(7,end=20,err=20) (y(i),i=1,n)
        read(7,end=20,err=20) (f(i),i=1,n)
        read(7,end=20,err=20) erp,sm,cs,bd,pts,tl
C        print*, x
c  flags for normal time-stepping operation:
        ltm = 0
c  maximum order of backward differencing set by  |bd|
        mbd = nint(bd)
      if (tl.le.t+1.d-6) then
c  test viability of time-limit and correct if necessary
        write(*,*) 'time-limit, ',tl,' is not after initial time,',t
        write(*,*) 'enter new initial time and time-limit'
        read(*,*) t,tl
        if (tl.le.t+1d-6) istop = 1
        write(*,*) 'printout time-step  (<0 => ',pts,' = previous):'
        read(*,*) dt
        if (dt.ge.0.d0) pts = dt
        write(*,*) 'precision, smoothing-type, correction-type,',
     *             ' order of backward differencing, depth :'
        write(*,*)  erp,', ',lsm,', ',lcs,', ',mbd,', ',h
        write(*,*) 'enter new values:'
        read(*,*)   erp,lsm,lcs,mbd,h
      end if
c  defaults on backward differencing:
        if (mbd.eq.0.or.mbd.gt.5) mbd = 5
        if (mbd.lt.-5) mbd = -5
        if (pts.ge.0.d0) then
          pt = 1
        else
c  if,  pts < 0,  flag for standing-wave behaviour, producing
c  additional output at times of max/min kinetic energy:
          pt = 2
        end if
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  post-processing initial entry point:
      entry inpost
      lsm = nint(sm)
      lcs = nint(cs)
      lpt = nint(pt)
      if (ltm.lt.0) then
c  post-processing options:
        if (lpt.eq.3) ltm = -3
        if (lpt.lt.0.and.lpt.ge.-3) ltm = lpt
        if (lpt.gt.3.and.lpt.lt.2000) ltm = -1
        if (lpt.lt.-1000) ltm = -1
      end if
      ak = abs(pt-lpt)
      istop = 0
      xtf = 1.d0
      if (abs(lpt).gt.3) xtf = 1.d0/py2
c  defaults for non-positive  wl
      if (wl.eq.0.d0) wl = py2
      if (wl.lt.0.d0) wl = -wl*py2
c  rescalings if spatial period, wl, is non-zero so that internally the program
c  can take the period as  2.py  while all output etc. is correctly scaled
        xtf = 1.d0/py2
        rx = py2/wl
        tx = sqrt(rx)
        fx = rx*tx
        if (ltm.ge.0) then
          do 7 i=1,n
c  scalings:
            x(i) = rx*x(i)
            y(i) = rx*y(i)
    7       f(i) = fx*f(i)
c  also preset higher derivatives to zero:
          do 8 j=1,5
            do 8 i=1,n
              ubd(i,j) = 0.d0
              vbd(i,j) = 0.d0
    8         fbd(i,j) = 0.d0
          t = tx*t
          h = rx*h
          pts = tx*pts
          tl = tx*tl
        end if
        xr = 1.d0/rx
        xt = 1.d0/tx
        xf = 1.d0/fx
        spd = wl
        xtf = xtf*xt
c  operate defaults on input data
      if (abs(abs(lsm)-11).ne.3) lsm = 1+2*nint((lsm-.9)*.5)
      if (abs(lsm).gt.15.or.abs(lsm).lt.5) lsm = 0
      if (lsm.gt.0) call smxyf
      e = erp
      if (e.eq.0.) e=-.001d0
      if (e.gt.0.and.wl.gt.0) e=rx*e
c  set non-input control variables
c  flag for non-symmetry:
      lsym = 0
      pt = t+abs(pts)
      if (pt.ge.tl) pt = tl
c  flag to trigger initial printout
      mpt = 1
c  lbp acts as a flag: -1  no normal potential gradients estimated,
c                       1  linear estimate possible for next gradient.
      lbp = -1
c  similarly,  lbd = 0  for no estimates for backward-differencing,
c                    1  sufficient previous data for linear estimate,
c         up  to     5  5th-order estimate possible.
c (correct startup with -1)
      lbd = -1
c  also, whether  dt  or  dtl  are non-zero serves a similar purpose
      dt = 0.d0
      dtl = 0.d0
      if (ltm.lt.0) go to 13
c - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  Normal time-stepping:
c  read initial surface-strain normalisation if available, otherwise this is
c  defined so that the initial strain is uniformly unity
      read(7,err=16,end=16) (strn0(i),i=1,n)
      do 9 i=1,n
    9   strn0(i) = rx*strn0(i)
      go to 18
   16 call ddi(x,n,py2,x1)
      call ddi(y,n,0.d0,y1)
      do 17 i=1,n
   17   strn0(i) = sqrt(x1(i)*x1(i)+y1(i)*y1(i))
      go to 12
c  read remaining input data if available
   18 do 19 j=1,3
   19   read(7,err=12,end=12) (pno(i,j),i=1,n)
   12 close(7,err=20)
      if (lpt.eq.1.or.lpt.eq.2) then
c  open post-processing data-file:
        open(7,file='ww.p',form='unformatted')
        rewind(7)
        write(7) t,-n,h,gty,wl,uv
        write(7) erp,sm,cs,pts,tl
        write(7) (real(strn0(i)),i=1,n)
      end if
      return
c==============================================================================
   13   lbp = 3
      entry begin
c  check surface data is non-trivial and reasonably sensible
c  or whether symmetry persists about  i = 1
      x(n+1) = x(1)+py2
      xc = x(n+1)+x(1)
      y(n+1) = y(1)
      f(n+1) = f(1)
      tpx = 0.d0
      tpy = 0.d0
      tpf = 0.d0
      x(0) = x(n)-py2
      n2 = n+2
      do 15 i=1,n+1
        j = n2-i
        tpx = tpx+abs(x(i)+x(j)-xc)
        tpy = tpy+abs(y(i)-y(j))
        tpf = tpf+abs(f(i)-f(j))
   15   if (x(i).eq.x(i-1).and.y(i).eq.y(i-1))
     *    stop '<no difference between successive points>'
c  identification of symmetry:
      if (tpx+tpy+tpf.lt.1.d-6.and.2*(n/2).eq.n) then
        write(*,*) 'even symmetry in input data'
        lsym = 1
        call setsym
      end if
c  examine range of data and check consistency etc.
      call hilow(y,n,yh,yl)
      call hilow(f,n,fh,fl)
      fh = .5d0*(fh-fl)
      yh = .5d0*(yh-yl)
      fh = sqrt(fh*fh+yh*yh)
      if (fh.eq.0.d0) stop '<all is totally calm>'
      yho = yh-yl
      if (yho.le.0.d0) yho = 2.d0*fh
c  For  e<0  to be an error estimate RELATIVE to the order of magnitude of the
c  overall variations in y and f, it is now normalised with respect to the
c  initial magnitudes of these variations
      if (e.lt.0)  e = -e*fh
      tol = xr*e
      if (lcs.lt.0) then
        ym = 0.d0
        call meany
        lcs = -lcs
      else
        call findym
      end if
      ymo = ym
      if (lcs.gt.3) lcs = 0
      if (h.gt.0..and.h-(ym-yl).le.e) stop '<surface touches bottom>'
c  if negative set  h  to its positive definition
      if (h.lt.0.) h = ym-(yl+h)
c  water is effectively "deep" if  exp(-2*(h+fh)) < e*e  (see "transh")
      if (exp(-(h+fh)).lt.e) then
        write(*,*) 'water is  effectively "deep"'
        h = 0.d0
      end if
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  set the transformation of the bottom in subroutine "kernel"
      call transh
      if (ltm.ne.-999) then
c  calculation of initial motion
        call kernel
        call prsure
        call times(ltm,e)
c  calculation of initial energy
        if (ltm.ge.0) call seteng
c  initialise "roughness" and set time-step counter
        if (ltm.ge.0.or.ltm.eq.-3)
     *    call rough(11,ud2,vd2,fd3,n,0.d0,err,erm)
        if (ltm.ge.0) ndtn = 0
      end if
      if (lpt.ne.1.and.lpt.ne.2) then
c  prepare for output of post-processing
       if(uflag) then
       open(12,file='ww.und',form='unformatted')
       rewind(12)
        write(12) n,real(h),real(gty),real(omg),real(wl)
        print*,'n,h,gty,omg,wl',n,h,gty,omg,wl
       else
        open(8,file='ww.out',form='formatted',err=21)
        rewind(8)
        if (abs(ig).eq.1) then
          open(10,file='ww.g',form='unformatted')
          rewind(10)
          write(10) real(spd)
        end if
c  heading for output
        if (h.gt.0.) then
          write(8,95) xr*h,gty,spd,xr*ym,xr*ent,n,tol,lsm,lcs
        else
          write(8,96) gty,spd,xr*ym,xr*ent,n,tol,lsm,lcs
        end if
        write(8,97)
       end if
      end if
      if (abs(lpt).gt.3) call fixtyp(ak)
      if (lpt.eq.1.or.lpt.eq.2) write(7) real(en0)
      return
c------------------------------------------------------------------------------
   20 if (n.gt.0) stop '<input data file incomplete>'
      if (n.lt.0) stop '<post-processing data incomplete>'
      return
   21 stop '<problem concerning output file - ww.out>'
c------------------------------------------------------------------------------
   90 format(f10.5,i6,f5.2,f8.5,f8.4,f8.5,f7.2,f6.2,f11.8,f9.5)
c   90 format(f9.5,'[',i6,f5.2,']',f8.5,f9.7,f7.5,'S',f7.2,f6.2,
c     *       '|',f10.7,f6.3,'%K')
   91 format(' ',130('-'))
   92 format('time =',f14.6,' (',i6,' ;',f6.3,')',
     *       ' <y>:',f9.6,', E:',f10.7,' [',f7.4,'% KE]  y:',
     *       f10.7,',',f10.7,'; slope (',f8.3,',',f8.3,')')
   95 format('depth ',f9.6,', gty ',f9.6,', x-period ',
     *  f9.6,', y-mean ',f9.6,' energy ',f13.9,'  [',i3,' pts, tol ',
     *  e8.2,', smooth:',i3,'.',i1,']')
   96 format(' infinite depth, gty ',f9.6,', x-period ',
     *  f9.6,', y-mean ',f9.6,', energy',f13.9,'  [',i3,' pts, tol ',
     *  e8.2,', smooth:',i3,'.',i1,']')
   97 format(' ',130('='))
      end
c==============================================================================
      subroutine kernel
c------------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (NN=4100, N7=NN+7)
      common /profil/ x(-6:N7),y(-6:N7),f(-6:N7),
     *                gty,h,t,ym,enk,enp,ent,en0,dt,dtl,
     *                n,ndtn,uv,lbd,mbd
      common /auxlry/ ft(-6:N7),sf(-6:N7),sft(-6:N7),sft2(-6:N7)
      common /scalng/ wl,rx,tx,fx,xr,xt,xf
      common /drvtvs/   u(-6:N7),  v(-6:N7), fd(-6:N7),
     *                 ud(-6:N7), vd(-6:N7),fd2(-6:N7),
     *                ud2(-6:N7),vd2(-6:N7),fd3(-6:N7),
     *                udb(NN,5), vdb(NN,5), fdb(NN,5),
     *                ubd(NN,5), vbd(NN,5), fbd(NN,5),
     *      ut(-6:N7), vt(-6:N7),ut2(-6:N7),vt2(-6:N7),
     *      ux(-6:N7), vx(-6:N7),utx(-6:N7),vtx(-6:N7),
     *     uxx(-6:N7),vxx(-6:N7)
      common /pressr/ p0(-6:N7),pd1(-6:N7),pd2(-6:N7)
      common /potntl/ p(-6:N7),ps(-6:N7),pn(-6:N7),
     *                pno(-6:N7,3),pnl(-6:N7,3),k(3),lbp
      common /matrix/ a(NN,NN),b(NN,NN),
     *                zwa(-6:N7),zwr(-6:N7),zwi(-6:N7)
      common /contrl/ e,erp,lsm,lcs,istop,lsym
      common /rufnes/ err(-6:N7),erm
      common /print/  pt,pts,tl,mpt,lpt,ig
      common /propty/ slp(-6:N7),curv(-6:N7),angvel(-6:N7),
     *                strn(-6:N7),rstr(-6:N7),astr(-6:N7),
     *                Prn(-6:N7),strn0(-6:N7),erl(-6:N7)
      common /utilty/ s(-6:N7),d(-6:N7)
      common /lbldrv/ x1(-6:N7),y1(-6:N7),r1(-6:N7),
     *                x2(-6:N7),y2(-6:N7),
     *                u1(-6:N7),v1(-6:N7),
     *                ut1(-6:N7),vt1(-6:N7),
     *                ux1(-6:N7),vx1(-6:N7)
      logical uflag
      common /under/ uflag
      save
c------------------------------------------------------------------------------
c  set up transformation  ( w  =  u + i.v = exp ( y - i.x ) )
      do 1 i=1,n
        ey = exp(y(i))
        u(i) =  ey*cos(x(i))
    1   v(i) = -ey*sin(x(i))
c  calculation of first and second derivatives w.r.t. integer surface
c  parameter:  u1 + i.v1,  and  s + i.d
      call ddi(u,n,0.d0,u1)
      call d2di(u,n,0.d0,s)
      call ddi(v,n,0.d0,v1)
      call d2di(v,n,0.d0,d)
      do 3 i=1,n
c  factors in transforming to real gradients of p
        cs = u1(i)*u1(i)+v1(i)*v1(i)
        zwr(i) = (v(i)*u1(i)-u(i)*v1(i))/cs
        zwi(i) = (u(i)*u1(i)+v(i)*v1(i))/cs
c  absolute value of transformation factor:
    3   zwa(i) = sqrt(zwr(i)*zwr(i)+zwi(i)*zwi(i))
c  setting up matrices a(i,j) & b(i,j),  (n x n),
c  embodying the integration kernels from Cauchy's theorem for operating on
c  periodic functions, ie gradients of periodic harmonic potential functions:
c (w, w1 and w2 represent complex vectors u+i.v, u1+i.v1 and s+i.d
c         in the explanations of the calculations below)
c  components of  w1(i)/(w(i)-w(j))  :  i.ne.j
      do 5 j=1,n
        do 4 i=1,j-1
c  division by  (w(i)-w(j))
          cr = u(i)-u(j)
          ci = v(i)-v(j)
          cs = cr*cr+ci*ci
          cr = cr/cs
          ci = ci/cs
c  real and imaginary parts of  w1(i)/(w(i)-w(j))  =  b(i,j) + i.a(i,j)
          a(i,j) = v1(i)*cr-u1(i)*ci
          b(i,j) = u1(i)*cr+v1(i)*ci
c  real and imaginary parts of  w1(j)/(w(j)-w(i))  =  b(j,i) + i.a(j,i)
          a(j,i) =  u1(j)*ci-v1(j)*cr
    4     b(j,i) = -u1(j)*cr-v1(j)*ci
c  diagonal terms after accounting for the singularity:  .5*w2(j)/w1(j)
        cs = u1(j)*u1(j)+v1(j)*v1(j)
        cs = cs+cs
        a(j,j) = (d(j)*u1(j)-s(j)*v1(j))/cs
    5   b(j,j) = (s(j)*u1(j)+d(j)*v1(j))/cs
c--------------------------------------------------------------------------
c  calculation of contribution from flat bottom at y = -h,
c              if water is not "deep":
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  Older form which is discontinuous with "deep" water formulation:
c (not used now because of poor convergence---in fact divergence---in
c  iterations, although `solution' is correct, ie. same as continuous form)
c      if (g.ne.0.d0) then
cc  components of  w1(i)/(w(i)-g/wc(j))   :  c => complex conjugate,
c        do 6 j=1,n
cc  evaluation of g/wc(j)
c          gs = g/(u(j)*u(j)+v(j)*v(j))
c          gr = gs*u(j)
c          gi = gs*v(j)
c          do 6 i=1,n
cc  (whole column, no singularities to consider)
cc  division by w(i)-g/wc(j)
c            cr = u(i)-gr
c            ci = v(i)-gi
c            cs = cr*cr+ci*ci
cc  real and im. parts of  w1(i)/(w(i)-g/wc(j))  =  -b(i,j) + i.a(i,j)
c            a(i,j) = a(i,j) + (v1(i)*cr-u1(i)*ci)/cs
c    6       b(i,j) = b(i,j) - (u1(i)*cr+v1(i)*ci)/cs
c      end if
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  More useful form which is continuous with "deep" water formulation:
c (from,   w1(i)/(w(i)-g/w(i)(j)) - w1(i)/w(i),  which  -> 0  as   g -> 0;
c  results are the same, at least for zero mean current,  ie. f, periodic)
      if (g.ne.0.d0) then
c                    w1(i) . g/wc(j)
c  components of   -------------------   :  c => complex conjugate,
c                  w(i).(w(i)-g/wc(j))
        do 6 j=1,n
c  evaluation of g/wc(j)
          gs = g/(u(j)*u(j)+v(j)*v(j))
          gr = gs*u(j)
          gi = gs*v(j)
          do 6 i=1,n
c  (whole column, no singularities to consider)
c  evaluation of  w1(i).g/wc(j)
            wgr = u1(i)*gr-v1(i)*gi
            wgi = v1(i)*gr+u1(i)*gi
c  division by  w(i).(w(i)-g/wc(j))
            ug = u(i)-gr
            vg = v(i)-gi
            cr = u(i)*ug-v(i)*vg
            ci = v(i)*ug+u(i)*vg
            cs = 1./(cr*cr+ci*ci)
c                            w1(i) . g/wc(j)
c  real and im. parts of   -------------------   =  -b(i,j) + i.a(i,j)
c                          w(i).(w(i)-g/wc(j))
            a(i,j) = a(i,j) + (wgi*cr-wgr*ci)*cs
    6       b(i,j) = b(i,j) - (wgr*cr+wgi*ci)*cs
      end if
c  symmetry :
      if (lsym.eq.1) then
c  case of even symmetry (surface and potentials symmetric about  i = 1 )
        n2 = n+2
        do 69 i=2,n/2
          n2i = n2-i
          do 69 j=1,n/2+1
             a(j,i) = a(j,i)+a(j,n2i)
   69        b(j,i) = b(j,i)-b(j,n2i)
      else if (lsym.eq.-1) then
c  case of uneven symmetry (surface symmetric and potentials anti-symmetric
c  about  i = 1 )
        n2 = n+2
        do 96 i=2,n/2
          n2i = n2-i
          do 96 j=1,n/2+1
             a(j,i) = a(j,i)-a(j,n2i)
   96        b(j,i) = b(j,i)+b(j,n2i)
      end if
      return
c==============================================================================
      entry transh
c  "deep" water (giving g = 0) is signified by h = 0
c  (bottom  :  u*u + v*v = g  at  y = ym-h  i.e. depth relative to mean level)
      if (h.gt.0.) then
        g = exp(2.d0*(ym-h))
      else
        g = 0.d0
      end if
c  label to signify no previous computation of normal potential gradients
      if (lbp.lt.0) then
        lpn = 0
      else
        lpn = 3
      end if
      return
c==============================================================================
      entry derive(l,er)
c------------------------------------------------------------------------
c  if it exists, set the projected estimate of the normal gradient term
      if (lpn.ge.l) then
        do 8 i=1,n
    8     pn(i) = pno(i,l)
      end if
c  solve for potential gradient using Cauchy's integral theorem
      call tangnt
      call normal(k(l),er)
      do 9 i=1,n
c  record normal gradient term for projection to next time-step
    9 pno(i,l) = pn(i)
      if (lpn.lt.l) lpn = l
c  velocity (l=1) or higher Eulerian derivatives of velocity (l>1)
      go to (10,20,30) l
   10 do 11 i=1,n
        u(i) = zwr(i)*ps(i)-zwi(i)*pn(i)
   11   v(i) = zwr(i)*pn(i)+zwi(i)*ps(i)
      return
   20 do 21 i=1,n
        ut(i) = zwr(i)*ps(i)-zwi(i)*pn(i)
   21   vt(i) = zwr(i)*pn(i)+zwi(i)*ps(i)
      return
   30 do 31 i=1,n
        ut2(i) = zwr(i)*ps(i)-zwi(i)*pn(i)
   31   vt2(i) = zwr(i)*pn(i)+zwi(i)*ps(i)
      end
c==============================================================================
      subroutine times(l,er)
c------------------------------------------------------------------------------
c  l = order of Lagrangian time derivative of x, y and f to be calculated
c               or all, up to third order, if  l = 0  or  l < -3
c               or, up to  |l|-th order,  if  -4 < l < 0
c  er = precision of calculation
c------------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (NN=4100, N7=NN+7, py2=.62831853071795864769d1)
      common /profil/ x(-6:N7),y(-6:N7),f(-6:N7),
     *                gty,h,t,ym,enk,enp,ent,en0,dt,dtl,
     *                n,ndtn,uv,lbd,mbd
      common /auxlry/ ft(-6:N7),sf(-6:N7),sft(-6:N7),sft2(-6:N7)
      common /scalng/ wl,rx,tx,fx,xr,xt,xf
      common /drvtvs/   u(-6:N7),  v(-6:N7), fd(-6:N7),
     *                 ud(-6:N7), vd(-6:N7),fd2(-6:N7),
     *                ud2(-6:N7),vd2(-6:N7),fd3(-6:N7),
     *                udb(NN,5), vdb(NN,5), fdb(NN,5),
     *                ubd(NN,5), vbd(NN,5), fbd(NN,5),
     *      ut(-6:N7), vt(-6:N7),ut2(-6:N7),vt2(-6:N7),
     *      ux(-6:N7), vx(-6:N7),utx(-6:N7),vtx(-6:N7),
     *     uxx(-6:N7),vxx(-6:N7)
      common /pressr/ p0(-6:N7),pd1(-6:N7),pd2(-6:N7)
      common /potntl/ p(-6:N7),ps(-6:N7),pn(-6:N7),
     *                pno(-6:N7,3),pnl(-6:N7,3),k(3),lbp
      common /matrix/ a(NN,NN),b(NN,NN),
     *                zwa(-6:N7),zwr(-6:N7),zwi(-6:N7)
      common /contrl/ e,erp,lsm,lcs,istop,lsym
      common /rufnes/ err(-6:N7),erm
      common /print/  pt,pts,tl,mpt,lpt,ig
      common /propty/ slp(-6:N7),curv(-6:N7),angvel(-6:N7),
     *                strn(-6:N7),rstr(-6:N7),astr(-6:N7),
     *                Prn(-6:N7),strn0(-6:N7),erl(-6:N7)
      common /utilty/ s(-6:N7),d(-6:N7)
      common /lbldrv/ x1(-6:N7),y1(-6:N7),r1(-6:N7),
     *                x2(-6:N7),y2(-6:N7),
     *                u1(-6:N7),v1(-6:N7),
     *                ut1(-6:N7),vt1(-6:N7),
     *                ux1(-6:N7),vx1(-6:N7)
      logical uflag
      common /under/ uflag
      dimension rzs(-6:N7),sv(-6:N7),rpi(-6:N7)
      save
c------------------------------------------------------------------------------
c  branch according to order of derivative of x, y and f  to be calculated
c  all derivatives if  l = 0,  or up to |l|-th if  l < 0, to maximum of three
      if (abs(l).eq.999) return
      if (l.gt.0.and.l.le.3) go to (1,2,3) l
c==============================================================================
c set potential function for, and calculate, velocity and first derivative of f
    1 l1 = 1
      do 11 i=1,n
   11   p(i) = f(i)
c  consistent with time-stepping accuracy to use precision:  er^(35/24)
      tp = er**1.458333d0
      call derive(1,tp)
c  derivative of surface profile
      call ddi(x,n,py2,x1)
      call ddi(y,n,0.d0,y1)
c  calculation of the stream function sf
      do 50 i=1,n
   50   rpi(i) = u(i)*y1(i)-v(i)*x1(i)
      call integr(rpi,n,0.d0,sf,0.d0)
c  calculation of first Lagrangian time derivative of f
      do 21 i=1,n
        sv(i) = u(i)*u(i)+v(i)*v(i)
   21   fd(i) = .5*sv(i)-(p0(i)+gty*y(i)-uv*sf(i))
      if (abs(l).eq.1) return
c==============================================================================
    2 if (l1.ne.1) stop '<last call was not for first derivative>'
      l1 = 2
c  set potential function for, and calculate, acceleration and
c  second derivative of f,  using   p  =  ft  =  - [ p0 + gty*y + q.q/2 ]
      do 12 i=1,n
        p(i) = fd(i)-sv(i)-uv*y(i)*u(i)
   12   ft(i) = p(i)
c  consistent with time-stepping accuracy to use precision:  2*er^(7/6)
      tp = 2.d0*er**1.166667d0
      call derive(2,tp)
c  point label derivatives of velocity components at surface
      call ddi(u,n,0.d0,u1)
      call ddi(v,n,0.d0,v1)
c  calculation of the first time partial derivative sft of the stream function
      do 51 i=1,n
   51   rpi(i) = ut(i)*y1(i)-vt(i)*x1(i)
      call integr(rpi,n,0.d0,sft,0.d0)
      uv2 = uv*uv
      do 32 i=1,n
        rzs(i) = 1./(x1(i)*x1(i)+y1(i)*y1(i))
c  derivatives of velocity components w.r.t. x, using formulae derived
c  from incompressible (ux + vy = 0) and irrotational (uy - vx = 0) conditions
c  and the formulae:  u1 = ux*x1 + uy*y1   and   v1 = vx*x1 + vy*y1.
        ux(i) = (u1(i)*x1(i)-v1(i)*y1(i))*rzs(i)
        vx(i) = (u1(i)*y1(i)+v1(i)*x1(i))*rzs(i)
c  conversion to Lagrangian time derivative of velocity using the formula:
c   qd = qt + (q.grad)q   (q = vector velocity; d => Lagr. t=> Eul. derivative)
        ud(i) = ut(i)+u(i)*ux(i)+v(i)*vx(i)+uv*y(i)*ux(i)
        vd(i) = vt(i)+u(i)*vx(i)-v(i)*ux(i)+uv*y(i)*vx(i)
c  calculation of second Lagrangian time derivative of f
        sv(i) =  u(i)*ud(i)+v(i)*vd(i)
   32   fd2(i) = sv(i)-gty*v(i)+uv*sft(i)-y(i)*v(i)*uv2
      if (abs(l).eq.2) return
c==============================================================================
    3 if (l1.ne.2) stop '<last call was not for second derivative>'
c  set potential function for, and calculate, second derivative of velocity and
c  third derivative of f, using   p  =  ftt  =  - [ pd1 + gty*v + q.(qt + qd) ]
      do 13 i=1,n
   13   p(i) = fd2(i)-(u(i)+uv*y(i))*(ut(i)+ud(i))
     *          -v(i)*(vt(i)+vd(i))-sv(i)-uv*u(i)*v(i)
c  consistent with time-stepping accuracy to use precision:  6*er^(7/8)
      tp = 6.d0*er**.875d0
      call derive(3,tp)
c  point label derivatives of Eulerian time derivative components of velocity
c  for calculation (as above) of derivatives of these w.r.t. x
      call ddi(ut,n,0.d0,ut1)
      call ddi(vt,n,0.d0,vt1)
c  point label derivatives of derivatives of velocity components w.r.t. x for
c  calculation of second derivatives of velocity components w.r.t.  x
      call ddi(ux,n,0.d0,ux1)
      call ddi(vx,n,0.d0,vx1)
c calculation of second partial time derivative sft2 of the stream function
      do 52 i=1,n
   52   rpi(i) = ut2(i)*y1(i)-vt2(i)*x1(i)
      call integr(rpi,n,0.d0,sft2,0.d0)
      do 33 i=1,n
        utx(i) = (ut1(i)*x1(i)-vt1(i)*y1(i))*rzs(i)
        vtx(i) = (ut1(i)*y1(i)+vt1(i)*x1(i))*rzs(i)
        uxx(i) = (ux1(i)*x1(i)-vx1(i)*y1(i))*rzs(i)
        vxx(i) = (ux1(i)*y1(i)+vx1(i)*x1(i))*rzs(i)
c conversion to second Lagrangian time derivative of velocity using the formula
c   qdd = qtt + 2(q.grad)qt + (qt.grad)q + (q.grad)(q.grad)q
        o2yy = uv2*y(i)*y(i)
        uuvv = u(i)*u(i)-v(i)*v(i)
      ud2(i)=ut2(i)+2*(u(i)*utx(i)+v(i)*vtx(i))+ud(i)*ux(i)+vd(i)*vx(i)
     *  +uxx(i)*uuvv+2*u(i)*v(i)*vxx(i)+
     *  uv*ux(i)*v(i)+
     *  2*uv*y(i)*(utx(i)+u(i)*uxx(i)+v(i)*vxx(i))+
     *  o2yy*uxx(i)
      vd2(i) = vt2(i)+2*(u(i)*vtx(i)-v(i)*utx(i))+ud(i)*vx(i)-
     *         vd(i)*ux(i)+vxx(i)*uuvv-2*u(i)*v(i)*uxx(i)
     *         +uv*v(i)*vx(i)
     *         +2.*uv*y(i)*(vtx(i)+u(i)*vxx(i)-v(i)*uxx(i))
     *         +o2yy*vxx(i)
c  calculation of third Lagrangian time derivative of f
   33 fd3(i) = ud(i)*ud(i)+vd(i)*vd(i)+u(i)*ud2(i)+v(i)*vd2(i)-gty*vd(i)
     *         +uv*(ut(i)*v(i)-u(i)*vt(i)+sft2(i))
     *         -(uv2)*(y(i)*vt(i)+v(i)*v(i)+y(i)*vd(i))
c==============================================================================
c  Fourth derivatives would be obtained using:
c
c   p = pd2 + gty*vd + q.[3qtt + 3(q.grad)qt + 2(qt.grad)q + (q.grad)(q.grad)q]
c           + qd.qd  + qt.qt
c
c qddd = qttt + 3(q.grad)qtt + 3(qt.grad)qt + (qtt.grad)q
c             + (qt.grad)(q.grad)q + 2(q.grad)(qt.grad)q + 3(q.grad)(q.grad)qt
c             + (q.grad)(q.grad)(q.grad)q
c
c------------------------------------------------------------------------------
      end
c==============================================================================
      subroutine tangnt
c------------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (NN=4100, N7=NN+7, rpy=3.1830988618379067153d-1)
      common /profil/ x(-6:N7),y(-6:N7),f(-6:N7),
     *                gty,h,t,ym,enk,enp,ent,en0,dt,dtl,
     *                n,ndtn,uv,lbd,mbd
      common /auxlry/ ft(-6:N7),sf(-6:N7),sft(-6:N7),sft2(-6:N7)
      common /scalng/ wl,rx,tx,fx,xr,xt,xf
      common /drvtvs/   u(-6:N7),  v(-6:N7), fd(-6:N7),
     *                 ud(-6:N7), vd(-6:N7),fd2(-6:N7),
     *                ud2(-6:N7),vd2(-6:N7),fd3(-6:N7),
     *                udb(NN,5), vdb(NN,5), fdb(NN,5),
     *                ubd(NN,5), vbd(NN,5), fbd(NN,5),
     *      ut(-6:N7), vt(-6:N7),ut2(-6:N7),vt2(-6:N7),
     *      ux(-6:N7), vx(-6:N7),utx(-6:N7),vtx(-6:N7),
     *     uxx(-6:N7),vxx(-6:N7)
      common /pressr/ p0(-6:N7),pd1(-6:N7),pd2(-6:N7)
      common /potntl/ p(-6:N7),ps(-6:N7),pn(-6:N7),
     *                pno(-6:N7,3),pnl(-6:N7,3),k(3),lbp
      common /matrix/ a(NN,NN),b(NN,NN),
     *                zwa(-6:N7),zwr(-6:N7),zwi(-6:N7)
      common /contrl/ e,erp,lsm,lcs,istop,lsym
      common /rufnes/ err(-6:N7),erm
      common /print/  pt,pts,tl,mpt,lpt,ig
      common /propty/ slp(-6:N7),curv(-6:N7),angvel(-6:N7),
     *                strn(-6:N7),rstr(-6:N7),astr(-6:N7),
     *                Prn(-6:N7),strn0(-6:N7),erl(-6:N7)
      common /utilty/ s(-6:N7),d(-6:N7)
      common /lbldrv/ x1(-6:N7),y1(-6:N7),r1(-6:N7),
     *                x2(-6:N7),y2(-6:N7),
     *                u1(-6:N7),v1(-6:N7),
     *                ut1(-6:N7),vt1(-6:N7),
     *                ux1(-6:N7),vx1(-6:N7)
      logical uflag
      common /under/ uflag
      dimension c1(-6:N7),c2(-6:N7)
      save
c------------------------------------------------------------------------------
      ns = n
c  (case of symmetry)
      if (lsym.ne.0) then
        ns = n/2+1
        n2 = n+2
      end if
c   [p  represents the potential function of the calculation]
c  first and second derivatives w.r.t. integer surface parameter (point-label)
      call ddi(p,n,0.d0,ps)
      call d2di(p,n,0.d0,s)
c  tangential contribution to Cauchy integral for normal gradient (pn):
c (the values at the singularities include minus the second tangential
c   derivative of p at those points)
      do 53 i=1,ns
        ts = -s(i)
        do 52 j=1,ns
   52     ts = ts+b(i,j)*ps(j)
c  sum (integral) is divided by py
   53   s(i) = rpy*ts
      return
c==============================================================================
      entry normal(m,er)
c------------------------------------------------------------------------------
      m = 1
c  existence of an estimate for normal component of p (stored as pn) ?
      if (lbp.ge.0) then
      else
c  otherwise set pn equal to contribution, s, from tangential gradient
      do 1 i=1,ns
    1   pn(i) = s(i)
      end if
c  ensure at least  3  successive good iterations or one very good one
      itg = 3
c  matrix iteration steps for pn  (temporarily, c1 and c2),  with accelerated
c  convergence (Shank's) formula conditionally applied every two iterations
    2 do 11 l=1,24
      cq1 = 0.d0
      do 4 i=1,ns
        tn = a(i,1)*pn(1)
        do 3 j=2,ns
    3     tn = tn+a(i,j)*pn(j)
        c1(i) = rpy*tn+s(i)
        tn = zwa(i)*abs(c1(i)-pn(i))
    4   if (tn.gt.cq1) cq1 = tn
      if (cq1.lt.er) then
        itg = itg-1
      else
        itg = 3
      end if
      if (itg.eq.0.or.100.d0*cq1.lt.er)  then
        if (lsym.eq.1) then
          do 20 i=2,ns-1
   20       c1(n2-i) = c1(i)
        else if (lsym.eq.-1) then
          do 21 i=2,ns-1
   21       c1(n2-i) = -c1(i)
        end if
        do 5 i=1,n
    5   pn(i) = c1(i)
          return
      end if
      m = m+1
      cq2 = 0.d0
      do 7 i=1,ns
        tn = a(i,1)*c1(1)
        do 6 j=2,ns
    6     tn = tn+a(i,j)*c1(j)
        c2(i) = rpy*tn+s(i)
        tn = zwa(i)*abs(c2(i)-c1(i))
    7   if (tn.gt.cq2) cq2 = tn
      if (cq2.lt.er) then
        itg = itg-1
      else
        itg = 3
      end if
      if (itg.eq.0.or.100.d0*cq2.lt.er)  then
        if (lsym.eq.1) then
          do 22 i=2,ns-1
   22       c2(n2-i) = c2(i)
        else if (lsym.eq.-1) then
          do 23 i=2,ns-1
   23       c2(n2-i) = -c2(i)
        end if
        do 8 i=1,n
    8   pn(i) = c2(i)
          return
      end if
      do 9 i=1,ns
    9   call Shank(9,c2(i),c1(i),pn(i))
   11 m = m+1
      m = -nint(cq2/er+.5)
      end
c==============================================================================
      subroutine miscel
c------------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (NN=4100, N7=NN+7, MF=1031,
     *  py=.31415926535897932384d1, py2=.62831853071795864769d1)
      common /profil/ x(-6:N7),y(-6:N7),f(-6:N7),
     *                gty,h,t,ym,enk,enp,ent,en0,dt,dtl,
     *                n,ndtn,uv,lbd,mbd
      common /auxlry/ ft(-6:N7),sf(-6:N7),sft(-6:N7),sft2(-6:N7)
      common /scalng/ wl,rx,tx,fx,xr,xt,xf
      common /drvtvs/   u(-6:N7),  v(-6:N7), fd(-6:N7),
     *                 ud(-6:N7), vd(-6:N7),fd2(-6:N7),
     *                ud2(-6:N7),vd2(-6:N7),fd3(-6:N7),
     *                udb(NN,5), vdb(NN,5), fdb(NN,5),
     *                ubd(NN,5), vbd(NN,5), fbd(NN,5),
     *      ut(-6:N7), vt(-6:N7),ut2(-6:N7),vt2(-6:N7),
     *      ux(-6:N7), vx(-6:N7),utx(-6:N7),vtx(-6:N7),
     *     uxx(-6:N7),vxx(-6:N7)
      common /pressr/ p0(-6:N7),pd1(-6:N7),pd2(-6:N7)
      common /potntl/ p(-6:N7),ps(-6:N7),pn(-6:N7),
     *                pno(-6:N7,3),pnl(-6:N7,3),k(3),lbp
      common /matrix/ a(NN,NN),b(NN,NN),
     *                zwa(-6:N7),zwr(-6:N7),zwi(-6:N7)
      common /contrl/ e,erp,lsm,lcs,istop,lsym
      common /rufnes/ err(-6:N7),erm
      common /print/  pt,pts,tl,mpt,lpt,ig
      common /propty/ slp(-6:N7),curv(-6:N7),angvel(-6:N7),
     *                strn(-6:N7),rstr(-6:N7),astr(-6:N7),
     *                Prn(-6:N7),strn0(-6:N7),erl(-6:N7)
      common /utilty/ s(-6:N7),d(-6:N7)
      common /lbldrv/ x1(-6:N7),y1(-6:N7),r1(-6:N7),
     *                x2(-6:N7),y2(-6:N7),
     *                u1(-6:N7),v1(-6:N7),
     *                ut1(-6:N7),vt1(-6:N7),
     *                ux1(-6:N7),vx1(-6:N7)
      logical uflag
      common /under/ uflag
      dimension xdv(-6:N7),ydv(-6:N7),fdv(-6:N7),
     *          smd(-6:N7),spd(-6:N7)
      dimension yfx(-6:MF),ffx(-6:MF),ufx(-6:MF),vfx(-6:MF),dta(5)
      real o1,o2,o3,xre(NN),yre(NN),fre(NN),pr1(NN),pr2(NN),pr3(NN)
      integer ii
      CHARACTER(LEN=20) :: fname
      save
c==============================================================================
      entry prsure
      do 3 i=1,n
      p0(i) = 0.
      pd1(i) = 0.
    3 pd2(i) = 0.
      return
c==============================================================================
      entry sizedt(nd,xdv,ydv,fdv,dm)
c  setting time-step (dm) so that the terms of order  nd  in the Taylor
c  series (5-point-"smoothed" and "spreadd") have a maximum magnitude of  e:
        call copy(xdv,n,smd)
        call smooth(5,smd,n,0.d0)
        call spreadd(smd,n,0.d0,spd)
        call mabsm(spd,n,0.d0,dxm)
        call copy(ydv,n,smd)
        call smooth(5,smd,n,0.d0)
        call spreadd(smd,n,0.d0,spd)
        call mabsm(spd,n,0.d0,dym)
        call copy(fdv,n,smd)
        call smooth(5,smd,n,0.d0)
        call spreadd(smd,n,0.d0,spd)
        call hilosm(spd,n,0.d0,dfm,dfl)
        dfm = .5*(dfm-dfl)
        dm = sqrt(dxm*dxm+dym*dym+dfm*dfm)
c  calculations for terms of order   1 =< nd =< 8:
        if (nd.eq.5) then
          dm = (60.d0*e/dm)**.2d0
        else if (nd.eq.4) then
          dm = (24.d0*e/dm)**.25d0
        else if (nd.eq.3) then
          dm = (6.d0*e/dm)**.333333d0
        else if (nd.eq.2) then
          dm = (2.d0*e/dm)**.5d0
        else if (nd.eq.6) then
          dm = (120.d0*e/dm)**.166667d0
        else if (nd.eq.7) then
          dm = (210.d0*e/dm)**.142857d0
        else if (nd.eq.8) then
          dm = (336.d0*e/dm)**.125d0
        else if (nd.eq.1) then
          dm = (e/dm)
        else
          dm = (e/dm)
        end if
      return
c==============================================================================
      entry nicedt(dtime)
c------------------------------------------------------------------------------
c  routine to adjust the size of the time-step to any prescribed form,
c  for instance to ensure that wave profiles at "specified" times are
c  obtained by time-stepping precisely to such times
c  or to give printout (lpt = 2) at max/min kinetic energies
c------------------------------------------------------------------------------
c  restricting the rate of time-step growth to quadruple in about 5 steps:
      if (dt.gt.0.d0.and.dtime.gt.1.32*dt) dtime = 1.32*dt
c  capturing time-step sizes to land precisely on an output time:
      if (t+1.1*dtime.ge.pt) then
        dtime = pt-t
        pt = pt+abs(pts)
        if (pt.ge.tl) pt = tl
        mpt = 1
c  homing in smoothly from as many as almost 4 steps distant:
      else if (t+1.9*dtime.ge.pt) then
        dtime = .5d0*(pt-t)
      else if (t+2.9*dtime.ge.pt) then
        dtime = .33333d0*(pt-t)
      else if (t+3.9*dtime.ge.pt) then
        dtime = .25d0*(pt-t)
      end if
      if (lpt.ne.2.or.en3.lt.0.d0) return
        tp = en2+enr-2.d0*en1
        if (abs(en2-enr-2*tp).gt.2.3*abs(tp)) return
c  estimate time of max/min kinetic energy:
          tp = .5d0*(en2-enr)/tp
          dtime = -dt+.5d0*tp*(dt+dtl+tp*(dt-dtl))
          en1 = -1.d0
          en2 = -1.d0
          en3 = -1.d0
          mpt = 1
      return
c==============================================================================
c  using backward-differencing on previous highest derivative values to
c  extend use of higher derivatives:
      entry bakxyf
c  first update "previous" times of derivative calculations:
c  note that the proper definition of  "dta(.)"  increases with  "lbd"
        dta(5) = dta(4)-dt
        dta(4) = dta(3)-dt
        dta(3) = dta(2)-dt
        dta(2) = dta(1)-dt
        dta(1) = -dt
c  and increment the level of backward differencing, if appropriate:
        if (lbd.lt.abs(mbd)) lbd = lbd+1
c  preset the routine "bakdif" with these times of calculation:
        call setbak(lbd,dta)
c  and calculate the appropriate higher derivatives:
        call bakdif(lbd,ud2,udb,ubd,n,NN)
        call bakdif(lbd,vd2,vdb,vbd,n,NN)
        call bakdif(lbd,fd3,fdb,fbd,n,NN)
      return
c==============================================================================
      entry nextpn(dtn)
c  linear extrapolation (if possible) to next estimate for  pn:
      if (lbp.eq.1) then
        rdt = dtn/dt
        do 16 j=1,3
          do 16 i=1,n
            tp = pno(i,j)+rdt*(pno(i,j)-pnl(i,j))
            pnl(i,j) = pno(i,j)
   16       pno(i,j) = tp
      else
        do 17 j = 1,3
          do 17 i=1,n
   17       pnl(i,j) = pno(i,j)
        lbp = 1
      end if
      return
c==============================================================================
c  Time-stepping:
      entry Taylor
      if (ndtn.ge.0) then
        ndtn = ndtn+1
      else
        ndtn = ndtn-1
      end if
      dt2 = dt/2.d0
      dt3 = dt/3.d0
      dt4 = dt/4.d0
      dt5 = dt/2.5d0
      dt6 = dt/2.d0
      dt7 = dt/1.75d0
      dt8 = dt/1.6d0
      do 20 i=1,n
c  stepping of profile by taylor series expansion
c  with up to 5 orders of backward differencing:
      x(i) = x(i)+dt*(u(i)+uv*y(i)+dt2*(ud(i)+uv*v(i)
     *         +dt3*(ud2(i)+uv*vd(i)+dt4*(ubd(i,1)+uv*vd2(i)
     *      +dt5*(ubd(i,2)+uv*vbd(i,1)+dt6*(ubd(i,3)+uv*vbd(i,2)
     *   +dt7*(ubd(i,4)+uv*vbd(i,3)+dt8*(ubd(i,5)+uv*vbd(i,4)))))))))
      y(i) = y(i)+dt*(v(i)+dt2*(vd(i)+dt3*(vd2(i)+dt4*(vbd(i,1)
     *   +dt5*(vbd(i,2)+dt6*(vbd(i,3)+dt7*(vbd(i,4)+dt8*vbd(i,5))))))))
   20 f(i) = f(i)+dt*(fd(i)+dt2*(fd2(i)+dt3*(fd3(i)+dt4*(fbd(i,1)
     *   +dt5*(fbd(i,2)+dt6*(fbd(i,3)+dt7*(fbd(i,4)+dt8*fbd(i,5))))))))
c  time corresponding to this profile
      t = t+dt
c - - - - - - - - - - - - - - - - - - - - - - - -
      entry setsym
c  ensure symmetry is maintained if appropriate
      if (lsym.gt.0) then
        n2 = n/2+1
        xc = .5d0*(x(1)+x(n2)-py)
        x(1) = xc
        x(n2) = xc+py
        n2 = n+2
        do 21 i=2,n/2
          j = n2-i
          tp = .5d0*(x(i)+py2-x(j))
          x(i) = xc+tp
          x(j) = xc+py2-tp
          y(i) = .5d0*(y(i)+y(j))
          y(j) = y(i)
          f(i) = .5d0*(f(i)+f(j))
   21     f(j) = f(i)
      end if
      return
c==============================================================================
      entry findym
        ym = -9.d9
      entry meany
        ymn = 0.d0
c  derivative of surface profile x-component
        call ddi(x,n,py2,x1)
        do 25 i=1,n
c  calculation of mean surface level
   25     ymn = ymn+y(i)*x1(i)
        ymn = ymn/py2
        if (lcs.eq.0.or.lcs.eq.2.or.ym.eq.-9.d9) then
          ym = ymn
        else
c  remove any shift in mean level
          ymn = ym-ymn
          do 26 i=1,n
   26       y(i) = y(i)+ymn
        end if
      return
c==============================================================================
      entry seteng
        en0 = -9.d9
        en2 = -1.d0
        en1 = -1.d0
        enr = -1.d0
      entry energy
        en3 = en2
        en2 = en1
        en1 = enr
c  calculation of kinetic and potential energies per period
        enk = 0.d0
        enp = 0.d0
        call ddi(x,n,py2,x1)
        call ddi(y,n,0.d0,y1)
        do 27 i=1,n
          enk = enk+f(i)*(v(i)*x1(i)-u(i)*y1(i))
     *          +uv*(uv*x1(i)*y(i)*y(i)/3-2*f(i)*y1(i))*y(i)
   27     enp = enp+x1(i)*(y(i)-ym)*(y(i)-ym)
      enk = .5d0*enk
      enp = gty*.5d0*enp
      ent = enk+enp
      enr = enk/ent
      if (en0.eq.-9.d9) then
        en0 = ent
      else if (lcs.gt.1) then
c  restore total energy (assuming proportionate drift in both
c  kinetic and potential energies):
        tp = sqrt(en0/ent)
        do 28 i=1,n
          f(i) = tp*f(i)
   28     y(i) = ym+tp*(y(i)-ym)
      end if
C     print*, enk  
       ii =ii +1
C       WRITE(fname,'(A,I4.4,A)') 'enkS25BW25_',ii,'.txt'
C       OPEN(11,file=fname,form='formatted')
C       WRITE(11,*) enk
C       CLOSE(UNIT=11)
C        WRITE(fname,'(A,I4.4,A)') 'enpS25BW25_',ii,'.txt'
C       OPEN(11,file=fname,form='formatted')
C       WRITE(11,*) enp
C       CLOSE(UNIT=11)      
       
C      print*, enp
C       open(unit=10,file='enk.txt', status = 'unknown')
C       write(10,*) (i*uu(i), i = 1,neq) 
C      WRITE(*,'(A,I4.4,A)') 'enp.txt'
C       OPEN(11,form='formatted')
C       WRITE(11,*) enp
C       close(unit =11)
      return
c==============================================================================
      entry smxyf
      if (lsm.gt.0) then
        call smooth(lsm,x,n,py2)
        call smooth(lsm,y,n,0.d0)
        call smooth(lsm,f,n,0.d0)
      else if (lsm.lt.0) then
        call smthwd(-lsm,x,n,py2,err)
        call smthwd(-lsm,y,n,0.d0,err)
        call smthwd(-lsm,f,n,0.d0,err)
      end if
      return
c==============================================================================
      entry prtxyf
      write(8,99)
      if (n.le.10)  then
       write(8,100)  'x : ', (x(i),i=1,n)
       write(8,100)  'y : ', (y(i),i=1,n)
       write(8,100)  'f : ', (f(i),i=1,n)
      else
       write(8,100)  'x : ',(x(i),i=1,10)
       write(8,101)  (x(i),i=11,n)
       write(8,100)  'y : ',(y(i),i=1,10)
       write(8,101)  (y(i),i=11,n)
       write(8,100)  'f : ',(f(i),i=1,10)
       write(8,101)  (f(i),i=11,n)
      end if
      return
c------------------------------------------------------------------------------
      entry prLagr
      write(8,99)
      if (n.le.10)  then
       write(8,100)  'u : ', (u(i),i=1,n)
       write(8,100)  'v : ', (v(i),i=1,n)
       write(8,100)  'f1: ', (fd(i),i=1,n)
       if (lpt.eq.-1) return
      write(8,*)
       write(8,100)  'u1: ', (ud(i),i=1,n)
       write(8,100)  'v1: ', (vd(i),i=1,n)
       write(8,100)  'f2: ', (fd2(i),i=1,n)
       if (lpt.eq.-2) return
      write(8,*)
       write(8,100)  'u2: ', (ud2(i),i=1,n)
       write(8,100)  'v2: ', (vd2(i),i=1,n)
       write(8,100)  'f3: ', (fd3(i),i=1,n)
      else
       write(8,100)  'u : ',(u(i),i=1,10)
       write(8,101)  (u(i),i=11,n)
       write(8,100)  'v : ',(v(i),i=1,10)
       write(8,101)  (v(i),i=11,n)
       write(8,100)  'f1: ',(fd(i),i=1,10)
       write(8,101)  (fd(i),i=11,n)
       if (lpt.eq.-1) return
      write(8,*)
       write(8,100)  'u1: ',(ud(i),i=1,10)
       write(8,101)  (ud(i),i=11,n)
       write(8,100)  'v1: ',(vd(i),i=1,10)
       write(8,101)  (vd(i),i=11,n)
       write(8,100)  'f2: ',(fd2(i),i=1,10)
       write(8,101)  (fd2(i),i=11,n)
       if (lpt.eq.-2) return
      write(8,*)
       write(8,100)  'u2: ',(ud2(i),i=1,10)
       write(8,101)  (ud2(i),i=11,n)
       write(8,100)  'v2: ',(vd2(i),i=1,10)
       write(8,101)  (vd2(i),i=11,n)
       write(8,100)  'f3: ',(fd3(i),i=1,10)
       write(8,101)  (fd3(i),i=11,n)
      end if
      return
c------------------------------------------------------------------------------
      entry prter
      write(8,99)
      if (n.le.10)  then
       write(8,102)  'er: ',(err(i),i=1,n)
      else
       write(8,102)  'er: ',(err(i),i=1,10)
       write(8,103)  (err(i),i=11,n)
      end if
      return
c------------------------------------------------------------------------------
      entry column
c  post-processing output of columns of information
      ns = n
      if (lsym.ne.0) ns = n/2+1
      kol = 0
      go to 30
      entry slope
      ns = n
      kol = 1
   30 deg = 180.d0/py
      slp(0) = deg*atan2(y1(1),x1(1))
      do 31 i=1,ns
        slp(i) = deg*atan2(y1(i),x1(i))
        tp = slp(i-1)-slp(i)
   31   if (abs(tp).gt.180.d0) slp(i) = slp(i)+sign(360.d0,tp)
      if (kol.eq.1) return
      call d2di(x,n,py2,x2)
      call d2di(y,n,0.d0,y2)
      do 32 i=1,ns
        r1(i) = x1(i)*x1(i)+y1(i)*y1(i)
        curv(i) = (x1(i)*y2(i)-y1(i)*x2(i))/r1(i)
        angvel(i) = (x1(i)*v1(i)-y1(i)*u1(i)-uv*y1(i)*y1(i))/r1(i)
        r1(i) = sqrt(r1(i))
        curv(i) = curv(i)/r1(i)
        Prn(i) = (x1(i)*(vd(i)+gty)-y1(i)*(ud(i)+uv*v(i)))/r1(i)
        strn(i) = r1(i)/strn0(i)
        rstr(i) = sqrt(ux(i)*ux(i)+vx(i)*vx(i))
        astr(i) = deg*atan2(-vx(i),ux(i))
  32    erl(i) = log10(xr*err(i)+1.02d-10)
      write(8,99)
      write(8,80)
      if (wl.gt.0) then
        do 33 i=1,ns
          tp = slp(i)*slp(i)/(deg*deg)+rx*rx*curv(i)*curv(i)
     *         +tx*tx*angvel(i)*angvel(i)
          if (kol*24.lt.n.and.tp.le.3*e) then
            kol = kol+1
          else
            if (kol.eq.1) write(8,81)
     *        i-1,xr*x(i-1),xr*y(i-1),xf*f(i-1),xt*u(i-1),xt*v(i-1),
     *        Prn(i-1),slp(i-1),rx*curv(i-1),tx*angvel(i-1),
     *        strn(i-1),tx*rstr(i-1),astr(i-1),xr*r1(i-1),
     *        erl(i-1),i-1
            if (kol.gt.1) write(8,*) '    * ',kol,' *'
            kol = 0
            write(8,81) i,xr*x(i),xr*y(i),xf*f(i),xt*u(i),xt*v(i),
     *                  Prn(i),slp(i),rx*curv(i),tx*angvel(i),
c     *                  strn(i),tx*rstr(i),astr(i),xr*r1(i),
     *                  ud(i),vd(i),xr*r1(i),         
     *                  erl(i),i
          end if
   33   continue
      else
        do 34 i=1,ns
          tp = slp(i)*slp(i)/(deg*deg)+curv(i)*curv(i)
     *         +angvel(i)*angvel(i)
          if (kol.lt.n/24.and.tp.le.3*e) then
            kol = kol+1
          else
            if (kol.eq.1) write(8,81)
     *        i-1,x(i-1),y(i-1),f(i-1),u(i-1),v(i-1),
     *        Prn(i-1),slp(i-1),curv(i-1),angvel(i-1),
     *        strn(i-1),rstr(i-1),astr(i-1),r1(i-1),erl(i-1),i-1
            if (kol.gt.1) write(8,*) '    * ',kol,' *'
            kol = 0
            write(8,81) i,x(i),y(i),f(i),u(i),v(i),
     *                   Prn(i),slp(i),curv(i),angvel(i),
     *                  strn(i),rstr(i),astr(i),r1(i),erl(i),i
          end if
   34   continue
      end if
      write(8,99)
      return
c------------------------------------------------------------------------------
      entry fixtyp(ak)
      if (abs(lpt).le.3) return
      if (lps.eq.0) lps = 1
      lfx = lpt
      if (lfx.gt.3) then
        nfx = 0
        cw = 1.d0+.5d0*ak*ak
        if (lfx.le.2000) then
          if (lfx.le.1000) then
            gv = 1.d0
          else
            lfx = lfx-1000
            if (gv.lt.0.d0) gv = cw
          end if
          write(8,90) abs(lps)*xt*pts,lfx,gv,ak
          gv = gv*.5d0*tx
        else
          nfx = -1
          if (gv.lt.0.d0) gv = 1.d0
          lfx = lfx-2000
          write(8,91) abs(lps)*xt*pts,lfx,gv
          gv = gv*.5d0*tx
        end if
        if (abs(ig).eq.1) write(10) lfx
      else if (lfx.ge.-1000) then
        nfx = 1
        lfx = -lfx
        write(8,92) abs(lps)*xt*pts,lfx
      else
        nfx = 2
        lfx = -1000-lfx
        write(8,93) abs(lps)*xt*pts,lfx
      end if
   90 format(' time-steps,',f13.7,': at ',i3,' values of x-',f7.5,
     *       't/2: amplitude, A : f = Re[A.exp(y+ix-ict)] : c = ',
     *       '1+ak.ak/2 : ak = ',f5.4)
   91 format(' time-steps,',f13.7,': at ',i3,' evenly-spaced values',
     *       ' of x-',f7.5,'t/2:   values of  y(x)')
   92 format(' time-steps,',f13.7,':    ',i3,' Fourier modes of  y')
   93 format(' time-steps,',f13.7,':    ',i3,' Fourier modes of',
     *       '  y  and  v')
      return
c - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      entry xfixed
c  post-processing examination of data using Fast Fourier Transform
c  on data interpolated to evenly-spaced points (in x):
      if (nfx.le.0) then
        npt = lfx
        x0 = gv*t
      else
        npt = 4
   40     npt = 2*npt
          if (npt.lt.lfx) go to 40
        if (lfx.gt.16.and.3*npt/4.ge.lfx) npt = 3*npt/4
        npt = 3*npt
        x0 = 0.d0
      end if
      if (npt.gt.MF-7) npt = MF-7
      call locate (x,n,x0,d,nz,py2)
      call shift (d,n,s,nz,py2)
      call shift (y,n,d,nz,0.d0)
      if (nfx.eq.0) call shift (u,n,smd,nz,0.d0)
      if (abs(nfx).ne.1) call shift (v,n,spd,nz,0.d0)
      do 41 i=1,5
        s(1-i) = s(n+1-i)-py2
        d(1-i) = d(n+1-i)
        smd(1-i) = smd(n+1-i)
        spd(1-i) = spd(n+1-i)
        s(n+i) = s(i)+py2
        d(n+i) = d(i)
        smd(n+i) = smd(i)
   41   spd(n+i) = spd(i)
      dx = py2/npt
      nx = 0
      do 46 i=1,npt
        xx = (i-1)*dx
        do 42 j=nx,n
          jnx = j
   42     if(j.gt.0.and.s(j+1).gt.s(j).and..5*(s(j+1)+s(j)).gt.xx)
     *      go to 43
   43   xj = 2.d0*(xx-s(jnx))/(s(jnx+1)-s(jnx-1))
        if (jnx.gt.nx) then
          nx = jnx
          call pol11(s,nx,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
          a0 = a0-xx
          ab = 2.d0*a2
          ac = 3.d0*a3
          ad = 4.d0*a4
          ae = 5.d0*a5
          af = 6.d0*a6
          ag = 7.d0*a7
          ah = 8.d0*a8
          ai = 9.d0*a9
          aj = 10.d0*a10
          call pol11(d,nx,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)
          if (nfx.eq.0)
     *      call pol11(smd,nx,c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)
          if (abs(nfx).ne.1)
     *      call pol11(spd,nx,d0,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)
        else
          a0 = a0-dx
        end if
          do 44 its=1,16
            tp = xj
            xj = tp-(a0+tp*(a1+tp*(a2+tp*(a3+tp*(a4+tp*(a5+tp*(a6
     *                 +tp*(a7+tp*(a8+tp*(a9+tp*a10))))))))))
     *            / (a1+tp*(ab+tp*(ac+tp*(ad+tp*(ae+tp*(af+tp*(ag
     *                 +tp*(ah+tp*(ai+tp*aj)))))))))
   44       if (abs(xj-tp).lt.1.d-12) go to 45
   45   yfx(i)  =  b0+xj*(b1+xj*(b2+xj*(b3+xj*(b4+xj*(b5+xj*(b6
     *               +xj*(b7+xj*(b8+xj*(b9+xj*b10)))))))))
        if (abs(nfx).ne.1)
     *    vfx(i) = d0+xj*(d1+xj*(d2+xj*(d3+xj*(d4+xj*(d5+xj*(d6
     *               +xj*(d7+xj*(d8+xj*(d9+xj*d10)))))))))
        if (nfx.eq.0) then
          ufx(i) = c0+xj*(c1+xj*(c2+xj*(c3+xj*(c4+xj*(c5+xj*(c6
     *               +xj*(c7+xj*(c8+xj*(c9+xj*c10)))))))))
          tp = xr*(xx+x0)-cw*xt*t
          csx = cos(tp)
          snx = sin(tp)
          tp = xt/exp(xr*yfx(i))
c  redefine  u  and  v  as real and imaginary parts of the complex amplitude
          tpu = tp*(ufx(i)*csx+vfx(i)*snx)
          vfx(i) = tp*(vfx(i)*csx-ufx(i)*snx)
          ufx(i) = tpu
        end if
   46   continue
      if (ig.eq.1) then
        write(10) real(xt*t/py2)
        write(10) (real(yfx(i)),i=1,npt)
      end if
      if (nfx.eq.0) then
        if (lfa.le.npt) then
c  removing higher than second harmonics
          call fFt(ufx,npt,0.d0,yfx)
          call fFt(vfx,npt,0.d0,ffx)
          do 47 i=lfa,npt
            yfx(i) = 0.d0
   47       ffx(i) = 0.d0
          call revfFt(yfx,npt,0.d0,ufx)
          call revfFt(ffx,npt,0.d0,vfx)
        end if
c  redefining  y  and  f  as amplitude and phase of complex amplitude
        do 48 i=1,npt
          yfx(i) = sqrt(ufx(i)*ufx(i)+vfx(i)*vfx(i))
   48     ffx(i) = 57.29577951d0*atan2(vfx(i),ufx(i))
        if (abs(ig).eq.1) then
          write(10) (real(yfx(i)),i=1,npt)
          write(10) (real(ffx(i)),i=1,npt)
        end if
        if (npt.le.10) then
          write(8,108) 'fA: ',(yfx(i),i=1,npt)
          write(8,106) 'fP: ',(ffx(i),i=1,npt)
        else
          write(8,108) 'fA: ',(yfx(i),i=1,10)
          write(8,109) (yfx(i),i=11,npt)
          write(8,106) 'fP: ',(ffx(i),i=1,10)
          write(8,107) (ffx(i),i=11,npt)
        end if
        return
      end if
      if (nfx.eq.-1) then
        if (npt.le.10) then
          write(8,108) 'Y : ',(xr*yfx(i),i=1,npt)
        else
          write(8,108) 'Y : ',(xr*yfx(i),i=1,10)
          write(8,109) (xr*yfx(i),i=11,npt)
        end if
        return
      end if
      call Fourir(0.d0,h,gty,yfx,npt,0.d0,ufx,ffx,.3d-9)
c      tp = 180.d0/py
c      tp = tp*tp
      if (lfx.le.10) then
        write(8,104) 'yA: ',(xr*ufx(i),i=1,lfx)
        write(8,106) 'yP:',(ffx(i),i=1,lfx)
cc  (for use in determining numerical dispersion relation:)
c        if(t.gt.0) then
c          if (h.le.0.) write(8,104) 'w/rk',
c     *     (1.-ffx(i)/(t*sqrt(tp*i)),i=1,lfx)
c          if (h.gt.0.) write(8,104) 'w/rk',
c     *     (1.-ffx(i)/(t*tanh(i*h)*sqrt(tp*i)),i=1,lfx)
c        end if
      else
        write(8,104) 'yA: ',(xr*ufx(i),i=1,10)
        write(8,105) (xr*ufx(i),i=11,lfx)
        write(8,106) 'yP:',(ffx(i),i=1,10)
        write(8,107) (ffx(i),i=11,lfx)
c        if(t.gt.0) then
c          if (h.le.0.) then
c            write(8,104) 'w/rk',(1.-ffx(i)/(t*sqrt(tp*i)),i=1,10)
c            write(8,105) (1.-ffx(i)/(t*sqrt(tp*i)),i=11,lfx)
c          else
c            write(8,104) 'w/rk',
c     *       (1.-ffx(i)/(t*tanh(i*h)*sqrt(tp*i)),i=1,10)
c            write(8,105)
c     *       (1.-ffx(i)/(t*tanh(i*h)*sqrt(tp*i)),i=11,lfx)
c          end if
c        end if
      end if
      if (nfx.lt.2) return
      call Fourir(0.d0,h,gty,vfx,npt,0.d0,ufx,ufx,.3d-9)
      if (lfx.le.10) then
        write(8,104) 'vA: ',(xt*ufx(i),i=1,lfx)
        write(8,106) 'vP:',(ufx(i),i=1,lfx)
      else
        write(8,104) 'vA: ',(xt*ufx(i),i=1,10)
        write(8,105) (xt*ufx(i),i=11,lfx)
        write(8,106) 'vP:',(ufx(i),i=1,10)
        write(8,107) (ufx(i),i=11,lfx)
      end if
      return
c------------------------------------------------------------------------------
      entry postps
   50 write(*,120)
  120 format(/,' Post-processing;  the following options are ',
     *  'available:',//,'ENTER:',/,'  1  : to print a column-wise ',
     *  'list of properties at surface points',/,'  2  : for a ',
     *  'Fourier analysis of the wave-surface',/,'  3  : for the ',
     *  'complex amplitude of a modulated wave-train',/,'  4  : for',
     *  ' y(x) at evenly spaced values of x-t/2',/,'  5  : to ',
     *  'ww.und file',/,' 0 : to ',
     *  'print profile and velocity potential only, at surface ',
     *  'points',/,' -1  : to print profile, potential and 1st ',
     *  'derivatives',/,' -2  : to print profile, potential, 1st and',
     *  ' 2nd derivatives',/,' -3  : to print profile, potential, ',
     *  '1st, 2nd and 3rd derivatives')
      read(*,*) lpt
      if (lpt.lt.-3.or.lpt.gt.5) then
        write(*,*) lpt,' out of range, please re-enter'
        go to 50
      end if
      if (lpt.le.0.and.lpt.ge.-3) then
        pt = lpt
      else if (lpt.eq.1) then
        pt = 3
      else if (lpt.eq.2) then
        write(*,*) ' How many Fourier modes ?'
        read(*,*) lpt
        lpt = abs(lpt)
        if (lpt.lt.8) lpt = 16
        if (lpt.gt.(MF-7)/2) lpt = (MF-7)/2
        pt = -lpt
   51   write(*,121)
  121   format('ENTER:',/,'  1  : for Fourier modes of y(x) only',/,
     *         '  2  : for Fourier modes of y(x) and v(x)')
        read(*,*) lpt
        if (lpt.lt.1.or.lpt.gt.2) then
          write(*,*) lpt,' out of range, please re-enter'
          go to 51
        end if
        if (lpt.eq.2) pt = pt-1000
       else if (lpt.eq.5) then
          uflag=.true.
          pt=3       
      else if (lpt.eq.3) then
        write(*,122)
  122   format('Complex amplitude of velocity potential at evenly-',
     *   'spaced points',/,'moving at the group velocity, 1/2:',/,
     *   ' How many points ?')
        read(*,*) lpt
        lpt = abs(lpt)
        if (lpt.lt.5) lpt = 12
        if (lpt.gt.MF-7) lpt = MF-7
   52   write(*,123)
  123   format('Carrier-wave phase-speed is estimated by  cw = ',
     *   '1+.5.ak^2',/,'with  |ak| < .5 ;  ENTER  ak',/,' [give ak<0',
     *   ' to define group velocity = cw/2, or a specified value]')
        read(*,*) pt
        if (abs(pt).ge..5d0) then
          write(*,*) pt,' out of range, please re-enter'
          go to 52
        end if
        if (pt.lt.0) then
          lpt = lpt+1000
          write(*,124) 1.d0+.5d0*pt*pt
  124     format('Group velocity assumed as',f9.6,'/2',/,'ENTER : ',
     *      ' <0 to accept this, or a preferred value  ( x 2 )')
          read(*,*) gv
        end if
        pt = lpt+abs(pt)
        write(*,125)
  125   format('Higher-order effects should be removed:',/,'ENTER ',
     *   'number of carrier-mode (or zero to retain all modes)')
        read(*,*) lfa
        if (lfa.lt.0) lfa = 0
        lfa = lfa*4-1
        write(*,126)
  126   format('ENTER: 1 for unformatted output to "ww.g" of  Y  ',
     *    'and Complex amplitude',/,'      -1 for Complex amplitude',
     *    ' only')
        read(*,*) ig
        if (abs(ig).ne.1) ig = 0
      else if (lpt.eq.4) then
        write(*,127)
  127   format('How many values of x-t/2 ?  ( <0 for ',
     *    'a different choice of group velocity)')
        read(*,*) lpt
        gv = -1
        if (lpt.lt.0) then
          write(*,124) 1.d0
          read(*,*) gv
        end if
        write(*,128)
  128   format('ENTER: 1 for unformatted output to "ww.g" of  Y')
        read(*,*) ig
        if (ig.ne.1) ig = 0
        lpt = abs(lpt)
        if (lpt.lt.5) lpt = 12
        if (lpt.gt.MF-7) lpt = MF-7
        pt = 2000+lpt
      end if
c - - - - - - - - - - - -
      n = abs(n)
      read(7,end=61,err=61) (xre(i),i=1,n)
      do 55 i=1,n
   55   strn0(i) = xre(i)
      read(7,end=61,err=61) o1
      en0 = o1
      xt = 1.d0
      if (wl.gt.0.d0) xt = sqrt(wl/py2)
      tp = xt/py2
      write(*,200) xt*t,xt*tl,xt*pts,tp*t,tp*tl,tp*pts
  200 format(' Output from time, ',f10.5,', to,',f13.5,' in steps of ',
     *  f13.5,/,' [ ( x 2.py ) :    ',f10.5,', to,',f13.5,' in steps ',
     *  'of ',f13.5,' ]',/,'ENTER:',/,'  N  : for printout at ',
     *  'every  N-th  output time',/,' -N  : to be queried for ',
     *  'printout at every  N-th  output time')
      read(*,*) lps
      if (lps.eq.0) lps = -1
      j = abs(lps)
      if (lps.lt.0) write(*,201)
  201   format('at times shown, ENTER: 0 to skip,',/,23x,'1 for ',
     *    'printout,',/,23x,'2 to also create restarting input ',
     *    'file,',/,23x,'3 to alter printout query',/,22x,
     *    '-1 to stop')
      go to 56
      entry nextps
      j = 0
      if (lres.eq.1) then
        open(9,file='ww.res',form='unformatted')
        go to 70
      end if
   56 j = j+1
      lres = 0
      if (j.lt.abs(lps)) then
        do 62 i=1,7
   62     read(7,end=61,err=61)
        go to 56
      end if
      read(7,end=61,err=61) t,n,ndtn,o1,o2,o3
      erm = o1
      enk = o2
      enp = o3
      ent = enp+enk
      read(7,end=61,err=61) (xre(i),i=1,n)
      read(7,end=61,err=61) (yre(i),i=1,n)
      read(7,end=61,err=61) (fre(i),i=1,n)
      read(7,end=61,err=61) (pr1(i),i=1,n)
      read(7,end=61,err=61) (pr2(i),i=1,n)
      read(7,end=61,err=61) (pr3(i),i=1,n)
      if (lps.le.0) then
        read(7,end=57,err=57) tp,i,j,o1,o2,o3
        backspace(7)
        go to 58
   57     write(*,*) 'Last entry in data:'
   58   write(*,202) xt*t,xt*t/py2,100*enk/ent
  202   format('time =',f10.5,' [ =',f9.5,' x 2.py ] :',f7.3,'% KE')
        read(*,*) j
        if (j.lt.0) stop '<post-processing stopped>'
        if (j.eq.0) go to 56
        if (j.eq.2) lres = 1
        if (j.eq.3) then
          write(*,203)
  203     format('ENTER  N  : to be queried for printout at every ',
     *      ' N-th  output time')
          read(*,*) lps
          lps = -abs(lps)
          if (lps.eq.0) lps = -1
          write(*,201)
          go to 58
        end if
      end if
      do 59 i=1,n
        x(i) = pr1(i)
        y(i) = pr2(i)
   59   f(i) = pr3(i)
      do 60 i=1,n
        pno(i,1) = x(i)
        pno(i,2) = y(i)
        pno(i,3) = f(i)
        x(i) = xre(i)
        y(i) = yre(i)
   60   f(i) = fre(i)
      return
   61 if(uflag) then
       close(11)
      else
       close(8)
      endif
      if (ig.eq.1) close(10)
      stop '<end of post-processing data;  output in "ww.out">'
c------------------------------------------------------------------------------
   80 format('pt :',5x,'x',9x,'y',9x,'f',5x,':',4x,'u',7x,'v',
     *    4x,':',3x,'Pn',3x,':   slope   curve   ang.vel:',
c     *    ' strain, rate at angle :  dr  "rough":  pt')
     *    ' Du/Dt   Dv/Dt   :  dr  "rough":  pt')    
   81 format(i3,':',3f10.6,' :',2f8.4,' :',f7.4,' :',
     *       f8.2,f9.4,f8.4,' :',2f7.4,f8.2,' :',f6.4,f6.2,' :',i4)
   99 format(' ',130('-'))
  100 format(a4,10f12.7)
  101 format(4x,10f12.7)
  102 format(a4,10(1x,f11.9))
  103 format(5x,f11.9,1x,f11.9,1x,f11.9,1x,f11.9,1x,f11.9,1x,
     * f11.9,1x,f11.9,1x,f11.9,1x,f11.9,1x,f11.9,1x,f11.9,1x,
     * f11.9,1x,f11.9,1x,f11.9,1x,f11.9,1x,f11.9)
  104 format(a4,10f12.9)
  105 format(4x,10f12.9)
  106 format(a3,f13.7,9f12.6)
  107 format(4x,10f12.6)
  108 format(a3,10f12.8)
  109 format(3x,10f12.8)
c------------------------------------------------------------------------------
      entry wrudat
      write(12) real(t)
      print*,'t=',t
      write(12) (real(x(i)),i=1,n)
      write(12) (real(y(i)),i=1,n)
      write(12) (real(x1(i)),i=1,n)
      write(12) (real(y1(i)),i=1,n)
      write(12) (real(u(i)),i=1,n)
      write(12) (real(v(i)),i=1,n)
      write(12) (real(ux(i)),i=1,n)
      write(12) (real(vx(i)),i=1,n)
      write(12) (real(ut(i)),i=1,n)
      write(12) (real(vt(i)),i=1,n)
      write(12) (real(uxx(i)),i=1,n)
      write(12) (real(vxx(i)),i=1,n)
      write(12) (real(utx(i)),i=1,n)
      write(12) (real(vtx(i)),i=1,n)
      write(12) (real(f(i)),i=1,n)
      write(12) (real(ft(i)),i=1,n)
      write(12) (real(sf(i)),i=1,n)
      write(12) (real(sft(i)),i=1,n)
      return
c---------------------------------------------------------------------------
      entry resump
   70 rewind(9)
      write(9) xt*t,n,xr*h,gty,wl,uv
      write(9) (xr*x(i),i=1,n)
      write(9) (xr*y(i),i=1,n)
      write(9) (xf*f(i),i=1,n)
      write(9) erp,dble(lsm),dble(lcs),dble(mbd),xt*pts,xt*t
      write(9) (xr*strn0(i),i=1,n)
      write(9) (pno(i,1),i=1,n)
      write(9) (pno(i,2),i=1,n)
      write(9) (pno(i,3),i=1,n)
      if (lres.eq.1) then
        close(9)
        go to 56
      end if
c - - - - - - - - - - - - -
      end
c==============================================================================
      subroutine soothe(f,n,pf)
      implicit double precision (a-h,o-z)
      parameter (NN=4100, N7=NN+7)
      dimension x(-6:N7),y(-6:N7),f(-6:N7),r(-6:N7),
     *          s(-6:N7),g(-6:N7),h(-6:N7)
c==============================================================================
c  routines with various possible ideas related to "smoothing" and "roughness"
c------------------------------------------------------------------------------
      call interp(f,n,pf,h)
      call interp(h,n,pf,f)
      tp = f(n)
      do 1 i=n,2,-1
    1 f(i) = f(i-1)
      f(1) = tp-pf
      return
c------------------------------------------------------------------------------
      entry smooth(m,f,n,pf)
        smo = 9.d9
        do 6 j=1,m/4
          call smcalc(m,f,n,pf,s)
          do 5 i=1,n
    5       f(i) = f(i)-s(i)
          call maxabs(s,n,sm,im)
          if (sm.gt..8*smo) return
    6     smo = sm
      return
c------------------------------------------------------------------------------
      entry smthwa(m,f,n,pf)
        call smcalc(m,f,n,pf,s)
        call absval(s,n,g)
        go to 9
c - - - - - - - - - - - - - - - - -
      entry smthwd(m,f,n,pf,x)
        call absval(x,n,g)
    9   call maxabs(g,n,gm,im)
        if (gm.lt..1d-12) return
        gm = 3./gm
        do 10 i=1,n
          g(i) = g(i)*gm
   10     if (g(i).gt.1.)  g(i) = 1.
          call spreadd(g,n,0.d0,h)
        smo = 9.d9
        do 12 j=1,m/4
          call smcalc(m,f,n,pf,s)
          do 11 i=1,n
   11       f(i) = f(i)-s(i)*h(i)
          call maxabs(s,n,sm,im)
          if (sm.gt..8*smo) return
   12     smo = sm
      return
c------------------------------------------------------------------------------
      entry rough(m,x,y,f,n,px,r,er)
c  estimate "roughness", r, on vector <x,y,f>, and maximum "roughness", er,
c  using an m-point smoothing formula
c--------------------------------------------
      call smcalc(m,x,n,px,g)
      call smcalc(m,y,n,0.d0,h)
      call smcalc(m,f,n,0.d0,s)
      do 20 i=1,n
   20 r(i) = sqrt(g(i)*g(i)+h(i)*h(i)+s(i)*s(i))
      call maxabs(r,n,er,im)
      r(0) = r(n)
      r(n+1) = r(1)
      er1 = r(im+1)-r(im-1)
      er2 = r(im)+r(im)-r(im+1)-r(im-1)
      if (abs(er2).gt..00001*r(im)) then
        er = r(im)+.125*er1*er1/er2
      end if
      return
c------------------------------------------------------------------------------
      entry hilow(f,n,fh,fl)
      ih = 1
      il = 1
      do 51 i=2,n
        if (f(ih).lt.f(i)) ih = i
   51   if (f(il).gt.f(i)) il = i
      do 52 i=1,5
        f(1-i) = f(n+1-i)
   52   f(n+i) = f(i)
      call pol11(f,ih,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
      xx = 0.d0
      fi = a1
      fii = 2.d0*a2
      do 53 i=1,16
        if (abs(fii).lt.1.d-9) go to 54
        xx = xx-fi/fii
        if (abs(fi).lt.1.d-7) go to 54
        fi = a1+xx*(a2*2.d0+xx*(a3*3.d0+xx*(a4*4.d0+xx*(a5*5.d0
     *         +xx*(a6*6.d0+xx*(a7*7.d0+xx*(a8*8.d0+xx*(a9*9.d0
     *         +xx*a10*10.d0))))))))
   53   fii= 2.d0*a2+xx*(6.d0*a3+xx*(12.d0*a4+xx*(20.d0*a5+xx*(30.d0*a6
     *     +xx*(42.d0*a7+xx*(56.d0*a8+xx*(72.d0*a9+xx*90.d0*a10)))))))
   54 fh = a0+xx*(a1+xx*(a2+xx*(a3+xx*(a4+xx*(a5+xx*(a6+xx*(a7+xx*(a8
     *       +xx*(a9+xx*a10)))))))))
      call pol11(f,il,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
      xx = 0.d0
      fi = a1
      fii = 2.d0*a2
      do 55 i=1,16
        if (abs(fii).lt.1.d-9) go to 56
        xx = xx-fi/fii
        if (abs(fi).lt.1.d-7) go to 56
        fi = a1+xx*(a2*2.d0+xx*(a3*3.d0+xx*(a4*4.d0+xx*(a5*5.d0
     *         +xx*(a6*6.d0+xx*(a7*7.d0+xx*(a8*8.d0+xx*(a9*9.d0
     *         +xx*a10*10.d0))))))))
   55   fii= 2.d0*a2+xx*(6.d0*a3+xx*(12.d0*a4+xx*(20.d0*a5+xx*(30.d0*a6
     *      +xx*(42.d0*a7+xx*(56.d0*a8+xx*(72.d0*a9+xx*90.d0*a10)))))))
   56 fl = a0+xx*(a1+xx*(a2+xx*(a3+xx*(a4+xx*(a5+xx*(a6+xx*(a7+xx*(a8
     *       +xx*(a9+xx*a10)))))))))
      end
c==============================================================================
      subroutine operat
      implicit double precision (a-h,o-z)
      parameter (NN=4100, N7=NN+7, MF=1031,
     *  py=.31415926535897932384d1, py2=.62831853071795864769d1)
      dimension w(-6:N7),z(-6:N7),x(MF),wf(-6:MF),zf(-6:MF),ph(-6:MF)
      save
c==============================================================================
c  routines using  11-point  polynomial-fitting based formulae (unless stated
c  otherwise)
c------------------------------------------------------------------------------
      entry ddi(w,n,pw,z)
        do 1 i=1,5
          w(1-i) = w(n+1-i)-pw
    1     w(n+i) = w(i)+pw
        do 2 i=1,n
    2     z(i) = (((((w(i+5)-w(i-5))-12.5d0*(w(i+4)-w(i-4)))
     *          +75.d0*(w(i+3)-w(i-3)))-300.d0*(w(i+2)-w(i-2)))
     *          +1050.d0*(w(i+1)-w(i-1)))/1260.d0
      return
c------------------------------------------------------------------------------
      entry d2di(w,n,pw,z)
        do 11 i=1,5
          w(1-i) = w(n+1-i)-pw
   11     w(n+i) = w(i)+pw
        do 12 i=1,n
   12     z(i) = (((((w(i+5)+w(i-5))-15.625d0*(w(i+4)+w(i-4)))
     *          +125.d0*(w(i+3)+w(i-3)))-750.d0*(w(i+2)+w(i-2)))
     *          +5250.d0*(w(i+1)+w(i-1))-9220.75d0*w(i))/3150.d0
      return
c------------------------------------------------------------------------------
      entry smcalc(m,w,n,pw,z)
c  m-point smoothing
        do 21 i=1,5
          w(1-i) = w(n+1-i)-pw
   21     w(n+i) = w(i)+pw
       if (m.eq.5) then
        do 22 i=1,n
   22     z(i) = (6.d0*w(i)-(4.d0*(w(i+1)+w(i-1))-(w(i+2)+w(i-2))))
     *           /16.d0
       else if (m.eq.7) then
        do 23 i=1,n
   23     z(i) = (20.d0*w(i)-(15.d0*(w(i+1)+w(i-1))-(6.d0*(w(i+2)
     *           +w(i-2))-(w(i+3)+w(i-3)))))/64.d0
       else if (m.eq.9) then
         do 24 i=1,n
   24     z(i) = (70.d0*w(i)-(56.d0*(w(i+1)+w(i-1))
     *          -(28.d0*(w(i+2)+w(i-2))-(8.d0*(w(i+3)+w(i-3))
     *          -(w(i+4)+w(i-4))))))/256.d0
       else if (m.eq.8) then
        do 25 i=1,n
   25     z(i) =(910.d0*w(i)-(553.d0*(w(i+1)+w(i-1))
     *          -(14.d0*(w(i+2)+w(i-2))+(121.d0*(w(i+3)+w(i-3))
     *           -37.d0*(w(i+4)+w(i-4))))))/1728.d0
       else if (m.eq.11) then
        do 26 i=1,n
   26     z(i) = (252.d0*w(i)-(210.d0*(w(i+1)+w(i-1))
     *          -(120.d0*(w(i+2)+w(i-2))-(45.d0*(w(i+3)+w(i-3))
     *          -(10.d0*(w(i+4)+w(i-4))-(w(i+5)+w(i-5)))))))/1024.d0
       else if (m.eq.13) then
        w(-5) = w(n-5)-pw
        w(n+6) = w(6)+pw
        do 27 i=1,n
   27     z(i) =  ( 924.d0*w(i) -(792.d0*(w(i+1)+w(i-1))
     *  -(495.d0*(w(i+2)+w(i-2))-(220.d0*(w(i+3)+w(i-3))
     *   -(66.d0*(w(i+4)+w(i-4)) -(12.d0*(w(i+5)+w(i-5))
     *          -(w(i+6)+w(i-6)) ))))))/4096.d0
      else
        w(-5) = w(n-5)-pw
        w(-6) = w(n-6)-pw
        w(n+6) = w(6)+pw
        w(n+7) = w(7)+pw
        if (m.eq.14) then
          do 28 i=1,n
   28       z(i) =(1562484.d0*w(i)-(1132923.d0*(w(i+1)+w(i-1))
     *  -(286781.d0*(w(i+2)+w(i-2))+(247027.d0*(w(i+3)+w(i-3))
     *  -(288586.d0*(w(i+4)+w(i-4))-(136033.d0*(w(i+5)+w(i-5))
     *   -(32941.d0*(w(i+6)+w(i-6)) -  3367.d0*(w(i+7)+w(i-7))
     *     )))))))/2985984.d0
        else
          do 29 i=1,n
   29       z(i) = (3432.d0*w(i) -(3003.d0*(w(i+1)+w(i-1))
     *  -(2002.d0*(w(i+2)+w(i-2))-(1001.d0*(w(i+3)+w(i-3))
     *   -(364.d0*(w(i+4)+w(i-4)) - (91.d0*(w(i+5)+w(i-5))
     *    -(14.d0*(w(i+6)+w(i-6)) - (w(i+7)+w(i-7)) )))))))/16384.d0
        end if
      end if
      return
c------------------------------------------------------------------------------
      entry double (w,n,pw)
        j = 1
        go to 30
      entry interp(w,n,pw,z)
        j = 0
   30  do 31 i=1,5
        w(1-i) = w(n+1-i)-pw
   31   w(n+i) = w(i)+pw
       if (j.eq.0) then
        do 32 i=1,n
   32     z(i) = (39690.d0*(w(i)+w(i+1))-(8820.d0*(w(i-1)+w(i+2))
     *          -(2268.d0*(w(i-2)+w(i+3))-(405.d0*(w(i-3)+w(i+4))
     *          -(35.d0*(w(i-4)+w(i+5)))))))/65536.d0
       else
        do 33 i=1,n
   33     x(i) = (39690.d0*(w(i)+w(i+1))-(8820.d0*(w(i-1)+w(i+2))
     *          -(2268.d0*(w(i-2)+w(i+3))-(405.d0*(w(i-3)+w(i+4))
     *          -(35.d0*(w(i-4)+w(i+5)))))))/65536.d0
        do 34 i=n,1,-1
          j = 2*i
          w(j-1) = w(i)
   34     w(j) = x(i)
       end if
      return
c-------------------------------------------------------------------
      entry divide(w,n,m)
        n = n/m
        do 41 i=1,n
   41     w(i) = w(1+(i-1)*m)
      return
c------------------------------------------------------------------------------
      entry hilo(w,n,wh,wl)
      wh = w(1)
      wl = wh
      do 51 i=2,n
      if (wh.lt.w(i)) wh = w(i)
   51 if (wl.gt.w(i)) wl = w(i)
      return
c------------------------------------------------------------------------------
      entry hilosm(w,n,pw,wh,wl)
      w(n+1) = w(1)+pw
      wh = w(1)+w(1)+w(2)+w(n)-pw
      wl = wh
      do 61 i=2,n
      tp = w(i)+w(i)+w(i+1)+w(i-1)
      if (wh.lt.tp) wh = tp
   61 if (wl.gt.tp) wl = tp
      wh = .25*wh
      wl = .25*wl
      return
c------------------------------------------------------------------------------
      entry maxabs(w,n,wm,im)
      im = 1
      do 71 i=2,n
   71 if (abs(w(im)).lt.abs(w(i))) im = i
      wm = abs(w(im))
      return
c------------------------------------------------------------------------------
      entry mabsm(w,n,pw,wm)
      w(n+1) = w(1)+pw
      wm = abs(w(1)+w(1)+w(2)+w(n)-pw)
      do 81 i=2,n
      tp = abs(w(i)+w(i)+w(i+1)+w(i-1))
   81 if (wm.lt.tp) wm = tp
      wm = .25*wm
      return
c------------------------------------------------------------------------------
      entry spreadd(w,n,pw,z)
      do 91 i=1,5
      w(1-i) = w(n+1-i)-pw
   91 w(n+i) = w(i)+pw
      do 92 i=1,n
   92 z(i) = (252.d0*w(i)+(210.d0*(w(i+1)+w(i-1))
     *      +(120.d0*(w(i+2)+w(i-2))+(45.d0*(w(i+3)+w(i-3))
     *      +(10.d0*(w(i+4)+w(i-4))+(w(i+5)+w(i-5)))))))/1024.d0
      return
c------------------------------------------------------------------------------
      entry copy(w,n,z)
      do 101 i=1,n
  101 z(i) = w(i)
      return
c------------------------------------------------------------------------------
      entry absval(w,n,z)
      do 111 i=1,n
  111 z(i) = abs(w(i))
      return
c------------------------------------------------------------------------------
      entry locate(z,n,z0,w,nz,wvl)
      nz = (z(n)-z0+9999.d0*wvl)/wvl
      nz = nz-9999
      wln = nz*wvl+z0
      w(1) = z(1)-wln
      nz = 1
      do 121 i=2,n
      w(i) = z(i)-wln
  121 if(w(i).ge.w(i-1).and.w(i)+w(i-1).lt.0.d0) nz = i
      return
c------------------------------------------------------------------------------
      entry shift(z,n,w,nz,wvl)
      nz1 = nz-1
      do 131 i=nz,n
  131   w(i-nz1) = z(i)
      nnz = n-nz1
      do 132 i=1,nz1
  132   w(i+nnz) = z(i)+wvl
      return
c------------------------------------------------------------------------------
      entry integr(w,n,pw,z,z0)
c  integral from one point to the next using a  10-point  polynomial-formula
      do 141 i=1,5
      w(1-i) = w(n+1-i)-pw
  141 w(n+i) = w(i)+pw
      z(1) = z0
      do 142 i=2,n+1
  142   z(i) = z(i-1)
     *       +(4134338.d0*(w(i-1)+w(i)) -641776.d0*(w(i-2)+w(i+1))
     *         +162680.d0*(w(i-3)+w(i+2))-28939.d0*(w(i-4)+w(i+3))
     *           +2497.d0*(w(i-5)+w(i+4)))/7257600.d0
      return
c------------------------------------------------------------------------------
      entry pol11(z,n,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
c  fitting an 11-point polynomial centred on the n-th value of z
c  NOTE:  where necessary  f  should be extended outside  [1,N]  BEFORE
c         calling this routine.
      a0 = z(n)
        d1 = z(n+1)-z(n-1)
        d2 = z(n+2)-z(n-2)
        d3 = z(n+3)-z(n-3)
        d4 = z(n+4)-z(n-4)
        d5 = z(n+5)-z(n-5)
      a1 = (1050.d0*d1-300.d0*d2+75.d0*d3-12.5d0*d4+d5)/1260.d0
      a3 = (-70098.d0*d1+52428.d0*d2-14607.d0*d3+2522.d0*d4
     *      -205.d0*d5)/181440.d0
      a5 = (1938.d0*d1-1872.d0*d2+783.d0*d3-152.d0*d4+13.d0*d5)
     *     /34560.d0
      a7 = (-378.d0*d1+408.d0*d2-207.d0*d3+52.d0*d4-5.d0*d5)
     *     /120960.d0
      a9 = (42.d0*d1-48.d0*d2+27.d0*d3-8.d0*d4+d5)/725760.d0
        d1 = z(n+1)+z(n-1)
        d2 = z(n+2)+z(n-2)
        d3 = z(n+3)+z(n-3)
        d4 = z(n+4)+z(n-4)
        d5 = z(n+5)+z(n-5)
      a2 = (-73766.d0*a0+42000.d0*d1-6000.d0*d2+1000.d0*d3
     *      -125.d0*d4+8.d0*d5)/50400.d0
      a4 = (192654.d0*a0-140196.d0*d1+52428.d0*d2-9738.d0*d3
     *      +1261.d0*d4-82.d0*d5)/362880.d0
      a6 = (-12276.d0*a0+9690.d0*d1-4680.d0*d2+1305.d0*d3
     *      -190.d0*d4+13.d0*d5)/172800.d0
      a8 = (462.d0*a0-378.d0*d1+204.d0*d2-69.d0*d3+13.d0*d4-d5)
     *     /120960.d0
      a10= (-252.d0*a0+210.d0*d1-120.d0*d2+45.d0*d3-10.d0*d4+d5)
     *     /3628800.d0
        return
c-----------------------------------------------------------------------------
      entry Shank(lg,z2,z1,zo)
      tp = zo+z2-2.d0*z1
      if (lg*abs(tp).gt.abs(z2-z1)) then
c  accelerated convergence (for behaviour not too close to linear growth)
        zo = z1+(zo-z1)*(z2-z1)/tp
      else
c  otherwise a linear growth limited by  "lg"
        zo = z1+.5*lg*(zo-z2)
      end if
      return
c-----------------------------------------------------------------------------
      entry fFt(wf,n,pw,zf)
          md = 1
          go to 170
      entry Fourir(t,h,gty,wf,n,pw,zf,ph,er)
          md = 0
  170   x(1) = wf(1)
        tp = pw/n
        do 171 i=2,n
  171     x(i) = wf(i)-(i-1)*tp
        ifail = 0
C        call c06faf(x,n,zf,ifail)
        tp = sqrt(4.d0/n)
        do 172 i=0,n/2
  172     zf(2*i) = tp*x(i+1)
        do 173 i=1,(n+1)/2-1
  173     zf(2*i-1) = -tp*x(n+1-i)
        zf(2*((n+1)/2-1)+1) = 0.d0
        zf(0) = .5d0*zf(0)
        zf(n) = .5d0*zf(n)
        if (md.eq.0) then
          tp = 180.d0/py
          do 174 i=1,n/2
            sn = zf(2*i-1)
            cs = zf(2*i)
            zf(i) = sqrt(sn*sn+cs*cs)
            if (zf(i).ge.er) then
c  if  t  and  gty  are provided non-zero then a dispersion-relation correction
c  is made on the phase-measurement (useful for checking dispersion-relation):
              gkt = gty*i
              if (h.gt.er) gkt = gkt*tanh(i*h)
              ph(i) = tp*(mod(atan2(-sn,cs)+t*sqrt(gkt)+py,py2)-py)
            else
              ph(i) = 0.d0
            end if
  174       continue
        end if
      return
c ---------------------------------------------
      entry fFtint(zf,m,n,pw,wf)
        do 181 i=m+1,n
  181     zf(i) = 0.d0
c - - - - - - - - - - - - - - - - - - - -
      entry revfFt(zf,n,pw,wf)
        tp = sqrt(.25d0*n)
        do 182 i=0,n/2
  182     x(i+1) = tp*zf(2*i)
        do 183 i=1,(n+1)/2-1
  183     x(n+1-i) = tp*zf(2*i-1)
        x(1) = 2.d0*x(1)
        if (2*(n/2).eq.n) x(n/2+1) = 2.d0*x(n/2+1)
        ifail = 0
C        call c06fbf (x,n,wf,ifail)
        tp = pw/n
        do 184 i=1,n
  184     wf(i) = x(i)+(i-1)*tp
c------------------------------------------------------------------------------
      end
c==============================================================================
      subroutine bakdif(k,do,db,bd,n,nn)
      implicit double precision (a-h,o-z)
      dimension do(-6:nn),db(nn,5),bd(nn,5),dd(5)
      save
c==============================================================================
        go to (11,21,31,41,51,61) k+1
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      entry setbak(k,dd)
        go to (10,20,30,40,50,60) k+1
        write(*,*) '"bakdif" paraketer k = ',k,' out of range'
      return
c==============================================================================
c  fitting to 6 points at  t = 0, a, b, c, d, e:
   60   m = 6
        a = dd(1)
        b = dd(2)
        c = dd(3)
        d = dd(4)
        e = dd(5)
c  column multipliers:  reciprocal product of differences from  a,  b,  etc.
        ax = 1.d0/((a-b)*(a-c)*(a-d)*(a-e)*a)
        bx = 1.d0/((a-b)*(b-c)*(b-d)*(b-e)*b)
        cx = 1.d0/((a-c)*(b-c)*(c-d)*(c-e)*c)
        dx = 1.d0/((a-d)*(b-d)*(c-d)*(d-e)*d)
        ex = 1.d0/((a-e)*(b-e)*(c-e)*(d-e)*e)
        ae = a*b*c*d*e
        ea = a+b+c+d+e
c  quartic combinations not containing  a,  b,  etc.
        x11 = ae/a*ax
        x21 =-ae/b*bx
        x31 = ae/c*cx
        x41 =-ae/d*dx
        x51 = ae/e*ex
        x01 =-(x11+x21+x31+x41+x51)
c  sums of cubic combinations not containing  a,  b,  etc.
        ab = a*b
        de = d*e
        ed = d+e
        ba = a+b
        x12 =-(b*c*(ed)+(b+c)*de)*ax
        x22 = (a*c*(ed)+(a+c)*de)*bx
        x32 =-(ab* (ed)+(ba) *de)*cx
        x42 = (ab*(c+e)+(ba)*c*e)*dx
        x52 =-(ab*(c+d)+(ba)*c*d)*ex
        x02 =-(x12+x22+x32+x42+x52)
c  sums of quadratic combinations not containing  a,  b,  etc.
        ae = a*(ea-a)+b*(c+ed)+c*(ed)+de
        x13 = (ae-a*(ea-a))*ax
        x23 =-(ae-b*(ea-b))*bx
        x33 = (ae-c*(ea-c))*cx
        x43 =-(ae-d*(ea-d))*dx
        x53 = (ae-e*(ea-e))*ex
        x03 =-(x13+x23+x33+x43+x53)
c  sums of linear combinations not containing  a,  b,  etc.
        x14 =-(ea-a)*ax
        x24 = (ea-b)*bx
        x34 =-(ea-c)*cx
        x44 = (ea-d)*dx
        x54 =-(ea-e)*ex
        x04 =-(x14+x24+x34+x44+x54)
c  unity:
        x15 = ax
        x25 =-bx
        x35 = cx
        x45 =-dx
        x55 = ex
        x05 =-(x15+x25+x35+x45+x55)
      return
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   61   if (m.ne.6) write(*,*) 'not set to fit to 6 points',m
        do 62 i=1,n
          bd(i,1) = x01*do(i)+x11*db(i,1)+x21*db(i,2)+x31*db(i,3)
     *         +x41*db(i,4)+x51*db(i,5)
          bd(i,2) = x02*do(i)+x12*db(i,1)+x22*db(i,2)+x32*db(i,3)
     *         +x42*db(i,4)+x52*db(i,5)
          bd(i,3) = x03*do(i)+x13*db(i,1)+x23*db(i,2)+x33*db(i,3)
     *         +x43*db(i,4)+x53*db(i,5)
          bd(i,4) = x04*do(i)+x14*db(i,1)+x24*db(i,2)+x34*db(i,3)
     *         +x44*db(i,4)+x54*db(i,5)
          bd(i,5) = x05*do(i)+x15*db(i,1)+x25*db(i,2)+x35*db(i,3)
     *         +x45*db(i,4)+x55*db(i,5)
c  shifting, ready for next differencing:
          db(i,5) = db(i,4)
          db(i,4) = db(i,3)
          db(i,3) = db(i,2)
          db(i,2) = db(i,1)
   62     db(i,1) = do(i)
      return
c---------------------------------------------------------------
c  fitting to 5 points at  t = 0, a, b, c, d:
   50   m = 5
        a = dd(1)
        b = dd(2)
        c = dd(3)
        d = dd(4)
        ax = 1.d0/((a-b)*(a-c)*(a-d)*a)
        bx = 1.d0/((a-b)*(b-c)*(b-d)*b)
        cx = 1.d0/((a-c)*(b-c)*(c-d)*c)
        dx = 1.d0/((a-d)*(b-d)*(c-d)*d)
        x11 =-b*c*d*ax
        x21 = a*c*d*bx
        x31 =-a*b*d*cx
        x41 = a*b*c*dx
        x01 =-(x11+x21+x31+x41)
        x12 = (b*(c+d)+c*d)*ax
        x22 =-(a*(c+d)+c*d)*bx
        x32 = (a*(b+d)+b*d)*cx
        x42 =-(a*(b+c)+b*c)*dx
        x02 =-(x12+x22+x32+x42)
        x13 =-(b+c+d)*ax
        x23 = (a+c+d)*bx
        x33 =-(a+b+d)*cx
        x43 = (a+b+c)*dx
        x03 =-(x13+x23+x33+x43)
        x14 = ax
        x24 =-bx
        x34 = cx
        x44 =-dx
        x04 =-(x14+x24+x34+x44)
      return
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   51   if (m.ne.5) write(*,*) 'not set to fit to 5 points',m
        do 52 i=1,n
          bd(i,1) = x01*do(i)+x11*db(i,1)+x21*db(i,2)
     *         +x31*db(i,3)+x41*db(i,4)
          bd(i,2) = x02*do(i)+x12*db(i,1)+x22*db(i,2)
     *         +x32*db(i,3)+x42*db(i,4)
          bd(i,3) = x03*do(i)+x13*db(i,1)+x23*db(i,2)
     *         +x33*db(i,3)+x43*db(i,4)
          bd(i,4) = x04*do(i)+x14*db(i,1)+x24*db(i,2)
     *         +x34*db(i,3)+x44*db(i,4)
c  shifting:
          db(i,5) = db(i,4)
          db(i,4) = db(i,3)
          db(i,3) = db(i,2)
          db(i,2) = db(i,1)
   52     db(i,1) = do(i)
      return
c---------------------------------------------------------------
c  fitting to 4 points at  t = 0, a, b, c:
   40   m = 4
        a = dd(1)
        b = dd(2)
        c = dd(3)
        ax = 1.d0/((a-b)*(a-c)*a)
        bx = 1.d0/((a-b)*(b-c)*b)
        cx = 1.d0/((a-c)*(b-c)*c)
        x11 = b*c*ax
        x21 =-a*c*bx
        x31 = a*b*cx
        x01 =-(x11+x21+x31)
        x12 =-(b+c)*ax
        x22 = (a+c)*bx
        x32 =-(a+b)*cx
        x02 =-(x12+x22+x32)
        x13 = ax
        x23 =-bx
        x33 = cx
        x03 =-(x13+x23+x33)
      return
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   41   if (m.ne.4) write(*,*) 'not set to fit to 4 points',m
        do 42 i=1,n
          bd(i,1) = x01*do(i)+x11*db(i,1)+x21*db(i,2)+x31*db(i,3)
          bd(i,2) = x02*do(i)+x12*db(i,1)+x22*db(i,2)+x32*db(i,3)
          bd(i,3) = x03*do(i)+x13*db(i,1)+x23*db(i,2)+x33*db(i,3)
c  shifting:
          db(i,4) = db(i,3)
          db(i,3) = db(i,2)
          db(i,2) = db(i,1)
   42     db(i,1) = do(i)
      return
c---------------------------------------------------------------
c  fitting to 3 points at  t = 0, a, b:
   30   m = 3
        a = dd(1)
        b = dd(2)
        ax = 1.d0/((a-b)*a)
        bx = 1.d0/((a-b)*b)
        x11 =-b*ax
        x21 = a*bx
        x01 =-(x11+x21)
        x12 = ax
        x22 =-bx
        x02 =-(x12+x22)
      return
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   31   if (m.ne.3) write(*,*) 'not set to fit to 3 points',m
        do 32 i=1,n
          bd(i,1) = x01*do(i)+x11*db(i,1)+x21*db(i,2)
          bd(i,2) = x02*do(i)+x12*db(i,1)+x22*db(i,2)
c  shifting:
          db(i,3) = db(i,2)
          db(i,2) = db(i,1)
   32     db(i,1) = do(i)
      return
c---------------------------------------------------------------
c  fitting to 2 points at  t = 0, a:
   20   m = 2
        a = dd(1)
        ax = 1.d0/a
        x11 = ax
        x01 =-ax
      return
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   21   if (m.ne.2) write(*,*) 'not set to fit to 2 points',m
        do 22 i=1,n
          bd(i,1) = x01*do(i)+x11*db(i,1)
c  shifting:
          db(i,2) = db(i,1)
   22     db(i,1) = do(i)
      return
c---------------------------------------------------------------
c  fitting to 1 point at  t = 0 (does nothing):
   10   m = 1
      return
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   11  if (m.ne.1) write(*,*) 'not set to fit to 1 point',m
        do 12 i=1,n
c  shifting:
   12     db(i,1) = do(i)
      return
c==============================================================================
      end

