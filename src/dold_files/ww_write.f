      program wwrite
      
      implicit double precision (a-h,o-z)
      integer niv
      parameter(niv = 32)
      common /profil/ x(niv),y(niv),f(niv),po(niv),po1(niv),
     *					po2(niv), 
     *                 gty,h,t,ym,enk,enp,ent,en0,dt,dtl,En,
     *                 n,ndtn,uv,lbd,mbd,BW,eps,BWo, C,eps0,k
      CHARACTER(LEN=20) :: fname
	  t=0.d0
	  n=niv
	  h=0d0
C	  gty=1d0
      gty=9.81d0 
	  uv=0.d0 
	  erp=1e-6
	  sm=0.d0
	  cs=4.d0
	  bd=0.d0
      pts =0.1d0
      open(7,file='S.txt',action='read', status='old',iostat=ierrorx)
      read(7,*) eps
      close(7)
      open(7,file='tl.txt',action='read', status='old',iostat=ierrorx)
      read(7,*) tl
      close(7)
      open(12,file='wl.txt',action='read',status='old',iostat=ierrorx)
      read(12,*) wl
      close(12)	     
      open(8,file='xc.txt',action='read',status='old',iostat=ierrorx)
      read(8,*) x
	   close(8)
	   open(9,file='yc.txt',action='read',status='old',iostat=ierrory)
	   read(9,*) y
	   close(9)
	  open(10,file='fc.txt',action='read',status='old',iostat=ierrory)
	   read(10,*) f
	   close(10)	   	   	   
	   BW=1d0
C  open(20,file='k.txt',action='read',status='old',iostat=ierrory)
C   read(20,*) k
C   close(20)	   
      open(11,file='ww.in',form='unformatted')
       rewind(11)
       write(11) t,n,h,gty,wl,uv,eps,BW
       write(11) (x(i),i=1,n)
       write(11) (y(i),i=1,n)
       write(11) (f(i),i=1,n)
       write(11) erp,sm,cs,bd,pts,tl
       close(11)    
C      WRITE(fname,'(A,F0.3,F0.3,A)') 'IC',eps, BW,'.txt'
C      OPEN(81,file=fname,form='formatted')
C      write(81,*) t,n,h,gty,wl,eps,BW,erp,sm,cs,bd,pts,tl
C      close(81) 
C      print*, BW*1d0
       end
