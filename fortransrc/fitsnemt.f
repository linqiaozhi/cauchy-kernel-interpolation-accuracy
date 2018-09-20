      implicit real *8(a-h,o-z)
C$      include 'omp_lib.h'      
      real *8 w(1 000 000)

      call prini(6,13)
      call prin2('Enter n*',n,0)
      read *, n

      call testequi(w)
      
      stop
      end

c-----------------------------------------------------

      subroutine testequi(w)
      implicit real *8 (a-h,o-z)
      real *8 w(*)

      parameter (nmax = 1 000 000)

      real *8 xs(nmax),xsp(nmax),xsort(nmax)
      real *8 ys(nmax),ysp(nmax),ysort(nmax)
      real *8 sigma1(nmax),sigma2(nmax),sigma3(nmax),sigma4(nmax)
      real *8 sigma1sort(nmax),sigma2sort(nmax),sigma3sort(nmax),
     1          sigma4sort(nmax)
      real *8 pot1(nmax),pot2(nmax),pot3(nmax),pot4(nmax)
      real *8 pot1sort(nmax),pot2sort(nmax),pot3sort(nmax),
     1     pot4sort(nmax)
      real *8 pot1ex(nmax),pot2ex(nmax),pot3ex(nmax),pot4ex(nmax)

      integer, allocatable :: icnt(:),ioffst(:),iboxsort(:)
      real *8, allocatable :: boxl(:,:),boxr(:,:),tvals(:,:)
      integer nterms,iarr(nmax)
      integer, allocatable :: irearr(:)
      real *8 xpts(10000),prods(10000),tinfo(10)
      real *8, allocatable :: xptsall(:),yptsall(:)
      complex *16, allocatable :: zkvals(:,:)
      complex *16, allocatable :: zfftin(:,:),zfftout(:,:)
      complex *16, allocatable :: zksave(:)
      real *8, allocatable :: workspace(:)
      integer *8  plan

      integer omp_get_num_threads

C$OMP PARALLEL
C$      nthreads = omp_get_num_threads()
C$OMP END PARALLEL
C$      call prinf('nthreads=*',nthreads,1)      


      bmin = -10.0d0
      bmax = 10.0d0

      nlat = 50

      nboxes = nlat*nlat
      nterms = 3
      rlbox = (bmax - bmin)/(nlat+0.0d0)
      n = 1000000

      xmin = -bmax
      xmax = bmax

      ymin = -bmax
      ymax = bmax

      do i=1,n
         xs(i) = xmin + hkrand(0)*(xmax-xmin)
         ys(i) = ymin + hkrand(0)*(ymax-ymin)

         sigma1(i) = hkrand(0)-0.5d0
         sigma2(i) = hkrand(0)-0.5d0
         sigma3(i) = hkrand(0)-0.5d0
         sigma4(i) = hkrand(0)-0.5d0

      enddo

      t1tot = second()
C$    t1tot = omp_get_wtime()      

      t1 = second()
C$       t1 = omp_get_wtime()
      allocate(boxl(2,nboxes),boxr(2,nboxes))
      allocate(icnt(nboxes+10),ioffst(nboxes+10),iboxsort(n))


      do i=1,nboxes
         icnt(i) = 0
      enddo

      ndim = 4

c
cc      create sources and binsort

      do 1000 i=1,n

      ixs = (xs(i)-bmin)/rlbox
      iys = (ys(i)-bmin)/rlbox
      if(ixs.ge.nboxes) ixs = nboxes-1
      if(iys.ge.nboxes) iys = nboxes-1
      if(ixs.lt.0) ixs = 0
      if(iys.lt.0) iys = 0

      icol = 1 + ixs
      irow = 1 + iys

      iadr = (irow-1)*nlat + icol

      icnt(iadr) = icnt(iadr)+1
      iboxsort(i) = iadr
      
 1000 continue


      ioffst(1) = 0
      do j=2,nboxes+1

      ioffst(j) = ioffst(j-1)+icnt(j-1)

      enddo

      do i=1,nboxes

      icnt(i) = 1

      enddo

      do i=1,n
        
         iadr = iboxsort(i)
         indx = ioffst(iadr)+icnt(iadr)
         iarr(indx) = i
         icnt(iadr) = icnt(iadr)+1

      enddo

      
      ii = 1
      do 1150 i=1,nlat
      do 1100 j=1,nlat

      boxl(1,ii) = bmin + (j-1)*rlbox
      boxr(1,ii) = bmin + j*rlbox

      boxl(2,ii) = bmin + (i-1)*rlbox
      boxr(2,ii) = bmin + i*rlbox

      ii = ii + 1

 1100 continue  
 1150 continue  



c
cc     sort everything
c

      call rsortsrcall(xs,xsort,ys,ysort,sigma1,sigma1sort,
     1      sigma2,sigma2sort,sigma3,sigma3sort,sigma4,
     2      sigma4sort,iarr,n)
      t2 = second()
C$     t2 = omp_get_wtime()


      call prin2('bin sorting time=*',t2-t1,1)

      t1 = second()
C$      t1 = omp_get_wtime()      
c
cc     precompute all precomutable operators
c

      h = 2.0d0/(nterms+0.0d0)
      xpts(1) = -1.0d0 + h/2
      do i=2,nterms

      xpts(i) = xpts(1) + (i-1)*h

      enddo

      call getprods(nterms,xpts,prods)

      nn = nterms*nterms*nboxes

      allocate(xptsall(nn),yptsall(nn),irearr(nn))
c
cc      get xptsall
c

 
      ii = 1

      nfourh = nterms*nlat

      nfour = 2*nterms*nlat


      h = h*rlbox/2.0d0

      xstart = bmin + h/2
      ystart = bmin + h/2
      ii = 1
      do i=1,nfourh
      do j=1,nfourh

      xptsall(ii) = xstart + (i-1)*h
      yptsall(ii) = ystart + (j-1)*h
      ii = ii +1


      enddo
      enddo

      ii = 1
      do ilat=1,nlat
      do jlat=1,nlat


      do i=1,nterms
      do j=1,nterms

      iy = (ilat-1)*nterms + j
      ix = (jlat-1)*nterms + i

      ipt = (ix-1)*nlat*nterms + iy

      irearr(ii) = ipt
      ii = ii + 1

      enddo
      enddo

      enddo
      enddo

      allocate(zkvals(nfour,nfour),zksave(10*nfour))
      allocate(workspace(10*nfour))

      do i=1,nfour
      do j=1,nfour

      zkvals(i,j) = 0

      enddo
      enddo


      ii = 1
      do i=1,nfourh
      do j=1,nfourh

      call fker(tmp,xptsall(ii),yptsall(ii),xptsall(1),yptsall(1))
      zkvals(i+nfourh,j+nfourh) = tmp
      zkvals(nfourh-i+2,j+nfourh) = tmp
      zkvals(i+nfourh,nfourh-j+2) = tmp
      zkvals(nfourh-i+2,nfourh-j+2) = tmp
      ii = ii+1
      
      enddo
      enddo

      call dcffti(nfour,zksave)
      call fourt2d_f_pre(zkvals,nfour,workspace,zksave)

      do i=1,nfour
      do j=1,nfour

      zkvals(i,j) = zkvals(i,j)/(nfour*nfour+0.0d0)

      enddo
      enddo


      t2 = second()
C$      t2 = omp_get_wtime()      
      call prin2('pre comp time=*',t2-t1,1)


      t1 = second()
C$      t1 = omp_get_wtime()
      call fasttsne2d(n,ndim,xsort,ysort,sigma1sort,sigma2sort,
     1     sigma3sort,sigma4sort,nboxes,nterms,xpts,
     1     boxl,boxr,ioffst,prods,nn,nfour,irearr,zkvals,zksave,
     2     workspace,pot1sort,pot2sort,pot3sort,pot4sort,tinfo)


      t2 = second()
C$      t2 = omp_get_wtime()     

      t1 = second()
C$      t1 = omp_get_wtime()      

      call rsortrevtargall(pot1sort,pot1,pot2sort,pot2,pot3sort,pot3,
     1        pot4sort,pot4,iarr,n)

      t2 = second()
C$      t2 = omp_get_wtime()     

      call prin2('unsorting time=*',t2-t1,1)

      t2tot = second()
C$      t2tot = omp_get_wtime()

      call prin2('total time=*',t2tot-t1tot,1)



c
cc     Compute error
c
      ntest = min(50,n)
      erra = 0
      ra = 0
      do i=1,ntest

      pot1ex(i) = 0
      pot2ex(i) = 0
      pot3ex(i) = 0
      pot4ex(i) = 0

      do j=1,n

      call fker(ff,xs(i),ys(i),xs(j),ys(j))
      pot1ex(i) = pot1ex(i) + sigma1(j)*ff
      pot2ex(i) = pot2ex(i) + sigma2(j)*ff
      pot3ex(i) = pot3ex(i) + sigma3(j)*ff
      pot4ex(i) = pot4ex(i) + sigma4(j)*ff

      enddo

      ra = ra + pot1ex(i)**2 + pot2ex(i)**2 + pot3ex(i)**2 + 
     1           pot4ex(i)**2
      erra = erra + (pot1(i)-pot1ex(i))**2 + 
     1              (pot2(i)-pot2ex(i))**2 +  
     1              (pot3(i)-pot3ex(i))**2 +  
     1              (pot4(i)-pot4ex(i))**2   
      enddo

      erra = sqrt(erra/ra)

      call prin2('erra=*',erra,1)


       return
       end
c------------------------------------------------------------------

      subroutine fasttsne2d(n,ndim,xs,ys,sigma1,sigma2,sigma3,sigma4,
     1  nboxes,nterms,xpts,
     1  boxl,boxr,ioffst,prods,nn,nfour,irearr,zkvals,zksave,
     2  workspace,pot1,pot2,pot3,pot4,tinfo)

      implicit real *8 (a-h,o-z)
      integer irearr(*)
      real *8 xs(*),ys(*),boxl(2,*),boxr(2,*),prods(*)
      real *8 pot1(*),pot2(*),pot3(*),pot4(*)
      real *8 sigma1(*),sigma2(*),sigma3(*),sigma4(*)

      real *8 xpts(*)
      real *8, allocatable :: mpol1(:),mpol2(:),mpol3(:),mpol4(:)
      real *8, allocatable :: loc1(:),loc2(:),loc3(:),loc4(:)

      real *8, allocatable :: mpol(:),loc(:,:),svalsx(:),
     1   svalsy(:),svalsxt(:),svalsyt(:)
      real *8, allocatable :: mpolsort(:)
      real *8 tinfo(*)
      complex *16 zkvals(nfour,*),zksave(*)
      real *8 workspace(*)
      complex *16, allocatable :: zmpol(:,:)

      real *8, allocatable :: xsp(:),ysp(:)
      integer ioffst(*),nmax,omp_get_max_threads,omp_get_num_procs

      allocate(svalsx(nterms*n),svalsy(nterms*n))
      allocate(svalsxt(nterms*n),svalsyt(nterms*n))
      allocate(mpol(nn*ndim),loc(nn,ndim),mpolsort(nn))
      allocate(xsp(n),ysp(n))
      

c
cc     transform and x coordinate to be in [-1,1]
c
      allocate(mpol1(nn),mpol2(nn),mpol3(nn),mpol4(nn))
      allocate(loc1(nn),loc2(nn),loc3(nn),loc4(nn))


      bsize = boxr(1,1) - boxl(1,1)

      t1 = second()
C$       t1 = omp_get_wtime() 

      bsizeinv = 2/bsize
      do 1050 ibox=1,nboxes
      do 1000 i=ioffst(ibox)+1,ioffst(ibox+1)

       xsp(i) = (xs(i)-boxl(1,ibox))*bsizeinv-1
       ysp(i) = (ys(i)-boxl(2,ibox))*bsizeinv-1


 1000 continue
 1050 continue

cc       call prin2('xsp=*',xsp,12)
cc       call prin2('prods=*',prods,nterms)
cc       call prin2('xpts=*',xpts,nterms)

      call getljvals2(nterms,xpts,prods,n,xsp,svalsx)
      call getljvals2(nterms,xpts,prods,n,ysp,svalsy)

cc      call prin2('svals=*',svals,16)
c
cc      compute mpol
c


      do 1300 i=1,nn

      mpol1(i) = 0
      mpol2(i) = 0
      mpol3(i) = 0
      mpol4(i) = 0
 1300 continue     
      t2 = second()
C$       t2 = omp_get_wtime()      


      tinfo(1) = t2-t1
cc      call prin2('compute interp time=*',t2-t1,1)


cc      call prinf('ioffst=*',ioffst,nboxes)


      t1 = second()
C$      t1 = omp_get_wtime()     


C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,istart,impx,impy,ii,i,rtmp)
      do ibox=1,nboxes

      istart = (ibox-1)*nterms*nterms

      do impx=1,nterms
      do impy=1,nterms 

      ii = istart + (impx-1)*nterms+impy

      do i=ioffst(ibox)+1,ioffst(ibox+1)

      rtmp = svalsy((impy-1)*n+i)*svalsx((impx-1)*n+i)


      mpol1(ii) = mpol1(ii) + rtmp*sigma1(i)
      mpol2(ii) = mpol2(ii) + rtmp*sigma2(i)
      mpol3(ii) = mpol3(ii) + rtmp*sigma3(i)
      mpol4(ii) = mpol4(ii) + rtmp*sigma4(i)

      enddo
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO    

      t2 = second()
C$      t2 = omp_get_wtime()      
cc      call prin2('time for form mpol=*',t2-t1,1)
      tinfo(2) = t2-t1

      do i=1,nn

      mpol(4*(i-1)+1) = mpol1(i)
      mpol(4*(i-1)+2) = mpol2(i)
      mpol(4*(i-1)+3) = mpol3(i)
      mpol(4*(i-1)+4) = mpol4(i)

      enddo

c
cc      mpol to loc
c


       allocate(zmpol(nfour,nfour))
       t1 = second()
C$       t1 = omp_get_wtime()      

       call mpoltoloc(nn,nfour,mpolsort,zmpol,irearr,zkvals,
     1  workspace,zksave,mpol1,loc1)
       call mpoltoloc(nn,nfour,mpolsort,zmpol,irearr,zkvals,
     1  workspace,zksave,mpol2,loc2)
       call mpoltoloc(nn,nfour,mpolsort,zmpol,irearr,zkvals,
     1  workspace,zksave,mpol3,loc3)
       call mpoltoloc(nn,nfour,mpolsort,zmpol,irearr,zkvals,
     1  workspace,zksave,mpol4,loc4)

       t2=second()
C$       t2 = omp_get_wtime()       

       tinfo(3) = t2-t1

c
cc   evaluate potential
c

      t1 = second()
C$      t1 = omp_get_wtime()     

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,n
      pot1(i) = 0
      pot2(i) = 0
      pot3(i) = 0
      pot4(i) = 0

      enddo
C$OMP END PARALLEL DO      

C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP$PRIVATE(ibox,istart,j,l,ii,rtmp1,rtmp2,rtmp3,rtmp4,i,rtmp)
      do ibox=1,nboxes
      istart = (ibox-1)*nterms*nterms

      do j=1,nterms
      do l=1,nterms
      ii = (j-1)*nterms + l

      rtmp1 = loc1(istart+ii)
      rtmp2 = loc2(istart+ii)
      rtmp3 = loc3(istart+ii)
      rtmp4 = loc4(istart+ii)


      do i=ioffst(ibox)+1,ioffst(ibox+1)

      rtmp = svalsx((j-1)*n+i)*svalsy((l-1)*n+i)

      pot1(i) = pot1(i) + rtmp*rtmp1
      pot2(i) = pot2(i) + rtmp*rtmp2
      pot3(i) = pot3(i) + rtmp*rtmp3
      pot4(i) = pot4(i) + rtmp*rtmp4

      enddo
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO      



      t2 = second()
C$       t2 = omp_get_wtime()      
      tinfo(4)  =t2-t1

      ttot = 0
      do i=1,4

      ttot = ttot + tinfo(i)

      enddo

      call prin2('tinfo=*',tinfo,4)
      call prin2('total time=*',ttot,1)



      return
      end
c----------------------------------------------------

       subroutine mpoltoloc(nn,nfour,mpolsort,zmpol,irearr,
     1  zkvals,workspace,zksave,mpol,loc)
       implicit real *8 (a-h,o-z)

       real *8 mpolsort(*),loc(*),mpol(*)
       complex *16 zkvals(nfour,*),zksave(*)
       complex *16 zmpol(nfour,*)
       integer irearr(*)
       real *8 workspace(*)


       call rsortrev(mpol,mpolsort,irearr,nn)

       do i=1,nfour
       do j=1,nfour

       zmpol(i,j) = 0

       enddo
       enddo

       nfourh = nfour/2
       
       do i=1,nfourh
       do j=1,nfourh

       ii = (i-1)*nfourh + j
       zmpol(i,j) = mpolsort(ii)

       enddo
       enddo

       call fourt2d_f_pre(zmpol,nfour,workspace,zksave)

       do i=1,nfour
       do j=1,nfour

       zmpol(i,j) = zmpol(i,j)*zkvals(i,j)

       enddo
       enddo

       call fourt2d_b_pre(zmpol,nfour,workspace,zksave)

       do i=1,nfourh
       do j=1,nfourh

       ii = (i-1)*nfourh+j

       mpolsort(ii) = real(zmpol(i,j))

       enddo
       enddo

       call rsort(mpolsort,loc,irearr,nn)



       return
       end
     

c------------------------------------------------------

      subroutine fker(f,xs,ys,xt,yt)
      implicit real *8 (a-h,o-z)
      

      rr = (xs-xt)**2 + (ys-yt)**2
      f = 1.0d0/(1.0d0+rr)**2

      return
      end

c--------------------------------------------------

      subroutine rsort(x,xsort,iarr,n)
      implicit real *8 (a-h,o-z)
      real *8 x(n),xsort(n)
      integer iarr(*)
 
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) 
      do i=1,n

      xsort(i) = x(iarr(i))

      enddo
C$OMP END PARALLEL DO      

      return
      end
c---------------------------------------------------
      
      subroutine rsortsrcall(x,xsort,y,ysort,sigma1,sigma1sort,
     1   sigma2,sigma2sort,sigma3,sigma3sort,sigma4,sigma4sort,
     2   iarr,n)
      implicit real *8 (a-h,o-z)
      real *8 x(*),xsort(*),y(*),ysort(*)
      real *8 sigma1(*),sigma1sort(*)
      real *8 sigma2(*),sigma2sort(*)
      real *8 sigma3(*),sigma3sort(*)
      real *8 sigma4(*),sigma4sort(*)

      integer iarr(*)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,n
         xsort(i) = x(iarr(i))
         ysort(i) = y(iarr(i))
         sigma1sort(i) = sigma1(iarr(i))
         sigma2sort(i) = sigma2(iarr(i))
         sigma3sort(i) = sigma3(iarr(i))
         sigma4sort(i) = sigma4(iarr(i))
      enddo
C$OMP END PARALLEL DO      


      return
      end

c----------------------------------------------------      
c
      subroutine rsortrev(xsort,x,iarr,n)
      implicit real *8 (a-h,o-z)
      real *8 x(*),xsort(*)
      integer iarr(*)
     
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)     
      do i=1,n

      x(iarr(i))= xsort(i) 

      enddo
C$OMP END PARALLEL DO      

      return
      end
c----------------------------------------------

      subroutine rsortrevtargall(x1sort,x1,x2sort,x2,x3sort,x3,x4sort,
     1   x4,iarr,n)
      implicit real *8 (a-h,o-z)
      real *8 x1(*),x1sort(*)
      real *8 x2(*),x2sort(*)
      real *8 x3(*),x3sort(*)
      real *8 x4(*),x4sort(*)
      integer iarr(*)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,n
         x1(iarr(i)) = x1sort(i)
         x2(iarr(i)) = x2sort(i)
         x3(iarr(i)) = x3sort(i)
         x4(iarr(i)) = x4sort(i)
      enddo

C$OMP END PARALLEL DO



      return
      end

c-----------------------------------------------

      subroutine getprods(n,x,p)
      implicit real *8 (a-h,o-z)

      real *8 x(*),p(*)

      do i=1,n
      p(i) = 1

      do j=1,n

      if(i.ne.j) then

      p(i) = p(i)*(x(i)-x(j))

      endif
      
      enddo

c
cc      note products are inverses of what they used to be
c

      p(i) = 1/p(i)
      enddo

      return
      end
c--------------------      
      subroutine getljvals2(np,x,prods,n,y,vals)
      implicit real *8 (a-h,o-z)
      real *8 x(*),y(*),vals(*),prods(*)
      real *8 ydiff(np,n),yprods(n)


      thresh=1.0d-6

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,i)
      do j=1,n
      yprods(j) = 1
      do i=1,np

      ydiff(i,j) = y(j) - x(i)
      yprods(j) = yprods(j)*ydiff(i,j)

      enddo
      enddo
C$OMP END PARALLEL DO      

      do i=1,np
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,k)
      do j=1,n

      if(abs(ydiff(i,j)).ge.thresh) vals((i-1)*n+j) =
     1    yprods(j)*prods(i)/ydiff(i,j)
      if(abs(ydiff(i,j)).le.thresh) then

      vals((i-1)*n+j) = prods(i)
      
      do k=1,np

      if(i.ne.k) vals((i-1)*n+j) = vals((i-1)*n+j)*ydiff(k,j)

      enddo

      endif

      enddo
C$OMP END PARALLEL DO      
      enddo


      return
      end
c---------------------------------     

      subroutine getkermat(n,x,y,amat)

      implicit real *8 (a-h,o-z)
      real *8 x(*),y(*),amat(n,*)

      do i=1,n

      do j=1,n

      call fker(amat(i,j),x(i),y(i),x(j),y(j))

      enddo
      enddo

      return
      end
