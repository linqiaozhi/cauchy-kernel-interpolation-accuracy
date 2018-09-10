c 
c       This file contains 4 user-callable subroutines: fourt2d_f,
c       fourt2d_b, fourt2d_evalders, fourt2d_eval. Following is a
c       brief description of the said four subroutines.
c
c   fourt2d_b - evaluates the two-dimensional inverse Fourier 
c        Transform. The transform is dimensioned n * n, and n 
c        (user-specified) must be even. Conceptually, this 
c        subroutine converts the Fourier series of the function 
c        F: [0, 2pi] * [0, 2pi] into its values on the said square.
c
c   fourt2d_f - inverse of the subroutine fourt2d_b
c
c   fourt2d_evalder - Given the Fourier transform of a function on 
c        the square  [-1,1] \times [-1,1] and a point on the said 
c        square, fourt2d_evalders evaluates the said function at the 
c        said point, together with the first derivatives with respecxt 
c        to x, y of the said function at the said point. 
c
c   fourt2d_eval - same as fourt2d_evalder, but does not evaluate the
c        derivatives, and is faster.
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c 
        subroutine fourt2d_b(f,n,w)
        implicit real *8 (a-h,o-z)
        save
        complex *16 f(n,n)
        real *8 w(1)
c
c        This entry evaluates the two-dimensional inverse
c        Fourier Transform. The transform is dimensioned 
c        n * n, and n (user-specified) must be even. Conceptually,
c        this subroutine converts the Fourier series of the 
c        function F: [0, 2pi] * [0, 2pi] into its values on the 
c        said square. Please note that the notation used is
c        conventional on input, and unconventional on output.
c
c     EXPLANATION: 
c
c        On input,
c    coefficient number (0,0) is located in the  element (1,1) of array f; 
c    coefficient number (1,0) is located in the  element (2,1) of array f; 
c    coefficient number (2,0) is located in the  element (3,1) of array f; 
c    . . . 
c    coefficient number (-1,0) is located in the  element (n,1) of array f; 
c    coefficient number (-2,0) is located in the  element (n-1,1) of array f; 
c    coefficient number (-3,0) is located in the  element (n-2,1) of array f; 
c    . . . 
c                        
c        On output, 
c    value F(0,0) is located in the element f(n/2+1,f/2+1),
c    value F(h,0) is located in the element f(n/2+2,f/2+1),  
c    value F(2*h,0) is located in the element f(n/2+3,f/2+1),
c    
c    and so on.
c
        iw1=1
        lw1=2*n+4
c
        iwsave=iw1+lw1
        lwsave=4*n+30
c
        call fourt2d_b0(f,n,w(iw1),w(iwsave))
c
        return
c
c
c
c
        entry fourt2d_f(f,n,w)
c
c        This entry evaluates the two-dimensional forward
c        Fourier Transform. The transform is dimensioned 
c        n * n, and n (user-specified) must be even. Conceptually,
c        this subroutine converts the values on the square 
c        F: [0, 2pi] * [0, 2pi] of the function F: [0, 2pi] * [0, 2pi] 
c        into its (two-dimensional) Fourier series of. Please note 
c        that the notation used is unconventional on input, and 
c        conventional on output.
c
c     EXPLANATION: 
c
c        On input,
c    value F(0,0) is located in the element f(n/2+1,f/2+1);
c    value F(h,0) is located in the element f(n/2+2,f/2+1);  
c    value F(2*h,0) is located in the element f(n/2+3,f/2+1);
c    
c    and so on.
c        On output, 
c    coefficient number (0,0) is located in the element (1,1) of array f; 
c    coefficient number (1,0) is located in the element (2,1) of array f; 
c    coefficient number (2,0) is located in the element (3,1) of array f; 
c    . . . 
c    coefficient number (-1,0) is located in the element (n,1) of array f; 
c    coefficient number (-2,0) is located in the element (n-1,1) of array f; 
c    coefficient number (-3,0) is located in the element (n-2,1) of array f; 
c    . . . 
c                        
        iw1=1
        lw1=2*n+4
c
        iwsave=iw1+lw1
        lwsave=4*n+30
c
        call fourt2d_f0(f,n,w(iw1),w(iwsave) )
c
        return
        end
c 
c 
c 
c 
c 
        subroutine fourt2d_f0(f,n,w1,wsave)
        implicit real *8 (a-h,o-z)
        save
        complex *16 f(n,n),w1(1),wsave(1)
c
c       initialize the FFT routine
c
        call DCFFTI(N,WSAVE)
c
c        transform along the y direction
c
        do 2000 i=1,n
c
        call fourt2d_setzero(w1,n*2)
c
        do 1200 j=1,n/2
        w1(j)=f(j+n/2,i)
c
        w1(n/2+j)=f(j,i)
 1200 continue
c
        call DCFFTf(n,w1,WSAVE)
c
        do 1400 j=1,n
c
        f(j,i)=w1(j)
 1400 continue
c
 2000 continue
c
c        transform along the x direction
c
        do 3000 i=1,n
c
        call fourt2d_setzero(w1,n*2)
c
        do 2200 j=1,n/2
        w1(j)=f(i,j+n/2)
c
        w1(n/2+j)=f(i,j)
 2200 continue
c
        call DCFFTf(n,w1,WSAVE)
c
        do 2400 j=1,n
c
        f(i,j)=w1(j)
 2400 continue
c
 3000 continue
c
        return
        end
c 
c 
c 
c 
c 
        subroutine fourt2d_b0(f,n,w1,wsave)
        implicit real *8 (a-h,o-z)
        save
        complex *16 f(n,n),w1(1),wsave(1)
c
c       initialize the FFT routine
c
        call DCFFTI(N,WSAVE)
c
c        transform along the y direction
c
        do 2000 i=1,n
c
        do 1200 j=1,n
        w1(j)=f(j,i)
 1200 continue
c
        call DCFFTb(n,w1,WSAVE)
c
        do 1600 j=1,n/2
c
        f(n/2+j,i)=w1(j)
        f(j,i)=w1(n/2+j)
 1600 continue
c
 2000 continue
c
c        transform along the x direction
c
        do 3000 i=1,n
c
        do 2200 j=1,n
        w1(j)=f(i,j)
 2200 continue
c
        call DCFFTb(n,w1,WSAVE)
c
        do 2600 j=1,n/2
c
        f(i,n/2+j)=w1(j)
        f(i,j)=w1(n/2+j)
 2600 continue
c
 3000 continue
c
        return
        end
c
c
c
c
c
        subroutine fourt2d_evalders(coefs,n,xx,yy,val,derxx,deryy)
        implicit real *8 (a-h,o-z)
        save
        complex *16 coefs(n,n),val,cd,ima,cx1,cx2,cy1,cy2,
     1      cxp,cxp_inv,cyp,cyp_inv,derxx,deryy
        data ima/(0.0d0,1.0d0)/
c
c        Given the Fourier transform of a function on the
c        square [-1,1] \times [-1,1] and a point on the 
c        said square, this subroutine evaluates the said
c        function at the said point, together with the 
c        first derivatives with respecxt to x, y of the 
c        said function at the said point.
c
c               Input parameters:
c
c  coefs - Fourier transform of the function to be evaluated, 
c        as produced by the subroutine fourt2d_f (see)
c  n - the dimensionality of the transform. MUST BE EVEN
c  (xx,yy) - the point in the square [-1,1]*[-1,1] at which 
c        the function and its derivatives are to be evaluated
c
c               Output parameters:
c
c  val - the function
c  derxx - the x-derivative
c  deryy - the y-derivative
c
c
c        . . . evaluate the user-supplied two-dimensional Fourier
c              series at the point (x,y) on the square [-1,1]^2
c
        done=1
        pi=atan(done)*4
        x=xx*pi
        y=yy*pi
c
        derxx=0
        deryy=0
        cd=0
        cxp=exp(ima*x)
        cxp_inv=1/cxp
        cx1=1
        cx2=cxp_inv
c
        do 1400 i=1,n/2
c
        cyp=exp(ima*y)
        cyp_inv=1/cyp
        cy1=1
        cy2=cyp_inv
c
        do 1200 j=1,n/2
c
cc        cd=cd+coefs(i,j)*exp(ima*(i-1)*x)*exp(ima*(j-1)*y)
cc        cd=cd+coefs(n-i+1,j)*exp(ima*(-i)*x)*exp(ima*(j-1)*y)
cc        cd=cd+coefs(i,n-j+1)*exp(ima*(i-1)*x)*exp(ima*(-j)*y)
cc        cd=cd+coefs(n-i+1,n-j+1)*exp(ima*(-i)*x)*exp(ima*(-j)*y)
c
        cd=cd+coefs(i,j)*cx1*cy1
        cd=cd+coefs(n-i+1,j)*cx2*cy1
        cd=cd+coefs(i,n-j+1)*cx1*cy2
        cd=cd+coefs(n-i+1,n-j+1)*cx2*cy2
c
        derxx=derxx+coefs(i,j)*cx1*cy1*ima*(i-1)
        derxx=derxx-coefs(n-i+1,j)*cx2*cy1*ima*i
        derxx=derxx+coefs(i,n-j+1)*cx1*cy2*ima*(i-1)
        derxx=derxx-coefs(n-i+1,n-j+1)*cx2*cy2*ima*i
c
        deryy=deryy+coefs(i,j)*cx1*cy1*ima*(j-1)
        deryy=deryy+coefs(n-i+1,j)*cx2*cy1*ima*(j-1)
        deryy=deryy-coefs(i,n-j+1)*cx1*cy2*ima*j
        deryy=deryy-coefs(n-i+1,n-j+1)*cx2*cy2*ima*j
c
        cy1=cy1*cyp
        cy2=cy2*cyp_inv
c
 1200 continue
c
        cx1=cx1*cxp
        cx2=cx2*cxp_inv
 1400 continue
c
        val=cd/n/n
        derxx=derxx/n/n*pi
        deryy=deryy/n/n*pi
c
        return
        end
c
c
c
c
c
        subroutine fourt2d_eval(coefs,n,xx,yy,val)
        implicit real *8 (a-h,o-z)
        save
        complex *16 coefs(n,n),val,cd,ima,cx1,cx2,cy1,cy2,
     1      cxp,cxp_inv,cyp,cyp_inv
        data ima/(0.0d0,1.0d0)/
c
c        Given the Fourier transform of a function on the
c        square [-1,1] \times [-1,1] and a point on the 
c        said square, this subroutine evaluates the said
c        function at the said point. If the gentle user 
c        wants also the partial derivatives of the function,
c        he is advised to use the subroutine fouer2d_evalders
c        (see).
c        
c
c               Input parameters:
c
c  coefs - Fourier transform of the function to be evaluated, 
c        as produced by the subroutine fourt2d_f (see)
c  n - the dimensionality of the transform. MUST BE EVEN
c  (xx,yy) - the point in the square [-1,1]*[-1,1] at which 
c        the function and its derivatives are to be evaluated
c
c               Output parameters:
c
c  val - the function
c
c
c        . . . evaluate the user-supplied two-dimensional Fourier
c              series at the point (x,y) on the square [-1,1]^2
c
        done=1
        pi=atan(done)*4
        x=xx*pi
        y=yy*pi
c
        cd=0
        cxp=exp(ima*x)
        cxp_inv=1/cxp
        cx1=1
        cx2=cxp_inv
c
        do 1400 i=1,n/2
c
        cyp=exp(ima*y)
        cyp_inv=1/cyp
        cy1=1
        cy2=cyp_inv
c
        do 1200 j=1,n/2
c
cc        cd=cd+coefs(i,j)*exp(ima*(i-1)*x)*exp(ima*(j-1)*y)
cc        cd=cd+coefs(n-i+1,j)*exp(ima*(-i)*x)*exp(ima*(j-1)*y)
cc        cd=cd+coefs(i,n-j+1)*exp(ima*(i-1)*x)*exp(ima*(-j)*y)
cc        cd=cd+coefs(n-i+1,n-j+1)*exp(ima*(-i)*x)*exp(ima*(-j)*y)
c
        cd=cd+coefs(i,j)*cx1*cy1
        cd=cd+coefs(n-i+1,j)*cx2*cy1
        cd=cd+coefs(i,n-j+1)*cx1*cy2
        cd=cd+coefs(n-i+1,n-j+1)*cx2*cy2
c
        cy1=cy1*cyp
        cy2=cy2*cyp_inv
c
 1200 continue
c
        cx1=cx1*cxp
        cx2=cx2*cxp_inv
 1400 continue
c
        val=cd/n/n
        return
        end
c 
c 
c 
c 
c 
        subroutine fourt2d_setzero(x,n)
        implicit real *8 (a-h,o-z)
        save
        real *8 x(1)
c
        do 1200 i=1,n
c
        x(i)=0
 1200 continue
c
        return
        end
c
c 
        subroutine fourt2d_b_pre(f,n,w,wsave)
        implicit real *8 (a-h,o-z)
        save
        complex *16 f(n,n)
        complex *16 wsave(*)
        real *8 w(1)
c
c        This entry evaluates the two-dimensional inverse
c        Fourier Transform. The transform is dimensioned 
c        n * n, and n (user-specified) must be even. Conceptually,
c        this subroutine converts the Fourier series of the 
c        function F: [0, 2pi] * [0, 2pi] into its values on the 
c        said square. Please note that the notation used is
c        conventional on input, and unconventional on output.
c
c     EXPLANATION: 
c
c        On input,
c    coefficient number (0,0) is located in the  element (1,1) of array f; 
c    coefficient number (1,0) is located in the  element (2,1) of array f; 
c    coefficient number (2,0) is located in the  element (3,1) of array f; 
c    . . . 
c    coefficient number (-1,0) is located in the  element (n,1) of array f; 
c    coefficient number (-2,0) is located in the  element (n-1,1) of array f; 
c    coefficient number (-3,0) is located in the  element (n-2,1) of array f; 
c    . . . 
c                        
c        On output, 
c    value F(0,0) is located in the element f(n/2+1,f/2+1),
c    value F(h,0) is located in the element f(n/2+2,f/2+1),  
c    value F(2*h,0) is located in the element f(n/2+3,f/2+1),
c    
c    and so on.
c
        iw1=1
        lw1=2*n+4
c
        iwsave=iw1+lw1
        lwsave=4*n+30
c
        call fourt2d_b0_pre(f,n,w,wsave)
c
        return
        end
c
c
c
c
        subroutine fourt2d_f_pre(f,n,w,wsave)
        implicit real *8 (a-h,o-z)
        complex *16 f(n,n)
        real *8 w(*)
        complex *16 wsave(*)
c
c        This entry evaluates the two-dimensional forward
c        Fourier Transform. The transform is dimensioned 
c        n * n, and n (user-specified) must be even. Conceptually,
c        this subroutine converts the values on the square 
c        F: [0, 2pi] * [0, 2pi] of the function F: [0, 2pi] * [0, 2pi] 
c        into its (two-dimensional) Fourier series of. Please note 
c        that the notation used is unconventional on input, and 
c        conventional on output.
c
c     EXPLANATION: 
c
c        On input,
c    value F(0,0) is located in the element f(n/2+1,f/2+1);
c    value F(h,0) is located in the element f(n/2+2,f/2+1);  
c    value F(2*h,0) is located in the element f(n/2+3,f/2+1);
c    
c    and so on.
c        On output, 
c    coefficient number (0,0) is located in the element (1,1) of array f; 
c    coefficient number (1,0) is located in the element (2,1) of array f; 
c    coefficient number (2,0) is located in the element (3,1) of array f; 
c    . . . 
c    coefficient number (-1,0) is located in the element (n,1) of array f; 
c    coefficient number (-2,0) is located in the element (n-1,1) of array f; 
c    coefficient number (-3,0) is located in the element (n-2,1) of array f; 
c    . . . 
c                        
        iw1=1
        lw1=2*n+4
c
        iwsave=iw1+lw1
        lwsave=4*n+30
c
        call fourt2d_f0_pre(f,n,w,wsave)
c
        return
        end
c 
c 
c 
c 
c 
        subroutine fourt2d_f0_pre(f,n,w1,wsave)
        implicit real *8 (a-h,o-z)
        save
        complex *16 f(n,n),w1(1),wsave(1)
c
c       initialize the FFT routine
c
c
c        transform along the y direction
c
        do 2000 i=1,n
c
        call fourt2d_setzero(w1,n*2)
c
        do 1200 j=1,n/2
        w1(j)=f(j+n/2,i)
c
        w1(n/2+j)=f(j,i)
 1200 continue
c
        call DCFFTf(n,w1,WSAVE)
c
        do 1400 j=1,n
c
        f(j,i)=w1(j)
 1400 continue
c
 2000 continue
c
c        transform along the x direction
c
        do 3000 i=1,n
c
        call fourt2d_setzero(w1,n*2)
c
        do 2200 j=1,n/2
        w1(j)=f(i,j+n/2)
c
        w1(n/2+j)=f(i,j)
 2200 continue
c
        call DCFFTf(n,w1,WSAVE)
c
        do 2400 j=1,n
c
        f(i,j)=w1(j)
 2400 continue
c
 3000 continue
c
        return
        end
c 
c 
c 
c 
c 
        subroutine fourt2d_b0_pre(f,n,w1,wsave)
        implicit real *8 (a-h,o-z)
        save
        complex *16 f(n,n),w1(1),wsave(1)
c
c       initialize the FFT routine
c
c        transform along the y direction
c
        do 2000 i=1,n
c
        do 1200 j=1,n
        w1(j)=f(j,i)
 1200 continue
c
        call DCFFTb(n,w1,WSAVE)
c
        do 1600 j=1,n/2
c
        f(n/2+j,i)=w1(j)
        f(j,i)=w1(n/2+j)
 1600 continue
c
 2000 continue
c
c        transform along the x direction
c
        do 3000 i=1,n
c
        do 2200 j=1,n
        w1(j)=f(i,j)
 2200 continue
c
        call DCFFTb(n,w1,WSAVE)
c
        do 2600 j=1,n/2
c
        f(i,n/2+j)=w1(j)
        f(i,j)=w1(n/2+j)
 2600 continue
c
 3000 continue
c
        return
        end
c
c
