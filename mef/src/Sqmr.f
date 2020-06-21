c *********************************************************************  
      subroutine rpsqmr(neq   ,nequ   ,nad   ,ia  ,ja
     1                 ,ad    ,au     ,al    ,m   ,b  ,x  
     2                 ,t     ,r      ,q     ,d   
     3                 ,tol   ,maxit
     4                 ,matvec,dot
     5                 ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6                 ,i_xfi ,i_rcvsi,i_dspli
     7                 ,fprint,flog   ,fhist  ,fnew
     8                 ,nprcs  ,mpi   ,fsqmr)
c **********************************************************************
c * Data de criacao    : 28/06/2016                                    *
c * Data de modificaco : 17/06/2020                                    * 
c * ------------------------------------------------------------------ *   
c * RPSQRM : Solucao de sistemas de equacoes pelo metodo QMR simetrico *
c * diagonal a direita                                                 *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq       - numero de equacoes                                     *
c * nequ      - numero de equacoes no bloco Kuu                        *
c * nad       - numero de termos nao nulos no bloco Kuu e Kpu  ou K    *
c * ia(*)     - ponteiro do formato CSR                                *
c * ja(*)     - ponteiro das colunas no formato CSR                    *
c * ad(neq)   - diagonal da matriz A                                   *
c * au(*)     - parte triangular superior de A                         *
c * al(*)     - parte triangular inferior de A                         *
c * m(*)      - precondicionador diagonal                              *
c * b(neqovlp)- vetor de forcas                                        *
c * x(neqovlp)- chute inicial                                          *
c * r(neq)    - arranjo local de trabalho                              *
c * q(neqovlp)- arranjo local de trabalho                              *
c * t(neq)    - arranjo local de trabalho                              *
c * d(neq)    - arranjo local de trabalho                              *
c * tol       - tolerancia de convergencia                             *
c * maxit     - numero maximo de iteracoes                             *
c * matvec    - nome da funcao externa para o produto matrix-vetor     *
c * dot       - nome da funcao externa para o produto escalar          *
c * my_id     - MPI                                                    *
c * neqf1i    - MPI                                                    *
c * neqf2i    - MPI                                                    *
c * neq_doti  - MPI                                                    *
c * i_fmap    - MPI                                                    *
c * i_xfi     - MPI                                                    *
c * i_rvcs    - MPI                                                    *
c * i_dspli   - MPI                                                    *
c * fprint    - saida na tela                                          *
c * flog      - log do arquivo de saida                                *
c * fhist     - log dos resuduos por iteracao                          *
c * fnew      - .true.  -> x0 igual a zero                             *
c *             .false. -> x0 dado                                     *
c * mpi       - true|false                                             *
c * nprcs     - numero de processos mpi                                *
c * fsqmr    - nao definido                                            * 
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq) - modificado                                                *
c * ad(*),al(*),au(*) - inalterados                                    *
c * fsqmr    - true para falha do sqmr                                 *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * Fonte: A New Krylov-subspace method for symmetric indefinite       * 
c * linear systems                                                     *
c *                                                                    *
c * vetor que sera multiplicado pela matriz necessitam ter a dimensao  *
c * neqovlp                                                            *
c * Ex: Ax=y                                                           *
c * x(neqolp) e u(neq)                                                 *
c * vetor que usa subrotina comunicate necessitam ter a dimensao       *
c * neqovlp                                                            * 
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
c ...
      real*8 get_time
c ......................................................................
c ... mpi
      logical mpi 
      integer neqf1i,neqf2i,neq_doti,nprcs,ierr
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................      
      integer neq,nequ,maxit,i,j,jj
      integer*8 ia(*),nad
      integer ja(*),my_id
      real*8 ad(*),au(*),al(*),b(*),m(*),x(*)
      real*8 r(*),t(*),q(*),d(*)
      real*8 dot,tol,conv,xkx,norm,alpha,beta,tmp1,tmp2,tau,ro,vn,v0
      real*8 sigma,cn,norm_b 
      real*8 norm_r
      real*8 time0,time
      real*8 dum1
      logical flog,fprint,fnew,fhist,fsqmr
c ...
      real*8 flop_sqrm
      real*8  mflops,vmean
c .....................................................................
      external matvec,dot
c ======================================================================
      time0 = get_time()
c ......................................................................
c
c ...
      do 5 i = 1, neq
        if(ad(i) .eq. 0.d0 ) then
          write(*,1000) i,ad(i)
          call stop_mef()
        endif 
   5  continue
c ......................................................................
c
c ... Chute inicial:
c
      if(fnew) then  
        do 10 i = 1, neq
          x(i)  = 0.d0
   10   continue
      endif     
c .......................................................................
c
c ... conv = tol * |b| 
      norm_b = dot(b,b,neq_doti)
      conv = tol*dsqrt(dabs(norm_b))
c .......................................................................
c  
c ... Ax0                                                            
      call matvec(neq,ia ,ja,ad ,al,x  ,t)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r(i) = b(i) - t(i)
c ... q = t
         q(i) = r(i) * m(i)
c ... d = 0.0
         d(i) = 0.d0
  100 continue
c ... ( r,r ) 
      tau = dsqrt(dot(r,r,neq_doti))
c ... ( r,q ) 
      ro  = dot(r,q,neq_doti)
c ......................................................................
c
c ...
      v0 = 0.0d0
c ......................................................................
      jj = 1
      do 230 j = 1, maxit
c ... t = Aq(j-1)
         call matvec(neq,ia,ja,ad,al,q,t)
c .....................................................................
c
c ... sigma = ( q(j-1),t)         
         sigma = dot(q,t,neq_doti)
         if( sigma .eq. 0.0d0) then
           print*,"RSQMR fail (sigma)!"
           call stop_mef()  
         endif  
c .....................................................................
c
c ... alpha(j-1) = ro(j-1)/sigma(j-1)
         alpha = ro/sigma        
c .....................................................................
c
c ...
         do 210 i = 1, neq
c ... r(j) = r(j-1) - alpha(j-1)*t
            r(i) = r(i) - alpha * t(i)
  210    continue
c .....................................................................
c
c ... v(j) = ||r||/tau(j-1)
         vn   = dsqrt(dot(r,r,neq_doti))/tau
c ... c(j) = 1/sqrt(1+v(j)*v(j) 
         cn   = 1.0d0/dsqrt(1.d0+vn*vn)
c ... tau(j) = tau(j-1)*v(j)*c(j)
         tau = tau*vn*cn
c .....................................................................
c
c ... tau(j) = (c(j)*c(j)*v(j-1)*v(j-1)) d(j-1)
c            + c(j)*c(j)*alpha(j-1) * q(j-1)
         tmp1 = cn*cn*v0*v0
         tmp2 = cn*cn*alpha 
         do 215 i = 1, neq
c ... d(j) = (cj*cj*vj*vj) d(j-1) + cj*cj*alpha(j-1) * q(j-1)
           d(i) = tmp1*d(i) + tmp2*q(i) 
c ... x(j) = x(j-1) + d(j)
           x(i) = x(i) + d(i) 
c ... u(j) = M(-1)r
           t(i) = r(i)*m(i)
  215    continue 
c .....................................................................
c
c ...
         v0 = vn
c ......................................................................
c
c ... 
         if(ro .eq. 0.d0) then
           print*,"RSQMR fail (ro)!"
           call stop_mef()  
         endif  
c .....................................................................
c
c ... (r,u) 
         tmp1 = dot(r,t,neq_doti) 
c ... beta = (r,u)/ro(j-1)
         beta = tmp1/ro
c ...
         ro = tmp1
c .....................................................................
c
c ...
         do 220 i = 1, neq
c ... q(j+1) = (M-1)r(j) + beta*q(j-1) = t + beta*q(j-1)
            q(i) = t(i) + beta * q(i)
  220    continue
c .....................................................................
c
c ...
         norm_r = dsqrt(dot(r,r,neq_doti))
         if(fhist) then  
           if(my_id .eq.0) write(18,1500),j,norm_r/norm_b,norm_r 
         endif 
c .....................................................................
c
c ...
         if (norm_r .lt. conv) goto 300
c ......................................................................
         if( jj .eq. 500) then
           jj = 0
           if(my_id .eq.0 .and. fprint ) write(*,1300),j,norm_r,conv 
         endif  
         jj = jj + 1
c ......................................................................
  230 continue
c ----------------------------------------------------------------------
      if(my_id .eq. 0 .and. fprint ) then 
        write(*,1200) maxit,norm_r
        write(*,*)' *** Switching to Gmres(l) !'
        if(flog) write(10,1200) maxit,norm_r
      endif
      fsqmr = .true.
      return
  300 continue
c
c ... Energy norm:  x*Kx
      call matvec(neq,ia,ja,ad,al,x,t)
      xkx = dot(x,t,neq_doti)
c ......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r =(b - Ax) (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - t(i)
  310 continue
      norm_r = dot(r,r,neq_doti)
      norm_r = dsqrt(norm_r)
      if( norm_r .gt. conv ) then
         if(my_id .eq.0 .and. fprint )then
           write(*,1400) norm_r,conv
         endif 
      endif
c ......................................................................
      time = get_time()
      time = time-time0
c ......................................................................
c 
c ...  
      mflops = (flop_sqrm(neq,nad,j,2,mpi)*1.d-06)/time  
c ......................................................................
c
c ...
      if(my_id .eq.0 .and. fprint )then
        if(mpi) then
          write(*,1110)tol,conv,j,xkx,norm,norm_r,mflops,time
        else
          write(*,1100)tol,conv,neq,nad,j,xkx,norm,norm_r,mflops,time
        endif
      endif
c ......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,3(a,d20.10),2(a,f20.2))')
     .       'RPSQMR: ',' it ',j, ' x * Kx ',xkx,' ||x|| ',norm
     .      ,' tol ',tol,' Mflops ',mflops,' time ',time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA RPSQMR:',/,5x,'Diagonal coefficient ' 
     . '- equation ',i9,d20.10)
 1100 format(' (RPSQMR) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||b||          = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'Mflops               = ',f20.2/
     . 5x,'CPU time (s)         = ',f20.2/)
 1110 format(' (RPSQMR_MPI) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||b||          = ',d20.6/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'Mflops               = ',f20.2/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,d20.10,
     .        ' iterations !',/)
 1300 format (' RPSQMR:',5x,'It',i7,5x,2d20.10)
 1400 format (' RPSQMR:',1x,'Explicit residual > tol * ||b|| :'
     .       ,1x,d20.10,1x,d20.10)
 1500 format ( 5x,i7,5x,2es20.10)
      end
c **********************************************************************