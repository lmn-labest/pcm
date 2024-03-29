c *********************************************************************
c * Metodos iterativos para solucao de sistemas lineares              *
c * ----------------------------------------------------------------- *
c * simetricos:                                                       *
c * ----------------------------------------------------------------- *
c * CG   - gradiente conjugados                                       *
c *                                                                   *
c * PCG  - gradiente conjugados com precondicionador diagonal         *
c *                                                                   *
c * BCG  - gradiente conjugados com precondicionador bloco diagonal   *
c *                                                                   *
c * ICCG - gradiente conjugados com precondicionador de fatoracoes    *   
c * incompletas  (LLT e LDLt)                                         *
c *                                                                   *
c *                                                                   *
c * SQRM - QRM simetrico                                              *
c *                                                                   *
c * RSQRM - QRM simetrico com precondicionador diagonal a direita     *
c *                                                                   *
c * LSQRM - QRM simetrico com precondicionador diagonal a esquerda    *
c *                                                                   *
c * ----------------------------------------------------------------- *
c * nao-simetricos:                                                   *
c * ----------------------------------------------------------------- *
c * bicgstab - gradiente bi-conjugados estabilizados                  *
c *                                                                   *
c * pbicgstab - gradiente bi-conjugados estabilizados  com            *
c * precondicionador diagonal                                         *
c *                                                                   *
c * icbicgstab - gradiente bi-conjugados estabilizados fatoracoes     *             
c * incompletas (LLT e LDLt)                                          *  
c *                                                                   *
c * bicgstabl2 - gradiente bi-conjugados estabilizados (l=2)          *
c *                                                                   *
c * pbicgstabl2- gradiente bi-conjugados estabilizados (l=2) com      *
c * precondicionador diagonal                                         *
c *                                                                   *
c * gmres(m) - GMRES com precondicionador diagonal                    *
c *                                                                   *
c * gmres2(m) - GMRES com precondicionador diagonal                   *
c *   (Ry = e) versao coluna                                          *
c *   (x = Vy) versao coluna                                          *
c *   versao com refined modified gram-schmidt                        * 
c *                                                                   *
c * block_it_pcg - resolucao iterativa do problema poro mecanico com  *
c * matriz blocada | Kuu Kup | onde kup = -kpu                        * 
c *                | kpu Kpp |                                        *
c * ----------------------------------------------------------------- *
c *********************************************************************  
      subroutine cg(neq   ,nequ  ,nad     ,ia      ,ja
     .             ,ad    ,au    ,al      ,b       ,x
     .             ,z     ,r     ,p       ,tol     ,maxit
     .             ,matvec,dot
     .             ,my_id ,neqf1i ,neqf2i ,neq_doti,i_fmapi
     .             ,i_xfi ,i_rcvsi,i_dspli
     .             ,fprint,flog   ,fhist  ,fnew)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 15/12/2016                                    * 
c * ------------------------------------------------------------------ *   
c * Subroutine CG : Solucao de sistemas de equacoes pelo metodo dos    *    
c * gradientes conjugados                                              *
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
c * b(neqovlp)- vetor de forcas                                        *
c * x(neqovlp)- chute inicial                                          *
c * z(neq)    - arranjo local de trabalho                              *
c * r(neq)    - arranjo local de trabalho                              *
c * p(neqovlp)- arranjo local de trabalho                              *
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
c * fnew      - .true.  x0 igual a zero                                *
c *            .false. x0 dado                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * ad(*),al(*),au(*) e b - inalterados                                *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
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
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................      
      integer neq,nequ,maxit,i,j,jj
      integer*8 ia(*),nad
      integer ja(*),my_id
      real*8  ad(*),au(*),al(*),x(*),b(*)
      real*8  r(*),z(*),p(*)
      real*8  dot,tol,conv,xkx,norm,d,di,alpha,beta,tmp,norm_b
      real*8  time0,time
      real*8 dum1
      logical flog,fprint,fnew,fhist
      external matvec,dot
c ======================================================================
      time0 = get_time()
c ......................................................................
c
c ...
      do 5 i = 1, neq
        if(ad(i) .eq. 0.d0 ) then
          write(*,1000) i
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
      d      = dot(b,b,neq_doti)
      norm_b = dsqrt(d)
      conv   = tol*dsqrt(d)
c .......................................................................
c  
c ... Ax0
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)  
     .           ,x,z 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r(i) = b(i) - z(i)
c ... p0 = r0
         p(i) = r(i)
  100 continue
c ... ( r(0),r(0) )
      d    = dot(r,r,neq_doti)
c ......................................................................
      jj = 1
      do 230 j = 1, maxit
c ... z = Ap(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .              ,p,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli 
     .              ,dum1)
c .....................................................................
c
c ... alpha = ( r(j),r(j) ) / ( Ap(j), p(j) ))
         alpha = d / dot(z,p,neq_doti)
c .....................................................................
c
c ...
         do 210 i = 1, neq
c ... x(j+1) = x(j) + alpha*p
            x(i) = x(i) + alpha * p(i)
c ... r(j+1) = r(j) - alpha*Ap
            r(i) = r(i) - alpha * z(i)
  210    continue
c .....................................................................
c
c ...    
         di   = dot(r,r,neq_doti) 
c ... beta = ( r(j+1),r(j+1) ) / ( r(j),r(j) )
         beta = di / d
c .....................................................................
c
c ...
         do 220 i = 1, neq
c ... p(j+1) = r(j+1) + beta*p(j)
            p(i) = r(i) + beta * p(i)
  220    continue
c .....................................................................
c
c ...
         if(fhist) write(18,1500),j,dsqrt(d)/norm_b 
c .....................................................................
c
c ...
         d =  di
         if (dsqrt(d) .lt. conv) goto 300
c ......................................................................
c
c ...
         if( jj .eq.500) then
           jj = 0
           write(*,1300),j,dsqrt(d),conv 
         endif  
         jj = jj + 1
c ......................................................................
  230 continue
c ......................................................................
      write(*,1200) maxit
      if(flog) write(10,1200) maxit
      call stop_mef()
  300 continue
c
c ... produto:  x*Kx
c
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .           ,x,z 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
       xkx = dot(x,z,neq_doti)
c ......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r = b - Ax (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - z(i)
  310 continue
      tmp  = dot(r,r,neq_doti)
      tmp = dsqrt(tmp)
      if( tmp .gt. 3.16d0*conv ) then
        write(*,1400) tmp,conv
      endif
c ......................................................................
      time = get_time()
      time = time-time0
c ......................................................................
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,tmp,time
      endif
c ......................................................................
c     Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       "CG: "," it ",j, " x * Kx ",xkx," ||x|| ",norm
     .      ," tol ",tol," time ",time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA CG:',/,5x,'Coeficiente da diagonal ' 
     . '- equacao ',i9)
 1100 format(' (CG) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * || b ||        = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' CG:',5x,'It',i7,5x,2d20.10)
 1400 format (' CG:',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
 1500 format (' CG: ',5x,i7,5x,2es20.10)
      end
c *********************************************************************  
c
c *********************************************************************  
      subroutine pcg(neq   ,nequ   ,nad   ,ia       ,ja
     1              ,ad    ,au     ,al    ,m        ,b      
     2              ,x     ,z      ,r     ,p     
     3              ,tol   ,maxit
     4              ,matvec,dot
     5              ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6              ,i_xfi ,i_rcvsi,i_dspli
     7              ,fprint,flog   ,fhist  ,fnew   
     8              ,nprcs ,mpi)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 15/12/2016                                    * 
c * ------------------------------------------------------------------ *   
c * Subroutine PCG : Solucao de sistemas de equacoes pelo metodo dos   *
c * gradientes conjugados com precondicionador diagonal para matrizes  *
c * simetricas.                                                        *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq      - numero de equacoes                                      *
c * nequ     - numero de equacoes no bloco Kuu                         *
c * nad      - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * ad(neq)  - diagonal da matriz A                                    *
c * au(*)    - parte triangular superior de A                          *
c * al(*)    - parte triangular inferior de A                          *
c * m(*)     - precondicionador diagonal                               *
c * b(neqovlp)- vetor de forcas                                        *
c * x(neqovlp)- chute inicial                                          *
c * z(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * p(neqovlp)- arranjo local de trabalho                              *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * matvec   - nome da funcao externa para o produto matrix-vetor      *
c * dot      - nome da funcao externa para o produto escalar           *
c * my_id    - MPI                                                     *
c * neqf1i   - MPI                                                     *
c * neqf2i   - MPI                                                     *
c * neq_doti - MPI                                                     *
c * i_fmap   - MPI                                                     *
c * i_xfi    - MPI                                                     *
c * i_rvcs   - MPI                                                     *
c * i_dspli  - MPI                                                     *
c * fprint   - saida na tela                                           *
c * flog     - log do arquivo de saida                                 *
c * fhist    - log dos resuduos por iteracao                           *
c * fnew     - .true.  -> x0 igual a zero                              *
c *            .false. -> x0 dado                                      *
c * mpi      - true|false                                              *
c * nprcs    - numero de processos mpi                                 *  
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq) - modificado                                                *
c * ad(*),al(*),au(*) - inalterados                                    *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
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
c
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
      real*8  ad(*),au(*),al(*),b(*),m(*),x(*)
      real*8  r(*),z(*),p(*)
      real*8  dot,tol,conv,xkx,norm,d,di,alpha,beta,tmp,norm_b
      real*8  norm_r,norm_m_r
      real*8  time0,time
      real*8 dum1
      logical flog,fprint,fnew,fhist
c ...
      real*8 flop_cg,memory_pcg,mem
      real*8  mflops,vmean
c .....................................................................
      external matvec,dot
c ======================================================================
      time0 = get_time()
c ......................................................................
c
c ...
      mem = memory_pcg(neq,nad)
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
c ... conv = tol * |(M-1)b|m = tol *(b,M-1b)
      do 15 i = 1, neq
         z(i) = b(i) * m(i)
   15 continue
      d      = dot(b,z,neq_doti)
      norm_b = dsqrt(dabs(d))  
      conv   = tol*dsqrt(dabs(d))
c .......................................................................
c  
c ... Ax0                                                            
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .           ,x,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .           ,dum1)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r(i) = b(i) - z(i)
c ... z0 = (M-1)r0
         z(i) = r(i) * m(i)
c ... p0 = r0
         p(i) = z(i)
  100 continue
c ... ( r(0),z(0) ) = ( r(0), (M-1)r0 )
      d    = dot(r,z,neq_doti)
c ......................................................................
      jj = 1
      do 230 j = 1, maxit
c ... z = Ap(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .              ,p,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .              ,dum1)
c .....................................................................
c
c ... alpha = ( r(j),z(j) ) / ( Ap(j), p(j) ))
         alpha = d / dot(z,p,neq_doti)
c .....................................................................
c
c ...
         do 210 i = 1, neq
c ... x(j+1) = x(j) + alpha*p
            x(i) = x(i) + alpha * p(i)
c ... r(j+1) = r(j) - alpha*Ap
            r(i) = r(i) - alpha * z(i)
c ... z  = (M-1)r0
            z(i) = r(i) * m(i)
  210    continue
c .....................................................................
c
c ... ( r(j+1),(M-1)r(j+1) ) = ( r(j+1),z )
         di   = dot(r,z,neq_doti) 
c ... beta = ( r(j+1),(M-1)r(j+1) ) / ( r(j),r(j) ) 
         beta = di / d
c .....................................................................
c
c ...
         do 220 i = 1, neq
c ... p(j+1) = (M-1)r(j+1) + beta*p(j) = z + beta*p(j)
            p(i) = z(i) + beta * p(i)
  220    continue
c .....................................................................
c
c ...
         if(fhist) then
           if(my_id .eq.0) write(18,1500),j,dsqrt(dabs(d))/norm_b
         endif  
c .....................................................................
c
c ...
         d =  di
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ......................................................................
         if(jj.eq.2000) then
           jj = 0
           if(my_id.eq.0) write(*,1300),j,dsqrt(dabs(d)),conv 
         endif  
         jj = jj + 1
c ......................................................................
  230 continue
c ----------------------------------------------------------------------
      write(*,1200) maxit,dsqrt(dabs(d))
      if(flog) write(10,1200) maxit,dsqrt(dabs(d))
      call stop_mef()
  300 continue
c
c ... Energy norm:  x*Kx
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .           ,x,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      xkx = dot(x,z,neq_doti)
c ......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r =M(-1)(b - Ax) (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - z(i)
        z(i) = r(i)*m(i)
  310 continue
      norm_m_r = dot(r,z,neq_doti)
      norm_m_r = dsqrt(dabs(norm_m_r))
      norm_r = dot(r,r,neq_doti)
      norm_r = dsqrt(norm_r)
      if( norm_m_r .gt. conv ) then
         if(my_id .eq.0 )then
           write(*,1400) norm_m_r,conv
         endif 
      endif
c ......................................................................
      time = get_time()
      time = time-time0
c ......................................................................
c 
c ...  
      mflops = (flop_cg(neq,nad,j,2,mpi)*1.d-06)/time  
c ......................................................................
c
c ...
      if(my_id .eq.0 .and. fprint )then
        if(mpi) then
          write(*,1110)tol,conv,j,xkx,norm,norm_r,norm_m_r,mflops,time
        else
          write(*,1100)tol,conv,neq,nad,j,xkx,norm,norm_r,norm_m_r
     .                ,mflops,time
        endif
      endif
c ......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
         write(10,'(a,a,i9,3(a,d20.10),3(a,f20.2))')
     .       'PCG: ',' it ',j, ' x * Kx ',xkx,' ||x|| ',norm
     .      ,' tol ',tol,' Mflops ',mflops,' memory(MB) ',mem
     .      , ' time ',time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA PCG:',/,5x,'Diagonal coefficient ' 
     . '- equation ',i9,d20.6)
 1100 format(' (PCG) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||b||m         = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'|| b - Ax ||m        = ',d20.10/
     . 5x,'Mflops               = ',f20.2/
     . 5x,'CPU time (s)         = ',f20.2/)
 1110 format(' (PCG_MPI) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||b||m         = ',d20.6/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'|| b - Ax ||m        = ',d20.10/
     . 5x,'Mflops               = ',f20.2/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,d20.10
     .        ' iterations !',/)
 1300 format (' PCG:',5x,'It',i7,5x,2d20.10)
 1400 format (' PCG:',1x,'Explicit residual > tol * ||b||| :'
     .       ,1x,d20.10,1x,d20.10)
 1500 format ( 5x,i7,5x,2es20.10)
      end
c *********************************************************************  
c
c *********************************************************************
      subroutine iccg(neq   ,nequ   ,nad     ,ia      ,ja
     .               ,ad    ,au     ,al      ,m       ,b       
     .               ,x     ,z      ,r       ,p   
     .               ,tol   ,maxit
     .               ,matvec,dot    ,triasolv
     .               ,my_id ,neqf1i ,neqf2i  ,neq_doti,i_fmapi
     .               ,i_xfi ,i_rcvsi,i_dspli
     .               ,fprint,flog   ,fhist   ,fnew)
c **********************************************************************
c * Data de criacao    : 10/04/2016                                    *
c * Data de modificaco : 20/04/2016                                    * 
c * ------------------------------------------------------------------ *   
c * IcCG : Solucao de sistemas de equacoes pelo metodo dos gradientes  *
c * conjugados com precondicionador incompleto                         *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq      - numero de equacoes                                      *
c * nequ     - numero de equacoes no bloco Kuu                         *
c * nad      - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * ad(neq)  - diagonal da matriz A                                    *
c * au(*)    - parte triangular superior de A                          *
c * al(*)    - parte triangular inferior de A                          *
c * m(*)     - precondicionador diagonal                               *
c * b(neq)   - vetor de forcas                                         *
c * x(neq)   - chute inicial                                           *
c * z(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * p(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * matvec   - nome da funcao externa para o produto matrix-vetor      *
c * dot      - nome da funcao externa para o produto escalar           *
c * triasolv - nome da funcao externa para o o solver triangular       *
c * my_id    -                                                         *
c * neqf1i   -                                                         *
c * neqf2i   -                                                         *
c * neq_doti -                                                         *
c * i_fmap   -                                                         *
c * i_xfi    -                                                         *
c * i_rvcs   -                                                         *
c * i_dspli  -                                                         *
c * fprint   - saida na tela                                           *
c * flog     - log do arquivo de saida                                 *
c * fhist    - log dos resuduos por iteracao                           *
c * fnew     - .true.  -> x0 igual a zero                              *
c *            .false. -> x0 dado                                      *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq) - modificado                                                *
c * ad(*),al(*),au(*) - inalterados                                    *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * Arranjos jat,iat e kat s�o utilizados na retrosubstituizao do      *
c * solver iLDLt                                                       * 
c **********************************************************************
      implicit none
c ...
      real*8 get_time
c ......................................................................
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................      
      integer neq,nequ,nad,maxit,i,j,jj
      integer ia(*),ja(*),my_id
      real*8  ad(*),au(*),al(*),x(*),b(*),r(*),z(*),p(*)
      real*8  m(*)
      real*8  dot,ddot,tol,conv,xkx,norm,d,di,alpha,beta,tmp
      real*8  norm_r,norm_m_r,norm_b
      real*8  time0,time
      real*8 dum1
      logical flog,fnew,fprint,fhist
      external matvec,dot,triasolv
c ======================================================================
      time0 = get_time()
c ......................................................................
c
c ...
      do 5 i = 1, neq
        if(ad(i) .eq. 0.d0 ) then
          write(*,1000) i
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
c ......................................................................
c
c ... conv = tol * |(M-1)b|m = tol *(b,M-1b)
      call triasolv(neq,ia,ja,m,m(neq+1),b,z)
      d      = dot(b,z,neq_doti)
      norm_b = dsqrt(dabs(d))
      conv   = tol*dsqrt(dabs(d))
c .......................................................................
c  
c ... Ax0
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)  
     .           ,x,z 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................      
c
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r(i) = b(i) - z(i)
  100 continue
c .......................................................................      
c
c ... Mz=r  
      call triasolv(neq,ia,ja,m,m(neq+1),r,z)
c .......................................................................      
c
c ...
      do 105 i = 1, neq
c ... p0 = r0
         p(i) = z(i)
  105 continue
c .......................................................................      
c
c ... ( z(0),z(0) )m = ( r(0), z0 ) (r(0), (M-1)r0 )
      d    = dot(r,z,neq_doti) 
c .......................................................................
      jj = 1
      do 230 j = 1, maxit
c ... z = Ap(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .              ,p,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli 
     .              ,dum1)
c .....................................................................
c
c ... alpha = ( r(j),z(j) ) / ( Ap(j), p(j) ))
         alpha = d / dot(z,p,neq_doti)
c .....................................................................
c
c ...
         do 210 i = 1, neq
c ... x(j+1) = x(j) + alpha*p
            x(i) = x(i) + alpha * p(i)
c ... r(j+1) = r(j) - alpha*Ap
            r(i) = r(i) - alpha * z(i)
  210    continue
c ......................................................................
c
c ... Mz=r  
         call triasolv(neq,ia,ja,m,m(neq+1),r,z)
c .......................................................................

c ... ( r(j+1),(M-1)r(j+1) ) = ( r(j+1),z )
         di   = dot(r,z,neq_doti) 
c ... beta = ( r(j+1),(M-1)r(j+1) ) / ( r(j),r(j) ) 
         beta = di / d
c .....................................................................
c
c ...         
         do 220 i = 1, neq
            p(i) = z(i) + beta * p(i)
  220    continue
c .....................................................................
c
c ...
         if(fhist) write(18,1500),j,dsqrt(dabs(d))/norm_b 
c .....................................................................
c
c ...
         d = di           
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ......................................................................
         if( jj .eq.500) then
           jj = 0
           write(*,1300),j,dsqrt(dabs(d)),conv 
         endif  
         jj = jj + 1
c ......................................................................
  230 continue
c ----------------------------------------------------------------------
      write(*,1200) maxit
      if(flog) write(10,1200) maxit
      call stop_mef()
  300 continue
c
c ... Energy norm: x*Kx
c
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .          ,x,z 
     .          ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      xkx    = dot(x,z,neq_doti)
c ......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r =(b - Ax) (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - z(i)
  310 continue
c ... Mz=r  
      call triasolv(neq,ia,ja,m,m(neq+1),r,z)
c .......................................................................
      norm_m_r = dot(r,z,neq_doti)
      norm_m_r = dsqrt(norm_m_r)
      norm_r = dot(r,r,neq_doti)
      norm_r = dsqrt(norm_r)
      if( norm_m_r .gt. 3.16d0*conv ) then
         if(my_id .eq.0 )then
           write(*,1400) norm_m_r,conv
         endif 
      endif
c ......................................................................
      time = get_time()
      time = time-time0
c ----------------------------------------------------------------------
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,norm_r,norm_m_r,time
      endif
c ......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       "ICCG: "," it ",j, " x * Kx ",xkx," ||x|| ",norm
     .      ," tol ",tol," time ",time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA ICCG:',/,5x,'Coeficiente da diagonal nulo
     . - equacao ',i9)
 1100 format(' (ICCG) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||b||m         = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'|| b - Ax ||m        = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' ICCG: ',5x,'It',i7,5x,2d20.10)
 1400 format (' ICCG:',1x,'Residuo exato > 3.16*conv ',1x,d20.10
     .       ,1x,d20.10)
 1500 format ( 'ICCG ',5x,i7,5x,2es20.10)
      end
c **********************************************************************
c
c *********************************************************************
      subroutine bpcg(neq   ,nequ   ,nad     ,ia      ,ja
     .               ,ad    ,au     ,al      ,m       ,b       
     .               ,x     ,z      ,r       ,p   
     .               ,tol   ,maxit  
     .               ,matvec,dot    
     .               ,my_id ,neqf1i ,neqf2i  ,neq_doti,i_fmapi
     .               ,i_xfi ,i_rcvsi,i_dspli
     .               ,fprint,flog   ,fhist   ,fnew)
c **********************************************************************
c * Data de criacao    : 17/06/2016                                    *
c * Data de modificaco : 18/06/2016                                    * 
c * ------------------------------------------------------------------ *   
c * BPCG : Solucao de sistemas de equacoes pelo metodo dos gradientes  *
c * conjugados com precondicionador bloco diagonal                     *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq      - numero de equacoes                                      *
c * nequ     - numero de equacoes no bloco Kuu                         *
c * nad      - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * ad(neq)  - diagonal da matriz A                                    *
c * au(*)    - parte triangular superior de A                          *
c * al(*)    - parte triangular inferior de A                          *
c * m(*)     - precondicionador bloco diagonal                         *
c * b(neq)   - vetor de forcas                                         *
c * x(neq)   - chute inicial                                           *
c * z(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * p(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * matvec   - nome da funcao externa para o produto matrix-vetor      *
c * dot      - nome da funcao externa para o produto escalar           *
c * my_id    -                                                         *
c * neqf1i   -                                                         *
c * neqf2i   -                                                         *
c * neq_doti -                                                         *
c * i_fmap   -                                                         *
c * i_xfi    -                                                         *
c * i_rvcs   -                                                         *
c * i_dspli  -                                                         *
c * fprint   - saida na tela                                           *
c * flog     - log do arquivo de saida                                 *
c * fhist    - log dos resuduos por iteracao                           *
c * fnew     - .true.  -> x0 igual a zero                              *
c *            .false. -> x0 dado                                      *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq) - modificado                                                *
c * ad(*),al(*),au(*) - inalterados                                    *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************
      implicit none
      include 'precond.fi'
      integer neqf1i,neqf2i,neq_doti
c ...
      real*8 get_time
c ......................................................................
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................      
      integer neq,nequ,nad,maxit,i,j,jj
      integer ia(*),ja(*),my_id,bsize
      real*8  ad(*),au(*),al(*),x(*),b(*),r(*),z(*),p(*)
      real*8  m(*),max_block_a(max_block*max_block)
      real*8  dot,ddot,tol,conv,xkx,norm,d,di,alpha,beta,tmp
      real*8  norm_r,norm_m_r,norm_b
      real*8  time0,time
      real*8 dum1
      logical flog,fnew,fprint,fhist
      external matvec,dot
c ======================================================================
      time0 = get_time()
c ......................................................................
c
c ...
      do 5 i = 1, neq
        if(ad(i) .eq. 0.d0 ) then
          write(*,1000) i
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
c ......................................................................
c
c ... conv = tol * |(M-1)b|m = tol *(b,M-1b)
      call op_block_precond(m,b,z,max_block_a,iparam)
      d      = dot(b,z,neq_doti)
      norm_b = dsqrt(dabs(d))
      conv   = tol*dsqrt(dabs(d))
c .......................................................................
c  
c ... Ax0
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)  
     .           ,x,z 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................      
c
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r(i) = b(i) - z(i)
  100 continue
c .......................................................................      
c
c ... z=M(-1)r  
      call op_block_precond(m,r,z,max_block_a,iparam)   
c .......................................................................      
c
c ...
      do 105 i = 1, neq
c ... p0 = r0
         p(i) = z(i)
  105 continue
c .......................................................................      
c
c ... ( z(0),z(0) )m = ( r(0), z0 ) (r(0), (M-1)r0 )
      d    = dot(r,z,neq_doti) 
c .......................................................................
      jj = 1
      do 230 j = 1, maxit
c ... z = Ap(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .              ,p,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli 
     .              ,dum1)
c .....................................................................
c
c ... alpha = ( r(j),z(j) ) / ( Ap(j), p(j) ))
         alpha = d / dot(z,p,neq_doti)
c .....................................................................
c
c ...
         do 210 i = 1, neq
c ... x(j+1) = x(j) + alpha*p
            x(i) = x(i) + alpha * p(i)
c ... r(j+1) = r(j) - alpha*Ap
            r(i) = r(i) - alpha * z(i)
  210    continue
c ......................................................................
c
c ... z=M(-1)r 
         call op_block_precond(m,r,z,max_block_a,iparam)    
c .......................................................................

c ... ( r(j+1),(M-1)r(j+1) ) = ( r(j+1),z )
         di   = dot(r,z,neq_doti) 
c ... beta = ( r(j+1),(M-1)r(j+1) ) / ( r(j),r(j) ) 
         beta = di / d
c .....................................................................
c
c ...         
         do 220 i = 1, neq
            p(i) = z(i) + beta * p(i)
  220    continue
c .....................................................................
c
c ...
         if(fhist) write(18,1500),j,dsqrt(dabs(d))/norm_b 
c .....................................................................
c
c ...
         d = di
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ......................................................................
         if( jj .eq.500)then
           jj = 0
           write(*,1300),j,dsqrt(dabs(d)),conv 
         endif  
         jj = jj + 1
c ......................................................................
  230 continue
c ----------------------------------------------------------------------
      write(*,1200) maxit
      if(flog) write(10,1200) maxit
      call stop_mef()
  300 continue
c
c ... Energy norm: x*Kx
c
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .          ,x,z 
     .          ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      xkx    = dot(x,z,neq_doti)
c ......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r =M(-1)(b - Ax) (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - z(i)
  310 continue
c ... z=M(-1)r
      call op_block_precond(m,r,z,max_block_a,iparam) 
c
      norm_m_r = dot(r,z,neq_doti)
      norm_m_r = dsqrt(dabs(norm_m_r))
      norm_r = dot(r,r,neq_doti)
      norm_r = dsqrt(norm_r)
      if( norm_m_r .gt. 3.16d0*conv ) then
         if(my_id .eq.0 )then
           write(*,1400) norm_m_r,conv
         endif 
      endif
c ......................................................................
      time = get_time()
      time = time-time0
c ----------------------------------------------------------------------
      bsize = iparam(3)*iparam(1)+ iparam(2) 
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,iparam(4),bsize
     .              ,xkx,norm,norm_r,norm_m_r,time
      endif
c ......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,i3,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       'BPCG: ',' it ',j,' block ',iparam(4), ' x * Kx '
     .      ,xkx,' ||x|| ',norm,' tol ',tol,' time ',time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA BPCG:',/,5x,'Coeficiente da diagonal nulo
     . - equacao ',i9)
 1100 format(' (BPCG) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||b||m         = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Size of blocks       = ',i20/
     . 5x,'Size of M            = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'|| b - Ax ||m        = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' BPCG: ',5x,'It',i7,5x,2d20.10)
 1400 format (' BPCG:',1x,'Residuo exato > 3.16*conv ',1x,d20.10
     .       ,1x,d20.10)
 1500 format ( 'BPCG: ',5x,i7,5x,2es20.10)
      end
c **********************************************************************
c

c **********************************************************************
      subroutine bicgstab(neq      ,nequ  ,nad,ia ,ja 
     .                    ,ad      ,au    ,al,b  ,x   
     .                    ,t       ,v     ,r ,p  ,r0
     .                    ,tol     ,maxit  
     .                    ,matvec  ,dot    
     .                    ,my_id   ,neqf1i,neqf2i 
     .                    ,neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .                    ,fprint  ,flog   ,fnew)
c **********************************************************************
c * Data de criacao    : 15/04/2016                                    *
c * Data de modificaco : 15/12/2016                                    * 
c * ------------------------------------------------------------------ *   
c * BICGSTAB  : Solucao de sistemas de equacoes pelo metodo dos        * 
c * gradientes biconjugados para matrizes nao-simetricas.              *                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *                                                                   *
c * neq      - numero de equacoes                                      *
c * nequ     - numero de equacoes no bloco Kuu                         *
c * nad      - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * ad(neq)  - diagonal da matriz A                                    *
c * au(*)    - parte triangular superior de A                          *
c * al(*)    - parte triangular inferior de A                          *
c * b(neq)   - vetor de forcas                                         *
c * x(neq)   - chute inicial                                           *
c * t(neq)   - arranjo local de trabalho                               *
c * v(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * p(neq)   - arranjo local de trabalho                               *
c * z(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * matvec   - nome da funcao externa para o produto matrix-vetor      *
c * dot      - nome da funcao externa para o produto escalar           *
c * my_id    -                                                         *
c * neqf1i   -                                                         *
c * neqf2i   -                                                         *
c * neq_doti -                                                         *
c * i_fmap   -                                                         *
c * i_xfi    -                                                         *
c * i_rvcs   -                                                         *
c * i_dspli  -                                                         *
c * fprint   - saida na tela                                           *
c * flog     - log do arquivo de saida                                 *
c * fnew     - .true.  x0 igual a zero                                 *
c *            .false. x0 dado                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq) - modificado                                                *
c * ad(*),al(*),au(*) e b - inalterados                                *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *                                                                      *
c **********************************************************************
      implicit none
c ...
      real*8 get_time
c ......................................................................
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................   
      integer neq,nequ
      integer maxit,i,j,jj,k
      integer*8 ia(*),nad
      integer ja(*),my_id
      real*8  ad(*),au(*),al(*),x(*),b(*)
      real*8  r(*),p(*),t(*),v(*),r0(*)
      real*8  dot,tol,conv,xkx,norm,d,alpha,beta,rr0,w,tmp,norm_r
      real*8  time0,time
      real*8  dum1 
      logical flog,fprint,fnew
      external matvec,dot
c ======================================================================
      time0 = get_time()
c ......................................................................
c     if(my_id.eq.0) print *, 'nad :',nad
c ......................................................................
c
c ...
      do 5 i = 1, neq
        if(ad(i) .eq. 0.d0 ) then
          write(*,1000) i
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
      d    = dot(b,b,neq_doti)
      conv = tol*dsqrt(dabs(d))
c .......................................................................
c
c ...
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1),
     .            x,p,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r0(i) = b(i) - p(i)
c ... r = r0
         r(i)  = r0(i)
c ... p = r0
         p(i)  = r0(i)
  100 continue
c .......................................................................
c
c ...
      jj = 1
      do 230 j = 1, maxit
c ... v = Ap(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .              ,p,v
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ... alpha = ( r(j),r0 ) / ( Ap(j), r0 ))
         rr0   = dot(r,r0,neq_doti)
         alpha = rr0/dot(v,r0,neq_doti)
c .......................................................................
c
c ...
         do 210 i = 1, neq
c ... x(j+1) = x(j) + alpha*p
            x(i) = x(i) + alpha * p(i)
c ... s(j)   = r(j) - alpha*Ap(j)
            r(i) = r(i) - alpha * v(i)
  210    continue
c ........................................................................
c
c ... (s,s)
         d = dot(r,r,neq_doti)
c ...
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ........................................................................
c
c ... t = Ar(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .              ,r,t 
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c ........................................................................
c
c ... w = ( Ar(j),r(j) ) / ( Ar(j), Ar(j) ))
         w = dot(t,r,neq_doti) / dot(t,t,neq_doti)
c ........................................................................
c
c ... 
         do 220 i = 1, neq
c ... x(j+1) = x(j) + w*r(j)
            x(i) = x(i) + w*r(i)
c ... r(j+1) = s(j) - w*As(j)
            r(i) = r(i) - w*t(i)
  220    continue
c .......................................................................
c
c ... (r,r)
         d = dot(r,r,neq_doti)
c ...
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c .......................................................................
c
c ... beta = ( r(j+1),r0 ) / ( r(j), r0 )) * (alpha/w) 
         beta = (dot(r,b,neq_doti) / rr0)*(alpha/w)
c .......................................................................
c
c ...
         do 225 i = 1, neq
             p(i) = r(i) + beta*(p(i)-w*v(i))
  225    continue
c .......................................................................
c
c ...
         if( jj .eq. 500) then
           jj = 0
           if(my_id.eq.0) write(*,1300),j,dsqrt(dabs(d)),conv 
         endif  
         jj = jj + 1
c ......................................................................
  230 continue
c ......................................................................
c
c ...
      if(my_id.eq.0) then
        write(*,1200) maxit
        if(fLog) write(10,1200) maxit
      endif
      call stop_mef()
  300 continue
c ......................................................................
c
c ... produto:  x*Kx
c
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .           ,x,p 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      xkx =  dot(x,p,neq_doti)
c .......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r = b - Ax (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - p(i)
  310 continue
      tmp  = dot(r,r,neq_doti)
      norm_r = dsqrt(tmp)
      if( tmp .gt. 3.16d0*conv ) then
         if(my_id .eq.0 )then
           write(*,1400) tmp,conv
         endif 
      endif
c ......................................................................
c
c ...
      time = get_time()
      time = time-time0
c .......................................................................
c
c ...
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,norm_r,time
      endif
c .......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       "BICGSTAB: "," it ",j, " x * Kx ",xkx," ||x|| ",norm
     .      ," tol ",tol," time ",time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (/,1x,'SUBROTINA BICCSTAB:',/,1x
     .       ,'Coeficiente da diagonal nulo ',i9)
 1100 format(' (BICGSTAB) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * || b ||        = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' BICGSTAB:',5x,'It',i7,5x,2d20.10)
 1400 format (' BICCSTAB:',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
      end
c *********************************************************************
c
c **********************************************************************
      subroutine bicgstabl2(neq     ,nequ  ,nad,ia ,ja 
     .                     ,ad      ,au    ,al ,b  ,x   
     .                     ,t       ,v     ,r  ,u  ,r0
     .                     ,w       ,s 
     .                     ,tol     ,maxit  
     .                     ,matvec  ,dot    
     .                     ,my_id   ,neqf1i,neqf2i 
     .                     ,neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .                     ,fprint  ,flog   ,fnew)
c **********************************************************************
c * Data de criacao    : 01/05/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * BICGSTABL2: Solucao de sistemas de equacoes pelo metodo dos        * 
c * gradientes biconjugados(2) para matrizes nao-simetricas.           *                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *                                                                   *
c * neq      - numero de equacoes                                      *
c * nequ     - numero de equacoes no bloco Kuu                         *
c * nad      - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * ad(neq)  - diagonal da matriz A                                    *
c * au(*)    - parte triangular superior de A                          *
c * al(*)    - parte triangular inferior de A                          *
c * b(neq)   - vetor de forcas                                         *
c * x(neq)   - chute inicial                                           *
c * t(neq)   - arranjo local de trabalho                               *
c * v(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * u(neq)   - arranjo local de trabalho                               *
c * r0(neq)  - arranjo local de trabalho                               *
c * w(neq)   - arranjo local de trabalho                               *
c * s(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * matvec   - nome da funcao externa para o produto matrix-vetor      *
c * dot      - nome da funcao externa para o produto escalar           *
c * my_id    -                                                         *
c * neqf1i   -                                                         *
c * neqf2i   -                                                         *
c * neq_doti -                                                         *
c * i_fmap   -                                                         *
c * i_xfi    -                                                         *
c * i_rvcs   -                                                         *
c * i_dspli  -                                                         *
c * fprint   - saida na tela                                           *
c * flog     - log do arquivo de saida                                 *
c * fnew     - .true.  x0 igual a zero                                 *
c *            .false. x0 dado                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * ad(*),al(*),au(*) e b - inalterados                                *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * Vers�o do livro Iterative krylov Method for large linear Systems   *                                                                      *
c **********************************************************************
      implicit none
      integer neqf1i,neqf2i,neq_doti
c ...
      real*8 get_time
c ......................................................................
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................   
      integer neq,nequ
      integer*8 ia(*),nad
      integer maxit,i,j,jj,k
      integer ja(*),my_id
      real*8  ad(*),au(*),al(*),x(*),b(*)
      real*8  r(*),u(*),t(*),v(*),r0(*),w(*),s(*)
      real*8  dot,tol,conv,xkx,norm,d
      real*8  alpha,beta,rr0,rr1,omega1,omega2,mi,ni,gamma,tau,tmp
      real*8  time0,time
      real*8  dum1 
      logical flog,fprint,fnew
      external matvec,dot
c ======================================================================
      time0 = get_time()
c ......................................................................
c
c ...
      do 5 i = 1, neq
        if(ad(i) .eq. 0.d0 ) then
          write(*,1000) i
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
      d    = dot(b,b,neq_doti)
      conv = tol*dsqrt(dabs(d))
c .......................................................................
c
c ... Ax
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1),
     .            x,t,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r0(i) = b(i) - t(i)
c ... r = r0
         r(i)  = r0(i)
c ... u = 0.0d0
         u(i)  = 0.d0 
  100 continue
c .......................................................................
c
c ...
      rr0     = 1.0d0
      alpha   = 0.d0
      omega2  = 1.d0
c .......................................................................
c
c ...
      jj = 1
      do 230 j = 1, maxit, 1
c ... rro = -w2*rr0
         rr0 = -omega2*rr0 
c ... even BiCG step:
c ... rr1 = (r,r0)
         rr1  = dot(r,r0,neq_doti)
c ... beta = alpha * rr1/ rr0
         beta = alpha*rr1/rr0
c ... rr0 = rr1
         rr0  = rr1
c .......................................................................
c
c ...
         do i = 1, neq
c ... u(j)   = r(j) - beta*u(i)
           u(i) = r(i) - beta * u(i)
         enddo        
c .......................................................................
c
c ... v = Au(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .              ,u,v
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...  gamma = (v,r0)
         gamma = dot(v,r0,neq_doti)
c ... alpha = rr0 / gamma
         alpha = rr0/gamma
c .......................................................................
c
c ...
         do i = 1, neq
c ... r(j)   = r(j) - alpha*v(i)
           r(i) = r(i) - alpha * v(i)
         enddo        
c .......................................................................
c
c ... s = Ar(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .              ,r,s
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...
         do i = 1, neq
c ... x(j) = x(j) + alpha*u
            x(i) = x(i) + alpha * u(i)
         enddo   
c ........................................................................
c
c ... odd BiCG step:
c ... rr1 = (s,r0)
         rr1  = dot(s,r0,neq_doti)
c ... beta = alpha*rr1/rr0
         beta = alpha*rr1/rr0
c ... rr0 = rr1
         rr0  = rr1
c .......................................................................
c
c ...
         do i = 1, neq
c ... v(j) = s(j) - beta*v
            v(i) = s(i) - beta * v(i)
         enddo   
c ........................................................................
c
c ... w = Av(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .              ,v,w
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ... gamma = (w,r0) 
         gamma = dot(w,r0,neq_doti)
c ... alpha = rr0/gamma
         alpha = rr0/gamma
c .......................................................................
c
c ...
         do i = 1, neq
c ... u(j)   = r(j) - beta*u(i)
           u(i) = r(i) - beta * u(i)
c ... r(j)   = r(j) - alpha*v(i)
           r(i) = r(i) - alpha * v(i)
c ... s(j)   = s(j) - alpha*w(i)
           s(i) = s(i) - alpha * w(i)
         enddo        
c .......................................................................
c
c ... t = As(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .              ,s,t
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ... GCR(2)-part
c ... w1 = (r,s)
         omega1 = dot(r,s,neq_doti)
c ... mi = (s,s)
         mi = dot(s,s,neq_doti)
c ... ni = (s,t)
         ni = dot(s,t,neq_doti)
c ... tau = (t,t)
         tau = dot(t,t,neq_doti)
c ... w2
         omega2 = dot(r,t,neq_doti)
c ... tau = tau - ni*ni/mi
         tau = tau - ni*ni/mi
c ... w2 = (w2 - ni*w1/mi)/tau
         omega2 = (omega2 - ni*omega1/mi)/tau
c ... w1 = (w1 - ni*w2)/mi
         omega1 =(omega1 - ni*omega2)/mi
c .......................................................................
c
c ...
         do i = 1, neq
c ... x(j+2) = x(i) + w1 * r(i) + w2 * s(i) + alpha * u(i)
           x(i) = x(i) + omega1 * r(i) + omega2 * s(i) + alpha * u(i)
c ... r(j+2) = r(j) - w1 * s(i) - w2 * t(i)
           r(i) = r(i) - omega1 * s(i) - omega2 * t(i)
         enddo        
c .......................................................................
c
c ...
         d  = dot(r,r,neq_doti)
c ...
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c .......................................................................
c
c ...
         do i = 1, neq
c ... u(j+1) = u(j) - w1 * v(i) - w2 * w(i)
           u(i) = u(i) - omega1 * v(i) - omega2 * w(i)
         enddo    
c .......................................................................
c
c ...
         if( jj .eq.500) then
           jj = 0
           if(my_id.eq.0) write(*,1300),j,dsqrt(dabs(d)),conv 
         endif  
         jj = jj + 1
c ......................................................................
  230 continue
c ......................................................................
c
c ...
      if(my_id.eq.0) then
        write(*,1200) maxit
        if(fLog) write(10,1200) maxit
      endif
      call stop_mef()
  300 continue
c ......................................................................
c
c ... produto:  x*Kx
c
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .           ,x,t 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      xkx =  dot(x,t,neq_doti)
c .......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r = b - Ax (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - t(i)
  310 continue
      tmp  = dot(r,r,neq_doti)
      tmp = dsqrt(tmp)
      if( tmp .gt. 3.16d0*conv ) then
         if(my_id .eq.0 )then
           write(*,1400) tmp,conv
         endif 
      endif
c ......................................................................
c
c ...
      time = get_time()
      time = time-time0
c .......................................................................
c
c ...
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,time
      endif
c .......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       "BICGSTAB(2): "," it ",j, " x * Kx ",xkx," ||x|| ",norm
     .      ," tol ",tol," time ",time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (/,1x,'SUBROTINA BICCSTAB(2):',/,1x
     .       ,'Coeficiente da diagonal nulo ',i9)
 1100 format(' (BICGSTAB(2)) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * || b ||        = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' BICGSTAB(2):',5x,'It',i7,5x,2d20.10)
 1400 format (' BICCSTAB(2):',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
      end
c **********************************************************************
c
c **********************************************************************
      subroutine pbicgstabl2(neq     ,nequ  ,nad,ia ,ja 
     1                     ,ad      ,au    ,al  ,m  ,b  ,x   
     2                     ,t       ,v     ,r  ,u   ,r0
     3                     ,w       ,s     ,h  ,p   ,z
     4                     ,tol     ,maxit  
     5                     ,matvec  ,dot    
     6                     ,my_id   ,neqf1i,neqf2i 
     7                     ,neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli
     8                     ,fprint  ,flog   ,fnew
     9                     ,nprcs ,mpi)
c **********************************************************************
c * Data de criacao    : 01/05/2016                                    *
c * Data de modificaco : 11/11/2016                                    * 
c * ------------------------------------------------------------------ *   
c * PBICGSTABl2  : Solucao de sistemas de equacoes pelo metodo dos     * 
c * gradientes biconjugados para matrizes nao-simetricas com           *
c * prendicionador diagonal                                            *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *                                                                   *
c * neq      - numero de equacoes                                      *
c * nequ     - numero de equacoes no bloco Kuu                         *
c * nad      - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * ad(neq)  - diagonal da matriz A                                    *
c * au(*)    - parte triangular superior de A                          *
c * al(*)    - parte triangular inferior de A                          *
c * b(neq)   - vetor de forcas                                         *
c * x(neq)   - chute inicial                                           *
c * t(neq)   - arranjo local de trabalho                               *
c * v(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * u(neq)   - arranjo local de trabalho                               *
c * r0(neq)  - arranjo local de trabalho                               *
c * w (neq)  - arranjo local de trabalho                               *
c * s (neq)  - arranjo local de trabalho                               *
c * p (neq)  - arranjo local de trabalho                               *
c * h (neq)  - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * matvec   - nome da funcao externa para o produto matrix-vetor      *
c * dot      - nome da funcao externa para o produto escalar           *
c * my_id    -                                                         *
c * neqf1i   -                                                         *
c * neqf2i   -                                                         *
c * neq_doti -                                                         *
c * i_fmap   -                                                         *
c * i_xfi    -                                                         *
c * i_rvcs   -                                                         *
c * i_dspli  -                                                         *
c * fprint   - saida na tela                                           *
c * flog     - log do arquivo de saida                                 *
c * fnew     - .true.  x0 igual a zero                                 *
c *            .false. x0 dado                                         *
c * mpi      - true|false                                              *
c * nprcs    - numero de processos mpi                                 * 
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * ad(*),al(*),au(*) e b - inalterados                                *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * Vers�o do livro Iterative krylov Method for large linear Systems   * 
c * A(M-1)y=b precondicionador a direita                               *                                                                      *
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
      integer neq,nequ
      integer maxit,i,j,jj,k
      integer*8 ia(*),nad
      integer ja(*),my_id
      real*8  ad(*),au(*),al(*),m(*),x(*),b(*)
      real*8  r(*),u(*),t(*),v(*),r0(*),w(*),s(*),p(*),z(*),h(*)
      real*8  dot,tol,conv,xkx,norm,d,morm_r
      real*8  alpha,beta,rr0,rr1,omega1,omega2,mi,ni,gamma,tau,tmp
      real*8  time0,time
      real*8  dum1 
      real*8  breaktol,btol
      parameter (btol = 1.d-32)
      logical flog,fprint,fnew
c ...
      real*8 flop_bicgstab2
      real*8 mflops,vmean
c .....................................................................
      external matvec,dot
c ======================================================================
      time0 = get_time()
c ......................................................................
c
c ...
      do 5 i = 1, neq
        if(ad(i) .eq. 0.d0 ) then
          write(*,1000) i
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
      d    = dot(b,b,neq_doti)
      conv = tol*dsqrt(dabs(d))
c .......................................................................
c
c ... Ax
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1),
     .            x,t,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r0(i) = b(i) - t(i)
c ... r = r0
         r(i)  = r0(i)
c ... u = 0.0d0
         u(i)  = 0.d0 
  100 continue
c .......................................................................
c
c ...
      rr0     = 1.0d0
      alpha   = 0.d0
      omega2  = 1.d0
c .......................................................................
c
c ...
      jj = 1
      do 230 j = 1, maxit, 1
c ... rro = -w2*rr0
         rr0 = -omega2*rr0 
         if( dabs(rr0) .lt. breaktol) then
           write(*,1510)dabs(rr0)
           call stop_mef() 
         endif  
c ... even BiCG step:
c ... rr1 = (r,r0)
         rr1  = dot(r,r0,neq_doti)
c ... beta = alpha * rr1/ rr0
         beta = alpha*rr1/rr0
c ... rr0 = rr1
         rr0  = rr1
c .......................................................................
c
c ...
         do i = 1, neq
c ... u(j) = r(j) - beta*u(i)
           u(i) = r(i) - beta * u(i)
c ... p    = M(-1)u
           z(i)= u(i)*m(i) 
         enddo        
c .......................................................................
c
c ... v = Ap(j) = AM(-1)u
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .              ,z,v
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...  gamma = (v,r0)
         gamma = dot(v,r0,neq_doti)
c ... alpha = rr0 / gamma
         alpha = rr0/gamma
c .......................................................................
c
c ...
         do i = 1, neq
c ... r(j)   = r(j) - alpha*v(i)
           r(i) = r(i) - alpha * v(i)
c ... p    = M(-1)r
           p(i) = r(i)*m(i) 
         enddo        
c .......................................................................
c
c ... s = Ap(j) = AM(-1)r
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .              ,p,s
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...
         do i = 1, neq
c ... x(j) = x(j) + alpha*(M-1)u
            x(i) = x(i) + alpha * z(i)
         enddo   
c ........................................................................
c
c ... odd BiCG step:
c ... rr1 = (s,r0)
         rr1  = dot(s,r0,neq_doti)
c ... beta = alpha*rr1/rr0
         beta = alpha*rr1/rr0
c ... rr0 = rr1
         rr0  = rr1
c .......................................................................
c
c ...
         do i = 1, neq
c ... v(j) = s(j) - beta*v
           v(i) = s(i) - beta * v(i)
c ... p    = M(-1)v
           p(i) = v(i)*m(i) 
         enddo   
c ........................................................................
c
c ... w = Av(j) = AM(-1)v
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .              ,p,w
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ... gamma = (w,r0) 
         gamma = dot(w,r0,neq_doti)
c ... alpha = rr0/gamma
         alpha = rr0/gamma
c .......................................................................
c
c ...
         do i = 1, neq
c ... u(j)   = r(j) - beta*u(i)
           u(i) = r(i) - beta * u(i)
c ... r(j)   = r(j) - alpha*v(i)
           r(i) = r(i) - alpha * v(i)
c ... s(j)   = s(j) - alpha*w(i)
           s(i) = s(i) - alpha * w(i)
c ... p    = M(-1)s
           p(i) = s(i)*m(i)
c ... z    = M(-1)u
           z(i) = u(i)*m(i) 
c ... h    = M(-1)r
           h(i) = r(i)*m(i)  
         enddo        
c .......................................................................
c
c ... t = Ap(j) = AM(-1)s
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1)
     .              ,p,t
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ... GCR(2)-part
c ... w1 = (r,s)
         omega1 = dot(r,s,neq_doti)
c ... mi = (s,s)
         mi = dot(s,s,neq_doti)
c ... ni = (s,t)
         ni = dot(s,t,neq_doti)
c ... tau = (t,t)
         tau = dot(t,t,neq_doti)
c ... w2
         omega2 = dot(r,t,neq_doti)
c ... tau = tau - ni*ni/mi
         tau = tau - ni*ni/mi
c ... w2 = (w2 - ni*w1/mi)/tau
         omega2 = (omega2 - ni*omega1/mi)/tau
c ... w1 = (w1 - ni*w2)/mi
         omega1 =(omega1 - ni*omega2)/mi
c .......................................................................
c
c ...
         do i = 1, neq
c ... x(j+2) = x(i) + w1 * r(i) + w2 * s(i) + alpha * u(i)
           x(i) = x(i) + omega1 * h(i) + omega2 * p(i) + alpha * z(i)
c ... r(j+2) = r(j) - w1 * s(i) - w2 * t(i)
           r(i) = r(i) - omega1 * s(i) - omega2 * t(i)
         enddo        
c .......................................................................
c
c ...
         d  = dot(r,r,neq_doti)
c ...
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c .......................................................................
c
c ...
         do i = 1, neq
c ... u(j+1) = u(j) - w1 * v(i) - w2 * w(i)
           u(i) = u(i) - omega1 * v(i) - omega2 * w(i)
         enddo    
c .......................................................................
c
c ...
         if(jj.eq.500) then
           jj = 0
           if(my_id.eq.0) write(*,1300),j,dsqrt(dabs(d)),conv 
         endif  
         jj = jj + 1
c ......................................................................
  230 continue
c ......................................................................
c
c ...
      if(my_id.eq.0) then
        write(*,1200) maxit
        if(fLog) write(10,1200) maxit
      endif
      call stop_mef()
  300 continue
c ......................................................................
c
c ... produto:  x*Kx
c
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .           ,x,t 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      xkx =  dot(x,t,neq_doti)
c .......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r = b - Ax (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - t(i)
  310 continue
      tmp  = dot(r,r,neq_doti)
      morm_r = dsqrt(tmp)
      if( tmp .gt. 3.16d0*conv ) then
         if(my_id .eq.0 )then
           write(*,1400) tmp,conv
         endif 
      endif
c ......................................................................
c
c ...
      time = get_time()
      time = time-time0
c .......................................................................
c 
c ...    
      mflops = (flop_bicgstab2(neq,nad,j,2,mpi)*1.d-06)/time  
c ......................................................................
c
c ...
      if(my_id .eq.0 .and. fprint )then
        if(mpi) then
          write(*,1110)tol,conv,j,xkx,norm,morm_r,mflops,time
        else
          write(*,1100)tol,conv,neq,nad,j,xkx,norm,morm_r,mflops,time
        endif
      endif
c .......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       "PBICGSTAB(2): "," it ",j, " x * Kx ",xkx," ||x|| ",norm
     .      ," tol ",tol," time ",time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (/,1x,'SUBROTINA BICCSTAB(2):',/,1x
     .       ,'Coeficiente da diagonal nulo ',i9)
 1100 format(' (PBICGSTAB(2)) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * || b ||        = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1110 format('  (PBICGSTAB_MPI(2)) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||b||m         = ',d20.6/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'Mflops               = ',f20.2/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' PBICGSTAB(2):',5x,'It',i7,5x,2d20.10)
 1400 format (' PBICCSTAB(2):',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
 1510 format (' PBICGSTAB:',1x,'Breakdown:',1x,'(r0)',1x,d20.10)
      end
c **********************************************************************
c
c **********************************************************************
      subroutine pbicgstab(neq     ,nequ  ,nad,ia ,ja 
     .                    ,ad      ,au    ,al ,m  ,b ,x   
     .                    ,t       ,v     ,r  ,p  ,z ,r0 
     .                    ,tol     ,maxit  
     .                    ,matvec  ,dot    
     .                    ,my_id   ,neqf1i,neqf2i 
     .                    ,neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .                    ,fprint  ,flog   ,fnew)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 20/04/2016                                    * 
c * ------------------------------------------------------------------ *   
c * PBICGSTAB : Solucao de sistemas de equacoes pelo metodo dos        * 
c * gradientes biconjugados com precondicionador diagonal para         *
c * matrizes nao-simetricas.                                           *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *                                                                   *
c * neq      - numero de equacoes                                      *
c * nequ     - numero de equacoes no bloco Kuu                         *
c * nad      - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * ad(neq)  - diagonal da matriz A                                    *
c * au(*)    - parte triangular superior de A                          *
c * al(*)    - parte triangular inferior de A                          *
c * m(*)     - precondicionador diagonal                               *
c * b(neq)   - vetor de forcas                                         *
c * x(neq)   - chute inicial                                           *
c * t(neq)   - arranjo local de trabalho                               *
c * v(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * p(neq)   - arranjo local de trabalho                               *
c * z(neq)   - arranjo local de trabalho                               *
c * r0(neq)  - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * matvec   - nome da funcao externa para o produto matrix-vetor      *
c * dot      - nome da funcao externa para o produto escalar           *
c * my_id    -                                                         *
c * neqf1i   -                                                         *
c * neqf2i   -                                                         *
c * neq_doti -                                                         *
c * i_fmap   -                                                         *
c * i_xfi    -                                                         *
c * i_rvcs   -                                                         *
c * i_dspli  -                                                         *
c * fprint   - saida na tela                                           *
c * flog     - log do arquivo de saida                                 *
c * fnew     - .true.  x0 igual a zero                                 *
c *            .false. x0 dado                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * ad(*),al(*),au(*) e b(*)  inalterados                              *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * A(M-1)y=b precondicionador a direita                               *
c * ------------------------------------------------------------------ *                                                                      *
c **********************************************************************
      implicit none
c ...
      real*8 get_time
c ......................................................................
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................   
      integer neq,nequ
      integer maxit,i,j,jj,k
      integer*8 ia(*),nad
      integer ja(*),my_id
      real*8  ad(*),au(*),al(*),m(*),x(*),b(*)
      real*8  r(*),p(*),t(*),v(*),z(*),r0(*)
      real*8  dot,tol,conv,d,alpha,beta,rr0,w,xkx,norm,tmp,norm_r
      real*8  time0,time
      real*8  breaktol,btol
      parameter (btol = 1.d-32)
      real*8  dum1 
      logical flog,fprint,fnew,f
      external matvec,dot
c ======================================================================
      time0 = get_time()
c ...................................................................... 
c
c ...
      do 5 i = 1, neq
        if(ad(i) .eq. 0.d0 ) then
          write(*,1000) i
          call stop_mef()
        endif 
   5  continue
c ......................................................................
c
c ... Chute inicial:
c 
      do 10 i = 1, neq
         x(i) = 0.d0
   10 continue
c .......................................................................
c
c ... conv = tol * |b|
      d    = dot(b,b,neq_doti)
      conv = tol*dsqrt(d)
      breaktol = btol*dsqrt(d)
c .......................................................................
c
c ... Ax0
   
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .           ,x,z 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r0(i) = b(i) - z(i)
c ... p = r0
         p(i)  = r0(i)
c ... r = r0
         r(i)  = r0(i)
c ... z = M(-1)p
         z(i)  = p(i)*m(i) 
  100 continue
c .......................................................................
c
c ...
      jj = 1
      do 230 j = 1, maxit
c ... v = Az(j) = AM(-1)p(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .              ,z,v 
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ... alpha = ( r(j),r0 ) / ( AM(-1)p(j), r0 ))
         rr0   = dot(r,r0,neq_doti)
         if( dsqrt(dabs(rr0)) .lt. breaktol) then
           write(*,1510)dabs(rr0)
           call stop_mef() 
         endif 
         alpha = rr0/dot(v,r0,neq_doti)
c .......................................................................
c
c ...
         do 210 i = 1, neq
c ... x(j+1) = x(j) + alpha*M(-1)p
            x(i) = x(i) + alpha * z(i)
c ... s(j)   = r(j) - alpha*AM(-1)p(j)
            r(i) = r(i) - alpha * v(i)
c ... z = M(-1)s
            z(i) = r(i) * m(i)
  210    continue
c ........................................................................
c
c ... (s,s)
         d = dot(r,r,neq_doti)
c ...
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ........................................................................
c
c ... t = Az = AM(-1)s(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1) 
     .              ,z,t 
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c ........................................................................
c
c ... w = ( AM(-1)s(j) ,s(j) ) / ( AM(-1)s(j), AM(-1)s(j) )
         w = dot(t,r,neq_doti) / dot(t,t,neq_doti)
         if( dabs(w) .lt. breaktol) then
           print*,breaktol
           write(*,1515),dabs(w)
           call stop_mef()
         endif 
c ........................................................................
c
c ... 
         do 220 i = 1, neq
c ... x(j+1) = x(j) + w*M(-1)s
            x(i) = x(i) + w*z(i)
c ... r(j+1) = s(j) - w*AM(-1)s
            r(i) = r(i) - w*t(i)
  220    continue
c ........................................................................
c
c ... (r,r)
         d = dot(r,r,neq_doti)
c ...
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ........................................................................
c
c ... beta = ( r(j+1),r0 ) / ( r(j), r0 )) * (alpha/w) 
         beta = (dot(r,r0,neq_doti) / rr0)*(alpha/w)
c .......................................................................
c
c ...
         do 225 i = 1, neq
c ... p(j+1) = r(i) + beta*(p(j)-w*v(i))
             p(i) = r(i) + beta*(p(i)-w*v(i))
c ... z = M(-1)p
             z(i) = p(i)*m(i)
  225    continue
c .......................................................................
c
c ...
         if( jj .eq. 500) then
           jj = 0
           if(my_id.eq.0) write(*,1300),j,dsqrt(dabs(d)),conv 
         endif  
         jj = jj + 1
c ......................................................................
  230 continue
c ......................................................................
c
c ...
      if(my_id.eq.0) then
        write(*,1200) maxit
        if(fLog) write(10,1200) maxit
      endif
      call stop_mef()
  300 continue
c ......................................................................
c
c ... produto:  x*Kx
c
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1),
     .            x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      xkx =  dot(x,z,neq_doti)
c .......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r = b - Ax (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - z(i)
  310 continue
      tmp  = dot(r,r,neq_doti)
      norm_r = dsqrt(tmp)
      if( tmp .gt. 3.16d0*conv ) then
        write(*,1400) tmp,conv
      endif
c ......................................................................
c
c ...
      time = get_time()
      time = time-time0
c .......................................................................
c
c ...
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,norm_r,time
      endif
c .......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       "PBICGSTAB: "," it ",j, " x * Kx ",xkx," ||x|| ",norm
     .      ," tol ",tol," time ",time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (/,1x,'SUBROTINA PBICCSTAB:',/,1x
     .       ,'Coeficiente da diagonal nulo ',i9)
 1100 format(' (PBICGSTAB) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * || b ||        = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' PBICGSTAB:',5x,'It',i7,5x,2d20.10)
 1400 format (' PBICGSTAB:',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
 1510 format (' PBICGSTAB:',1x,'Breakdown:',1x,'(r,r0)',1x,d20.10)
 1515 format (' PBICGSTAB:',1x,'Breakdown:',1x,'w',1x,d20.10)
      end
c *********************************************************************
c
c **********************************************************************
      subroutine icbicgstab(neq     ,nequ  ,nad,ia ,ja 
     1                     ,ad      ,au    ,al ,m  ,b     ,x   
     2                     ,t       ,v     ,r  ,p  ,z     ,r0
     3                     ,tol     ,maxit  
     4                     ,matvec  ,dot    
     5                     ,my_id   ,neqf1i,neqf2i 
     6                     ,neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli
     7                     ,fprint  ,flog   ,fnew)
c **********************************************************************
c * Data de criacao    : 15/04/2016                                    *
c * Data de modificaco : 15/12/2016                                    * 
c * ------------------------------------------------------------------ *   
c * PBICGSTAB : Solucao de sistemas de equacoes pelo metodo dos        * 
c * gradientes biconjugados com precondicionador diagonal para         *
c * matrizes nao-simetricas.                                           *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *                                                                   *
c * neq      - numero de equacoes                                      *
c * nequ     - numero de equacoes no bloco Kuu                         *
c * nad      - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * ad(neq)  - diagonal da matriz A                                    *
c * au(*)    - parte triangular superior de A                          *
c * al(*)    - parte triangular inferior de A                          *
c * m(*)     - precondicionador diagonal                               *
c * b(neq)   - vetor de forcas                                         *
c * x(neq)   - chute inicial                                           *
c * t(neq)   - arranjo local de trabalho                               *
c * v(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * p(neq)   - arranjo local de trabalho                               *
c * z(neq)   - arranjo local de trabalho                               *
c * r0(neq)  - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * matvec   - nome da funcao externa para o produto matrix-vetor      *
c * dot      - nome da funcao externa para o produto escalar           *
c * my_id    -                                                         *
c * neqf1i   -                                                         *
c * neqf2i   -                                                         *
c * neq_doti -                                                         *
c * i_fmap   -                                                         *
c * i_xfi    -                                                         *
c * i_rvcs   -                                                         *
c * i_dspli  -                                                         *
c * fprint   - saida na tela                                           *
c * flog     - log do arquivo de saida                                 *
c * fnew     - .true.  x0 igual a zero                                 *
c *            .false. x0 dado                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq) - modificado                                                *
c * ad(*),al(*),au(*) - inalterados                                    *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * A(M-1)y=b precondicionador direita                                 *
c * ------------------------------------------------------------------ *                                                                      *
c **********************************************************************
      implicit none
c ...
      real*8 get_time
c ......................................................................
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................   
      integer neq,nequ
      integer maxit,i,j,jj,k
      integer*8 ia(*),nad
      integer ja(*),my_id
      real*8  ad(*),au(*),al(*),m(*),x(*),b(*)
      real*8  r(*),p(*),t(*),v(*),z(*),r0(*)
      real*8  dot,tol,conv,xkx,norm,d,alpha,beta,rr0,w,tmp
      real*8  time0,time
      real*8  dum1 
      logical flog,fprint,fnew
      external matvec,dot
c ======================================================================
      time0 = get_time()
c ......................................................................
c
c ...
      do 5 i = 1, neq
        if(ad(i) .eq. 0.d0 ) then
          write(*,1000) i
          call stop_mef()
        endif 
   5  continue
c ......................................................................
c
c ... Chute inicial:
c 
      do 10 i = 1, neq
         x(i) = 0.d0
   10 continue
c .......................................................................
c
c ... conv = tol * |b|
      d    = dot(b,b,neq_doti)
      conv = tol*dsqrt(dabs(d))
c .......................................................................
c
c ... Ax0
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1),
     .            x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r0(i) = b(i) - z(i)
c ... r = r0
         p(i) = r0(i)
c ... p = r0
         r(i) = r0(i)
  100 continue
c .......................................................................
c
c ... Mz=p  
      call ildlt_solv(neq,ia,ja,m,m(neq+1),p,z)
c .......................................................................      
c
c ...
      jj = 1
      do 230 j = 1, maxit
c ... v = Az(j) = AM(-1)p(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1),
     .               z,v,
     .               neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ... alpha = ( r(j),r0 ) / ( AM(-1)p(j), r0 ))
         rr0   = dot(r,r0,neq_doti)
         alpha = rr0/dot(v,r0,neq_doti)
c .......................................................................
c
c ...
         do 210 i = 1, neq
c ... x(j+1) = x(j) + alpha*M(-1)p
            x(i) = x(i) + alpha * z(i)
c ... s(j)   = r(j) - alpha*AM(-1)p(j)
            r(i) = r(i) - alpha * v(i)
  210    continue
c ........................................................................
c
c
c ... (s,s)
         d = dot(r,r,neq_doti)
c ...
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ........................................................................
c
c ... Mz=r  
         call ildlt_solv(neq,ia,ja,m,m(neq+1),r,z)
c .......................................................................
c
c ... (r,r)
         d = dot(r,r,neq_doti)
c ...
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ........................................................................
c
c ...  t = Az = AM(-1)s(j)
         call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1),
     .               z,t,
     .               neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c ........................................................................
c
c ... w = ( AM(-1)s(j) ,s(j) ) / ( AM(-1)s(j), AM(-1)s(j) )
         w = dot(t,r,neq_doti) / dot(t,t,neq_doti)
c ........................................................................
c
c ... 
         do 220 i = 1, neq
c ... x(j+1) = x(j) + w*M(-1)s
            x(i) = x(i) + w*z(i)
c ... r(j+1) = s(j) - w*AM(-1)s
            r(i) = r(i) - w*t(i)
  220    continue
c ........................................................................
c
c ... (r,r)
         d = dot(r,r,neq_doti)
c ...
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ........................................................................
c
c ... beta = ( r(j+1),r0 ) / ( r(j), r0 )) * (alpha/w) 
         beta = (dot(r,r0,neq_doti) / rr0)*(alpha/w)
c .......................................................................
c
c ...
         do 225 i = 1, neq
c ... p(j+1) = r(i) + beta*(p(j)-w*v(i))
             p(i) = r(i) + beta*(p(i)-w*v(i))
  225    continue
c .......................................................................
c
c ... Mz=p  
         call ildlt_solv(neq,ia,ja,m,m(neq+1),p,z)
c .......................................................................
c
c ...
         if( jj .eq. 500) then
           jj = 0
           if(my_id.eq.0) write(*,1300),j,dsqrt(dabs(d)),conv 
         endif  
         jj = jj + 1
c ......................................................................
  230 continue
c ......................................................................
      if(my_id.eq.0) then
        write(*,1200) maxit
        if(fLog) write(10,1200) maxit
      endif
      call stop_mef()
  300 continue
c
c ... produto:  x*Kx
c
      call matvec(neq,nequ,ia,ja,ia(neq+2),ja(nad+1),ad,al,al(nad+1),
     .            x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
       xkx = dot(x,z,neq_doti)
c .......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r = b - Ax (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - z(i)
  310 continue
      tmp  = dot(r,r,neq_doti)
      tmp = dsqrt(tmp)
      if( tmp .gt. 3.16d0*conv ) then
         if(my_id .eq.0 )then
           write(*,1400) tmp,conv
         endif 
      endif
c ......................................................................
c
c ...
      time = get_time()
      time = time-time0
c .......................................................................
c
c ...
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,time
      endif
c .......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       "ICBICGSTAB: "," it ",j, " x * Kx ",xkx," ||x|| ",norm
     .      ," tol ",tol," time ",time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (/,1x,'SUBROTINA ICBICCSTAB:',/,1x
     .       ,'Coeficiente da diagonal nulo ',i9)
 1100 format(' (ICBICGSTAB) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * || b ||        = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' ICBICCSTAB:',5x,'It',i7,5x,2d20.10)
 1400 format (' ICBICCSTAB:',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
      end
c *********************************************************************
c
c *********************************************************************
      subroutine pcg_block_it(neq   ,nequ  ,neqp
     .                       ,nad   ,naduu ,nadpp
     .                       ,iau   ,jau
     .                       ,iap   ,jap
     .                       ,iapu  ,japu
     .                       ,adu   ,adp
     .                       ,alu   ,alp  ,alpu
     .                       ,mu    ,mp    ,b     ,x
     .                       ,z     ,r     ,s
     .                       ,bu    ,bp    ,bu0   ,bp0
     .                       ,u     ,p
     .                       ,tol   ,ctol  ,maxit ,cmaxit,alfap ,alfau
     .                       ,fnew  ,istep
     .                       ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .                       ,i_xfi ,i_rcvsi,i_dspli)
c **********************************************************************
c * Data de criacao    : 12/12/2015                                    *
c * Data de modificaco : 13/12/2016                                    *
c * ------------------------------------------------------------------ *
c * Subroutine PCG_BLOCK_IT                                            *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * neq      - numero de equacoes                                      *
c * nequ     - numero de equacoes no bloco Kuu                         *
c * neqp     - numero de equacoes no bloco Kpp                         *
c * nad      - numero de termos nao nulos no bloco Kuu, kpp e kpu      *
c * naduu    - numero de termos nao nulos no bloco Kuu                 *
c * nadpp    - numero de termos nao nulos no bloco Kpp                 *
c * iau(*)   - ponteiro do formato CSR (bloco Kuu)                     *
c * jau(*)   - ponteiro das colunas no formato CSR (bloco Kuu)         *
c * iap(*)   - ponteiro do formato CSR (bloco Kp)                      *
c * jap(*)   - ponteiro das colunas no formato CSR (bloco Kpp)         *
c * iapu(*)  - ponteiro do formato CSR (bloco Kpu )                    *
c * japu(*)  - ponteiro das colunas no formato CSR (bloco Kpu)         *
c * adu(*)   - diagonal da matriz A  (bloco Kuu)                       *
c * alu(*)   - parte triangular inferior de A (bloco Kuu)              *
c * adp(*)   - diagonal da matriz A  (bloco Kpp)                       *
c * alp(*)   - parte triangular inferior de A (bloco Kpp)              *
c * alpu(*)  - parte triangular inferior de A (bloco Kpu)              *
c * mu(*)    - precondicionador diagonal (bloco Kuu)                   *
c * mp(*)    - precondicionador diagonal (bloco Kuu)                   *
c * b(neq)   - vetor de forcas                                         *
c * x(neq)   - valores do passo anterior                               *
c * z(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * s(neq)   - arranjo local de trabalho                               *
c * bu(nequ) - arranjo local de trabalho                               *
c * bp(neqp) - arranjo local de trabalho                               *
c * bu0(nequ)- arranjo local de trabalho                               *
c * bp0(neqp)- arranjo local de trabalho                               *
c * u(nequ)  - arranjo local de trabalho                               *
c * p(neqp)  - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * ctol     - tolerancia de convergencia do ciclo externo             *
c * cmaxit   - numero maximo de iteracoes do ciclo externo             *
c * fnew     - .true. chute inicial nulo                               *
c *          .false. passo anterior                                    *
c * istep    - numero do passo de tempo                                *
c * my_id    -                                                         *
c * neqf1i   -                                                         *
c * neqf2i   -                                                         *
c * neq_doti -                                                         *
c * i_fmap   -                                                         *
c * i_xfi    -                                                         *
c * i_rvcs   -                                                         *
c * i_dspli  -                                                         *
c * matvec   - nome da funcao externa para o produto matrix-vetor      *
c * dot      - nome da funcao externa para o produto escalar           *
c * energy   - nao definido                                            *
c * ------------------------------------------------------------------ *
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * ------------------------------------------------------------------ *
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************
      implicit none
c ...
      real*8 get_time
c ......................................................................
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................
      integer neq,nequ,neqp,nad,naduu,nadpp,maxit,i,j,jj,istep
      integer cmaxit
      integer iau(*),jau(*),iap(*),jap(*),iapu(*),japu(*)
      integer my_id,idum
      real*8 adu(*),adp(*),alu(*),alp(*),alpu(*)
      real*8 mu(*),mp(*),x(*),r(*),z(*),b(*),s(*)
      real*8 bp(*),bu(*),bp0(*),bu0(*)
      real*8 u(*),p(*),xkx,norm
      real*8 dot,dot_par,tol,ctol,d,alpha,beta
      real*8 resid_u,resid_p,p_conv,u_conv,alfap,alfau
      real*8 time0,time,time_csr
      real*8 dum1
      logical l_u_conv,l_p_conv,fnew
      external matvec_csrc_sym_pm
      external dot_par,dot
c ======================================================================
      time0    = get_time()
      time_csr = 0.d0
c ......................................................................
c
c ...
      time0 = get_time()
c ... 
c     alfap = 0.1d0
c     alfau = 0.2d0
c.......................................................................
c 
c ...
      if(fnew) then
        do j = 1, nequ 
          u(j) = 0.d0 
        enddo
        do j = 1, neqp
          p(j) = 0.d0        
        enddo
c.......................................................................
c 
c ...
      else
        do j = 1, nequ 
          u(j) = x(j) 
        enddo
        do j = 1, neqp
          p(j) = x(nequ+j) 
        enddo
      endif
c.......................................................................
c 
c ... 
      l_u_conv = .false.
      l_p_conv = .false.
c.......................................................................
c 
c ... bu 
      call aequalb(bu0,b,nequ)
c ... bp
      call aequalb(bp0,b(nequ+1),neqp)
c.......................................................................
c
c ... 
c     print*,ctol,cmaxit,maxit,tol
      jj = 0
      do i = 1, cmaxit
c
c ...  rp= Fp-kpu*U
        time_csr = get_time() - time_csr
        call matvec_csr_pm(neqp,nequ,iapu,japu,alpu,u,r,.false.)
        call aminusb(bp0,r,bp,neqp) 
        time_csr = get_time() - time_csr
c ......................................................................
c
c ... P = inv(Kpp)*(Fp - kpu*U)= inv(Kpp)*rp
        call pcg(neqp      ,neqp       ,nadpp,iap       ,jap
     1          ,adp       ,alp        ,alp  ,mp        ,bp  
     2          ,x(nequ+1) ,z          ,r    ,s     
     3          ,tol       ,maxit
     4          ,matvec_csrc_sym_pm,dot_par 
     5          ,my_id     ,neqf1i     ,neqf2i,neqp    ,i_fmapi
     6          ,i_xfi     ,i_rcvsi    ,i_dspli
     7          ,.false.   ,.false.    ,.false.,.false.
     8          ,1         ,.false.)
c ......................................................................
c
c ... x - > p
       do j = 1, neqp
         p(j) = (1.d0-alfap)*p(j) + alfap*x(nequ+j)
       enddo
c ......................................................................
c
c ...  ru= Fu-kup*P
        time_csr = get_time() - time_csr
        call matvec_csr_pm(neqp,nequ,iapu,japu,alpu,p,r,.true.)
        call aminusb(bu0,r,bu,nequ)
        time_csr = get_time() - time_csr
c ......................................................................
c
c ... U = inv(Kuu)*(Fu - kup*P)=inv(Kuu)*ru
        call pcg(nequ   ,nequ   ,naduu,iau    ,jau
     1          ,adu    ,alu    ,alu  ,mu     ,bu   
     2          ,x      ,z      ,r    ,s      
     3          ,tol    ,maxit
     4          ,matvec_csrc_sym_pm,dot_par 
     5          ,my_id  ,neqf1i ,neqf2i,nequ    ,i_fmapi
     6          ,i_xfi  ,i_rcvsi,i_dspli
     7          ,.false.   ,.false.    ,.false.,.false.
     8          ,1         ,.false.)
c ......................................................................
c
c ... x - > u
       do j = 1, nequ 
         u(j) = (1.d0-alfau)*u(j) + alfau*x(j) 
       enddo
c ......................................................................
c           
c ... Kpp*P
        time_csr = get_time() - time_csr
        call matvec_csrc_sym_pm(neqp      ,neqp
     .                        ,iap       ,jap ,idum     ,idum        
     .                        ,adp       ,alp ,alp       
     .                        ,p         ,z
     .                        ,neqf1i,neqf2i
     .                        ,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c ... Kpu*U
        call matvec_csr_pm(neqp,nequ,iapu,japu,alpu,u,r,.false.)
c ...
        call aminusb(bp0,r,bp,neqp)
        call aminusb(bp ,z,bp,neqp)
        time_csr = get_time() - time_csr
        resid_p = dsqrt(dot(bp,bp,neqp))
c .....................................................................
c            
c ... Kuu*U
        time_csr = get_time() - time_csr
        call matvec_csrc_sym_pm(nequ      ,nequ
     .                        ,iau       ,jau ,idum     ,idum         
     .                        ,adu       ,alu ,alu              
     .                        ,u         ,z
     .                        ,neqf1i,neqf2i
     .                        ,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c ... Kpu*P
        call matvec_csr_pm(neqp,nequ,iapu,japu,alpu,p,r,.true.)
c ...
        call aminusb(bu0,r,bu,nequ)
        call aminusb(bu,z,bu,nequ)
        time_csr = get_time() - time_csr
        resid_u = dsqrt(dot(bu,bu,nequ))
c ... 
        if( i .eq. 1) then
          u_conv = dsqrt(dot(bu,bu,nequ))*ctol
          p_conv = dsqrt(dot(bp,bp,neqp))*ctol
c ......................................................................
c 
c ... 
        else 
          if( resid_u .lt. u_conv ) l_u_conv = .true.
          if( resid_p .lt. p_conv ) l_p_conv = .true.
        endif
c ......................................................................
        jj = jj + 1
        if( jj . eq. 10 ) then
          write(*,200),i,resid_p ,p_conv,resid_u,u_conv
          jj = 0
        endif
        if( l_u_conv .and. l_p_conv ) goto 100
        if( i .eq. 500) then
          print*, 'MAXIT'
          stop
        endif
      enddo
  100 continue   
c ......................................................................
      time = get_time()
c ......................................................................
c
c ... calculo de x*(kx)
c
c ... x*kx = u*bu + p*bp
      xkx = dot(u,bu0,nequ) + dot(p,bp0,neqp)
c ......................................................................
c
c ... norm-2 = || x || = || u || + || p ||
      norm = dsqrt(dot(u,u,nequ)) + dsqrt(dot(p,p,neqp))
c ......................................................................
c
c ...
      write(*,300)ctol,neq,i,xkx,norm,time-time0
c ......................................................................
c
c ...
      write(16,*) istep,i,time-time0,time_csr,xkx,norm
c ......................................................................
c
c ...
      do j = 1, nequ 
        x(j) = u(j) 
      enddo
      do j = 1, neqp
        x(nequ+j) = p(j) 
      enddo
c ......................................................................
      return
 200  format (1x,'It',1x,i4,1x,'rp',1x,2d10.2,1x,'ru',1x,2d10.2)
 300  format(' (BLOCKPCG) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
      end
c **********************************************************************
c
c **********************************************************************
      real*8 function smachn()
c **********************************************************************
c *                                                                    *
c *   SMACHN: calcula a precisao da maquina para real*8                *
c *                                                                    *
c **********************************************************************
      implicit none
      smachn = 1.0d0
100   smachn = smachn*0.5d0
      if(smachn+1.d0 .ne. 1.d0) go to 100
c     smachn = 2.0d0*smachn
      return
      end
c ***********************************************************************
c     subroutine bicgstab(neq,ia,ja,ad,au,al,m,b,x,y,z,p,r,s,tol,maxit,
c    .                    matvec,dot,my_id,neqf1i,neqf2i,neq_doti,
c    .                    i_fmapi,i_xfi,i_rcvsi,i_dspli)
c **********************************************************************
c *                                                                    *
c *   Subroutine BICGSTAB                                              *
c *                                                                    *
c *   Solucao de sistemas de equacoes pelo metodo dos gradientes       *
c *   conjugados com precondicionador diagonal para matrizes           *
c *   simetricas.                                                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq    - numero de equacoes                                      *
c *   ia(*)  - ponteiro do formato CSR                                 *
c *   ja(*)  - ponteiro das colunas no formato CSR                     *
c *   ad(neq)- diagonal da matriz A                                    *
c *   au(*)  - parte triangular superior de A                          *
c *   al(*)  - parte triangular inferior de A                          *
c *   m(*)   - precondicionador diagonal                               *
c *   b(neq) - vetor de forcas                                         *
c *   x(neq) - chute inicial                                           *
c *   z(neq) - arranjo local de trabalho                               *
c *   r(neq) - arranjo local de trabalho                               *
c *   tol    - tolerancia de convergencia                              *
c *   maxit  - numero maximo de iteracoes                              *
c *   matvec - nome da funcao externa para o produto matrix-vetor      *
c *   dot    - nome da funcao externa para o produto escalar           *
c *   energy - nao definido                                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   x(neq) - vetor solucao                                           *
c *   b(neq) - modificado                                              *
c *   ad(*),al(*),au(*) - modificado                                   *
c *                                                                    *
c **********************************************************************
c     implicit none
c     include 'mpif.h'
c     integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
c     integer*8 i_fmapi,i_xfi
c     integer*8 i_rcvsi,i_dspli
c .....................................................................
c     integer neq,maxit,nad,i,j,k
c     integer ia(*),ja(*),my_id
c     real*8  ad(*),au(*),al(*),m(*),x(*),r(*),p(*),b(*),y(*),z(*),s(*)
c     real*8  dot,ddot,tol,conv,energy,d,alpha,beta,rr0,w
c     real*8  time0,time
c     real*8 dum1
c     external matvec,dot
c ======================================================================
c     time0 = get_time() 
c ......................................................................
c     nad = ia(neq+1)-1
c     if(my_id.eq.0)print *, 'nad :',nad
c ......................................................................
c
c ... Chute inicial:
c
c     do 10 i = 1, neq
c        x(i) = 0.d0
c ...... pre-condicionador diagonal:         
c        b(i)  = b(i)/m(i)
c        ad(i) = ad(i)/m(i)
c        do 5 k = ia(i), ia(i+1)-1
c           j = ja(k)
c           al(k) = al(k) / m(i)
c           au(k) = au(k) / m(j)
c  5     continue      
c 10  continue
c ----------------------------------------------------------------------
c     call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,r,
c    .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)    
c     do 100 i = 1, neq
c        r(i) = b(i) - r(i)
c        p(i) = r(i)
c        b(i) = r(i) 
c 100 continue
c     d    = dot(r(1),r(1),neq)
c     conv = tol*dsqrt(dabs(d))
c ----------------------------------------------------------------------
c     do 230 j = 1, maxit
c        call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
c    .               p,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,
c    .               dum1)
c        rr0 = dot(r,b,neq)
c        alpha = rr0/dot(z,b,neq)
c         do 210 i = 1, neq
c           s(i) = r(i) - alpha * z(i)
c 210    continue
c        call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
c    .               s,y, neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,
c    .               dum1)
c        w = dot(y,s,neq) / dot(y,y,neq)
c        do 220 i = 1, neq
c           x(i) = x(i) + alpha*p(i) + w * s(i)
c           r(i) = s(i) - w*y(i)
c 220    continue
c        beta = (dot(r,b,neq) / rr0)*(alpha/w)
c        do 225 i = 1, neq
c            p(i) = r(i) + beta*(p(i)-w*z(i))
c 225    continue
c        d = dot(r,r,neq)  
c        if (dsqrt(dabs(d)) .lt. conv) goto 300
c 230 continue
c ----------------------------------------------------------------------
c     write(*,1200) maxit
c     call stop_mef()
c 300 continue
c
c ... Energy norm:
c
c      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,z,
c     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli)
c     b(1:neq) = b(1:neq)*m(1:neq)  
c     energy   = dot(x(1),b(1),neq)
c ......................................................................
c     time = get_time()
c     time = time-time0
c ----------------------------------------------------------------------
c     if(my_id.eq.0)write(*,1100) neq,j,energy,time
c ......................................................................
c     Controle de flops
c     if(my_id.eq.0)write(10,'(a,a,i9,a,d20.10,a,d20.10,f20.2)')
c    .              "BICGSTAB: ", "it",j, " energy norm ",energy,
c    .              " tol ",tol,"time",time
c ......................................................................
c     return
c ======================================================================
c1000 format (//,5x,'SUBROTINA BICGSTAB:',/,5x,'Coeficiente da diagonal
c    . nulo ou negativo - equacao ',i7)
c1100 format(' (BICGSTAB) solver:'/
c    . 5x,'Number of equations  = ',i20/
c    . 5x,'Number of iterations = ',i20/
c    . 5x,'Energy norm          = ',d20.6/
c    . 5x,'CPU time (s)         = ',f20.2/)
c1200 format (' *** WARNING: No convergence reached after ',i4,
c    .        ' iterations !',/)
c     end
c **********************************************************************
      subroutine sym_ortho(a,b,c,s,r)
c **********************************************************************
c * Data de criacao    : 18/05/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *
c * SYM_ORTHO: Givens rotation                                        * 
c * (versa com melhor comportamento numerico)                          *      
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * a   - paramentro                                                   *
c * b   - paramentro s da diagonal principal                           *
c * c   - nao definido                                                 *
c * s   - nao definido                                                 *
c * r   - nao definido                                                 *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * c   - cos(teta)                                                    *
c * s   - seno(teta)                                                   *
c * r   - raiz(a^2 + b^2)                                              *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * | c  s | | a |    | raiz(a^2 + b^2) |   | r |                      *
c   |      | |   |  = |                 | = |   |                      *
c * | s -c | | b |    |        0        |   | 0 |                      *
c * ------------------------------------------------------------------ *
c **********************************************************************
      implicit none
      real*8 a,b,c,s,r,t,ma,mb,sa,sb
      real*8 sign1
c ...      
      ma = dabs(a)
      mb = dabs(b)
      sa = sign1(a)
      sb = sign1(b)
c .....................................................................
c
c ...        
      if(b .eq. 0.d0) then
        s = 0.d0
        r = ma
        if ( a .eq. 0.d0) then
          c = 1.0d0
        else
          c = sa
        endif
      else if( a .eq. 0.d0) then
        c = 0.d0
        s = sb
        r = mb
      else if( mb .gt. a) then
        t = a/b
        s = sb/dsqrt(1.d0 + t*t)
        c = s*t
        r = b/s
      else if( ma .gt. mb ) then
        t = b/a
        c = sa/dsqrt(1.d0+t*t)
        s = c*t
        r = a/c
      endif
      return
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 18/05/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *
c * SIGN1 : retorna o sinal de a                                       * 
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * a   - paramentro                                                   *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * 0.0d0, 1.d0 ou -1.d0
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************
      real*8 function sign1(a)
      implicit none
      real*8 a
      sign1 = 0.0d0
      if( a .gt. 0.d0) then
        sign1 = 1.d0
      else if( a .lt. 0.d0) then
        sign1 = -1.d0
      endif
      return
      end
c ***********************************************************************
c
c ***********************************************************************
c * Data de criacao    : 30/05/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *
c * SYM_ORTHO: Givens rotation                                         * 
c * (versa com melhor comportamento numerico)                          *      
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * a   - paramentro                                                   *
c * b   - paramentro s da diagonal principal                           *
c * c   - nao definido                                                 *
c * s   - nao definido                                                 *
c * r   - nao definido                                                 *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * c   - cos(teta)                                                    *
c * s   - seno(teta)                                                   *
c * r   - raiz(a^2 + b^2)                                              *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * | c  s | | a |    | raiz(a^2 + b^2) |   | r |                      *
c   |      | |   |  = |                 | = |   |                      *
c * |-s  c | | b |    |        0        |   | 0 |                      *
c * ------------------------------------------------------------------ *
c **********************************************************************
      subroutine sym_ortho2(a,b,c,s,r)
      implicit none
      real*8 a,b,c,s,r,t,ma,mb
      real*8 sign1
c ...      
      ma = dabs(a)
      mb = dabs(b)
c .....................................................................
c
c ...     
      c = 1.d0
      s = 0.d0   
      r = a
      if(b .ne. 0.d0) then
        if( mb .gt. ma) then
          t = a/b
          s = 1.d0/dsqrt(1.d0 + t*t)
          c = s*t
          r = b/s
        else 
          t = b/a
          c = 1.d0/dsqrt(1.d0+t*t)
          s = c*t
          r = a/c
        endif
      endif
      return
      end
c ********************************************************************** 
c
c ***********************************************************************
c * Data de criacao    : 05/04/2019                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *
c * MEMORY_PCG : calcula a memoria necessaria no CG                    *    
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * neq - numero de equacoes                                           *
c * nad - numero de termos nao nulos fora da diagonal principal        *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************
      real*8 function memory_pcg(neq,nad)
      implicit none
      integer neq
      integer*8 nad
      real*8 memory_csrc,solver,matvec
c ...      
      solver = 8.0*4.0*neq
      matvec = 12.0*neq + 12*nad + 4.0
c .....................................................................
c
c ...     
      memory_pcg = (solver + matvec)/(1024.d0**2)
      return
      end
c ********************************************************************** 
