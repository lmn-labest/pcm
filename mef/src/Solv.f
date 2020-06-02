c **********************************************************************
c *                                                                    *
c *   SOLV_PCM                                             06/12/2015  *
c *                                                        01/06/2020  *
c *   Metodos iterativos de solucao:                                   *
c *                                                                    *
c *   cg                                                               *
c *   pcg                                                              *
c *   bcg                                                              *
c *   iccg                                                             *
c *   sqrm                                                             *
c *   rsqmr                                                            *
c *   lsqmr                                                            *
c *   gmres2(m)                                                        *
c **********************************************************************
      subroutine solv_pm(neq    ,nequ    ,neqp
     1                ,nad    ,naduu   ,nadpp
     2                ,ip     ,ja      ,ad         ,au     ,al 
     3                ,m      ,b       ,x          ,tol    ,maxit
     4                ,ngram  ,block_pu,n_blocks_up,solver ,istep
     5                ,cmaxit ,ctol    ,alfap      ,alfau  ,precond
     6                ,fmec   ,fporomec,fterm      ,fhist_solv ,fprint)
      use Malloc
      implicit none
      include 'precond.fi'
      include 'openmp.fi'
      include 'time.fi'
c ... ponteiros      
      integer*8 i_a,i_b,i_c,i_d,i_g,i_h,i_y,i_z,i_r,i_s
c ......................................................................
      integer*8 nad,ip(*)
      integer ja(*),neq,nequ,neqp,naduu,nadpp
      integer maxit,solver,ngram,istep,n_blocks_up
      real*8  ad(*),au(*),al(*),m(*),x(*),b(*),tol,energy
c ... pcg duplo
      integer cmaxit
      real*8  ctol,alfap,alfau
c ......................................................................
      logical fmec,fporomec,fterm,fhist_solv,fprint
      logical fsqmr,block_pu
c ... precondicionador
      integer precond
      real*8  max_block_a(max_block*max_block)
c ...................................................................... 
c
c ...
      if (omp_solv) then
c ... matriz aramazenada em csrc blocado (Kuu,Kpp,Kpu)
         if(block_pu) then
           pmatrixtime = get_time() - pmatrixtime 
           i_threads_y = alloc_8('buffer_y',nth_solv,neq)
           call partition_matrix(ip,ja,neq,nequ,nad,.false.,block_pu)
           pmatrixtime = get_time() - pmatrixtime
c .......................................................................
c
c ... matriz aramazenada no csrc simetrico (Kuu,-Kpp,-Kpu)
         else   
           pmatrixtime = get_time() - pmatrixtime 
           i_threads_y = alloc_8('buffer_y',nth_solv,neq)
           call partition_matrix(ip,ja,neq,nequ,nad,.false.,block_pu)
           pmatrixtime = get_time() - pmatrixtime
         endif
      endif
c ......................................................................
c
c ... usa o gmres caso o sqmr falhar
      fsqmr = .false.
  100 continue
c ......................................................................
c  
c ... PCG
      if (solver .eq. 1) then
c ...
        i_z = alloc_8('zsolver ',1,neq)
        i_r = alloc_8('rsolver ',1,neq)
        i_s = alloc_8('psolver ',1,neq)
c ......................................................................
c
c ... calculo do precondicionador
        call cal_precond(ip     ,ja    ,m
     1                  ,ad     ,al    ,ia(i_z)
     2                  ,precond,neq   ,nequ
     3                  ,0      ,0     ,0      
     4                  ,0      ,0     ,0      
     5                  ,.false.,0    )
c .....................................................................
c
c ... matriz aramazena em csrc blocado (Kuu,Kpp,Kpu)
        if(block_pu) then
           print*,"PCG nao disponivel para a matriz", 
     .            " blocada nao simetrica !!"
          stop   
c .....................................................................
c
c ... matriz aramazenada no csrc simetrico (Kuu,-Kpp,-Kpu)
        else
c ... omp
          if(omp_solv) then
            call call_cg_omp(neq   ,nequ    ,nad   ,ip      ,ja
     1                   ,ad       ,al     ,m      ,b       ,x   
     2                   ,ia(i_z)  ,ia(i_r),ia(i_s) 
     3                   ,tol      ,maxit  ,precond,iparam ,fhist_solv 
     4                   ,0        ,0      ,0     ,neq     ,0      
     5                   ,0       ,0       ,0      ,ia(i_threads_y)
     6                   ,1       ,.false. ,.false.)  
c .....................................................................
c
c ... sequencial (cg, pcg e iccg)
          else
            call call_cg(neq      ,nequ   ,nad    ,ip      ,ja
     1                   ,ad      ,al     ,m      ,b       ,x   
     2                   ,ia(i_z) ,ia(i_r),ia(i_s) 
     3                   ,tol     ,maxit  ,precond,iparam ,fhist_solv 
     4                   ,0       ,0      ,0      ,neq     ,0      
     5                   ,0       ,0      ,0      
     6                   ,1       ,.false.,.false.)  
c .....................................................................
c
          endif      
c .....................................................................
        endif
c .....................................................................
c
c ...
        i_s = dealloc('psolver ')
        i_r = dealloc('rsolver ')
        i_z = dealloc('zsolver ')
c ......................................................................
c
c ... GMRES :
      elseif(solver .eq. 2 .or. fsqmr) then
c ...
        if(block_pu) then
          if(n_blocks_up .eq. 1 ) then
            print*,"GMRES nao disponivel para a matriz",
     .             " blocada [Kuu Kpp]."
            stop   
          else if( n_blocks_up .eq. 3 ) then
            print*,"GMRES nao disponivel para a matriz",
     .           " blocada [Kuu] [Kpp] [Kpu]."
            stop   
          endif
        endif
c ......................................................................
c
c ...
         i_g = alloc_8('gsolver ',neq,ngram+1)         
         i_h = alloc_8('hsolver ',ngram+1,ngram)
         i_y = alloc_8('ysolver ',1,ngram)
         i_c = alloc_8('csolver ',1,ngram)
         i_s = alloc_8('ssolver ',1,ngram)
         i_r = alloc_8('rsolver ',1,ngram+1)
c .....................................................................
c
c ... calculo do precondicionador
        call cal_precond(ip     ,ja    ,m
     1                  ,ad     ,al    ,ia(i_g) 
     2                  ,precond,neq   ,nequ
     3                  ,0      ,0     ,0      
     4                  ,0     ,0      ,0      
     5                  ,.false.,0)
c .....................................................................
c
c ... matriz aramazenada no csrc simetrico (Kuu,-Kpp,-Kpu)
c ... omp
         if(omp_solv) then
           call call_gmres_omp(neq   ,neq,nad  
     1                     ,ip     ,ja     ,ad     ,al 
     2                     ,m      ,b      ,x      ,ia(i_g)
     3                     ,ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r)  
     4                     ,tol    ,maxit  ,precond,fhist_solv,fprint
     5                     ,ngram  
     5                     ,0      ,0      ,0     ,neq     ,0      
     6                     ,0      ,0      ,0      ,ia(i_threads_y)
     7                     ,1      ,.false. ,.false.) 
c .....................................................................
c
c ... (sequencial+mpi)
         else 
           call call_gmres(neq   ,neq,nad  
     1                     ,ip     ,ja     ,ad     ,au     ,al 
     2                     ,m      ,b      ,x      ,ia(i_g)
     3                     ,ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r)  
     4                     ,tol    ,maxit  ,precond,fhist_solv,fprint
     5                     ,ngram
     6                     ,0      ,0      ,0     ,neq     ,0      
     7                     ,0      ,0      ,0      
     7                     ,1      ,.false. ,.false.) 
         endif 
c .....................................................................
c
c ......................................................................
         i_r = dealloc('rsolver ')
         i_s = dealloc('ssolver ')
         i_c = dealloc('csolver ')
         i_y = dealloc('ysolver ')
         i_h = dealloc('hsolver ')
         i_g = dealloc('gsolver ')
c ......................................................................
c
c ...                                         
      elseif(solver .eq. 3) then
        print*,"GAUSS: Not available !!"
c ......................................................................
c
c ... BICGSTAB :
      else if (solver .eq. 4) then
c ...
        print*,"BICGSTAB: Not available !!"
c ......................................................................  
c
c ... PCG_BLOCK_IT                                                                   
      else if (solver .eq. 5) then
c ...
        print*,"PCG_BLOCK_ITBICGSTAB(2): Not available !!"
c ......................................................................
c
c ... BICGSTAB(2):
      else if(solver .eq. 6 ) then
c ...
        print*,"BICGSTAB(2): Not available !!"
c ...................................................................... 
c
c ... MINRES:
      else if(solver .eq. 7 ) then
c ...
        print*,"MINRES: Not available !!"   
c .....................................................................
c
c ... CR - Conjugate Residual
      else if(solver .eq. 8 ) then
c ...
        print*,"CR: Not available !!"
c .....................................................................
c
c ... SYMMLQ:
      else if(solver .eq. 9 ) then
c ...
        print*,"SYMMLQ: Not available !!" 
c .....................................................................
c
c ... mkl_pardiso
      else if(solver .eq. 10 ) then
        i_z = alloc_8('zsolver ',1,neq)
        i_y = alloc_4('ysolver ',1,neq+1)
        if(fmec .or. fterm) then
          call call_mkl_pardiso(neq,nad,ip,ja,ad,b,x,ia(i_z),ia(i_y),2)
        else if(fporomec) then
          call call_mkl_pardiso(neq,nad,ip,ja,ad,b,x,ia(i_z),ia(i_y),-2)
        endif
        i_y  = dealloc('ysolver ')
        i_z  = dealloc('zsolver ') 
c ......................................................................
c
c ... SQRM - QRM simetrico
      else if(solver .eq. 11) then
c ... matriz aramazena em csrc blocado nao simentrico (Kuu,Kpp,Kpu)
        if(block_pu) then
          print*,"SQRM nao disponivel para a matriz", 
     .           " blocada nao simetrica !!"
          stop 
        endif    
c .....................................................................
c
c ...
        i_z = alloc_8('zsolver ',1,neq)
        i_h = alloc_8('hsolver ',1,neq)
        i_r = alloc_8('rsolver ',1,neq)
        i_s = alloc_8('psolver ',1,neq)
c .....................................................................
c
c ... calculo do precondicionador
        call cal_precond(ip     ,ja    ,m
     1                  ,ad     ,al    ,ia(i_z)
     2                  ,precond,neq   ,nequ
     3                  ,0      ,0     ,0      
     4                  ,0     ,0      ,0      
     5                  ,.false.,0)
c .....................................................................
c
c ...
        if(omp_solv)then
          call call_sqmr_omp(neq  ,nequ  ,nad   ,ip      ,ja
     1                    ,ad     ,al    ,m     ,b       ,x    
     2                    ,ia(i_z),ia(i_h),ia(i_r) ,ia(i_s)     
     3                    ,tol    ,maxit ,precond,iparam ,fhist_solv
     4                    ,fprint
     5                    ,0      ,0     ,0      ,neq     ,0      
     6                    ,0      ,0     ,0      ,ia(i_threads_y)
     7                    ,1      ,.false.,.false.,fsqmr) 
c .....................................................................
c
c ...
        else
          call call_sqmr(neq    ,nequ  ,nad     ,ip      ,ja
     1                ,ad       ,al    ,m       ,b       ,x    
     2                ,ia(i_z),ia(i_h),ia(i_r)  ,ia(i_s) 
     3                ,tol      ,maxit ,precond ,iparam  ,fhist_solv
     4                ,fprint
     5                ,0        ,0     ,0       ,neq     ,0      
     6                ,0        ,0     ,0      
     7                ,1        ,.false.,.false.,fsqmr) 
       endif
c .....................................................................
c
c ...
        i_s = dealloc('psolver ')
        i_r = dealloc('rsolver ')
        i_h = dealloc('hsolver ')
        i_z = dealloc('zsolver ')
c ......................................................................
c
c ...
        if(fsqmr) goto 100
c ......................................................................
      endif
c ......................................................................           
c
c ...                                                                     
      if (omp_solv) then
        pmatrixtime = get_time() - pmatrixtime 
        i_threads_y = dealloc('buffer_y')
        pmatrixtime = get_time() - pmatrixtime
      endif
c ......................................................................         
c
c ...   
      return
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 11/04/2016                                    *
c * Data de modificaco : 21/05/2020                                    * 
c * ------------------------------------------------------------------ *   
c * CALL_CG : chama a versao do gradiente conjudado desejada           *    
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq      - numero de equacoes                                      *
c * nequ     - numero de equacoes no bloco Kuu                         *
c * nad      - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * ad(neq)  - diagonal da matriz A                                    *
c * al(*)    - parte triangular inferior de A                          *
c * b(neq)   - vetor de forcas                                         *
c * m(*)     - precondicionador                                        *
c * x(neq)   - chute inicial                                           *
c * z(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * s(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * pc       - precondicionador                                        *
c *            1 - nenhum                                              *
c *            2 - diaggonal                                           *
c *            3 - iLDLt(0)                                            *
c *            4 - iLLt(0)                                             *
c *            5 - modulo da diagonal                                  *
c *            6 - bloco diagonal                                      *
c *            7 - bloco diagonal com proximacao de schur              *
c * iparam   - parametros do bloco diagonal                            *
c *          - iparam(1) - numero de sub matriz em blocos              *
c *          - iparam(2) - numero de inversos da diagonal simples      *
c *          - iparam(3) - numero de termos nos bloco                  *
c *          - iparam(4) - tamanho do bloco                            *
c * fhist_log- log do residuo por iteracao                             *
c * my_id    - MPI                                                     *
c * neqf1i   - MPI                                                     *
c * neqf2i   - MPI                                                     *
c * neq_doti - MPI                                                     *
c * i_fmap   - MPI                                                     *
c * i_xfi    - MPI                                                     *
c * i_rvcs   - MPI                                                     *
c * i_dspli  - MPI                                                     *
c * mpi      - true|false                                              *
c * ovlp     - overllaping                                             *
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
c **********************************************************************  
      subroutine call_cg(neq      ,nequ  ,nad   ,ia      ,ja
     1                  ,ad       ,al    ,m     ,b       ,x    
     2                  ,z        ,r     ,s    
     3                  ,tol      ,maxit ,precond,iparam ,fhist_log
     4                  ,my_id    ,neqf1i,neqf2i ,neq_doti,i_fmapi
     5                  ,i_xfi    ,i_rcvsi,i_dspli
     6                  ,nprcs    ,ovlp   ,mpi    )
      implicit none
      include 'time.fi'
c ... mpi
      integer my_id
      integer neqf1i,neqf2i,nprcs
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
      logical ovlp,mpi
c .....................................................................
      integer neq,nequ,neq_doti 
      integer*8 ia(*),nad
      integer ja(*)
      real*8  ad(*),al(*),x(*),b(*)
c ... arranjos auxiliares
      real*8 z(*),r(*),s(*)
c ...
      real*8  tol
      integer maxit  
      logical fhist_log
c ... precondicionador
      integer precond,iparam(*)
      real*8 m(*)
      real*8 dot_par
c ...
      external dot_par
      external matvec_csrc_sym_pm,matvec_csrcr_sym_pm  
      external ildlt_solv,illt_solv  
c ......................................................................

c ... cg
      if(precond .eq. 1) then
c ...  
        call cg(neq    ,nequ   ,nad,ia   ,ja
     1         ,ad     ,al     ,al ,b    ,x
     2         ,z      ,r      ,s  ,tol  ,maxit
c ... matvec comum:
     3         ,matvec_csrc_sym_pm,dot_par 
     4         ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     5         ,i_xfi ,i_rcvsi,i_dspli
     6         ,.true.,.true. ,fhist_log ,.true.)
c .....................................................................
c
c ... pcg - cg com precondicionador diagonal
      else if(precond .eq. 2 .or. precond .eq. 5 ) then
c ... overllaping
         if(ovlp) then  
           call pcg(neq    ,nequ   ,nad,ia   ,ja
     1             ,ad     ,al     ,al ,m    ,b
     2             ,x      ,z      ,r  ,s
     3             ,tol    ,maxit
c ... matvec comum:
     4             ,matvec_csrcr_sym_pm,dot_par 
c
     5             ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6             ,i_xfi ,i_rcvsi,i_dspli
     7             ,.true.,.true. ,fhist_log,.true.
     8             ,nprcs ,mpi)
c .....................................................................
c
c ...non-overllaping
         else 
           call pcg(neq    ,nequ   ,nad,ia   ,ja
     1             ,ad     ,al     ,al ,m    ,b
     2             ,x      ,z      ,r  ,s
     3             ,tol    ,maxit
c ... matvec comum:
     4             ,matvec_csrc_sym_pm,dot_par 
c
     5             ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6             ,i_xfi ,i_rcvsi,i_dspli
     7             ,.true.,.true. ,fhist_log,.true.
     8             ,nprcs ,mpi)
         endif
c .....................................................................
c
c ... iccg - cg com precondicionador LDLT(0) imcompleto
      elseif(precond .eq. 3 ) then
        call iccg(neq      ,nequ  ,nad   ,ia      ,ja
     1           ,ad       ,al    ,al    ,m       ,b    
     2           ,x        ,z     ,r     ,s   
     3           ,tol      ,maxit
c ... matvec comum:
     4           ,matvec_csrc_sym_pm,dot_par,ildlt_solv 
     5           ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6           ,i_xfi ,i_rcvsi,i_dspli
     7           ,.true.,.true. ,fhist_log,.true.)
c .....................................................................
c
c ... iccg - cg com precondicionador LLT(0) imcompleto
      elseif(precond .eq. 4 ) then
        call iccg(neq      ,nequ  ,nad   ,ia      ,ja
     1           ,ad       ,al    ,al    ,m       ,b    
     2           ,x        ,z     ,r     ,s   
     3           ,tol      ,maxit
c ... matvec comum:
     4           ,matvec_csrc_sym_pm,dot_par,illt_solv
     5           ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6           ,i_xfi ,i_rcvsi,i_dspli
     7           ,.true.,.true. ,fhist_log,.true.)
c .....................................................................
c
c ... bpcg - cg com bloco diagonal 
      elseif(precond .eq. 6 ) then
        call bpcg(neq      ,nequ  ,nad   ,ia      ,ja
     1           ,ad       ,al    ,al    ,m       ,b    
     2           ,x        ,z     ,r     ,s   
     3           ,tol      ,maxit 
c ... matvec comum:
     4           ,matvec_csrc_sym_pm,dot_par
     5           ,my_id ,neqf1i ,neqf2i ,neq_doti,i_fmapi
     6           ,i_xfi ,i_rcvsi,i_dspli
     7           ,.true.,.true. ,fhist_log ,.true.)
      endif  
c .....................................................................
      return
      end    
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 23/09/2016                                    *
c * Data de modificaco : 15/12/2016                                    * 
c * ------------------------------------------------------------------ *   
c * CALL_CG_OMP : chama a versao do gradiente conjudado desejada       *    
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq      - numero de equacoes                                      *
c * nequ     - numero de equacoes no bloco Kuu                         *
c * nad      - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * ad(neq)  - diagonal da matriz A                                    *
c * al(*)    - parte triangular inferior de A                          *
c * b(neq)   - vetor de forcas                                         *
c * m(*)     - precondicionador                                        *
c * x(neq)   - chute inicial                                           *
c * z(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * s(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * precond  - precondicionador                                        *
c *            1 - nenhum                                              *
c *            2 - diaggonal                                           *
c *            3 - iLDLt(0)                                            *
c *            4 - iLLt(0)                                             *
c *            5 - modulo da diagonal                                  *
c *            6 - bloco diagonal                                      *
c *            7 - bloco diagonal com proximacao de schur              *
c * iparam   - parametros do bloco diagonal                            *
c *          - iparam(1) - numero de sub matriz em blocos              *
c *          - iparam(2) - numero de inversos da diagonal simples      *
c *          - iparam(3) - numero de termos nos bloco                  *
c *          - iparam(4) - tamanho do bloco                            *
c * fhist_log- log do residuo por iteracao                             *
c * my_id    - MPI                                                     *
c * neqf1i   - MPI                                                     *
c * neqf2i   - MPI                                                     *
c * neq_doti - MPI                                                     *
c * i_fmap   - MPI                                                     *
c * i_xfi    - MPI                                                     *
c * i_rvcs   - MPI                                                     *
c * i_dspli  - MPI                                                     *
c * thread_y - buffer de equacoes para o vetor y (openmp)              *  
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
      subroutine call_cg_omp(neq      ,nequ  ,nad   ,ia      ,ja
     1                  ,ad       ,al    ,m     ,b       ,x    
     2                  ,z        ,r     ,s    
     3                  ,tol      ,maxit ,precond,iparam ,fhist_log
     4                  ,my_id    ,neqf1i,neqf2i ,neq_doti,i_fmapi
     5                  ,i_xfi    ,i_rcvsi,i_dspli,thread_y
     6                  ,nprcs    ,ovlp   ,mpi    )
      implicit none
      include 'time.fi'
c ... mpi
      integer my_id
      integer neqf1i,neqf2i,nprcs
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
      logical ovlp,mpi
c .....................................................................
      integer neq,nequ,neq_doti 
      integer*8 ia(*),nad
      integer ja(*)
      real*8  ad(*),al(*),x(*),b(*)
c ... arranjos auxiliares
      real*8 z(*),r(*),s(*)
c ... buffer do CSR omp
      real*8 thread_y(*)
c ...
      real*8  tol
      integer maxit  
      logical fhist_log
c ... precondicionador
      integer precond,iparam(*)
      real*8 m(*)
c ...
      external dot_par_omp
      external matvec_csrc_sym_pm_omp,matvec_csrcr_sym_pm_omp  
      external ildlt_solv,illt_solv  
c ......................................................................

c ... cg
      if(precond .eq. 1) then
        print*,"CG: Openmp not available !!"  
c
c ... pcg - cg com precondicionador diagonal
      else if(precond .eq. 2 .or. precond .eq. 5 ) then
c ... overllaping
        if(ovlp) then  
          call pcg_omp(neq    ,nequ   ,nad,ia   ,ja
     1                ,ad     ,al     ,al ,m    ,b 
     2                ,x      ,z      ,r  ,s
     3                ,tol    ,maxit
c ... matvec comum:
     4                ,matvec_csrcr_sym_pm_omp,dot_par_omp
c             
     5                ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6                ,i_xfi ,i_rcvsi,i_dspli,thread_y
     7                ,.true.,.true. ,fhist_log,.true.
     8                ,nprcs ,mpi)
c .....................................................................
c
c ...non-overllaping
        else 
          call pcg_omp(neq    ,nequ   ,nad,ia   ,ja
     1                ,ad     ,al     ,al ,m    ,b 
     2                ,x      ,z      ,r  ,s
     3                ,tol    ,maxit
c ... matvec comum:
     4                ,matvec_csrc_sym_pm_omp,dot_par_omp
c             
     5                ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6                ,i_xfi ,i_rcvsi,i_dspli,thread_y
     7                ,.true.,.true. ,fhist_log,.true.
     8                ,nprcs ,mpi)
        endif
c .....................................................................
c
c ... iccg - cg com precondicionador LDLT(0) imcompleto
      elseif(precond .eq. 3 ) then
        print*,"ICCG: Openmp not available !!" 
c .....................................................................
c
c ... iccg - cg com precondicionador LLT(0) imcompleto
      elseif(precond .eq. 4 ) then
        print*,"ICCG: Openmp not available !!" 
c .....................................................................
c
c ... bpcg - cg com bloco diagonal 
      elseif(precond .eq. 6 ) then
        print*,"BCCG: Openmp not available !!" 
      endif  
c .....................................................................
      return
      end    
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 28/06/2016                                    *
c * Data de modificaco : 02/02/2017                                    * 
c * ------------------------------------------------------------------ *   
c * CALL_SQMR:chama a versao do metodo QMR simetrico                   *    
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq      - numero de equacoes                                      *
c * nequ     - numero de equacoes no bloco Kuu                         *
c * nad      - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * ad(neq)  - diagonal da matriz A                                    *
c * al(*)    - parte triangular inferior de A                          *
c * b(neq)   - vetor de forcas                                         *
c * m(*)     - precondicionador                                        *
c * x(neq)   - chute inicial                                           *
c * c(neq)   - arranjo local de trabalho                               *
c * h(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * s(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * pc       - precondicionador                                        *
c *            1 - nenhum                                              *
c *            2 - diaggonal                                           *
c *            3 - iLDLt(0)                                            *
c *            4 - iLLt(0)                                             *
c *            5 - modulo da diagonal                                  *
c *            6 - bloco diagonal                                      *
c *            7 - bloco diagonal com proximacao de schur              *
c * iparam   - parametros do bloco diagonal                            *
c *          - iparam(1) - numero de sub matriz em blocos              *
c *          - iparam(2) - numero de inversos da diagonal simples      *
c *          - iparam(3) - numero de termos nos bloco                  *
c *          - iparam(4) - tamanho do bloco                            *
c * fhist_log- log do residuo por iteracao                             *
c * fprint   -                                                         *
c * my_id    - MPI                                                     *
c * neqf1i   - MPI                                                     *
c * neqf2i   - MPI                                                     *
c * neq_doti - MPI                                                     *
c * i_fmap   - MPI                                                     *
c * i_xfi    - MPI                                                     *
c * i_rvcs   - MPI                                                     *
c * i_dspli  - MPI                                                     *
c * mpi      - true|false                                              *
c * ovlp     - overllaping                                             *
c * nprcs    - numero de processos mpi                                 *
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
c **********************************************************************  
      subroutine call_sqmr(neq      ,nequ  ,nad      ,ia      ,ja
     1                    ,ad       ,al    ,m        ,b       ,x    
     2                    ,c        ,h     ,r        ,s    
     3                    ,tol      ,maxit ,pc       ,iparam ,fhist_log
     4                    ,fprint
     5                    ,my_id    ,neqf1i,neqf2i   ,neq_doti,i_fmapi
     6                    ,i_xfi    ,i_rcvsi,i_dspli
     7                    ,nprcs    ,ovlp   ,mpi     ,fsqmr)
      implicit none
c ... mpi
      integer my_id
      integer neqf1i,neqf2i,nprcs
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
      logical ovlp,mpi
c .....................................................................
      integer neq,nequ,neq_doti
      integer*8 ia(*),nad 
      integer ja(*)
      real*8  ad(*),al(*),x(*),b(*)
c ... arranjos auxiliares
      real*8 c(*),h(*),r(*),s(*)
c ...
      real*8  tol
      integer maxit  
      logical fhist_log,fprint,fsqmr
c ... precondicionador
      integer pc,iparam(*)
      real*8 m(*)
c ...
      external dot_par
      external matvec_csrc_sym_pm,matvec_csrcr_sym_pm  
      external ildlt_solv,illt_solv  
c ......................................................................
c
      if(pc .eq. 1 ) then
        call sqmr(neq    ,nequ   ,nad    ,ia    ,ja
     1           ,ad     ,al     ,al     ,m     ,b   ,x 
     2           ,c      ,h      ,r      ,s 
     3           ,tol   ,maxit
     4           ,matvec_csrc_sym_pm,dot_par
     5           ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6           ,i_xfi ,i_rcvsi,i_dspli
     7           ,.true.,.true. ,fhist_log,.true.)
c .....................................................................
c
c ...
       else if(pc .eq. 2 .or. pc .eq. 5 .or.  pc .eq. 7) then
c ... overllaping
         if(ovlp) then
           call rpsqmr(neq    ,nequ   ,nad    ,ia    ,ja
     1                ,ad     ,al     ,al      ,m    ,b   ,x 
     2                ,c      ,h      ,r       ,s 
     3                ,tol   ,maxit
     4                ,matvec_csrcr_sym_pm,dot_par
     5                ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6                ,i_xfi ,i_rcvsi,i_dspli
     7                ,fprint,.true. ,fhist_log,.true.
     8                ,nprcs ,mpi,fsqmr)
c .....................................................................
c
c ...non-overllaping
         else 
           call rpsqmr(neq    ,nequ   ,nad    ,ia    ,ja
     1                ,ad     ,al     ,al      ,m     ,b   ,x 
     2                ,c      ,h      ,r       ,s 
     3                ,tol   ,maxit
     4                ,matvec_csrc_sym_pm,dot_par
     5                ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6                ,i_xfi ,i_rcvsi,i_dspli
     7                ,fprint,.true. ,fhist_log,.true.
     8                ,nprcs ,mpi,fsqmr)
          endif
c .....................................................................
        endif
c .....................................................................
      return
      end    
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 22/09/2016                                    *
c * Data de modificaco : 02/02/2017                                    * 
c * ------------------------------------------------------------------ *   
c * CALL_SQMR_OMP:chama a versao do metodo QMR simetrico               *    
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq      - numero de equacoes                                      *
c * nequ     - numero de equacoes no bloco Kuu                         *
c * nad      - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * ad(neq)  - diagonal da matriz A                                    *
c * al(*)    - parte triangular inferior de A                          *
c * b(neq)   - vetor de forcas                                         *
c * m(*)     - precondicionador                                        *
c * x(neq)   - chute inicial                                           *
c * c(neq)   - arranjo local de trabalho                               *
c * h(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * s(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * pc       - precondicionador                                        *
c *            1 - nenhum                                              *
c *            2 - diaggonal                                           *
c *            3 - iLDLt(0)                                            *
c *            4 - iLLt(0)                                             *
c *            5 - modulo da diagonal                                  *
c *            6 - bloco diagonal                                      *
c *            7 - bloco diagonal com proximacao de schur              *
c * iparam   - parametros do bloco diagonal                            *
c *          - iparam(1) - numero de sub matriz em blocos              *
c *          - iparam(2) - numero de inversos da diagonal simples      *
c *          - iparam(3) - numero de termos nos bloco                  *
c *          - iparam(4) - tamanho do bloco                            *
c * fhist_log- log do residuo por iteracao                             *
c * fprint   -                                                         *
c * my_id    -                                                         *
c * neqf1i   -                                                         *
c * neqf2i   -                                                         *
c * neq_doti -                                                         *
c * i_fmap   -                                                         *
c * i_xfi    -                                                         *
c * i_rvcs   -                                                         *
c * i_dspli  -                                                         *
c * thread_y - buffer de equacoes para o vetor y (openmp)              * 
c * mpi      - true|false                                              *
c * ovlp     - overllaping                                             *
c * nprcs    - numero de processos mpi                                 *
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
c **********************************************************************  
      subroutine call_sqmr_omp(neq      ,nequ  ,nad   ,ia      ,ja
     1                    ,ad       ,al    ,m     ,b       ,x    
     2                    ,c        ,h     ,r     ,s    
     3                    ,tol      ,maxit ,pc    ,iparam ,fhist_log
     4                    ,fprint
     5                    ,my_id    ,neqf1i,neqf2i ,neq_doti,i_fmapi
     6                    ,i_xfi    ,i_rcvsi,i_dspli,thread_y
     7                    ,nprcs    ,ovlp   ,mpi    ,fsqmr )
      implicit none
c ... mpi
      integer my_id
      integer neqf1i,neqf2i,nprcs
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
      logical ovlp,mpi
c .....................................................................
      integer neq,nequ,neq_doti 
      integer*8 ia(*),nad 
      integer ja(*)
      real*8  ad(*),al(*),x(*),b(*)
c ... arranjos auxiliares
      real*8 c(*),h(*),r(*),s(*)
c ... buffer do CSR omp
      real*8 thread_y(*)
c ...
      real*8  tol
      integer maxit  
      logical fhist_log,fprint,fsqmr
c ... precondicionador
      integer pc,iparam(*)
      real*8 m(*)
c ...
      external dot_par_omp
      external matvec_csrc_sym_pm_omp,matvec_csrcr_sym_pm_omp
      external ildlt_solv,illt_solv  
c ......................................................................
c
      if(pc .eq. 1 ) then
        print*,"SQMR: Openmp not available !!"  
c .....................................................................
c
c ...
      else if(pc .eq. 2 .or. pc .eq. 5 .or.  pc .eq. 7) then
c ... overllaping
        if(ovlp) then
          call rpsqmr_omp(neq    ,nequ   ,nad    ,ia    ,ja
     1             ,ad     ,al     ,al   ,m     ,b   ,x 
     2             ,c      ,h      ,r    ,s 
     3             ,tol   ,maxit
     4             ,matvec_csrcr_sym_pm_omp,dot_par_omp
     5             ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6             ,i_xfi ,i_rcvsi,i_dspli,thread_y
     7             ,fprint,.true. ,fhist_log,.true.
     8             ,nprcs ,mpi,fsqmr)
c .....................................................................
c
c ...non-overllaping
        else 
          call rpsqmr_omp(neq    ,nequ   ,nad    ,ia    ,ja
     1             ,ad     ,al     ,al      ,m     ,b   ,x 
     2             ,c      ,h      ,r       ,s 
     3             ,tol   ,maxit
     4             ,matvec_csrc_sym_pm_omp,dot_par_omp
     5             ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     6             ,i_xfi ,i_rcvsi,i_dspli,thread_y
     7             ,fprint,.true. ,fhist_log,.true.
     8             ,nprcs ,mpi,fsqmr)
        endif
c .....................................................................
      endif
c .....................................................................
      return
      end    
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 28/11/2016                                    *
c * Data de modificaco : 29/01/2017                                    * 
c * ------------------------------------------------------------------ *   
c * CALL_GMRES :chama a versao do GMRES desejada                       *    
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq      - numero de equacoes                                      *
c * neqovlp  - numero de eqquacoes em overlapping
c * nad      - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * ad(neq)  - diagonal da matriz A                                    *
c * al(*)    - parte triangular inferior de A                          *
c * m(*)     - precondicionador                                        *
c * b(*)     - vetor de forcas                                         *
c * x(*)     - chute inicial                                           *
c * g(*)     - arranjo local de trabalho                               *
c * h(*)     - arranjo local de trabalho                               *
c * y(*)     - arranjo local de trabalho                               *
c * v(*)     - arranjo local de trabalho                               *
c * s(*)     - arranjo local de trabalho                               *
c * r(*)     - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * precond  - precondicionador                                        *
c *            1 - nenhum                                              *
c *            2 - diag                                                *
c *            3 - iLDLt                                               *
c *            4 -                                                     *
c *            5 - modulo da diagonal                                  *
c * fhist_log- log do residuo por iteracao                             *
c * fprint   - saida na tela                                           *
c * nkrylov  - numero de bases krylov                                  *
c * my_id    -                                                         *
c * neqf1i   -                                                         *
c * neqf2i   -                                                         *
c * neq_doti -                                                         *
c * i_fmap   -                                                         *
c * i_xfi    -                                                         *
c * i_rvcs   -                                                         *
c * i_dspli  -                                                         *
c * ovlp     - overllaping                                             *
c * nprcs    - numero de processos mpi                                 *
c * mpi      - true|false                                              *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq) - modificado                                                *
c * ad(*),al(*) - inalterados                                          *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************  
      subroutine call_gmres(neq    ,neqovlp,nad  
     1                     ,ia     ,ja     ,ad     ,au       ,al 
     2                     ,m      ,b      ,x      ,g
     3                     ,h      ,y      ,c      ,s        ,r     
     4                     ,tol    ,maxit  ,precond,fhist_log,fprint
     5                     ,nkrylov
     6                     ,my_id  ,neqf1i ,neqf2i,neq_doti,i_fmapi
     7                     ,i_xfi  ,i_rcvsi,i_dspli
     8                     ,nprcs  ,ovlp   ,mpi) 
      implicit none
      include 'time.fi'
c ... mpi
      integer my_id
      integer neqf1i,neqf2i,nprcs
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
      logical ovlp,mpi
c .....................................................................
      integer neq,neqovlp,neq_doti,nkrylov 
      integer*8 ia(*),nad
      integer ja(*)
      real*8  ad(*),au(*),al(*),x(*),b(*)
c ... arranjos auxiliares
      real*8 g(*),h(*),y(*),c(*),r(*),s(*)
c ...
      real*8  tol
      integer maxit  
      logical fhist_log,fprint
c ... precondicionador
      integer precond
      real*8  m(*)
c ...
      external dot_par
      external matvec_csrc_pm,matvec_csrcr_sym_pm        
c ......................................................................
c
c ... gmres sem precondicionador
      if(precond .eq. 1) then
        print*,'Not implemented GMRES algorithm without',
     .        ' preconditioner !!!'
        call stop_mef()
c ......................................................................
c
c ... pbicgstabl2 - bicgstabl2 com precondicionador diagonal
      else if(precond .eq. 2 .or. precond .eq. 5 ) then
c ... overllaping
         if(ovlp) then
           call gmres2(neq,neq,nad,ia,ja 
     1                ,ad ,au ,al ,m ,b ,x,nkrylov
     2                ,g  ,h  ,y  ,c ,s ,r       
     3                ,tol    ,maxit   
c ... matvec comum:
     4                ,matvec_csrcr_sym_pm,dot_par 
     5                ,neqovlp
     6                ,my_id  ,neqf1i ,neqf2i ,neq_doti,i_fmapi
     7                ,i_xfi  ,i_rcvsi,i_dspli 
     8                ,fprint ,.true. ,fhist_log,.true.
     8                ,nprcs  , mpi) 
c .....................................................................
c
c ...non-overllaping
         else 
           call gmres2(neq,neq,nad,ia,ja 
     1                ,ad ,au ,al ,m ,b ,x,nkrylov
     2                ,g  ,h  ,y  ,c ,s ,r       
     3                ,tol    ,maxit   
c ... matvec comum:
     4                ,matvec_csrc_pm ,dot_par 
     5                ,neqovlp
     6                ,my_id  ,neqf1i ,neqf2i ,neq_doti,i_fmapi
     7                ,i_xfi  ,i_rcvsi,i_dspli 
     8                ,fprint ,.true. ,fhist_log,.true.
     8                ,nprcs  , mpi) 
         endif
c .....................................................................
c
c ...
      elseif(precond .eq. 3 ) then
        print*,'iLDLt not implemented for GMRES!!!'
        stop
c ...
      elseif(precond .eq. 4 ) then
        print*,'iLLt not implemented for GMRES!!!'
        stop
      endif  
c .....................................................................
      return
      end    
c *********************************************************************
c
c **********************************************************************
c * Data de criacao    : 29/12/2016                                    *
c * Data de modificaco : 29/01/2017                                    * 
c * ------------------------------------------------------------------ *   
c * CALL_GMRES_OMP : chama a versao do GMRES OMP desejada              *    
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq      - numero de equacoes                                      *
c * neqovlp  - numero de eqquacoes em overlapping                      * 
c * nad      - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * ad(neq)  - diagonal da matriz A                                    *
c * al(*)    - parte triangular inferior de A                          *
c * m(*)     - precondicionador                                        *
c * b(*)     - vetor de forcas                                         *
c * x(*)     - chute inicial                                           *
c * g(*)     - arranjo local de trabalho                               *
c * h(*)     - arranjo local de trabalho                               *
c * y(*)     - arranjo local de trabalho                               *
c * v(*)     - arranjo local de trabalho                               *
c * s(*)     - arranjo local de trabalho                               *
c * r(*)     - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * precond  - precondicionador                                        *
c *            1 - nenhum                                              *
c *            2 - diag                                                *
c *            3 - iLDLt                                               *
c *            4 -                                                     *
c *            5 - modulo da diagonal                                  *
c * fhist_log- log do residuo por iteracao                             *
c * fprint   - saida na tela                                           *
c * nkrylov  - numero de bases krylov                                  *
c * my_id    -                                                         *
c * neqf1i   -                                                         *
c * neqf2i   -                                                         *
c * neq_doti -                                                         *
c * i_fmap   -                                                         *
c * i_xfi    -                                                         *
c * i_rvcs   -                                                         *
c * i_dspli  -                                                         *
c * thread_y - buffer de equacoes do omp                               *
c * ovlp     - overllaping                                             *
c * nprcs    - numero de processos mpi                                 *
c * mpi      - true|false                                              *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq) - modificado                                                *
c * ad(*),al(*) - inalterados                                          *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************  
      subroutine call_gmres_omp(neq    ,neqovlp,nad  
     1                     ,ia     ,ja     ,ad     ,al 
     2                     ,m      ,b      ,x      ,g
     3                     ,h      ,y      ,c      ,s        ,r     
     4                     ,tol    ,maxit  ,precond,fhist_log,fprint
     5                     ,nkrylov
     6                     ,my_id  ,neqf1i ,neqf2i,neq_doti,i_fmapi
     7                     ,i_xfi  ,i_rcvsi,i_dspli,thread_y
     8                     ,nprcs  ,ovlp   ,mpi) 
      implicit none
      include 'time.fi'
c ... mpi
      integer my_id
      integer neqf1i,neqf2i,nprcs
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
      logical ovlp,mpi
c .....................................................................
      integer neq,neqovlp,neq_doti,nkrylov 
      integer*8 ia(*),nad
      integer ja(*)
      real*8  ad(*),al(*),x(*),b(*)
c ... arranjos auxiliares
      real*8 g(*),h(*),y(*),c(*),r(*),s(*)
c ... buffer do CSR omp
      real*8 thread_y(*)
c ...
      real*8  tol
      integer maxit  
      logical fhist_log,fprint
c ... precondicionador
      integer precond
      real*8  m(*)
c ...
      external dot_par_omp
      external matvec_csrc_sym_pm_omp,matvec_csrcr_sym_pm_omp        
c ......................................................................
c
c ... gmres sem precondicionador
      if(precond .eq. 1) then
        print*,'Not implemented GMRES algorithm without',
     .        ' preconditioner !!!'
        call stop_mef()
c ......................................................................
c
c ... pbicgstabl2 - bicgstabl2 com precondicionador diagonal
      else if(precond .eq. 2 .or. precond .eq. 5 ) then
c ... overllaping
         if(ovlp) then
           call gmres2_omp(neq,neq,nad,ia,ja 
     1                ,ad ,al ,al ,m ,b ,x,nkrylov
     2                ,g  ,h  ,y  ,c ,s ,r       
     3                ,tol    ,maxit   
c ... matvec comum:
     4                ,matvec_csrcr_sym_pm_omp,dot_par_omp 
     5                ,neqovlp
     6                ,my_id  ,neqf1i ,neqf2i ,neq_doti,i_fmapi
     7                ,i_xfi  ,i_rcvsi,i_dspli,thread_y 
     8                ,fprint ,.true. ,fhist_log,.true.
     8                ,nprcs  , mpi) 
c .....................................................................
c
c ...non-overllaping
         else 
           call gmres2_omp(neq,neq,nad,ia,ja 
     1                ,ad ,al ,al ,m ,b ,x,nkrylov
     2                ,g  ,h  ,y  ,c ,s ,r       
     3                ,tol    ,maxit   
c ... matvec comum:
     4                ,matvec_csrc_sym_pm_omp,dot_par_omp 
     5                ,neqovlp
     6                ,my_id  ,neqf1i ,neqf2i ,neq_doti,i_fmapi
     7                ,i_xfi  ,i_rcvsi,i_dspli,thread_y 
     8                ,fprint ,.true. ,fhist_log,.true.
     8                ,nprcs  , mpi) 
         endif
c .....................................................................
c
c ...
      elseif(precond .eq. 3 ) then
        print*,'iLDLt not implemented for GMRES!!!'
        stop
c ...
      elseif(precond .eq. 4 ) then
        print*,'iLLt not implemented for GMRES!!!'
        stop
      endif  
c .....................................................................
      return
      end    
c *********************************************************************
c
c *********************************************************************
      subroutine get_res(u,x,id,fnno,nnode,ndf)
      implicit none
      integer nnode,ndf
      integer id(ndf,*),fnno(*),i,j,k
      real*8 u(ndf,*),x(*)
c ... loop nos
      do i = 1, nnode
        do j = 1, ndf - 1
          k    = id(j,i)
          if( k .gt. 0 ) then
            x(k) = u(j,i)
          endif
        enddo
        if(fnno(i) .eq. 1) then
          k    = id(ndf,i)
          if( k .gt. 0 ) then
            x(k) = u(ndf,i)
          endif
        endif
      enddo
c .....................................................................
c
c ...
      return
      end
c *************************************************************************
c
c **********************************************************************
c * Data de criacao    : 18/04/2016                                    *
c * Data de modificaco : 07/11/2018                                    * 
c * ------------------------------------------------------------------ *  
c * SET_PRECOND : escolhe o precondicionandor                          *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * macro   - precondicionador escolhido                               *
c * solver  - nao definido                                             *
c * nin     - aqruivo de entrada                                       *
c * my_id   - id do processo do mpi                                    *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * solver  - solver escolhido                                         *
c *         1 - CG                                                     *
c *         2 - GMRES                                                  *
c *         3 -                                                        *
c *         4 - BICGSTAB                                               *
c *         5 -                                                        *
c *         6 - BICGSTABL2                                             *
c *         7 - MINRES                                                 *
c *         8 - PCR                                                    *
c *         9 - SYMMLQ                                                 *
c *        10 - pardiso                                                *
c *        11 - SQRM                                                   *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c ********************************************************************** 
      subroutine set_solver(solver,nin,my_id)
      implicit none
      include 'string.fi'
      include 'precond.fi'
      character*8 macros(12),string
      integer solver,nin,my_id
      integer i,nmc 
      data macros/'cg      ','sqmr    ','symmlq  '
     .           ,'cr      ','minres  ','bicgsl2 '
     .           ,'bicgs   ','gmres   ','pardiso '
     .           ,'block_it','        ','        '/
      data nmc /12/
c ...
      write(string,'(8a)') (word(i),i=1,8)
c ... CG
      if( string .eq. macros(1)) then
        solver = 1
        if(my_id.eq.0) then
          write(*,'(a10,1x,a8)')'Solver :',adjustr(macros(1))
        endif
c .....................................................................
c
c ... SQRM
      else if( string .eq. macros(2)) then
        solver = 11
        if(my_id.eq.0) then
          write(*,'(a10,1x,a8)')'Solver :',adjustr(macros(2))
        endif
c .....................................................................
c
c ... symmlq
      else if( string .eq. macros(3)) then
        solver =  9
        if(my_id.eq.0) then
          write(*,'(a10,1x,a8)')'Solver :',adjustr(macros(3))
        endif
c .....................................................................
c
c ... cr     
      else if( string .eq. macros(4)) then
        solver =  8
        if(my_id.eq.0) then
          write(*,'(a10,1x,a8)')'Solver :',adjustr(macros(4))
        endif
c .....................................................................
c
c ... minres
      else if( string .eq. macros(5)) then
        solver =  7
        if(my_id.eq.0) then
          write(*,'(a10,1x,a8)')'Solver :',adjustr(macros(5))
        endif
c .....................................................................
c
c ... BICGSTABL2
      else if( string .eq. macros(6)) then
        solver =  6
        if(my_id.eq.0) then
          write(*,'(a10,1x,a8)')'Solver :',adjustr(macros(6))
        endif
c .....................................................................
c
c ... BICGSTAB
      else if( string .eq. macros(7)) then
        solver =  4
        if(my_id.eq.0) then
          write(*,'(a10,1x,a8)')'Solver :',adjustr(macros(7))
        endif
c .....................................................................
c
c ... GMRES   
      else if( string .eq. macros(8)) then
        solver =  2
        if(my_id.eq.0) then
          write(*,'(a10,1x,a8)')'Solver :', adjustr(macros(8))
        endif
c .....................................................................
c
c ... PARDISO 
      else if( string .eq. macros(9)) then
        solver =  10
        if(my_id.eq.0) then
          write(*,'(a10,1x,a8)')'Solver :', adjustr(macros(9))
        endif
c .....................................................................
c
c ... BLOK_IT_PCG
      else if( string .eq. macros(10)) then
        solver =  5
        if(my_id.eq.0) then
          write(*,'(a10,1x,a8)')'Solver :', adjustr(macros(9))
        endif
c .....................................................................
c
c ...                         
      else
        do i = 1, nmc
          if(my_id.eq.0) then
            print*,'Error reading macro (SOLVER) !'
            print*,'Solver available:'
            print*,'Solver: ',macros(i)
          endif
        enddo
        call stop_mef()
      endif 
c .....................................................................
c
c ...
      return
      end
c ********************************************************************** 