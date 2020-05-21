c **********************************************************************
c * Data de criacao    : 30/04/2016                                    *
c * Data de modificaco : 20/05/2020                                    * 
c * ------------------------------------------------------------------ *   
c * CALL_PARDISO : chama o sover pardiso mkl                           *    
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq      - numero de equacoes                                      *
c * nnz      - numero total de nao zeros                               *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * a(*)     - matriz A                                                *
c * b(neq)   - vetor de forcas                                         *
c * x(neq)   -                                                         *
c * mtype    - tipo da matriz                                          *
c *           -2 -> simetrico indefinido                               *
c *            2 -> simetrico definido positivo                        *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq)      - vetor solucao                                        *
c * a(*),b(neq) - inalterados                                          *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************  
      subroutine call_mkl_pardiso(neq,nnz,ia8,ja,a,b,x,z,ia4,mtype)
      implicit none
      include 'time.fi'
c ......................................................................
      integer*8 ia8(*),nnz
      integer ja(*),neq,ia4(*)
      real*8 a(*),b(*),x(*),z(*)
      real*8 norm,xkx,mem
c ... variavel interna do mkl( 64 btis)
      integer*8 pt(64) 
c ...
      integer iparm(64),msglvl,error
c ...
      integer maxfct,mnum,mtype,phase,nrhs
      integer idum
      real*8 ddum
c ...
      real*8 dot
c ...
      call ia8_to_ia4(ia4,ia8,neq)
c ...
      error  = 0
      maxfct = 1
      mnum   = 1
      nrhs   = 1
c ...
      msglvl = 1
c .....................................................................
c
c ...
      pt(1:64)    = 0
      iparm(1:64) = 0
c .....................................................................
c
c ... simetrico indefinido
      if( mtype .eq. -2) then
        iparm(1)  = 1 ! no solver default
        iparm(2)  = 2 ! fill-in reordering from METIS
        iparm(7)  = 2 ! numbers of iterative refinement steps
        iparm(10) = 8 ! perturbe the pivot elements with 1E-08
        iparm(21) = 1 ! Pivoting for symmetric indefinite matrices.                                   
        iparm(24) = 0 ! Parallel factorization control.
        iparm(27) = 0 ! 1 - checa a estrutura de dados  
c .....................................................................
c
c ... simetrico definido positivo
      else if( mtype .eq. 2) then
        iparm(1)  = 1 ! no solver default
        iparm(2)  = 2 ! fill-in reordering from METIS
        iparm(7)  = 2 ! numbers of iterative refinement steps 
        iparm(24) = 0 ! Parallel factorization control.
        iparm(27) = 0 ! 1 - checa a estrutura de dados  
      endif         
c .....................................................................
c
c ...
      time = get_time()  
c .....................................................................
c
c ... 
      phase = 13        
      msglvl = 0
#if _MKL_
      call pardiso (pt  , maxfct, mnum, mtype, phase, neq, a, ia4
     .            , ja  , idum, nrhs  , iparm, msglvl, b, x, error)
#endif
      time = get_time() - time 
c .....................................................................
c
c ...
      mem = (max(iparm(15),iparm(16) + iparm(17)))/1024.d0
c .....................................................................
c 
c ... produto:  x*Kx
      call matvec_csr_sym_v3(neq,ia4,ja,a,x,z,.true.)
      xkx = dot(x,z,neq)
c .......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq))
c ......................................................................
c
c ... Termination and release of memory
      phase = -1 ! release internal memory
      msglvl = 0
#if _MKL_
      call pardiso (pt, maxfct, mnum, mtype, phase,neq, ddum, idum,idum,
     .              idum, nrhs, iparm, msglvl, ddum, ddum, error)
#endif
c .....................................................................
c
c ...
      write(*,1100)neq,mem,xkx,norm,time
      write(10,'(a,a,d20.10,a,d20.10,a,f20.2,a,f20.2)')
     .       "PARDISO: "," x * Kx ",xkx," ||x|| ",norm
     .      ," Memory (MB)  ",mem," time ",time
c .....................................................................
c
c ...
      return
c ======================================================================
 1100 format(' (PARDISO) solver:'/
     1 5x,'Number of equations  = ',i20/
     2 5x,'Memory (MB)          = ',f20.2/
     4 5x,'x * Kx               = ',d20.10/
     5 5x,'|| x ||              = ',d20.10/
     6 5x,'CPU time (s)         = ',f20.2/)
      end
c ***********************************************************************  
c
c **********************************************************************
c * Data de criacao    : 05/04/2019                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * CALL_PARDISO : chama o sover pardiso mkl                           *    
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * ia4      - nao defindo                                             *
c * ia8      - vetor ia do CSR                                         *
c * nnz      - numero total de nao zeros                               *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * ia4      - vetor ia do CSR                                         *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************    
      subroutine ia8_to_ia4(ia4,ia8,neq)
      implicit none
      integer ia4(*),neq
      integer*8 ia8(*), i     
      do i = 1, neq+1        
        ia4(i) = ia8(i)
      enddo
      return
      end