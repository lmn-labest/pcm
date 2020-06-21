c **********************************************************************
c *                                                                    *
c *   MATVEC.F                                           31/08/2005    *
c *                                                                    *
c *   Este arquivo contem subrotinas para produto matriz-vetor e       *
c *   produto escalar de vetores                                       *
c *                                                                    *
c *   matvec_csr                                                       *
c *   matvec_csr_sym_v3                                                *
c *   matvec_csrsym                                                    *
c *   matvec_csrsym1*         (loops aninhados)                        *
c *   matvec_csrc                                                      *
c *   matvec_csrc1                                                     *
c *   matvec_csrc2                                                     *
c *   matvec_csrcsym                                                   *
c *   matvec_csrcsym1                                                  *
c *   matvec_csrcsym2                                                  *
c *   matvec_csrcr                                                     *
c *   matvec_csrcr1                                                    *
c *   matvec_csrcrsym                                                  *
c *   matvec_csrcrsym1                                                 *
c * --------------------------Poromecanico --------------------------- *
c *   matvec_csrc_block_pm                                                   *
c *   matvec_csrcr_sym_pm                                              *
c *   matvec_csrc_sym_pm                                               *
c *   matvec_csr_sym_pm                                                *
c * ------------------------------------------------------------------ *
c *   saxpb                                                            *
c *   dot                                                              *
c *   dot_par                                                          *
c *   ddot                                                             *
c *   aequalb                                                          *
c *   matvec_sym                                                       *
c *                                                                    *
c *   1 = loop interno desenrolado                                     *
c *   2 = loops desenrolados                                           *
c *                                                                    *
c **********************************************************************
      subroutine matvec_csr(neq,ia,ja,ad,a,al,x,y,neqf1i,neqf2i,
     .                      i_fmapi,i_xfi,i_rcvsi,i_dspli,dum4)
c **********************************************************************
c *                                                                    *
c *   MATVEC_CSR: produto matriz-vetor y = Ax  (A nao-simetrica),      *
c *   ----------                      coef. de A no formato CSR        *
c *                                    e grafo nao-simetrico.          *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq   - numero de equacoes                                       *
c *   ia(neq+1) - ia(i) informa a posicao no vetor a do primeiro       *
c *                     coeficiente nao-nulo da equacao i              *
c *   ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa      *
c *               a posicao k no vetor a                               *
c *   ad(neq) - diagonal da matriz A                                   *
c *   a(nad)  - coef. de A, sem a diagonal                             *
c *   al(*)   - nao utilizado                                          *
c *   x(neq+1)- vetor a ser multiplicado                               *
c *   y(neq+1)- nao definido                                           *
c *   neqf1i - numero de equacoes no buffer de recebimento (MPI)       *
c *   neqf2i - numero de equacoes no buffer de envio (MPI)             *
c *   i_fmapi- ponteiro para o mapa de comunicacao  (MPI)              *
c *   i_xfi  - ponteiro para o buffer de valores    (MPI)              *
c *   i_rcvsi- ponteiro extrutura da comunicacao    (MPI)              *
c *   i_dspli- ponteiro extrutura da comunicacao    (MPI)              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   y(neq) - vetor contendo o resultado do produto y=Ax              *
c *                                                                    *
c **********************************************************************      
      implicit none
      integer neq,ia(*),ja(*),i,k
      integer neqf1i,neqf2i
c ... ponteiros        
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
      real*8  ad(*),a(*),al(*),x(*),y(*),t
      real*8 dum4
c ......................................................................
      do 110 i = 1, neq
c
c ...    Produto da diagonal de A por x:
c
         t = ad(i)*x(i)
c
c ...    Produto da linha i pelo vetor x:
c
         do 100 k = ia(i), ia(i+1)-1 
            t = t + a(k)*x(ja(k))
  100    continue
c
c ...    Armazena o resultado em y(i):
c
         y(i) = t
  110 continue
c ......................................................................
      return
      end                
c **********************************************************************
c
c **********************************************************************
      subroutine matvec_csr_sym_v3(neq,ia  ,ja,a  ,x   ,y,flag)
c **********************************************************************
c * Data de criacao    : 01/05/2016                                    *
c * Data de modificaco : 00/00/0000                                    *
c * ------------------------------------------------------------------ * 
c * MATVEC_CSR_SYM_V3 : produto matriz-vetor CSR padrao simetrico      *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * neq       - numero de equacoes                                     *
c * ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro        *
c *                   coeficiente nao-nulo da linha   i                *
c * ja(*)     - ja(k) informa a coluna do coeficiente que ocupa        *
c *             a posicao k no vetor au                                *
c * a(*)   - coeficientes                                              *
c * x(neq) - vetor a ser multiplicado                                  *
c * y(neq) - nao definido                                              *
c * flag   - .true.  triangular superior                               *
c *        - .false. triangular inferior                               *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *                                                                    *
c * ------------------------------------------------------------------ *
c * y(neq) - vetor contendo o resultado do produto y = Ax              *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * CSR padrao com tem vetores (ia,ja,a)                               *
c * ja em ordem crescente                                              *
c **********************************************************************
      implicit none
      include 'time.fi'
      integer neq,ia(*),ja(*),i,k,kk,jak
      real*8  a(*),x(*),y(*),xi,t,s
      logical flag
c ......................................................................
c
c ...
      if(flag) then
        y(1:neq) = 0.0d0
        do i = 1, neq
          kk   = ia(i)
          xi   = x(i)
          y(i) = y(i) + a(kk)*xi
c ...
          do k = kk+1, ia(i+1)-1
            jak   = ja(k)
            s     = a(k)
c ... parte superior
            y(i)   =  y(i)  + s*x(jak)
c ... parte inferior
            y(jak) = y(jak) + s*xi
          enddo
c .....................................................................
        enddo
c .....................................................................
c
c ...
      else
       do i = 1, neq
          xi   = x(i)
c ...
          do k = ia(i), ia(i+1)-2
            jak   = ja(k)
            s     = a(k)
c ... parte superior
            t     =  t      + s*x(jak)
c ... parte inferior
            y(jak) = y(jak) + s*xi
          enddo
c .....................................................................
          y(i) = t + a(k)*xi
        enddo
c .....................................................................
      endif
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine matvec_csrc_block_pm(neq,nequ,ia,ja,iapu,japu,ad,al,
     1                        apul,x,y,
     2                        neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
     3                        i_dspli)
c **********************************************************************
c *                                                                    *
c *   matvec_csrc_block_pm:produto matriz-vetor y = Ax                 *
c *   (A Kuu, Kpp e Kpu )                                              *     
c *                   coef. de A no formato CSRC.                      *
c *       | Kuu  -Kpu |                                                *
c *   A = |           |                                                *
c *       | Kpu   Kpp |                                                *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq   - numero de equacoes                                       *
c *   nequ  - numero de equacoes no bloco Kuu                          *
c *   ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro      *
c *                     coeficiente nao-nulo da linha   i              *
c *   ja(nad  ) - ja(k) informa a coluna do coeficiente que ocupa      *
c *               a posicao k no vetor al                              *
c *   iapu(neqp+1)-ia(i) informa a posicao no vetor au do primeiro     *
c *                     coeficiente nao-nulo da linha   i              *
c *   japu(nadpu) - ja(k) informa a coluna do coeficiente que ocupa    *
c *               a posicao k no vetor apul                            *
c *   ad(neq)- diagonal da matriz A                                    *
c *   al(nad)- parte triangular inferior de Kuu e Kpp no formato CSR   *
c *   apul(*)- Kpu no formato CSR                                      *
c *   x(neq) - vetor a ser multiplicado                                *
c *   y(neq) - nao definido                                            *
c *   neqf1i - numero de equacoes no buffer de recebimento (MPI)       *
c *   neqf2i - numero de equacoes no buffer de envio (MPI)             *
c *   i_fmapi- ponteiro para o mapa de comunicacao  (MPI)              *
c *   i_xfi  - ponteiro para o buffer de valores    (MPI)              *
c *   i_rcvsi- ponteiro extrutura da comunicacao    (MPI)              *
c *   i_dspli- ponteiro extrutura da comunicacao    (MPI)              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   y(neq) - vetor contendo o resultado do produto y = Ax            *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'time.fi'
      integer neq,nequ,neqp,ia(*),ja(*),iapu(*),japu(*),i,ii,k,jak
      real*8  ad(*),al(*),apul(*),x(*),y(*),s,t,xi
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
c ......................................................................
      time0 = get_time()
c
c ... Loop nas linhas: Kuu e Kpp
c
      do 110 i = 1, neq
c
c ...    Produto da diagonal de A por x:
c
         xi = x(i)
         t  = ad(i)*xi
c
c ...    Loop nos coeficientes nao nulos da linha i:
c
         do 100 k = ia(i), ia(i+1)-1
            jak = ja(k)
            s   = al(k)
c
c ...       Produto da linha i pelo vetor x (parte triangular inferior):
c
            t   = t + s*x(jak)
c
c ...       Produto dos coef. da parte triangular superior por x(i):
c
            y(jak) = y(jak) + s*xi
  100    continue
c
c ...    Armazena o resultado em y(i):
c
         y(i) = t
  110 continue
c .....................................................................
c
c ... loop nas linha Kpu
      neqp = neq - nequ
      do 120 i = 1, neqp
        ii = nequ+i
        xi = x(ii)
        do 130 k = iapu(i), iapu(i+1)-1
          jak   = japu(k)
          s     = apul(k)
          
c ... Kpu
          y(ii)  =  y(ii)  + s*x(jak)
c ... Kup
          y(jak) =  y(jak) - s*xi
  130   continue
  120 continue
c .....................................................................
c
c ...
      matvectime = matvectime + get_time() - time0
c ......................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine matvec_csrcr_sym_pm(neq    ,dum0  ,ia
     1                          ,ja         ,iar   ,jar
     2                          ,ad         ,al    ,ar   
     3                          ,x          ,y    
     4                          ,neqf1i     ,neqf2i
     5                          ,i_fmapi    ,i_xfi ,i_rcvsi
     6                          ,i_dspli    ,dum4)
c **********************************************************************
c * Data de criacao    : 27/10/2016                                    *
c * Data de modificaco : 00/00/0000                                    *
c * ------------------------------------------------------------------ * 
c * MATVEC_CSRC_csr_sym_pm: produto matriz-vetor y = Ax  (A simetrica),*
c *                    coef. de A no formato CSRC e grafo simetrico.   *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq        - numero de equacoes                                    *
c * ia(neq+1)  - ia(i) informa a posicao no vetor au do primeiro       *
c *                   coeficiente nao-nulo da equacao i                *
c * ja(neq+1)  - ja(k) informa a coluna do coeficiente que ocupa       *
c *             a posicao k no vetor al                                *
c * iar(neq+1) - ia(i) informa a posicao no vetor ar do primeiro       *
c *             coeficiente nao-nulo da equacao i da parte retangular  *
c * jar(nadr)  - ja(k) informa a coluna do coeficiente que ocupa       *
c *             a posicao k no vetor ar                                *
c * ad(neq)    - diagonal da matriz A                                  *
c * al(nad)    - parte triangular inferior de A no formato CSR         *
c * ar(nad)    - parte retangular de A no formatp CSR                  *
c * x(neqovlp) - vetor a ser multiplicado                              *
c * y(neq)     - nao definido                                          *
c * neqf1i     - numero de equacoes no buffer de recebimento (MPI)     *
c * neqf2i     - numero de equacoes no buffer de envio (MPI)           *
c * i_fmapi    - ponteiro para o mapa de comunicacao  (MPI)            *
c * i_xfi      - ponteiro para o buffer de valores    (MPI)            *
c * i_rcvsi    - ponteiro extrutura da comunicacao    (MPI)            *
c * i_dspli    - ponteiro extrutura da comunicacao    (MPI)            *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * y(neq) - vetor contendo o resultado do produto y = Ax              *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'time.fi'
      integer*8 ia(*),iar(*),k
      integer neq,ja(*),jar(*),i,jak,dum0
      real*8  ad(*),al(*),ar(*),x(*),y(*),t,xi,s
      real*8 dum4   
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
c ......................................................................
      time0 = get_time()
c
c ... Loop nas linhas:
c
      do 200 i = 1, neq
c
c ...    Produto da diagonal de A por x:
c
         xi = x(i)
         t  = ad(i)*xi
c
c ...    Loop nos coeficientes nao nulos da linha i:
c
         do 100 k = ia(i), ia(i+1)-1
            jak = ja(k)
c
c ...       Produto da linha i pelo vetor x (parte triangular inferior):
c
            s = al(k)
            t   = t + s*x(jak)
c
c ...       Produto dos coef. da parte triangular superior por x(i):
c
            y(jak) = y(jak) + s*xi
  100    continue
c
         do 110 k = iar(i), iar(i+1)-1
            jak = jar(k)
c
c ...       Produto da linha i pelo vetor x (retangulo a direita):
c
            t   = t + ar(k)*x(jak)
  110    continue
 
c ...    Armazena o resultado em y(i):
 
         y(i) = t
  200 continue
      matvectime = matvectime + get_time() - time0
c ......................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine matvec_csrc_sym_pm(neq, ia, ja, ad, al, x, y)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 17/06/2020                                    *
c * ------------------------------------------------------------------ * 
c * MATVEC_CSRC_SYM_PM: produto matriz-vetor y = Ax  (A simetrica),    *
c *                   coef. de A no formato CSRC.                      *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * neq   - numero de equacoes                                         *
c * ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro        *
c *                   coeficiente nao-nulo da linha   i                *
c * ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa        *
c *             a posicao k no vetor au                                *
c * ad(neq)- diagonal da matriz A                                      *
c * al(nad)- parte triangular inferior de A, no formato CSR, ou        *
c *          parte triangular superior de A, no formato CSC            *
c * x(neq) - vetor a ser multiplicado                                  *
c * y(neq) - nao definido                                              *
c * ------------------------------------------------------------------ *
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * y(neq) - vetor contendo o resultado do produto y = Ax              *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'time.fi'
      integer*8 ia(*),k
      integer neq,ja(*),i,jak
      real*8  ad(*),al(*),x(*),y(*),s,t,xi
c ......................................................................
      time0 = get_time()
c
c ... Loop nas linhas:
c
      do 110 i = 1, neq
c
c ...    Produto da diagonal de A por x:
c
         xi = x(i)
         t  = ad(i)*xi
c
c ...    Loop nos coeficientes nao nulos da linha i:
c
         do 100 k = ia(i), ia(i+1)-1
            jak = ja(k)
            s   = al(k)
c
c ...       Produto da linha i pelo vetor x (parte triangular inferior):
c
            t   = t + s*x(jak)
c
c ...       Produto dos coef. da parte triangular superior por x(i):
c
            y(jak) = y(jak) + s*xi
  100    continue
c
c ...    Armazena o resultado em y(i):
c
         y(i) = t
  110 continue
      matvectime = matvectime + get_time() - time0
c ......................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine matvec_csrc_pm(neq      ,dum0  ,ia
     1                         ,ja       ,dum1  ,dum2
     2                         ,ad       ,au    ,al   
     3                         ,x        ,y   
     4                         ,neqf1i   ,neqf2i
     5                         ,i_fmapi  ,i_xfi ,i_rcvsi
     6                         ,i_dspli  ,dum4)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 15/12/2016                                    *
c * ------------------------------------------------------------------ * 
c * MATVEC_CSRC_PM: produto matriz-vetor y = Ax  (A nao-simetrica),    *
c *                   coef. de A no formato CSRC.                      *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * neq   - numero de equacoes                                         *
c * ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro        *
c *                   coeficiente nao-nulo da linha   i                *
c * ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa        *
c *             a posicao k no vetor au                                *
c * ad(neq)- diagonal da matriz A                                      *
c * al(nad)- parte triangular inferior de A, no formato CSR, ou        *
c *          parte triangular superior de A, no formato CSC            *
c * au(nad)- parte triangular superior de A, no formato CSR, ou        *
c *          parte triangular inferior de A, no formato CSC            *
c * x(neq) - vetor a ser multiplicado                                  *
c * y(neq) - nao definido                                              *
c * neqf1i - numero de equacoes no buffer de recebimento (MPI)         *
c * neqf2i - numero de equacoes no buffer de envio (MPI)               *
c * i_fmapi- ponteiro para o mapa de comunicacao  (MPI)                *
c * i_xfi  - ponteiro para o buffer de valores    (MPI)                *
c * i_rcvsi- ponteiro extrutura da comunicacao    (MPI)                *
c * i_dspli- ponteiro extrutura da comunicacao    (MPI)                *
c * ------------------------------------------------------------------ *
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * y(neq) - vetor contendo o resultado do produto y = Ax              *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'time.fi'
      integer*8 ia(*),k
      integer neq,ja(*),dum0,dum1,dum2,dum3,i,jak
      real*8  ad(*),al(*),au(*),x(*),y(*),s,t,xi
      real*8 dum4
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
c ......................................................................
      time0 = get_time()
c
c ... Loop nas linhas:
c
      do 110 i = 1, neq
c
c ...    Produto da diagonal de A por x:
c
         xi = x(i)
         t  = ad(i)*xi
c
c ...    Loop nos coeficientes nao nulos da linha i:
c
         do 100 k = ia(i), ia(i+1)-1
            jak = ja(k)
c
c ...       Produto da linha i pelo vetor x (parte triangular inferior):
c
            t   = t + al(k)*x(jak)
c
c ...       Produto dos coef. da parte triangular superior por x(i):
c
            y(jak) = y(jak) + au(k)*xi
  100    continue
c
c ...    Armazena o resultado em y(i):
c
         y(i) = t
  110 continue
      matvectime = matvectime + get_time() - time0
c ......................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine matvec_csr_pm(ni,nj,ia  ,ja,al  ,x   ,y,flag)
c **********************************************************************
c * Data de criacao    : 27/03/2016                                    *
c * Data de modificaco : 30/04/2016                                    *
c * ------------------------------------------------------------------ * 
c * MATVEC_CSR_PM: produto matriz-vetor retangular y = Ax              *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * neqi      - numero de linhas da matriz A                           *
c * neqj      - numero de colunas da matriz A                          *
c * ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro        *
c *                   coeficiente nao-nulo da linha   i                *
c * ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa        *
c *             a posicao k no vetor au                                *
c * al(nad)- coeficientes                                              *
c * x(neq) - vetor a ser multiplicado                                  *
c * y(neq) - nao definido                                              *
c * flag   - protudo transposto                                        *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *                                                                    *
c * ------------------------------------------------------------------ *
c * y(*)   - vetor contendo o resultado do produto                     * 
c *        .true.  - y(nj) = (At)x                                     *
c *        .false. - y(ni) = Ax                                        *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'time.fi'
      integer ni,nj,ia(*),ja(*),i,k,jak
      real*8  al(*),x(*),y(*),xi
c ... ponteiros      
      logical flag
c ......................................................................
c
c ... loop nas linha Kup
      if(flag) then
        y(1:nj) = 0.d0
        do i = 1, ni
          xi = x(i)
          do k = ia(i), ia(i+1)-1
            jak   = ja(k)
c ... Kup
            y(jak) =  y(jak) - al(k)*xi
          enddo      
        enddo   
c ......................................................................
c
c ... loop nas linha Kpu
      else 
        do i = 1, ni
          y(i) = 0.0d0
c ......................................................................
c
          do k = ia(i), ia(i+1)-1
            jak   = ja(k)
c ... Kpu
            y(i)  =  y(i)  + al(k)*x(jak)
          enddo
        enddo
      endif
c .....................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine saxpb(a,b,x,n,c)
c **********************************************************************
c *                                                                    *
c *   SAXPB: escalar . vetor + vetor                                   *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    a(n) - vetor de dimensao n                                      *
c *    b(n) - vetor de dimensao n                                      *
c *    x    - escalar                                                  *
c *    n    - dimensao                                                 *
c *    c(n) - nao definido                                             *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   c - resultado: a.x + b                                           *
c *                                                                    *
c **********************************************************************
      implicit none
      integer n,i
      real*8 a(*),b(*),c(*),x
c ......................................................................
      do 100 i = 1, n
         c(i) = a(i) * x + b(i)
  100 continue
      return
      end
      real*8 function dot(a,b,n)
c **********************************************************************
c *                                                                    *
c *   DOT: Produto escalar a.b                                         *
c *                                                                    *
c **********************************************************************
      implicit none
      integer n,i
      real*8  a(*),b(*)
c ......................................................................
      dot = 0.d0
      do 100 i = 1, n
         dot = dot + a(i)*b(i)
  100 continue
c ......................................................................    
      return
      end
      real*8 function dot_par(a,b,neq_doti)
c **********************************************************************
c *                                                                    *
c *   DOT: Produto escalar a.b                                         *
c *                                                                    *
c **********************************************************************
      implicit none        
      include 'time.fi'
      integer neq_doti
      integer n,i,k
      real*8  a(*),b(*),tmp
c ......................................................................
      time0 = get_time()
      tmp = 0.d0
      do 100 i = 1, neq_doti
         tmp = tmp + a(i)*b(i)
  100 continue
      dottime = dottime + get_time() - time0
c ......................................................................
      dot_par = tmp
c ......................................................................    
      return
      end      
      real*8 function ddot(x,y,n)
c **********************************************************************
c *                                                                    *
c *   DOT: Produto escalar x.y                                         *
c *                                                                    *
c **********************************************************************
      implicit none    
      real*8 x(*),y(*),temp
      integer i,m,mp1,n,inc
c ......................................................................
      ddot = 0.d0
      temp = 0.d0
      inc  = 10
c ......................................................................
      m = mod(n,inc)
      if(m .eq. 0) goto 40
      do 30 i = 1,m
         temp = temp + x(i)*y(i)
   30 continue
      if(n .lt. inc) goto 60
   40 mp1 = m + 1
      do 50 i = mp1, n, inc
         temp = temp + x(i  )*y(  i) + x(i+1)*y(i+1) + x(i+2)*y(i+2) +
     .                 x(i+3)*y(i+3) + x(i+4)*y(i+4) + x(i+5)*y(i+5) +
     .                 x(i+6)*y(i+6) + x(i+7)*y(i+7) + x(i+8)*y(i+8) +
     .                 x(i+9)*y(i+9)
   50 continue
   60 ddot = temp
      return
      end
c **********************************************************************
      subroutine matvec_sym(a,x,y,nl,nc)
c **********************************************************************
c * Data de criacao    : 18/01/2017                                    *
c * Data de modificaco :                                               * 
c * ------------------------------------------------------------------ *
c * MATVEC_SYM: produro matriz vetor simetrico                         *  
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * a    - matriz                                                      *
c * x    -                                                             *
c * y    -                                                             *
c * nl   - numero de linhas                                            *
c * nc   - numero de colunas                                           *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * y = ax                                                             *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * Matriz a com armazenamento da parte triangular inferior            *     
c **********************************************************************
      implicit none
      integer nc,nl,i,j
      real*8  a(nl,*),x(*),y(*),xj,aij
c ......................................................................
      y(1:nl) = 0.d0
      do 200 j = 1, nc 
         xj   = x(j) 
         y(j) = y(j) + a(j,j)*xj  
         do 100 i = j+1, nl 
            aij  = a(i,j) 
            y(i) = y(i) + aij*xj 
            y(j) = y(j) + aij*x(i)
  100    continue
  200 continue
c ......................................................................  
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine matvec(a,x,y,nl,nc)
c **********************************************************************
c * Data de criacao    : 18/01/2017                                    *
c * Data de modificaco :                                               * 
c * ------------------------------------------------------------------ *
c * MATVEC_SYM: produro matriz vetor simetrico                         *  
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * a    - matriz                                                      *
c * x    -                                                             *
c * y    -                                                             *
c * nl   - numero de linhas                                            *
c * nc   - numero de colunas                                           *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * y = ax                                                             *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      integer nc,nl,i,j
      real*8  a(nl,*),x(*),y(*),xj,aij
c ......................................................................
      y(1:nl) = 0.d0
      do 200 j = 1, nc 
         xj   = x(j) 
         do 100 i = 1, nl  
           y(i) = y(i) + a(i,j)*xj
  100    continue
  200 continue
c ......................................................................  
      return
      end
      subroutine aequalb(a,b,n)
c **********************************************************************
c *                                                                    *
c *   AEQUALB:                                                         *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    a(n) - vetor de dimensao n                                      *
c *    b(n) - vetor de dimensao n                                      *
c *    n    - dimensao                                                 *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   a = b                                                            *
c *                                                                    *
c **********************************************************************
      implicit none
      integer n,i
      real*8 a(*),b(*)
c ......................................................................
      do 100 i = 1, n
         a(i) = b(i)
  100 continue
      return
      end
      subroutine aminusb(a,b,c,n)
c **********************************************************************
c *                                                                    *
c *   AMINUSB:                                                         *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    a(n) - vetor de dimensao n                                      *
c *    b(n) - vetor de dimensao n                                      *
c *    n    - dimensao                                                 *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   c = a - b                                                        *
c *                                                                    *
c **********************************************************************
      implicit none
      integer n,i
      real*8 a(*),b(*),c(*)
c ......................................................................
      do 100 i = 1, n
         c(i) = a(i) - b(i)
  100 continue
      return
      end
      subroutine vsum(a,b,n,c)
c **********************************************************************
c *                                                                    *
c *   VSUM:                                                            *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    a(n) - vetor de dimensao n                                      *
c *    b(n) - vetor de dimensao n                                      *
c *    n    - dimensao                                                 *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   c = a + b                                                        *
c *                                                                    *
c **********************************************************************
      implicit none
      integer n,i
      real*8 a(*),b(*),c(*)
c ......................................................................
      do 100 i = 1, n
         c(i) = a(i) + b(i)
  100 continue
      return
      end      
      subroutine vsmul(a,x,n,c)
c **********************************************************************
c *                                                                    *
c *   VSUM:                                                            *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    a(n) - vetor de dimensao n                                      *
c *    x    - escalar                                                  *
c *    n    - dimensao                                                 *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   c = x*a                                                          *
c *                                                                    *
c **********************************************************************
      implicit none
      integer n,i
      real*8 a(*),x,c(*)
c ......................................................................
      do 100 i = 1, n
         c(i) = a(i)*x
  100 continue
      return
      end 
c **********************************************************************
      subroutine vet(a,b,c)
c **********************************************************************
c *                                                                    *
c *   VET: Produto vetorial axb                                        *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  c(3),a(3),b(3)
c ......................................................................
      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
c ......................................................................    
      return
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 23/10/2016                                    *
c * Data de modificaco : 31/10/2016                                    * 
c * ------------------------------------------------------------------ *   
c * FLOP_DOT : calcula o numero de operacoes de pontos flutuantes      *    
c * do produto interno                                                 *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * n        - dimensao dos vetores                                    *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      real*8 function flop_dot(n)
      implicit none
      integer n
      flop_dot = 2.d0*n
      return
      end
c ********************************************************************** 
c
c **********************************************************************
c * Data de criacao    : 23/10/2016                                    *
c * Data de modificaco : 31/10/2016                                    * 
c * ------------------------------------------------------------------ *   
c * FLOP_CSRC: calcula o numero de operacoes de pontos flutuantes      *    
c * da op matvec CSRC                                                  *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * nl       - numero de linhas                                        *
c * nad      - numero de termos fora da diagonal                       *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      real*8 function flop_csrc(nl,nad)
      implicit none
      integer nl
      integer*8 nad
      flop_csrc = nl + 4.d0*nad
      return
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 23/10/2016                                    *
c * Data de modificaco : 29/12/2016                                    * 
c * ------------------------------------------------------------------ *   
c * FLOP_CG : calcula o numero de operacoes de pontos flutuantes       *    
c * do CG                                                              *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq      - numero de equacores                                     *
c * nad      - numero de termos fora da diagonal                       *
c * it       - numero de iteracoes                                     *
c * icod     - 1 - CG                                                  *
c *            2 - PCG                                                 *      
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      real*8 function flop_cg(neq,nad,it,icod,mpi)
      implicit none
      integer neq,it,icod,ierr
      integer*8 nad
      real*8 flop_csrc,flop_dot,fmatvec,fdot,flops,gflops
      logical mpi
c
      fmatvec = flop_csrc(neq,nad)
      fdot    = flop_dot(neq)
c ... CG
      if(icod .eq. 1) then
        flops = (fmatvec + 2.d0*fdot + 6.d0*neq + 2.d0)*it     
c ... PCG
      elseif(icod .eq. 2) then
        flops = (fmatvec + 2.d0*fdot + 7.d0*neq + 2.d0)*it 
      endif
c .....................................................................
c
c ...
      flop_cg = flops
c .....................................................................
c
c ...      
      return
      end
c ********************************************************************** 
c
c **********************************************************************
c * Data de criacao    : 31/10/2016                                    *
c * Data de modificaco : 29/12/2016                                    * 
c * ------------------------------------------------------------------ *   
c * FLOP_SQRM: calcula o numero de operacoes de pontos flutuantes      *   
c * do SQRM                                                            *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq      - numero de equacores                                     *
c * nad      - numero de termos fora da diagonal                       *
c * it       - numero de iteracoes                                     *
c * icod     - 1 - SQRM                                                *
c *            2 - RSQRM                                               *      
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      real*8 function flop_sqrm(neq,nad,it,icod,mpi)
      implicit none
      integer neq,it,icod,ierr
      integer*8 nad
      real*8 flop_csrc,flop_dot,fmatvec,fdot,flops,gflops
      logical mpi
c
      fmatvec = flop_csrc(neq,nad)
      fdot    = flop_dot(neq)
c ... SQRM
      if(icod .eq. 1) then
        flops = (fmatvec + 5.d0*fdot + 6.d0*neq + 2.d0)*it     
c ... RSQRM
      elseif(icod .eq. 2) then
        flops = (fmatvec + 5.d0*fdot + 6.d0*neq + 14.d0)*it 
      endif
c .....................................................................
c
c ...
      flop_sqrm = flops
c .....................................................................
c
c ...      
      return
      end
c ********************************************************************** 
c
c **********************************************************************
c * Data de criacao    : 15/11/2016                                    *
c * Data de modificaco : 29/12/2016                                    * 
c * ------------------------------------------------------------------ *   
c * FLOP_BICGSTAB(2) : calcula o numero de operacoes de pontos         *
c * flutuantes do BICGSTAB(2)                                          *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq      - numero de equacores                                     *
c * nad      - numero de termos fora da diagonal                       *
c * it       - numero de iteracoes                                     *
c * icod     - 1 - BICGSTAB2                                           *
c *            2 - PBICGSTAB2                                          *      
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      real*8 function flop_bicgstab2(neq,nad,it,icod,mpi)
      implicit none
      integer neq,it,icod,ierr
      integer*8 nad
      real*8 flop_csrc,flop_dot,fmatvec,fdot,flops,gflops
      logical mpi
c
      fmatvec = flop_csrc(neq,nad)
      fdot    = flop_dot(neq)
c ... BICGSTAB2
      if(icod .eq. 1) then
        flops = 0.d0    
c ... PBICGSTAB2
      elseif(icod .eq. 2) then
        flops = (4.d0*fmatvec + 10.d0*fdot + 23.d0*neq + 17.d0)*it 
      endif
c .....................................................................
c
c ...
      flop_bicgstab2 = flops
c .....................................................................
c
c ...      
      return
      end
c ********************************************************************** 
 
