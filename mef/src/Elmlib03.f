c *********************************************************************
c * Biblioteca de elementos poro-mecanicos                            *
c * ----------------------------------------------------------------- *
c * Elastico:                                                         *
c * ----------------------------------------------------------------- *
c * ------------------ Elementos lineares --------------------------- *
c *                                                                   *
c * ------------------ Elementos quadraticos ------------------------ *
c *                                                                   *
c * ELMT12_PM - barra de 3 nos para o problema poro mecanico elastico *
c *                                                                   *
c * ----------------------------------------------------------------- *
c * Plastico-estatico:                                                *
c * ----------------------------------------------------------------- *
c * ------------------ Elementos lineares --------------------------- *
c *                                                                   *
c * ------------------ Elementos quadraticos ------------------------ *
c *********************************************************************
c
c **********************************************************************
      subroutine elmt12_pm(e,iqm,iqh,x,u,dp,p,s,txn,ndm,nst,nel,isw
     .                    ,block_pu)
c **********************************************************************
c * Data de criacao    : 26/11/2019                                    *
c * Data de modificaco : 21/05/2020                                    * 
c * ------------------------------------------------------------------ *       
c * ELMT12_PM: Elemento hexaedricos de 20 nos para problemas           *  
c * poromecanico elasticos                                             *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * e(10) - constantes fisicas                                         *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *           e(3) = coeficiente de Darcy                              *
c *           e(4) = modulo de Biot                                    *
c *           e(5) = coeficiente de Biot                               *
c *           e(6) = massa especifica homogenizada do meio poroso      *
c *           e(7) = massa especifica do fluido                        *
c * iqm(7) - cargas mecanicas por elemento                             *
c * iqh(7) - cargas hidraulicas por elemento                           *
c * x(ndm,nem) - coordenadas nodais locais                             *
c * u(nst)     - graus de liberade por elemento (u + p)                *
c * dp(*)      - delta p ( p(n  ,0  ) - p(0) )                         *
c * p(nst)     - nao definido                                          *
c * s(nst,nst) - nao definido                                          *
c * txn(6,nen) - tensoes nodais                                        *      
c * ndm - dimensao                                                     *
c * nst - numero de graus de liberdade por elemento (3*20 + 1*8)       *
c * nel - numero do elemento                                           *
c * isw - codigo de instrucao                                          *
c *  1 = delta t critico                                               *
c *  2 = matriz K e forcas internas Ku                                 *
c *  3 = tesnsoes e fluxos                                             *
c *  4 = forcas de volume, superficies e integrais do passo            *
c *    de tempo anterior                                               *
c *  5 = Tensoes iniciais                                              *
c *  6 =                                                               *
c *  7 = variacao da porosidade                                        *
c * block_pu - true - armazenamento em blocos Kuu,Kpp e kpu            *
c *            false- aramzenamento em unico bloco                     *      
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * e - constantes fisicas                                             *
c * s - matriz de elemento                                             *
c * p - isw = 2  residuo                                               *
c *     isw = 3  tensao e fluxo                                        *
c *     isw = 4  cargas de superfice, volume e integras do passo       *
c *     de tempo anterior                                              *
c *     isw = 7  porosidade                                            *
c *     isw = 8  porosidade                                            *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * u(1:3) - isw = 2 diferenca de deslocamentos                        *
c * u(1:3) - isw = 3 e isw = 4 deslocamentos                           *
c * u(4:5) - isw = 2 diferenca de pressao                              *
c * u(4:5) - isw = 3 e isw = 4 pressao                                 *
c **********************************************************************
      implicit none
      include 'transiente.fi'
      include 'gravity.fi'
      include 'load.fi'
      common /gauss/ pg, wg
      integer ndm,nst,nel,isw
      integer i,j,l,k,tpm,tph,tp,tp1
      integer nen,nint,nint_face,lx,ly,lz
      integer iqm(*),iqh(*)
c ...
      integer cargm,cargh
c ...
      real*8 u(*),dp(*)
      real*8 p(*),s(nst,*)
c ... funcoes de interpolacao
      real*8 hu(3),hux(3)
      real*8 hp(2),hpx(2)
c ...
      real*8 ri,w,wt1,wt2,wt3,det
      real*8 rn(3)
      real*8 pi,epsi,txi,txn,dpm,pm
c ... integracao numerica de hexaedros           
      real*8 pg(10,10),wg(10,10)
c ...
      real*8 dt_c
      real*8 e(*),x(ndm,*)
c ...
      real*8 l_c
      real*8 perm,a,b,c,ym,ps
      real*8 fluid_d,dt_perm,pm_d
      real*8 dt_fluid,dt_fluid_perm,lambda,mi
      real*8 imod_biot,coef_biot
      real*8 fluid_sw
      real*8 a1,a2,a3,tmp
c ...
      real*8 scale
      parameter (scale = 1.d-06)
c ...
      logical block_pu 
c
      data rn / -1.d0,
     1           1.d0,              
     2           0.0d0/                
c
      data nen/3/
      parameter (nint = 4)
c ......................................................................
c
c ...
      pm_d      = e(6)*scale            
      fluid_d   = e(7)*scale
c ......................................................................
c
c ... fluid specific weight
      fluid_sw  = fluid_d*gravity_mod
c ......................................................................
c
c ... matriz constitutiva:
      ym       = e(1)
      ps       = e(2)
c ......................................................................
c
c ... 
      perm     = e(3)/fluid_sw
c ......................................................................
c
c ...
      imod_biot= 1.d0/e(4)
      coef_biot= e(5)
c ......................................................................
c
c ... 
      dt_perm  = perm*dt
      dt_fluid_perm = fluid_d*dt_perm
c ......................................................................
c
c ...
      a1       = 1.d0 - ps
      a2       = 1.d0 - 2.d0*ps
      a3       = 1.d0 + ps
c ......................................................................
c
c ...
c     a        = (ym*a1)/(a3*a2)
c     b        = ps/a1
c     c        = 0.5d0*(a2/a1) 
      a        = ym
      b        = 0.d0
      c        = 0.d0
c .....................................................................
c
c ...
      if(isw .eq. 2 .or. isw .eq. 3 .or. isw .eq. 4 .or. isw .eq. 8)then
        det = 0.d0
        do i = 1, ndm
          det = det + (x(i,2)-x(i,1))*(x(i,2)-x(i,1))
        enddo   
        if (det .le. 0.d0) then 
          print*,'*** Subrotina ELMT: determinante <= 0 ',nel
          stop
        endif      
        det = 0.5d0*dsqrt(det)
      endif
c .....................................................................     
c
c ===
      goto (100, 200, 300, 400, 500, 600, 700, 800) isw
c ======================================================================
c
c.... calculo do delta t critico               
c
c ......................................................................
  100 continue
c ...
c     volum = hexa_vol(x)
c     l_c   = volum**(1.0d0/3.d0)
c ...
c     dt_c  = ((l_c*l_c)/perm) 
c    .      * ( imod_biot+( coef_biot*coef_biot*a3*a2 )/( ym*a1 ) )
c     p(1)  = dt_c
c .....................................................................
      return
c ======================================================================
c
c ... Matriz de rigidez:
c
c ......................................................................
  200 continue
c ... Matriz de rigidez:
      do i = 1, nst
        do j = 1, nst
          s(j,i) = 0.d0
        enddo
        p(i) = 0.d0
      enddo
c .....................................................................
c       
c ... 
      do 215 lx = 1, nint
        ri = pg(lx,nint)
c ...                                               
        call sfbar2_m(hp,hpx,ri,.true.,.true.)
        hpx(1:2)= hpx(1:2)/det
c .....................................................................
c
c ...
        w   = wg(lx,nint)*det
c .....................................................................
c
c ...
        call sfbar3_m(hu,hux,ri,.true.,.true.)
        hux(1:3)= hux(1:3)/det
c .....................................................................
c
c ... Kuu ( Int((Bt)*C*(B)*dV) )
        wt1 = w*a
        do i = 1, 3
          do j = 1, 3
            s(i,j) = s(i,j) + hux(i)*hux(j)*wt1 
          enddo
        enddo    
c .....................................................................
c
c ...-Kpu ( Int((b*(N~t)*(1t)(B)*dV) )  
        wt2 = w*coef_biot
        do  i = 1, 2
          l   = i + 3
          wt1 = hp(i)*wt2
          do  j = 1, 3
            s(l,j) = s(l,j) - hux(j)*wt1
          enddo
c .....................................................................
        enddo    
c .....................................................................
c
c ...    
        wt1 = w*imod_biot
        wt2 = w*dt_perm
c ...................................................................
c
c ... Kpp ( Int((1/M)*(N~t)*(N~)*dV) )    
        do i = 1, 2
          l = i + 3
          do j = 1, 2
            k = j + 3
            s(l,k) = s(l,k) - wt1*hp(i)*hp(j) - wt2*hpx(i)*hpx(j)
          enddo
c .....................................................................
        enddo    
c .....................................................................
  215 continue
c .....................................................................
c
c ...Kup = kpu  
      do i = 1, 2 
        l   = i + 3  
        do  j = 1, 3 
          s(j,l) = s(l,j)
        enddo 
      enddo   
c ................................................................ 
c
c ...
      if(block_pu) then
c ... -Kpp
        do  i = 1, 2
          l = i + 3
          do j = 1, 2
            k = j + 3
            s(l,k) = -s(l,k)
          enddo    
        enddo    
c ....................................................................
c
c ... Kpu = -Kpu
        do i = 1, 2
          l   = i + 3 
          do j = 1, 3
            s(l,j) = -s(l,j)
          enddo    
        enddo   
c .....................................................................
      endif
c .....................................................................
c
c ... Forcas internas:
      call lku_m(s,u,p,nst)
c     if(block_pu) then
c       call lku_m(s,u,p,nst)
c     else
c       call lku_sym(s,u,p,nst)    
c     endif    
c .....................................................................
c
c .....................................................................      
      return  
c ======================================================================
c
c ... Tensoes nodais e fluxo nodais:
c
c ......................................................................
  300 continue
c ...
c     p(1:2) = sigma
c     p(3:4) = sigmaBiot
c     p(5:6) = darcyVel
c 
c ... tensao nodal total
      do 310 i = 1, 2
c       tp = (i-1)*1 + 1
        tp  = i
c ... tensao efetiva de biot (1x1+1=2)         
c       tp1 = (i-1)*1 + 1   
        tp1 = i + 1   
c ...   p(1...8) = u(4 5) 
        l   = i  + 3 
c ... calculo do determinante
        call sfbar2_m(hp,hpx,0.d0,.false.,.true.)
        hpx(1:2)= hpx(1:2)/det
c ... calculo da derivadas das funcoes de interpolacao
        call sfbar3_m(hu,hux,rn(i),.false.,.true.)
        hux(1:3)= hux(1:3)/det
c .....................................................................
c ... deformacoes nos pontos de integracao
        epsi = hux(1)*u(1) + hux(2)*u(2) + hux(3)*u(3)
c ... tensao nos pontos de integracao
        txi  = ym*epsi
c .....................................................................            
c
c ... Tensao = D*deformacao elastica - biot*dp
        dpm     = coef_biot*dp(i)
         pm     = coef_biot*u(l)
c ... tensao total
        p(tp)   = txi   - dpm
c ... tensao efetiva de biot
        p(tp1)  = txi   + pm
c .....................................................................
  310 continue
c .....................................................................
c
c ... fluxo nodal 
      do 320 i = 1, 2
c ... fuxo (2x1 + 2x1 + 1 = 5)                 
        tp = i + 4
c ...
        call sfbar2_m(hp,hpx,0.d0,.false.,.true.)
        hpx(1:2)= hpx(1:2)/det
c ... p nos pontos de integracao
        call darcy_vel_1D(perm,hpx,u(4),2,p(tp))  
c .....................................................................
  320 enddo
c ......................................................................
      return
c ======================================================================
c
c ... Cargas distribuidas no volume, no contorno e variaveis no passo
c     anterior:
c
c ......................................................................
 400  continue
c ...                            
      do 415 lx = 1, nint
        ri = pg(lx,nint)
c ...                                               
        call sfbar2_m(hp,hpx,ri,.true.,.true.)
        hpx(1:2)= hpx(1:2)/det     
c .....................................................................
c
c ...
        w   = wg(lx,nint)*det
c .....................................................................
c
c ...
        call sfbar3_m(hu,hux,ri,.true.,.true.)
        hux(1:3)= hux(1:3)/det
c .....................................................................
c
c ... deformacoes nos pontos de integracao
        epsi = hux(1)*u(1) + hux(2)*u(2) + hux(3)*u(3)
c ... tensao nos pontos de integracao
        txi  = ym*epsi
c .....................................................................            
c
c ... delta p nos pontos de integracao
        dpm = hp(1)*dp(1) + hp(2)*dp(2)
c .....................................................................
c
c ...
        dpm = coef_biot*dpm
        txi = txi - dpm
c .....................................................................
c
c ...
        wt1 = w  
        wt2 = w*pm_d 
c .....................................................................
c
c ... Fu = int(BeT*sigma*dv)
        do i = 1, 3
          p(i) = p(i) + wt1*hux(i)*txi 
        enddo
c .....................................................................        
c
c ...
        wt1 = w*dt_perm 
c .....................................................................
c
c ... Fp = int(dt*k*(B~T)*(B~)*pe*dv)
        do i = 1, 2 
          l = i + 3
          tmp = 0.0d0
          do j = 1, 2
            k    = j + 3
            tmp  = tmp + wt2*u(k)*(hpx(i)*hpx(j))
          enddo
          p(l) = p(l) - tmp
c .....................................................................
        enddo       
  415 continue
c .....................................................................
c
c .....................................................................
c
c     
c ... bloco Fp = -Fp
      if(block_pu) then
        p(3:nst) =-p(3:nst)
      endif
c .....................................................................
      return
c ======================================================================
c
c ... Tensoes iniciais:                  
c
c ......................................................................
  500 continue
c ....
      return
c ======================================================================
c
c ... Matriz geometrica:
c
c ......................................................................
  600 continue
      stop
c ======================================================================
c
c ... Variacao da porosidade             
c
c ......................................................................
c 
c ...  
  700 continue                      
c ... variacao de porosidade
c     do 710 i = 1, 8
c ... calculo do terminante
c       call sfhexa8_m(hp,hpx,hpy,hpz,rn(i),sn(i),tn(i),.false.,.true.)
c       call jacob3d_m(hpx,hpy,hpz,xj,xji,x,det,8,nel ,.true.)
c .....................................................................
c
c ... calculo da derivadas das funcoes de interpolacao
c       call sfhexa20_m(hu,hux,huy,huz,rn(i),sn(i),tn(i),.false.,.true.)
c       call jacob3d_m(hux,huy,huz,xj,xji,x,det,20,nel ,.false.)
c .....................................................................
c
c ... variacao da deformacao
c       call deform3d(hux,huy,huz,u,epsi,20)
c ......................................................................
c
c ... variacao da porosidade = biot (eps11 + eps22 + eps33) + dp/M
c ...   p(1...8) = u(61...68) 
c       l    = i  + 60 
c       pi   = imod_biot*u(l)
c .....................................................................
c
c ... variacao da porosidade = biot*tr(deps)total + dp/M 
c       p(i) = coef_biot*(epsi(1) + epsi(2) + epsi(3) ) + pi
c ..................................................................... 
  710 continue
c .....................................................................
      return
c
c ... velocidade nodais:
c
c ......................................................................
  800 continue
c ...
c     p(1:2) = darcyVel
c 
c ... fluxo nodal 
      do 810 i = 1, 2
c ... fuxo (1)                 
        tp = i
c ...
        call sfbar2_m(hp,hpx,0.d0,.false.,.true.)
        hpx(1:2)= hpx(1:2)/det
c ... p nos pontos de integracao
        call darcy_vel_1D(perm,hpx,u(4),2,p(tp))  
c .....................................................................
  810 enddo
c ......................................................................
      return
c ======================================================================
      return
      end
c *********************************************************************    


