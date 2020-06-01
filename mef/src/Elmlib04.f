c *********************************************************************
c * Biblioteca de elementos de transporte                             * 
c *********************************************************************
c * ------------------ Elementos lineares --------------------------- *
c *                                                                   *
c * ELMT01_trans - elemento unidimensional (Garlekin)                 *
c *                                                                   *
c * ELMT02_trans - elemento unidimensional (conveccao-difusao, SUPG)  *
c *                                                                   *
c * ------------------ Elementos quadraticos ------------------------ *
c *                                                                   *
c * ----------------------------------------------------------------- *
c *                                                                   *
c ********************************************************************* 
c
c *********************************************************************
      subroutine elmt01_trans(e,iq,x,u,v,vel,p,s,ndm,nst,nel,isw)
c **********************************************************************
c * Data de criacao    : 23/05/2020                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *       
c * ELMT01_trans - elemento unidimensional (Garlekin)                  *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * e(10) - constantes fisicas                                         *
c *           e(14) = coeficiente de difusao                           *
c * iq(7)  - cargas mecanicas por elemento                             *
c * x(ndm,nem) - coordenadas nodais locais                             *
c * u(nst)     - graus de liberade por elemento (u)                    *
c * v(nst)     - derivada temporal                                     *
c * vel(nst)   - campo de velocidade de infiltracao(Darcy)             *
c * p(nst)     - nao definido                                          *
c * s(nst,nst) - nao definido                                          *  
c * ndm - dimensao                                                     *
c * nst - numero de graus de liberdade por elemento (2*1)              *
c * nel - numero do elemento                                           *
c * isw - codigo de instrucao                                          *
c *  1 =                                                               *
c *  2 = matriz K e forcas internas Ku                                 *
c *  3 = fluxos                                                        *
c *  4 =                                                               *
c *  5 =                                                               *
c *  6 =                                                               *
c *  7 =                                                               *     
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * e - constantes fisicas                                             *
c * s - matriz de elemento                                             *
c * p - isw = 2  residuo                                               *
c *     isw = 3                                                        *
c *     isw = 4                                                        *  
c *     isw = 7  porosidade                                            *
c *     isw = 8  porosidade                                            *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c *                                                                    *
c * vel(ksi) = vel1*h1(ksi)  + vel2*h2(ksi)                            * 
c * div_vel  = vel1*hx1(ksi) + vel2*hx2(ksi) = (vel2 - vel1)/l         *
c *                                                                    *
c * hx1 = -1.0/l                                                       *   
c * hx2 =  1.0/l                                                       *
c *                                                                    *
c * int(NiNi) = 2.d0/3.d0                                              *
c * int(NiNj) = 1.d0/3.d0                                              *
c *                                                                    *
c * int(Ni) = 1.0                                                      *
c **********************************************************************
      implicit none
      include 'transiente.fi'
      include 'gravity.fi'
      include 'load.fi'
      integer i,j,nel,isw
      integer tp
      integer nen,ndm,nst,iq(*)
c ...
      real*8 u(*),v(*),vel(*),vel_mod,vx1,vx2
      real*8 p(*),s(nst,*)
c ...
      real*8 m11,m12
c ...
      real*8 wt,det,l
      real*8 rn(3)      
c ...
      real*8 e(*),x(ndm,*)
c ...
      real*8 difus,poro 
c ...
      data rn / -1.d0,
     1           1.d0,              
     2           0.0d0/                
c
      data nen/3/
c ......................................................................
c
c ...
      difus = e(14)
      poro  = 1.d0            
c ......................................................................
c
c ...
      if(isw .eq. 2 .or. isw .eq. 3)then
        call jacob1d(x,det,ndm,nel)
c ... comprimento do elmento
        l = 2.d0*det
      endif      
c .....................................................................     
c
c ===
      goto (100, 200, 300, 400) isw
c ======================================================================
c
c.... calculo do delta t critico               
c
c ......................................................................
  100 continue
      return
c ======================================================================
c
c ... Matriz de rigidez:
c
c ......................................................................
  200 continue
c ... Matriz de rigidez:
      s(1,1) = 0.0d0
      s(1,2) = 0.0d0
      s(2,2) = 0.0d0
      s(2,1) = 0.0d0
      p(1)   = 0.0d0
      p(2)   = 0.0d0
c .....................................................................
c
c ...
      vx1 = vel(1)
      vx2 = vel(2)
c .....................................................................
c
c ... velocidade media
      vel_mod = dabs(vx1 + vx2)*0.5d0      
c .......................................................................
c
c ... Matriz de massa consistente:
      m11 = poro*(1.d0/3.0)*l
      m12 = m11*0.5d0
c .....................................................................
c
c ... Matriz de rigidez:     
c ... Kd
      wt     = difus/l
      s(1,1) = wt             ! hx(1)*hx(1)*difus*det
      s(1,2) =-wt             ! hx(1)*hx(2)*difus*det      
      s(2,1) =-wt             ! hx(2)*hx(1)*difus*det
      s(2,2) = wt             ! hx(2)*hx(2)*difus*det
c ... Kc
      wt = (vx2 - vx1)/3.d0 
      s(1,1) = s(1,1) + wt           ! div_u*h(1)*h(2)
      s(1,2) = s(1,2) + wt*0.5d0     ! div_u*h(1)*h(2)
      s(2,1) = s(2,1) + wt*0.5d0     ! div_u*h(2)*h(1)
      s(2,2) = s(2,2) + wt           ! div_u*h(2)*h(2)
c ... Ka      
      s(1,1) = s(1,1) - (2.d0*vx1 + vx2)/6.d0  ! h(1)*vel*hx(1)*det
      s(1,2) = s(1,2) + (2.d0*vx1 + vx2)/6.d0  ! h(1)*vel*hx(2)*det
      s(2,1) = s(2,1) - (2.d0*vx2 + vx1)/6.d0  ! h(2)*vel*hx(1)*det
      s(2,2) = s(2,2) + (2.d0*vx2 + vx1)/6.d0  ! h(2)*vel*hx(2)*det
c .....................................................................
c
c ... Forcas internas:
      p(1) = s(1,1)*u(1) + s(1,2)*u(2) + m11*v(1) + m12*v(2)
      p(2) = s(2,1)*u(1) + s(2,2)*u(2) + m12*v(1) + m11*v(2)       
c .....................................................................
c
c ...
      wt = alfa*dt  
      s(1,1) = m11 + s(1,1)*wt
      s(1,2) = m12 + s(1,2)*wt
      s(2,1) = m12 + s(2,1)*wt 
      s(2,2) = m11 + s(2,2)*wt
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
      return
c ======================================================================
c
c ... Cargas distribuidas no volume, no contorno e variaveis no passo
c     anterior:
c
c ......................................................................
 400  continue
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
      return
c
c ......................................................................
c
c ... velocidade nodais:
c
c ......................................................................
  800 continue
c ......................................................................
      return
c ======================================================================
      return
      end
c ********************************************************************* 
c
c *********************************************************************
      subroutine elmt02_trans(e,iq,x,u,v,vel,p,s,ndm,nst,nel,isw)
c **********************************************************************
c * Data de criacao    : 31/05/2020                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *       
c * ELMT02_trans - elemento unidimensional (conveccao-difusao, SUPG)   *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * e(10) - constantes fisicas                                         *
c *           e(14) = coeficiente de difusao                           *
c * iq(7)  - cargas mecanicas por elemento                             *
c * x(ndm,nem) - coordenadas nodais locais                             *
c * u(nst)     - graus de liberade por elemento (u)                    *
c * v(nst)     - derivada temporal                                     *
c * vel(nst)   - campo de velocidade de infiltracao(Darcy)             *
c * p(nst)     - nao definido                                          *
c * s(nst,nst) - nao definido                                          *  
c * ndm - dimensao                                                     *
c * nst - numero de graus de liberdade por elemento (2*1)              *
c * nel - numero do elemento                                           *
c * isw - codigo de instrucao                                          *
c *  1 =                                                               *
c *  2 = matriz K e forcas internas Ku                                 *
c *  3 = fluxos                                                        *
c *  4 =                                                               *
c *  5 =                                                               *
c *  6 =                                                               *
c *  7 =                                                               *     
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * e - constantes fisicas                                             *
c * s - matriz de elemento                                             *
c * p - isw = 2  residuo                                               *
c *     isw = 3                                                        *
c *     isw = 4                                                        *  
c *     isw = 7  porosidade                                            *
c *     isw = 8  porosidade                                            *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c *                                                                    *
c * vel(ksi) = vel1*h1(ksi)  + vel2*h2(ksi)                            * 
c * div_vel  = vel1*hx1(ksi) + vel2*hx2(ksi) = (vel2 - vel1)/l         *
c *                                                                    *
c * hx1 = -1.0/l                                                       *   
c * hx2 =  1.0/l                                                       *
c *                                                                    *
c * int(NiNi) = 2.d0/3.d0                                              *
c * int(NiNj) = 1.d0/3.d0                                              *
c *                                                                    *
c * int(Ni) = 1.0                                                      *
c **********************************************************************
      implicit none
      include 'transiente.fi'
      include 'gravity.fi'
      include 'load.fi'
      integer i,j,nel,isw
      integer tp
      integer nen,ndm,nst,iq(*)
c ...
      real*8 u(*),v(*),vel(*),vel_mod,vx1,vx2
      real*8 p(*),s(nst,*)
c ...
      real*8 m(2,2)
c ...
      real*8 wt,det,tau,upwind,l
      real*8 rn(3)      
c ...
      real*8 e(*),x(ndm,*)
c ...
      real*8 difus 
c ...
      real*8 peclet
c
      data rn / -1.d0,
     1           1.d0,              
     2           0.0d0/                
c
      data nen/3/
c ......................................................................
c
c ...
      difus = e(14)            
c ......................................................................
c
c ...
      if(isw .eq. 2 .or. isw .eq. 3)then
        call jacob1d(x,det,ndm,nel)
c ... comprimento do elmento
        l = 2.d0*det
      endif      
c .....................................................................     
c
c ===
      goto (100, 200, 300, 400) isw
c ======================================================================
c
c.... calculo do delta t critico               
c
c ......................................................................
  100 continue
      return
c ======================================================================
c
c ... Matriz de rigidez:
c
c ......................................................................
  200 continue
c ... Matriz de rigidez:
      s(1,1) = 0.0d0
      s(1,2) = 0.0d0
      s(2,2) = 0.0d0
      s(2,1) = 0.0d0
      p(1)   = 0.0d0
      p(2)   = 0.0d0
c .....................................................................
c
c ...
      vx1 = vel(1)
      vx2 = vel(2)
c .....................................................................
c
c ... velocidade media
      vel_mod = dabs(vx1 + vx2)*0.5d0      
c .......................................................................
c
c ...
      tau = upwind(vel_mod,l,difus)
c .......................................................................
c
c ... Matriz de massa consistente:
      m(1,1) = (1.d0/3.d0)*l - tau*(2.d0*vx1 + vx2)/6.d0
      m(1,2) = (1.d0/6.d0)*l - tau*(2.d0*vx2 + vx1)/6.d0
      m(2,1) = (1.d0/6.d0)*l + tau*(2.d0*vx1 + vx2)/6.d0
      m(2,2) = (1.d0/3.d0)*l + tau*(2.d0*vx2 + vx1)/6.d0
c .....................................................................
c
c ... Matriz de rigidez:     
c ... Kd
      wt     = difus/l
      s(1,1) = wt             ! hx(1)*hx(1)*difus*det
      s(1,2) =-wt             ! hx(1)*hx(2)*difus*det      
      s(2,1) =-wt             ! hx(2)*hx(1)*difus*det
      s(2,2) = wt             ! hx(2)*hx(2)*difus*det
c ... Kc
      wt = (vx2 - vx1)/3.d0 
      s(1,1) = s(1,1) + wt           ! div_u*h(1)*h(2)
      s(1,2) = s(1,2) + wt*0.5d0     ! div_u*h(1)*h(2)
      s(2,1) = s(1,2)                ! div_u*h(2)*h(1)
      s(2,2) = s(2,2) + wt           ! div_u*h(2)*h(2)
c ... Ka      
      s(1,1) = s(1,1) - (2.d0*vx1 + vx2)/6.d0  ! h(1)*vel*hx(1)*det
      s(1,2) = s(1,2) + (2.d0*vx1 + vx2)/6.d0  ! h(1)*vel*hx(2)*det
      s(2,1) = s(2,1) - (2.d0*vx2 + vx1)/6.d0  ! h(2)*vel*hx(1)*det
      s(2,2) = s(2,2) + (2.d0*vx2 + vx1)/6.d0  ! h(2)*vel*hx(2)*det
c .....................................................................
c
c ... Matriz de rigidez (SUPG):     
c ... sKd = 0.d0
c ... sKc
c     wt = 0.d0                
c     s(1,1) = s(1,1) + wt                              
c     s(1,2) = s(1,2) + wt                              
c     s(2,1) = s(1,2)                                   
c     s(2,2) = s(2,2) + wt                             
c ... sKa
      wt = tau*(vx1*vx1 + vx1*vx2 + vx2*vx2)/(3.d0*l)      
      s(1,1) = s(1,1) + wt  ! (tau*vel*hx(1))*(vel*hx(1))*det
      s(1,2) = s(1,2) - wt  ! (tau*vel*hx(1))*(vel*hx(2))*det
      s(2,1) = s(2,1) - wt  ! (tau*vel*hx(2))*(vel*hx(1))*det
      s(2,2) = s(2,2) + wt  ! (tau*vel*hx(2))*(vel*hx(2))*det
c .....................................................................

c
c ... Forcas internas:
      p(1) = s(1,1)*u(1)+s(1,2)*u(2) + m(1,1)*v(1) + m(1,2)*v(2)
      p(2) = s(2,1)*u(1)+s(2,2)*u(2) + m(2,1)*v(1) + m(2,2)*v(2)
c .....................................................................
c
c ...
      wt = alfa*dt  
      s(1,1) = m(1,1) + s(1,1)*wt
      s(1,2) = m(1,2) + s(1,2)*wt
      s(2,1) = m(2,1) + s(2,1)*wt 
      s(2,2) = m(2,2) + s(2,2)*wt
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
      return
c ======================================================================
c
c ... Cargas distribuidas no volume, no contorno e variaveis no passo
c     anterior:
c
c ......................................................................
 400  continue
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
      return
c
c ......................................................................
c
c ... velocidade nodais:
c
c ......................................................................
  800 continue
c ......................................................................
      return
c ======================================================================
      return
      end
c ********************************************************************* 