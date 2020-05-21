c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * PRINTNODE: imprime uma grandeza do no no tempo                     *
c * ------------------------------------------------------------------ *
c * Parametro de Entrada :                                             *
c * ------------------------------------------------------------------ *
c * u       - vetor com a grandeza a ser escrita                       *
c * no      - numero do no                                             *
c * ndf     - graus de liberdade da grandeza                           *
c * istep   - passo de tempo                                           *
c * t       - tempo                                                    *
c * nameres - nome da gradeza a ser escrita                            *
c * prename - nome do arquivo de saida                                 *
c * nout    - numero do arquivo a ser escrito                          *
c * open    - abertura de um novo arquivo .true.                       *
c * ------------------------------------------------------------------ *
c * Parametro de Saida :                                               *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************
      subroutine printnode(u,no,ndf,istep,t,nameres,prename,nout,code
     .                    ,open)
      implicit none
c ===        
      character*30 nameres
      real*8 u(ndf,*)
      real*8 t
      integer no,istep,nout,ndf,i,code
      character*80 fileout,prename,name
      logical open
c =====================================================================        
c
c ===
c ... abre um novo arquivo
      if(open) then
        fileout = name(prename,no,0,code)
        open(unit=nout,file=fileout)
        write(nout,'(a,a,a,i9)')
     .  '# Variacao ',trim(nameres),' no tempo no ponto',no
      else
c ... reabre o arquivo e add uma nova linha      
        fileout = name(prename,no,0,code)
        open(unit=nout,file=fileout,status ='old',access='append') 
      endif
c =====================================================================        
c
c === 
      write(nout,'(i9,10e18.10)')istep,t,(u(i,no), i=1, ndf)
c      write(nout,*)istep,',',istep*dt,',',u(no)
      close(nout)
      return
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 29/01/2017                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * PRINTPI: imprime uma grandeza do ponto de integracao no tempo      *
c * ------------------------------------------------------------------ *
c * Parametro de Entrada :                                             *
c * ------------------------------------------------------------------ *
c * u       - vetor com a grandeza a ser escrita                       *
c *elo      - numero do elemento                                       *
c * ndf     - graus de liberdade da grandeza                           *
c * npi     - pontos de integracao                                     *
c * istep   - passo de tempo                                           *
c * t       - tempo                                                    *
c * nameres - nome da gradeza a ser escrita                            *
c * prename - nome do arquivo de saida                                 *
c * nout    - numero do arquivo a ser escrito                          *
c * open    - abertura de um novo arquivo .true.                       *
c * ------------------------------------------------------------------ *
c * Parametro de Saida :                                               *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************
      subroutine printpi(u,el,ndf,npi,istep,t
     1                  ,nameres,prename,nout
     2                  ,code   ,open)
      implicit none
c ===        
      character*30 nameres
      integer el,istep,nout,ndf,npi,i,j,code
      real*8 u(ndf,npi,*)
      real*8 t
      character*80 fileout,prename,name
      logical open
c .....................................................................
c
c ...
c ... abre um novo arquivo
      if(open) then
        fileout = name(prename,el,0,code)
        open(unit=nout,file=fileout)
        write(nout,'(a,a,a,i9)')
     .  '# Variacao ',trim(nameres)
     .   ,' no tempo nos pontos de integracao no elemento',el
      else
c ... reabre o arquivo e add uma nova linha      
        fileout = name(prename,el,0,code)
        open(unit=nout,file=fileout,status ='old',access='append') 
      endif
c .....................................................................
c
c ...
      write(nout,'(i9,e18.10)',advance='no')istep,t
      do i = 1 , npi
        write(nout,'(i3)',advance='no')i
        do j = 1, ndf
          write(nout,'(1x,e18.10)',advance='no')u(j,i,el)
        enddo
      enddo
c .....................................................................
c
c ...
      close(nout)
      return
      end
c *********************************************************************