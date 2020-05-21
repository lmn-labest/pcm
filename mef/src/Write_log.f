c *********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 20/05/2020                                    * 
c * ------------------------------------------------------------------ *   
c * WRITE_LOG : Escrever o arquivo de log de excuacao                 *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * nnode     - numero de nos                                         *
c * numel     - numero de elementos                                   *
c * numel_nov - numero de elementos non-overllaping                   *  
c * numel_ov  - numero de elementos overllaping                       *
c * ndf       - graus de liberdade mecanico                           *
c * ndft      - graus de liberdade termico                            *
c * neq       - numero de equacoes total                              *
c * nequ      - numero de equacoes u                                  *
c * neqp      - numero de equacoes p                                  *
c * nad       - numero de elementos nao nulos fora da diag principal  *
c * nadu      - numero de coeficientes nao nulos do bloco u           *
c * nadu      - numero de coeficientes nao nulos do bloco p           *
c * nadup     - numero de coeficientes nao nulos do bloco up          *
c * nad1      - numero de elementos nao nulos di csrcr(overlaping)    *
c * neqt      - numero de equacoes do termico                         *
c * omp_elmt  - flag do openmp na fase de elemento                    *
c * nth_elmt  - numero de threads usado na fase de elemento           *
c * omp_solv  - flag do openmp na fase do solver                      *
c * nth_solv  - numero de threads usado do solver                     *
c * num_colors- numero de cores usado para colorir a malha            *
c * prename   - prefixo do nome do arquivo de saida                   *
c * nlog      - numero do arquivo de saida                            *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c * OBS:                                                              *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine write_log_file(nnode   ,numel,numel_nov ,numel_ov,ndf 
     1                         ,neq     ,nequ ,neqp 
     3                         ,nad     ,nadu ,nadp ,nadpu     ,nad1
     4                         ,omp_elmt,nth_elmt   ,omp_solv ,nth_solv
     5                         ,fporomec,num_colors ,prename  ,nlog)
      use Malloc
      implicit none
      include 'time.fi'
c ... malha
      integer nnode,numel,numel_nov,numel_ov
c ... informacoes do sistema      
      integer neq,nequ,neqp
      integer nadu,nadp,nadpu,nad1
      integer*8 nad      
c ...
      integer ndf
      logical fporomec
c ... openmp
      integer nth_elmt,nth_solv,num_colors
      logical omp_elmt,omp_solv
c ... variaveis de arquivos      
      character*80 fname,name,prename
      integer nlog
c ...
      real*8 use_work_vector,get_buffer_size
      integer buf
c .....................................................................
c
c ... abre o arquivo de logs
      fname = name(prename,0,0,14)
      open(nlog, file= fname)
      write(nlog,'(a)')"# Arquivo de log do poro mecanico"
      write(nlog,'(a)')"Tempos (seg):"
c .....................................................................
c
c ... Tempo levado na reord
      call twrite('REORD ',reordtime,1,nlog)
c
c
c ... Tempo levado na numeq
      call twrite('NUMEQ ',numeqtime,1,nlog)
c
c ... Tempo levado na struc
c
      call twrite('STRUC ',dstime,1,nlog)
c
c ...                        
c
      call twrite('VECTR ',vectime,1,nlog)
c
c ... Tempo levado na pform
c
      call twrite('ELMT  ',elmtime,1,nlog)
c
c ... Tempo levado na tform
c
      call twrite('TFORM ',tformtime,1,nlog)
c
c ... Tempo levado no precondicionador
c
      call twrite('PCOND ',precondtime,1,nlog)
c
c ... Tempo levado no solver triagonal do solver fatorado
c
      call twrite('IFSOLV',ifatsolvtime,1,nlog)
c
c ... Tempo levado no bloco diagonal solver                 
c
      call twrite('BDSOLV',prebdiagtime,1,nlog)
c
c ... Tempo levado no solver
c
      call twrite('SOLVR ',soltime,1,nlog)
c
c ... Tempo levado no matvec
c
      call twrite('MATVC ',matvectime,1,nlog)
c
c ... Tempo levado no produto interno (dot)
c
      call twrite('DOT   ',dottime,1,nlog)
c
c ... Tempo levado no envia e recebimento de dados (MPI)
c
      call twrite('SNDRC ',sendtime,1,nlog)
c
c ... Tempo levado na geracao/atualizacao do buffer de dados
c     rebidos e enviados (MPI)
c
      call twrite('OVERH ',ovhtime,1,nlog)
c
c ... Tempo levado no colormesh
c
      call twrite('COLOR ',colortime,1,nlog)
c
c ... omp
      if(omp_solv) then
c
c ... Tempo levado no matix partition
c
        call twrite('PMATRI ',pmatrixtime,1,nlog)
c
c ... Tempo levado no inicializacao do buffer
c
        call twrite('INITBUF',tinitbuffer,1,nlog)
c
c ... Tempo levado do acumolo do buffer
c
        call twrite('ACBUF  ',tacbuffer,1,nlog)
      endif   
c
c ... Tempo levado na escrita dos resultados
c
      call twrite('WRES   ',writetime,1,nlog)
c
c ... Tempo levado atualizacao da propriedades poromecanicas
c
      if(fporomec)then
        call twrite('UPPROP ',upproptime,1,nlog)
      endif
c
c ... Tempo Total
      call twrite('TOTAL ',totaltime,1,nlog)
c
c ... 
c
      write(nlog,'(a)')"Memoria:"
c
c ... Total de memoria usado no vetor de trabalho
      call twrite('ia MB    ',use_work_vector('MB'),1,nlog)
c
c ...
c
      write(nlog,'(a)')"Malha e sistema linear:"
c ... poromecanico
      if(fporomec)then
          write(nlog,'(a)')"Poromecanico:"
          call  itwrite('neq   ',neq  ,1,nlog)
          call  itwrite('nequ  ',nequ ,1,nlog)
          call  itwrite('neqp  ',neqp ,1,nlog)
          call ditwrite('nad   ',nad  ,1,nlog)
          call  itwrite('nadu  ',nadu ,1,nlog)
          call  itwrite('nadp  ',nadp ,1,nlog)
          call  itwrite('nadpu ',nadpu,1,nlog)
        endif
c .....................................................................
      call itwrite('nnode ',nnode,1,nlog)
      call itwrite('numel ',numel,1,nlog)
c .....................................................................
c
c ...
        if( ndf .gt. 0) then
c
c ... numero de coeficientes nao nulos    
c
           call ditwrite('nad   ',nad,1,nlog)
c
c ... numero de coeficieno nao nulos overlapping
c
           call itwrite('nadpu ',nadpu,1,nlog)
c
c ... numero de coeficieno nao nulos overlapping
c
           call itwrite('nadr  ',nad1,1,nlog)
        endif
c .....................................................................
c
c ... numero de nos      
c
      call itwrite('nnode ',nnode,1,nlog)
c
c ... numero de elementos
 
      call itwrite('nel   ',numel,1,nlog)
c .....................................................................
c
c ... openmp
c
      if(omp_elmt .or. omp_solv) then
        write(nlog,'(a)')"Openmp:"
c
c ... numero de cores usado         
c
        call itwrite('ncolor',num_colors,1,nlog)
c 
c ... numero de nthreads na fase do elemento
c 
        if(omp_elmt) then
          call itwrite('nth_elmt',nth_elmt,1,nlog)
        endif
c ... numero de nthreads na fase do elemento
c 
        if(omp_solv)then
          call itwrite('nth_solv',nth_solv,1,nlog)
        endif
c
c ... poro-mecanico      
        if(ndf .gt. 0 .and. omp_solv) then
          buf = neq
          call twrite('bfm MB',get_buffer_size('MB',buf),1,nlog)
        endif
      endif  
c .....................................................................
c
      close(nlog) 
      return
      end
c *********************************************************************
c
c c *********************************************************************
c ......................................................................
      subroutine twrite(text,ts,nprcs,nout)
      implicit none
      character*6 text
      real*8 ts(*)
      integer i,nprcs,nout
      write(nout,99) text,(ts(i),i=1,nprcs)
      return
   99 format(a6,128f18.6)      
      end
c ......................................................................
      subroutine itwrite(text,ts,nprcs,nout)
      implicit none
      character*6 text
      integer ts(*)
      integer i,nprcs,nout
      write(nout,99) text,(ts(i),i=1,nprcs)
      return
   99 format(a6,128i16)
      end
c ......................................................................
      subroutine ditwrite(text,ts,nprcs,nout)
      implicit none
      character*6 text
      integer*8 ts(*)
      integer i,nprcs,nout
      write(nout,99) text,(ts(i),i=1,nprcs)
      return
   99 format(a6,128i16)
      end
c **********************************************************************
