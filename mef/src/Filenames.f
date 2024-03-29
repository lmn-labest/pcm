      character*80 function name(NomeArqDados,NumArq1,NumArq2,code)
c *********************************************************************
c * Data de criacao    : 12/12/2015                                   *
c * Data de modificaco : 20/03/2019                                   * 
c * ----------------------------------------------------------------- *
c * NAME: nomes de aquivos                                            *
c * ----------------------------------------------------------------- *
c * Parametros de entrada:                                            *
c * ------------------------------------------------------------------*
c * NumArq       - numero do arquivo                                  *
c * code         - codigo de instrucao                                *
c * ----------------------------------------------------------------- *    
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * NomeArqDados - nome do arquivo                                    *
c * ----------------------------------------------------------------- * 
c * OBS:                                                              *
c * ----------------------------------------------------------------- *  
c *  code:                                                            *
c *  0 -> geometria vtk                                               *
c *  1 -> geometrica vtu                                              *
c *  2 -> resultado vtk                                               *
c *  3 -> resultado vtu                                               *
c *  4 ->                                                             *
c *  5 -> condicao de contorno .vtk                                   *
c *  6 -> condicao de contorno .vtu                                   *
c *  7 -> condicao de contorno face mecanico .vtk                     *
c *  8 -> condicao de contorno face mecanico .vtu                     *
c *  9 -> condicao de contorno face hidraulico .vtk                   *
c * 10 -> condicao de contorno face hidraulico .vtu                   *
c * 12 -> arquivo do tempos medios(mpi)                               *
c * 13 -> arquivo de entrada do mpi                                   *
c * 14 -> arquivo de tempos                                           *
c * 15 -> arquivo do log do solver                                    *
c * 16 -> arquivo do log nao linear                                   *
c * 17 -> arquivo do log nao linear                                   *
c * 20 -> cor da malha .vtk                                           *
c * 21 -> cor da malha .vtu                                           *
c * 30 -> arquivos de resultados por nos                              *
c * 31 -> arquivos de tensoes total por nos                           *
c * 32 -> arquivos de tensoes efetivas de terzaghi por nos            *
c * 33 -> arquivos de tensoes efetivas de biot por nos                *
c * 34 -> arquivos de fluxo de darcy por nos                          *
c * 35 -> arquivos de tensoes totais nos pontos de integracao         *
c * 36 -> arquivos de variaveis da plasticidade nos pontos de integra *
c *       cao                                                         *
c * 37 -> arquivos de porosidada por nos                              *
c *********************************************************************
      implicit none      
      character*80 NomeArqDados,NomeArqGerado
      character*30 StrExtensao,Str1,Str2
      integer      iPonto,TamanhoNome,NumArq1,NumArq2,code
c ......................................................................
c
c ... geom .vtk
      if(code .eq. 0) then
        StrExtensao='_geo.vtk'
c ... geom .vtu
      elseif(code .eq. 1) then
        StrExtensao='_geo.vtu'
c ... res .vtk
      elseif(code .eq. 2) then
        write(StrExtensao,'( I6 )') NumArq1
        write(StrExtensao,'( A  )') adjustl(StrExtensao)  
        StrExtensao='_istep_'//trim(StrExtensao)//'.vtk'
c ... res .vtu
      elseif(code .eq. 3) then
        write(StrExtensao,'( I6 )') NumArq1
        write(StrExtensao,'( A  )') adjustl(StrExtensao)  
        StrExtensao='_istep_'//trim(StrExtensao)//'.vtu'
c ... _tetra.dat
      elseif(code .eq. 4) then
        StrExtensao='_tetra.dat'
c ... geom .vtk
      else if(code .eq. 5) then
        StrExtensao='_bc.vtk'
c ... geom .vtu
      else if(code .eq. 6) then
        StrExtensao='_bc.vtu'
c ... geom .vtk
      else if(code .eq. 7) then
        StrExtensao='_face_mec_bc.vtk'
c ... geom .vtu
      else if(code .eq. 8) then
        StrExtensao='_face_mec_bc.vtu'
c ... geom .vtk
      else if(code .eq. 9) then
        StrExtensao='_face_pres_bc.vtk'
c ... geom .vtu
      else if(code .eq. 10) then
        StrExtensao='_face_pres_bc.vtu'
c ... log de tempos
      elseif(code .eq. 12) then
        write(StrExtensao,'( I6 )') NumArq1
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_t_mean_'//trim(StrExtensao)//'.txt'
c ...
      elseif(code .eq. 13) then
        write(StrExtensao,'( I6 )') NumArq1
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_'//trim(StrExtensao)//'.dat'
c ... log de tempos
      elseif(code .eq. 14) then
        write(StrExtensao,'( I6 )') NumArq1
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_t_'//trim(StrExtensao)//'.txt'
c ... log solver
      elseif(code .eq. 15) then
        write(StrExtensao,'( I6 )') NumArq1
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_log_'//trim(StrExtensao)//'.txt'
c ... log nao linear
      elseif(code .eq. 16) then
        write(StrExtensao,'( I6 )') NumArq1
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_log_nonlinear_'//trim(StrExtensao)//'.txt'
c ... log hist solver
      elseif(code .eq. 17) then
        write(StrExtensao,'( I6 )') NumArq1
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_log_hist_'//trim(StrExtensao)//'.txt'
c ... cores            
      elseif(code .eq. 20) then
        StrExtensao='_color.vtk'
c ... cores                  
      elseif(code .eq. 21) then
        StrExtensao='_color.vtu'
c ... pu node file 
      elseif(code .eq. 30) then
        write(StrExtensao,'( I6 )') NumArq1
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_up_node_'//trim(StrExtensao)//'.txt'
c ... stress node file 
      elseif(code .eq. 31) then
        write(StrExtensao,'( I6 )') NumArq1
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_stress_node_'//trim(StrExtensao)//'.txt'
c ... terzaghi stress node file 
      elseif(code .eq. 32) then
        write(StrExtensao,'( I6 )') NumArq1
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_stressE_node_'//trim(StrExtensao)//'.txt'
c ... biot stress node file 
      elseif(code .eq. 33) then
        write(StrExtensao,'( I6 )') NumArq1
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_stressB_node_'//trim(StrExtensao)//'.txt'
c ... flux node file 
      elseif(code .eq. 34) then
        write(StrExtensao,'( I6 )') NumArq1
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_flux_node_'//trim(StrExtensao)//'.txt'
c ... stress pi file 
      elseif(code .eq. 35) then
        write(StrExtensao,'( I6 )') NumArq1
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_stress_pi_'//trim(StrExtensao)//'.txt'
c ... plastic pi file 
      elseif(code .eq. 36) then
        write(StrExtensao,'( I6 )') NumArq1
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_plastic_pi_'//trim(StrExtensao)//'.txt'
c ... porosity node file 
      elseif(code .eq. 37) then
        write(StrExtensao,'( I6 )') NumArq1
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_poro_node_'//trim(StrExtensao)//'.txt'
c ... mtx 
      elseif(code .eq. 50) then
        write(Str1,'( I6 )') NumArq1
        write(Str2,'( I6 )') NumArq2
        Str1 = adjustl(Str1)
        Str2 = adjustl(Str2)
        StrExtensao='_'//trim(Str1)//'_'//trim(Str2)
c ... kuu mtx 
      elseif(code .eq. 51) then
        write(StrExtensao,'( I6 )') NumArq1
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_kuu_'//trim(StrExtensao)
c ... kpp mtx 
      elseif(code .eq. 52) then
        write(StrExtensao,'( I6 )') NumArq1
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_kpp_'//trim(StrExtensao)
c ... pcg duplo
      elseif(code .eq. 53) then
        write(StrExtensao,'( I6 )') NumArq1
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_dpcg_'//trim(StrExtensao)//'.txt'
      endif
c ......................................................................      
      TamanhoNome = INDEX( NomeArqDados, ' '  )
      iPonto = INDEX( NomeArqdados, '.' )
      if( iPonto .EQ. 0 ) then
          NomeArqGerado = NomeArqDados(1:TamanhoNome-1) // StrExtensao
      else
          NomeArqGerado = NomeArqDados(1:iPonto-1) // StrExtensao
      endif
      name = NomeArqGerado
c ......................................................................      
      return
      end
