C     PROYECTO D3

C     EN ESTE PROGRAMA GENERO LA DISTRIBUCION DE EMBRIONES, TENIENDO EN CUENTA
C     QUE ELLOS SE COMEN TODO EL MATERIAL DE SUS ZONAS DE ALIMENTACION. ESTA
C     SUPOSICION ES LOGICA PARA LOS EMBRIONES QUE EVOLUCION EN LA PARTE INTERNA
C     DE UN SISTEMA, YA QUE SUS ZONAS DE ALIMENTACION SON MUY ANGOSTAS.
C     PARA SIMULAR LA DISTRIBUCION DE EMBRIONES EN LA PARTE EXTERNA ASUMIMOS
C     QUE LOS MISMOS HAN LOGRADO ACRETAR EL 50 % DEL MATERIAL EXISTENTE EN
C     SUS ZONAS DE ALIMENTACION. eSTO ESTA ACORDE CON LOS RESULTADOS OBTENIDOS
C     EN PAPERS ANTERIORES USANDO EL MODELO SEMIANALTICO.


      implicit real*8(a-h,o-z)
      real*4 ran2
      integer semilla
      character(5) argumento

      
C     CONSTANTES Y SEMILLA PARA LA RAN2
     
      call getarg(1,argumento)     ! Para las semillas  
      read(argumento,'(i5)') iarg  ! agregar estas
      semilla = iarg               ! tres lineas 
      
      pi = 4.d0*datan(1.d0)
      ua = 1.4959787D+13          ! 1 UA en cm
      am_sol = 1.9891D+33         ! MASA DEL SOL en gr
      am_tierra = 5.9722D+27      ! MASA DE LA TIERRA en gr
      am_jupiter = 3.d0*9.54791938424326609E-04 ! MASA DE JUP en Masas Sol.
      hielo1 = 1.d0               ! 
      hielo2 = 2.d0               ! salto en la linea del hielo
!     hielo2 = 1.05d0             ! salto en la linea del hielo 5% (SIN SALTO)
      dens1 = 3.0d0               ! densidades de los embriones y plan.
      dens2 = 1.5d0               ! antes y despues de aice en gr/cm3
!     dens2 = 3.0d0               ! antes y despues de aice en gr/cm3 (SIN SALTO)
      aice = 2.7d0              ! posicion de la linea de hielo en UA
    
      z0 = 0.0153d0               ! factor de Lodders (2009). Abun.solar
      emet = 0.                   ! METALICIDAD ASUMIDA PARA LA ESTRELLA
      AMSTAR = 1.D0               ! MASA ASUMIDA PARA LA ESTRELLA EN MS

C     COSAS QUE PUEDEN IR CAMBIANDO      
      gama = 1.0d0                  ! exponente del perfil de densidad
      am = 0.01d0                   ! masa total del disco en masas sol.
      am_gr = am*am_sol            ! masa total del disco en gr
      ac = 25.d0                   ! radio caract. del disco en UA
      ac_cm = ac*ua                ! radio caract. del disco en cm 
      
C     LIMITES DEL DISCO EN UA
      a1 = 0.5d0
      aborde1 = 2.7d0
      aborde2 = 2.7d0
      a2 = 5.d0
      
C     CALCULO SIGMA_GAS_0
      anum = (2.-gama)*am_gr
      den = 2*pi*(ac_cm**2.d0)
      sigma_g_0 = anum/den  ! en gr/cm**2

C     CALCULO SIGMA_S_0
      sigma_s_0 = sigma_g_0*z0*(10.d0**(emet))  !en gr/cm*2
      
      open(1,file='embriones')
      open(2,file='big.in')
      
      write(2,'(")O+_06 Big-body initial data  (WARNING: Do not delete
     *this line!!)")')
      write(2,'(") Lines beginning with `)` are ignored.")')
      write(2,'(")--------------------------------------------------")')
      write(2,'("style (Cartesian, Asteroidal, Cometary) = Asteroidal")'
     *)
      write(2,'("epoch (in days) = 0.")')
      write(2,'(")--------------------------------------------------")')
     
100   format(a80)   

c     --------------------------------------------
c     AHORA PONGO LOS EMBRIONES
C     --------------------------------------------
      a = a1
      i = 0
      amamint = 0.d0
      amamext = 0.d0
      
C     PONEMOS LOS EMBRIONES      
      do while(a.le.aborde1)
      i = i + 1
      a_cm = a*ua     ! semieje en cm
C     CALCULAMOS LA DENSIDAD SUPERFICIAL EN GR/CM2
      coc = a/ac
      term1 = coc**(-gama)
      term2 = exp(-(coc**(2.d0-gama)))
c     p va a representar cuanta masa va a los embriones respecto
c     de la masa total del sistema
      p = 1.d0
      delta = 5.d0 + ((ran2(semilla)*5.d0))  !separacion entre embriones
      sigma = sigma_s_0*hielo1*term1*term2 !sigma_solidos en "a"
      
c     A partir de aqui defino la masa del embrion para un dado "a"       
       uno = 2.d0*pi*delta*(a_cm**2.d0)*sigma*p
       uno = uno**3.d0
       dos = 1.5d0*am_sol
       tres = uno/dos
       tres = sqrt(tres) ! esta es la masa en gr
       
       amamint = amamint + (tres/am_tierra)
       
       write(1,*)a,tres/am_tierra ! para un a (en UA), Masa_embrion
       
       am2 = tres/am_sol !escribimos la masa del embrion en masas 
                         !solares ya que el Mercury trabaja con estas 
                         !unidades
       
       write(2,*)i,"m=",am2,"r=",1.d0,"d=",dens1
       e = 0.02d0*ran2(semilla)   !excentricidades menores a 0.02
       ang1 = 0.5d0*ran2(semilla) ! inclinaciones menores a 0.5 grados
       ang2 = 360.d0*ran2(semilla)
       ang3 = 360.d0*ran2(semilla)
       ang4 = 360.d0*ran2(semilla)
       write(2,*)a,e,ang1
       write(2,*)ang2,ang3,ang4
       write(2,*)0.d0,0.d0,0.d0
       
C     CALCULAMOS EL RADIO DE HILL MUTUO
       rh = (2.d0*tres)/(3.d0*(am_sol))
       rh = a*(rh**(1.d0/3.d0))
       
       a = a + (delta*rh)
      enddo      
           
      
      a = 2.7d0
      do while(a.le.a2)
      i = i + 1
      a_cm = a*ua     ! semieje en cm
C     CALCULAMOS LA DENSIDAD SUPERFICIAL EN GR/CM2
      coc = a/ac
      term1 = coc**(-gama)
      term2 = exp(-(coc**(2.d0-gama)))
c     p va a representar cuanta masa va a los embriones respecto
c     de la masa total del sistema
c      pend = (1.-0.1)/(3.5-10.)
c      ord = 1 - pend*3.5
c      p = pend*a + ord
       p = 1
      delta = 5.d0 + ((ran2(semilla)*5.d0))  !separacion entre embriones
      sigma = sigma_s_0*hielo2*term1*term2 !sigma_solidos en "a"
      
       
c     A partir de aqui defino la masa del embrion para un dado "a"       
       uno = 2.d0*pi*delta*(a_cm**2.d0)*sigma*p
       uno = uno**3.d0
       dos = 1.5d0*(am_sol)
       tres = uno/dos
       tres = sqrt(tres) ! esta es la masa en gr
       
       amamext = amamext + (tres/am_tierra)
       
       write(1,*)a,tres/am_tierra ! para un a (en UA), Masa_embrion
       
       am2 = tres/am_sol !escribimos la masa del embrion en masas 
                         !solares ya que el Mercury trabaja con estas 
                         !unidades
       
       write(2,*)i,"m=",am2,"r=",1.d0,"d=",dens2
       e = 0.02d0*ran2(semilla)   !excentricidades menores a 0.02
       ang1 = 0.5d0*ran2(semilla) ! inclinaciones menores a 0.5 grados
       ang2 = 360.d0*ran2(semilla)
       ang3 = 360.d0*ran2(semilla)
       ang4 = 360.d0*ran2(semilla)
       write(2,*)a,e,ang1
       write(2,*)ang2,ang3,ang4
       write(2,*)0.d0,0.d0,0.d0
       
C     CALCULAMOS EL RADIO DE HILL MUTUO
       rh = (2.d0*tres)/(3.d0*(am_sol))
       rh = a*(rh**(1.d0/3.d0))
       
       a = a + (delta*rh)
      enddo
      write(*,*)amamint,amamext
      
      end

      FUNCTION ran2(idum)
C     ----------------------------------------------------------
C     Generador de números aleatorios de largo período (>2x10^18)
C     El valor devuelto es un valor de una muestra aleatoria
C     uniforme en el rango [0,1]
C     El valor de inicialización de idum debe ser un entero.
C     negativo.
C     -----------------------------------------------------------<----->
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *     IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *     NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
      idum=max(-idum,1)
      idum2=idum
      do 11 j=NTAB+8,1,-1
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      if (j.le.NTAB) iv(j)=idum
11      continue
      iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
C   -----------------------------------------------------
C  (C) Copr. 1986-92 Numerical Recipes Software w45r1&w.#.
