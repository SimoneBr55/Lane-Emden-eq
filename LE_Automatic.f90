PROGRAM LaneEmden
  IMPLICIT NONE
  !!---VARIABILI DI INPUT---!!
  REAL*8, ALLOCATABLE :: xi(:), theta(:), eta(:)
  INTEGER :: n, iterMax, pts
  REAL*8 :: h
  !!------------------------!!

  !!---VARIABILI DI LAVORO---!!
  INTEGER :: i, j
  REAL*8 :: xi_temp, theta_temp, eta_temp
  REAL*8 :: toll
  INTEGER, DIMENSION(:) :: poli(3)
  INTEGER :: pol, passo
  !REAL, DIMENSION(:,:) :: profPos(1,3), profNeg(1,3)
  REAL*8 :: xiL, xiR, thetaL, thetaR, etaL, etaR
  !REAL*8 :: f1, f2
  !!-------------------------!!

  !!--------------------------!!
  REAL*8, EXTERNAL :: pressure, density, analyticTheta, nPhi, aPhi, R0, massaTotale, poliConst, massaTotaleParam

  !!---VARIABILI DI OUTPUT---!!
  REAL*8 :: xi_0, theta_0, eta_0, R_0, K, M_Ch, M_Ch_Sol
  REAL*8, ALLOCATABLE :: prPrN(:,:), prPrA(:,:), prDenN(:,:), prDenA(:,:)
  REAL*8, ALLOCATABLE :: prMassAn(:,:), prMassNum(:,:)
  REAL*8, ALLOCATABLE :: diffPress(:,:), diffDen(:,:), diffMass(:,:)
  REAL*8 :: massa, massaPar
  !!-------------------------!!

  !!-----AZZERAMENTO-----!!
  !!---------------------!!

  !!-----ALLOCAMENTO-----!!
  !!---------------------!!

  PRINT*, "QUESTA È LA VERSIONE AUTOMATICA DELLO SCRIPT DELL'ESERCIZIO 2"
  PRINT*, "IL CODICE PROCEDE AUTOMATICAMENTE AD ANALIZZARE I SINGOLI PUNTI DELL'ESERCIZIO"
  PRINT*, "E A STAMPARE SU TERMINALE E SU FILE I VALORI NECESSARI ALLA RELAZIONE"
  PRINT*, ""

  poli(1) = 0
  poli(2) = 1
  poli(3) = 3

  DO pol=1,3

     n = poli(pol)


  DO passo=1,2
     iterMax = 1000000   ! Esagero con il numero di iterazioni massime, tanto comunque il codice interessato si ferma al compimento

     ALLOCATE(xi(iterMax), theta(iterMax), eta(iterMax))

     IF(passo == 1) THEN
        h = 0.001
     ELSEIF(passo == 2) THEN
        h = 0.0001
     ENDIF

  ! Scelgo punti iniziali

  xi(1) = TINY(h) ! Scelgo un valore di partenza appena superiore allo zero per non dover effettuare una divisione impossibile
  theta(1) = 1.
  eta(1) = 0.

  !Effettuo integrazione numerica con RungeKutta al quarto ordine
  i = 1
  DO WHILE ((theta(i) .GT. 0.) .AND. (i .LE. iterMax))
     CALL RungeKutta4(xi(i), theta(i), eta(i), n, h, xi(i+1), theta(i+1), eta(i+1))
     i = i + 1
  ENDDO
  pts = i-1

  !Parte rimossa per interazione con utente
  !IF((i .GE. iterMax) .OR. (theta(pts)*theta(pts+1) .GE. 0)) THEN
  !   PRINT*, "L'algoritmo di Runge-Kutta non ha raggiunto il profilo negativo."
  !   PRINT*, "Riprovare con piÃ¹ iterazioni o con meno risoluzione."
  !   DEALLOCATE(xi, theta, eta)
  !   GOTO 555
  !ENDIF
  
  !Scrivo su file i valori interpolati
  IF(pol == 1) THEN
     IF(passo == 1) THEN
        OPEN(77,file="puntiInterpolati_n0_h1.txt")
     ELSEIF(passo == 2) THEN
        OPEN(77,file="puntiInterpolati_n0_h2.txt")
     ENDIF
  ELSEIF(pol == 2) THEN
     IF(passo == 1) THEN
        OPEN(77,file="puntiInterpolati_n1_h1.txt")
     ELSEIF(passo == 2) THEN
        OPEN(77,file="puntiInterpolati_n1_h2.txt")
     ENDIF
  ELSEIF(pol == 3) THEN
     IF(passo == 1) THEN
        OPEN(77,file="puntiInterpolati_n3_h1.txt")
     ELSEIF(passo == 2) THEN
        OPEN(77,file="puntiInterpolati_n3_h2.txt")
     ENDIF
  ENDIF

  DO j=1,pts+1
     WRITE(77,*) j, xi(j), theta(j), eta(j)
  ENDDO

  CLOSE(77, status= 'keep')
  
  ! A questo punto, avrà i valori di passaggio tra profilo positivo e profilo negativo. Questi valori sono immagazzinati nelle righe di indice pts e pts+1.

  xiL = xi(pts)
  xiR = xi(pts+1)
  thetaL = theta(pts)
  thetaR = theta(pts+1)
  etaL = eta(pts)
  etaR = eta(pts+1)

 ! "ORA CERCO IL VALORE XI_0"
 ! "UTILIZZO UN ALGORITMO DI BISEZIONE CHE SFRUTTI IL METODO DI RUNGE-KUTTA"
 ! "CERCO LO ZERO CON UNA TOLLERANZA DI 8 CIFRE SIGNIFICATIVE"
  toll = 10**(-8.)

  CALL bisezioneRK(xiL, xiR, thetaL, thetaR, etaL, etaR, toll,n, xi_0, theta_0, eta_0)

  PRINT*, " n = ", n
  PRINT*, " h = ", h
  PRINT*, " xi_0 = ", xi_0
  IF(n == 0) THEN
     PRINT*, "La soluzione analitica per n =", n, " è ", SQRT(DBLE(6))
     PRINT*, "La differenza è ", ABS(xi_0 - SQRT(DBLE(6)))
  ELSEIF(n == 1) THEN
     PRINT*, "La soluzione analitica per n =", n, " è ", ACOS(DBLE(-1))
     PRINT*, "La differenza è ", ABS(xi_0 - ACOS((DBLE(-1))))
  ELSE
     PRINT*, "Per n = ", n, " non è possibile una soluzione analitica"
  ENDIF
  PRINT*, ""
  !PRINT*, xi_0, theta_0, Acos(-1.), SQRT(6.)

  !!___CALCOLA L'ERRORE CON SOLUZIONE ANALITICA DISPONIBILE___!!
  !! Userò solamente i risultati con passo h_2
  IF(pol == 1) THEN
     IF(passo == 1) THEN
        !Spreco tempo, ma lascio la retrocompatibilità
        !OPEN(86,file="prPressione_n0_h1.txt")
        !OPEN(88,file="prDensity_n0_h1.txt")
        !OPEN(87, file="prPrAnalitico_n0_h1.txt")
        !OPEN(89,file="prDenAnalitico_n0_h1.txt")
        !OPEN(90,file="diffPress_n0_h1.txt")
        !OPEN(91,file="diffDen_n0_h1.txt")
     ELSEIF(passo == 2) THEN
        OPEN(86,file="prPressione_n0_h2.txt")
        OPEN(88,file="prDensity_n0_h2.txt")
        OPEN(87,file="prPrAnalitico_n0_h2.txt")
        OPEN(89,file="prDenAnalitico_n0_h2.txt")
        OPEN(90,file="diffPress_n0_h2.txt")
        OPEN(91,file="diffDen_n0_h2.txt")
     ENDIF
  ELSEIF(pol == 2) THEN
     IF(passo == 1) THEN
        !Spreco tempo, ma lascio la retrocompatibilità
        !OPEN(86,file="prPressione_n1_h1.txt")
        !OPEN(88,file="prDensity_n1_h1.txt")
        !OPEN(87,file="prPrAnalitico_n1_h1.txt")
        !OPEN(89,file="prDenAnalitico_n1_h1.txt")
     ELSEIF(passo == 2) THEN
        OPEN(86,file="prPressione_n1_h2.txt")
        OPEN(88,file="prDensity_n1_h2.txt")
        OPEN(87,file="prPrAnalitico_n1_h2.txt")
        OPEN(89,file="prDenAnalitico_n1_h2.txt")
        OPEN(90,file="diffPress_n1_h2.txt")
        OPEN(91,file="diffDen_n1_h2.txt")
     ENDIF
  ELSEIF(pol == 3) THEN
     IF(passo == 1) THEN
        !Spreco tempo, ma lascio la retrocompatibilità
        !OPEN(86,file="prPressione_n3_h1.txt")
        !OPEN(88,file="prDensity_n3_h1.txt")
     ELSEIF(passo == 2) THEN
        OPEN(86,file="prPressione_n3_h2.txt")
        OPEN(88,file="prDensity_n3_h2.txt")
     ENDIF
  ENDIF
  
  IF(passo == 1) THEN
     GOTO 333 ! Vado a termine ciclo
  ENDIF
  ALLOCATE(prPrN(pts, 2))
  IF(n .NE. 3) THEN
     ALLOCATE(prPrA(pts,2), diffPress(pts,2)) !Il caso analitico posso farlo solo se n non è 3
  ENDIF
  
  DO i=1,pts-1
     prPrN(i,1) = xi(i)
     prPrN(i,2) = pressure(theta(i), n) !Calcolo tramite funzione
     WRITE(86,*) prPrN(i,1), prPrN(i,2)
     IF(n .NE. 3) THEN
        prPrA(i,1) = xi(i)
        prPrA(i,2) = pressure(analyticTheta(xi(i), n),n) ! Calcolo tramite funzione di funzione
        WRITE(87,*)  prPrA(i,1), prPrA(i,2)
        diffPress(i,1) = xi(i)
        diffPress(i,2) = ABS(prPrN(i,2) - prPrA(i,2))
        WRITE(90,*) diffPress(i,1), diffPress(i,2)
     ENDIF
  ENDDO
  prPrN(pts,1) = xi_0
  prPrN(pts,2) = pressure(theta_0, n)
  WRITE(86,*) prPrN(pts,1), prPrN(pts,2)
  PRINT*, "Profilo di Pressione in Integrazione Numerica Scritto"
  IF(n .NE. 3) THEN
     prPrA(pts,1) = xi_0
     prPrA(pts,2) = pressure(analyticTheta(xi_0, n), n)
     WRITE(87,*) prPrA(pts,1), prPrA(pts,2)
     diffPress(pts,1) = xi_0
     diffPress(pts,2) = ABS(prPrN(pts,2) - prPrA(pts,2))
     WRITE(90,*) diffPress(pts,1), diffPress(pts,2)
     PRINT*, "Profilo di Pressione in Integrazione Analitica Scritto"
  ELSE
     PRINT*, "No Profilo Analitico di Pressione"
  ENDIF

  ALLOCATE(prDenN(pts,2))
  IF(n .NE. 3) THEN
     ALLOCATE(prDenA(pts,2), diffDen(pts,2))
  ENDIF

  !Ripeto per la densità
  DO i=1, pts-1
     prDenN(i,1) = xi(i)
     prDenN(i,2) = density(theta(i), n)
     WRITE(88,*) prDenN(i,1), prDenN(i,2)
     IF(n .NE. 3) THEN
        prDenA(i,1) = xi(i)
        prDenA(i,2) = density(analyticTheta(xi(i), n), n)
        WRITE(89,*) prDenA(i,1), prDenA(i,2)
        diffDen(i,1) = xi(i)
        diffDen(i,2) = ABS(prDenN(i,2) - prDenA(i,2))
        WRITE(91,*) diffDen(i,1), diffDen(i,2)
     ENDIF
  ENDDO
  prDenN(pts,1) = xi_0
  prDenN(pts,2) = density(theta_0, n)
  WRITE(88,*) prDenN(pts,1), prDenN(pts,2)
  PRINT*, "Profilo di Densità in Integrazione Numerica Scritto"
  IF(n .NE. 3) THEN
     prDenA(pts,1) = xi_0
     prDenA(pts,2) = density(analyticTheta(xi_0, n), n)
     WRITE(89,*) prDenA(pts,1), prDenA(pts,2)
     diffDen(pts,1) = xi_0
     diffDen(pts,2) = ABS(prDenN(pts,2) - prDenA(pts,2))
     WRITE(91,*) diffDen(pts,1), diffDen(pts,2)
     PRINT*, "Profilo di Densità in integrazione Analitica Scritto"
  ELSE
     PRINT*, "No Profilo Analitico di Densità"
  ENDIF

  !PROVO A CALCOLARE IL PROFILO DI MASSA ANALITICO

  ALLOCATE(prMassNum(pts,2))
  IF(n .NE. 3) THEN
     ALLOCATE(prMassAn(pts,2), diffMass(pts,2))
  ENDIF
  IF(n == 0) THEN
     OPEN(92,file="prMassUniNum_n0.txt")
     OPEN(93,file="prMassUniAn_n0.txt")
     OPEN(94,file="diffMassUni_n0.txt")
  ELSEIF(n == 1) THEN
     OPEN(92,file="prMassUniNum_n1.txt")
     OPEN(93,file="prMassUniAn_n1.txt")
     OPEN(94,file="diffMassUni_n1.txt")
  ELSE
     OPEN(92,file="rMassUniNum_n3.txt") !caso n=3
  ENDIF

  DO i=1,pts-1
     prMassNum(i,1) = xi(i)
     prMassNum(i,2) = nPhi(xi(i), eta(i))
     WRITE(92,*) prMassNum(i,1), prMassNum(i,2)
     IF(n .NE. 3) THEN
        prMassAn(i,1) = xi(i)
        prMassAn(i,2) = aPhi(xi(i), n)
        WRITE(93,*) prMassAn(i,1), prMassAn(i,2)
        diffMass(i,1) = xi(i)
        diffMass(i,2) = ABS(prMassNum(i,2) - prMassAn(i,2))
        WRITE(94,*) diffMass(i,1), diffMass(i,2)
     ENDIF
  ENDDO
  prMassNum(pts,1) = xi_0
  prMassNum(pts,2) = nPhi(xi_0, eta_0)
  IF(n .NE. 3) THEN
     prMassAn(pts,1) = xi_0
     prMassAn(pts,2) = aPhi(xi_0, n)
     WRITE(93,*) prMassAn(pts,1), prMassAn(pts,2)
     diffMass(pts,1) = xi_0
     diffMass(pts,2) = ABS(prMassNum(pts,2) - prMassAn(pts,2))
     WRITE(94,*) diffMass(pts,1), diffMass(pts,2)
  ENDIF

  ! Calcolo R_0

  R_0 = R0(xi_0,n)
  PRINT*, R_0

  IF(n == 3) THEN
     PRINT*, "MASSA TOTALE IN FUNZIONE DI rho_c e K"
     PRINT*, "vedere relazione"
     massaPar = massaTotaleParam(n, xi_0, eta_0)
     PRINT*, massaPar ! Gli va moltiplicato K^(3/2), che determino come passaggio successivo
     PRINT*, "COSTANTE DELL'EQUAZIONE DI STATO POLITROPICA"
     K = poliConst()
     PRINT*, K
     M_ch = massaTotale(xi_0, eta_0, K)
     PRINT*, M_ch
     M_ch_Sol = M_ch/(DBLE(1.99)*(DBLE(10**33.)))
     PRINT*, M_ch_Sol
  ENDIF
  
  

  

  
  
  
  

  DO i=86,94
     CLOSE(i, status='keep')
  ENDDO
  IF(n .NE. 3) THEN
     DEALLOCATE(prDenA, prPrA, diffPress, diffDen)
     DEALLOCATE(diffMass, prMassAn)
  ENDIF
  DEALLOCATE(prDenN, prPrN, prMassNum)
333 CONTINUE
  DEALLOCATE(xi, theta, eta)
ENDDO
ENDDO

  
  
  


  !!-----DEALLOCAMENTO-----!!
  !!-----------------------!!

END PROGRAM LaneEmden

REAL*8 FUNCTION poliConst()
  IMPLICIT NONE
  REAL*8, PARAMETER :: h = 6.6260755E-27, c = 2.99792458E+10, pi = ACOS(DBLE(-1)), m_p = 1.6726231E-24
  poliConst = ((h*c)/DBLE(8))*(DBLE(3)/pi)**(DBLE(1)/DBLE(3))*(DBLE(1)/(DBLE(2)*m_p))**(DBLE(4)/DBLE(3))
END FUNCTION poliConst

REAL*8 FUNCTION massaTotaleParam(n, x, z)
  IMPLICIT NONE
  REAL*8, intent(in) :: x, z
  INTEGER, intent(in) :: n
  REAL*8, PARAMETER :: pi = ACOS(DBLE(-1)), G = 6.67408E-8
  massaTotaleParam = DBLE(4)*pi*(DBLE(n+1)/(DBLE(4)*pi*G))**(DBLE(3)/DBLE(2))*(-x**2 * z)
END FUNCTION massaTotaleParam

REAL*8 FUNCTION massaTotale(x,z,K)
  IMPLICIT NONE
  REAL*8, intent(in) :: x, z, K
  REAL*8, PARAMETER :: pi = ACOS(DBLE(-1)), G = 6.67408E-8
  massaTotale = K**(DBLE(3)/DBLE(2)) * DBLE(4)*pi*(DBLE(1)/(pi*G))**(DBLE(3)/DBLE(2))*(-x**2 * z)
END FUNCTION massaTotale
  
  

REAL*8 FUNCTION R0(x,n)
  IMPLICIT NONE
  REAL*8, intent(in) :: x
  INTEGER, intent(in) :: n
  REAL*8, PARAMETER :: pi = ACOS(DBLE(-1)), G = 6.67408E-8
  PRINT*, "G = ", G, "pi=", pi
  R0 = SQRT((DBLE(n+1))/(DBLE(4)*pi*G))*x !Va moltiplicata una f(rho_0, K)| su n
END FUNCTION R0
REAL*8 FUNCTION nPhi(x,z)
  IMPLICIT NONE
  REAL*8, intent(in) :: x, z
  nPhi = -z*(x**2)
END FUNCTION nPhi

REAL*8 FUNCTION aPhi(x,n)
  IMPLICIT NONE
  REAL*8, intent(in) :: x
  INTEGER, intent(in) :: n
  REAL*8 :: dy
  IF(n == 0) THEN
     dy = -(1./3.)*x
  ELSEIF(n == 1) THEN
     dy = (x*COS(x) - SIN(x))/(x**2)
  ENDIF
  aPhi = -(x**2)*dy
END FUNCTION aPhi

REAL*8 FUNCTION analyticTheta(x, n)
  IMPLICIT NONE
  REAL*8, intent(in) :: x
  INTEGER, intent(in) :: n
  IF(n == 0) THEN
     analyticTheta = 1. -(x**2)/6.
  ELSEIF(n == 1) THEN
     analyticTheta = (sin(x))/(x)
  ENDIF
END FUNCTION analyticTheta

REAL*8 FUNCTION density(y,n)
  IMPLICIT NONE
  REAL*8, intent(in) :: y
  INTEGER, intent(in) :: n
  density = y**n
END FUNCTION density

REAL*8 FUNCTION pressure(y,n)
  IMPLICIT NONE
  REAL*8, intent(in) :: y
  INTEGER, intent(in) :: n
  pressure = y**(n+1)
END FUNCTION pressure
  
REAL*8 FUNCTION f1(x,y,z)
  IMPLICIT NONE
  REAL*8, intent(in) :: x,y,z
  !REAL(DP), intent(out) :: df1
  f1 = z
END FUNCTION f1
  
REAL*8 FUNCTION f2(x,y,z,n)
  IMPLICIT NONE
  REAL*8, intent(in) :: x,y,z
  !REAL*8, intent(out) :: df2
  INTEGER, intent(in) :: n
  IF(y == 1.) THEN
     f2 = -(1.-(x**2./6.)+(n/120.)*x**4.)**n - (2./x)*z
  ELSE
     f2 = -y**REAL(n) - (2./x)*z
  ENDIF
END FUNCTION f2

SUBROUTINE RungeKutta4(x_0, y_0, z_0, n, t_step, x_1, y_1, z_1)
  IMPLICIT NONE
  REAL*8, intent(in) :: x_0, y_0, z_0
  REAL*8, intent(in) :: t_step
  INTEGER, intent(in) :: n
  REAL*8, intent(out) :: x_1, y_1, z_1
  REAL*8, EXTERNAL :: f1, f2
  !REAL*8 :: f1,f2
  ! k_* sono le funzioni relative alla funzione f1
  ! l_* sono le funzioni relative alla funzione f2
  REAL*8 :: k_1, l_1, k_2, l_2, k_3, l_3, k_4, l_4
  
  !!FUNZIONI AUSILIARIE
  k_1 = f1(x_0, y_0, z_0)
  l_1 = f2(x_0, y_0, z_0, n)
  
  k_2 = f1(x_0+(0.5)*t_step, y_0+(0.5)*t_step*k_1, z_0+(0.5)*t_step*l_1)
  l_2 = f2(x_0+(0.5)*t_step, y_0+(0.5)*t_step*k_1, z_0+(0.5)*t_step*l_1, n)
  
  k_3 = f1(x_0+(0.5)*t_step, y_0+(0.5)*t_step*k_2, z_0+(0.5)*t_step*l_2)
  l_3 = f2(x_0+(0.5)*t_step, y_0+(0.5)*t_step*k_2, z_0+(0.5)*t_step*l_2, n)
  
  k_4 = f1(x_0+t_step, y_0+t_step*k_3, z_0+t_step*l_3)
  l_4 = f2(x_0+t_step, y_0+t_step*k_3, z_0+t_step*l_3, n)
  
  !!CALCOLO PUNTI
  x_1 = x_0 + t_step
  y_1 = y_0 + (t_step/6.)*(k_1 + (2.)*k_2 + (2.)*k_3 + k_4)
  z_1 = z_0 + (t_step/6.)*(l_1 + (2.)*l_2 + (2.)*l_3 + l_4)
END SUBROUTINE RungeKutta4

SUBROUTINE bisezioneRK(xl_c, xr_c, yl_c, yr_c, zl_c, zr_c, toll, n, x_0, y_0, z_0)
  IMPLICIT NONE
  REAL*8, intent(in) :: xl_c, xr_c, yl_c, yr_c, zl_c, zr_c
  REAL*8 :: xl, xr, yl, yr, zl, zr
  REAL*8 :: temp_x, temp_y, temp_z
  REAL*8, intent(in) :: toll
  REAL*8 :: intervallo
  REAL*8 :: t_med
  INTEGER, intent(in) :: n
  REAL*8, intent(out) :: x_0, y_0, z_0
  
  xl = xl_c
  xr = xr_c
  yl = yl_c
  yr = yr_c
  zl = zl_c
  zr = zr_c
  
111 CONTINUE
  t_med = (xr - xl)/2.
  intervallo = ABS(xr-xl)
  
  CALL RungeKutta4(xl, yl, zl, n, t_med, temp_x, temp_y, temp_z) ! Mi sposto con uno step che mi porta al punto medio tra le due rilevazioni.
  
  IF(intervallo .LE. toll) THEN ! SE il valore che ho della y Ã¨ uguale a zero entro la tolleranza ALLORA
     GOTO 112 ! ESCO dal circolo
  ENDIF
  IF(temp_y .GT. 0.) THEN !Se il nuovo valore Ã¨ maggiore di zero
     !! I valori trovati sono ancora alla sinistra dello zero
     xl = temp_x
     yl = temp_y
     zl = temp_z
     GOTO 111
  ELSEIF(temp_y .LT. 0.) THEN !Se il nuovo valore Ã¨ minore di zero
     !! Ho superato lo zero, con il mio t_med e devo tornare indietro
     xr = temp_x
     yr = temp_y
     zr = temp_z
     GOTO 111
  ENDIF
112 CONTINUE ! Se sono arrivato qui, allora il valore che ho trovato per la y Ã¨ uguale a zero, entro la tolleranza.
  ! Di conseguenza, i temp_* trovati corrispondono a quelli dello x_0 ricercato.
  PRINT*, "XI_0 trovato (entro la tolleranza", toll,")"
  x_0 = temp_x
  y_0 = temp_y
  z_0 = temp_z
END SUBROUTINE bisezioneRK



SUBROUTINE trapez_custom(a, b, polInd, pts, int)
  IMPLICIT NONE
  INTEGER, intent(in) :: polInd, pts
  REAL*8, intent(in) :: a, b
  REAL*8, intent(out) :: int
  REAL*8, EXTERNAL :: integranda_0, integranda_1
  REAL*8 :: delta, area, x1, x2
  INTEGER :: i, j, n
 
  n=pts
  delta = (b-a)/n
  int = 0.
  DO i=1,n
     x1 = a+(i-1)*delta
     x2 = a+i*delta
     IF(polInd == 0) THEN
        area = delta*(integranda_0(x1) + integranda_0(x2))*0.5
     ELSEIF(polInd == 1) THEN
        area = delta*(integranda_1(x1) + integranda_1(x2))*0.5
     ENDIF
     int = int + area
  ENDDO
END SUBROUTINE trapez_custom
  

  
REAL*8 FUNCTION integranda_0(x)
  IMPLICIT NONE
  REAL*8 :: x
  REAL*8, EXTERNAL :: analyticTheta
  integranda_0=(x**2) 
END FUNCTION integranda_0

REAL*8 FUNCTION integranda_1(x)
  IMPLICIT NONE
  REAL*8 :: x
  REAL*8, EXTERNAL :: analyticTheta
  integranda_1=analyticTheta(x,1)*(x**2)
END FUNCTION integranda_1

REAL*8 FUNCTION int_try(x)
  REAL*8 :: x
  int_try = x**2
END FUNCTION int_try
