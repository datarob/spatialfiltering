'----------------------------------------------------------------------------------
'Globales und lokales Moran's I (inkl. Sattelpunktapproximation)
'EViews 5
'----------------------------------------------------------------------------------

'VERWENDUNG:
'run moransad [%0] [%1] [%2] [%3] [%4] [%5]

'Eingabeparameter:
'%0 = Reihe mit endogener Variable
'%1 = Reihe/Gruppe mit exogener(n) Variable(n) (const...Regression auf Konstante)
'%2 = Räumliche Gewichtungsmatrix
'%3 = Kodierungsverfahren: c...global standardisiert, 
'                          w...Zeilensummen standardisiert,
'                          s...Varianz stabilisierend) 
'%4 =	g... nur globales MI, 
'	   gl... globales u. lokales MI
'%5 = j... Sattelpunkapproximation durchführen
'		n... Sattelpunktapproximation nicht durchführen

mode quit	setzt Ausführungsmodus auf "quit"

'----------------------------------------------------------------------------------
'Fehlerüberprüfung:
'----------------------------------------------------------------------------------

!rG = @rows({%2})
!cG = @columns({%2})
if !rG<>!cG then
	statusline Programmabbruch: Gewichtungsmatrix nicht quadratisch
	stop
endif

'----------------------------------------------------------------------------------
'Variablen erzeugen:
'----------------------------------------------------------------------------------

tic                                'Startet Timer

stom({%0},mi_y)                    'Vektor der endogenen Variable  

series const = 1                   'Matrix der exogenen Variablen (inkl. Konstante)
if {%1} = const then               'bei Eingabe von "const" Regression auf Konstante
    group exo const
else
    group exo const {%1}
endif

stom(exo,mi_x)                     'Gruppe in Matrix unwandeln


matrix mi_g = {%2}                          'räumliche Gewichtungsmatrix    
matrix mi_g = (@transpose(mi_g)+mi_g)/2     'Symmetrie sicherstellen

%code = %3                          'Kodierungsverfahren
%GlobMI = %4                        '(nur) globales Moran's I berechnen
%Sad = %5                           'Sattelpunkapproximation (nicht) durchführen

statusline

delete const exo

'----------------------------------------------------------------------------------
'S U B R O U T I N E   Z U R   B E R E C H N U N G:
'----------------------------------------------------------------------------------

vector mi_mi                        'Objekte für Ergebnisse erzeugen
vector mi_my
vector mi_moments_exp
vector mi_moments_var
vector mi_moments_skw
vector mi_moments_kur
vector mi_z_mi
vector mi_prob_nv   
vector mi_sad
vector mi_sad_r
vector mi_sad_u
vector mi_prob_sad
table(41,4) mi_info

subroutine local calc_moran(vector mi_y,matrix mi_x,matrix mi_g,string %code, _
string %GlobMI,string %Sad,vector mi_mi,vector mi_my,vector mi_moments_exp, _
vector mi_moments_var,vector mi_moments_skw,vector mi_moments_kur,vector mi_z_mi, _
vector mi_prob_nv,vector mi_sad, vector mi_sad_r, vector mi_sad_u, _
vector mi_prob_sad,table mi_info)


!nRows = @rows(mi_y)            'bestimmt Zeilenanzahl

!nLoop = !nRows + 1             'Fallunterscheidung:
    
    if %GlobMI = "g" then	     'Globales MI		 
        !nLoop = 1
    else			     'Globales und lokales MI
            !nLoop = !nRows + 1
    endif

!k = @columns(mi_x)             'Spaltenanzahl der Matrix X
!df = !nRows - !k               'Freiheitsgrade

'----------------------------------------------------------------------------------
'Regression durchführen:
'----------------------------------------------------------------------------------

matrix mi_m = @identity(!nRows) - mi_x*@inverse((@transpose(mi_x)*mi_x)) _
					*@transpose(mi_x)			'Residualmaker
mi_my = mi_m*mi_y                           'Residuen
!denom = (@transpose(mi_my)*mi_my)(1,1)     'Nenner von Moran's I



'----------------------------------------------------------------------------------
'Schleife über globales und lokales Moran's I:
'----------------------------------------------------------------------------------

for !i = 1 to !nLoop
    statusline Iteration !i von !nLoop      'aktuelle Iteration anzeigen

    '-------------------------------------------------------------------------------
    'Gewichtungsmatrix kodieren für globales Moran's I:
    '-------------------------------------------------------------------------------
    
    if !i = 1 then
       if %code = "s" then                      's...Varianz stabilisierend
           vector(!nRows) a
               for !j = 1 to !nRows
                   a(!j) = @sqrt(@sumsq(@rowextract(mi_g,!j)))
               next
           matrix de = (@inverse(@makediagonal(a)))
           matrix vi = de*mi_g
               for !j = 1 to !nRows
                   a(!j) = @sum(@rowextract(vi,!j))
               next
           vi = !nRows/@sum(a)*vi
           vi = 0.5*(vi + @transpose(vi))       'Symmetrie sicherstellen
       endif
       
       if %code = "w" then                      'w...Zeilensummen standardisiert
           vector(!nRows) a
               for !j = 1 to !nRows
                   a(!j) = @sum(@rowextract(mi_g,!j))
               next
           matrix de = (@inverse(@makediagonal(a)))
           matrix vi = de*mi_g
           vi = 0.5*(vi + @transpose(vi))       'Symmetrie sicherstellen
       endif
       
       if %code = "c" then                      'c...global standardisiert
           vector(!nRows) a
               for !j = 1 to !nRows
                   a(!j) = @sum(@rowextract(mi_g,!j))
               next
           matrix vi = (!nRows/@sum(a))*mi_g
           vi = 0.5*(vi + @transpose(vi))       'Symmetrie sicherstellen
       endif
        
    '-------------------------------------------------------------------------------
    'Gewichtungsmatrix kodieren für lokales Moran's I:
    '-------------------------------------------------------------------------------

    else
      delete vi   
      matrix(!nRows,!nRows) vi                      'sternförmige
      rowplace(vi,@rowextract(mi_g,!i-1),!i-1)      'Gewichtungsmatrix
      colplace(vi,@columnextract(mi_g,!i-1),!i-1)                                    
        
      if %code = "s" then                           's...Varianz stabilisierend
          vi = (!nRows^2)*vi*de(!i-1,!i-1)/(2*@sum(a))
      endif
        
      if %code = "w" then                           'w...Zeilensummen standardisiert
          vi = !nRows*vi*de(!i-1,!i-1)/2
      endif
        
      if %code = "c" then                           'c...global standardisiert
          vi = (!nRows^2)*vi/(2*@sum(a))
      endif
    endif 

    '-------------------------------------------------------------------------------
    'Berechnung von Moran's I:
    '-------------------------------------------------------------------------------
    
    vector(!nLoop) mi_mi
    mi_mi(!i) = (@transpose(mi_my)*vi*mi_my)(1,1)/!denom 

    '-------------------------------------------------------------------------------
    'Eigenwerte berechnen und k Eigenwerte = 0 entfernen:
    '-------------------------------------------------------------------------------
   
    'Symmetrie sicherstellen
    sym mvm = 0.5*(@transpose((mi_m*vi*mi_m))+(mi_m*vi*mi_m)) 
    
    vector evalue =  @eigenvalues(mvm)          'Eigenwerte von MVM
    vector(!nLoop) mi_evalue_min
    mi_evalue_min(!i) = evalue(1)						'minimaler Eigenwert
    vector(!nLoop) mi_evalue_max
    mi_evalue_max(!i) = evalue(!nRows)    			'maximaler Eigenwert

    for !j = 1 to !nRows                        'sortiertes Eigenwertspektrum
        if evalue(!j) > 0.0000001 then
            !idxpos = !nRows - !j + 1
                exitloop
        endif
    next
   
    vector(!df) tau
    for !j = 1 to !idxpos
        tau(!j) = evalue(!nRows - !j + 1)
    next

    for !j = !idxpos + 1 to !df
        tau(!j) = evalue(!nRows - !k - !j + 1)
    next

    '-------------------------------------------------------------------------------
    'Momente der Verteilung:
    '-------------------------------------------------------------------------------
    
    vector(!nLoop) mi_moments_exp                        'Erwartungswert
    mi_moments_exp(!i) = @sum(tau)/!df

    vector(!nLoop) mi_moments_var                        'Varianz
    vector(!df) tau_e
    for !j = 1 to !df
        tau_e(!j) = (tau(!j) - mi_moments_exp(!i))^2
    next
    mi_moments_var(!i) = (2*@sum(tau_e))/(!df*(!df + 2))
    
    vector(!nLoop) mi_moments_skw                        'Schiefe
    for !j = 1 to !df
        tau_e(!j) = (tau(!j) - mi_moments_exp(!i))^3
    next
    mi_moments_skw(!i) = (8*@sum(tau_e))/(!df*(!df + 2)*(!df + 4))
    mi_moments_skw(!i) = mi_moments_skw(!i)/(mi_moments_var(!i)^1.5)
    
    vector(!nLoop) mi_moments_kur                        'Kurtosis
    for !j = 1 to !df
        tau_e(!j) = (tau(!j) - mi_moments_exp(!i))^4
    next    
    mi_moments_kur(!i) = (48*@sum(tau_e)+12*(@sumsq(tau - mi_moments_exp(!i))^2))/ _
								  (!df*(!df+2)*(!df+4)*(!df+6))
    mi_moments_kur(!i) = mi_moments_kur(!i)/(mi_moments_var(!i)^2)       


    '-------------------------------------------------------------------------------
    'Wahrscheinlichkeit bei Normalverteilung:
    '-------------------------------------------------------------------------------    
    
    vector(!nLoop) mi_z_mi               'z-standardisiertes MI
    mi_z_mi(!i) = (mi_mi(!i) - mi_moments_exp(!i))/@sqrt(mi_moments_var(!i))   

    vector(!nLoop) mi_prob_nv            'Wahrscheinlichkeit bei Normalverteilung
    mi_prob_nv(!i) = @cnorm((mi_mi(!i) - mi_moments_exp(!i))/ _
						    @sqrt(mi_moments_var(!i)))

    
    '-------------------------------------------------------------------------------
    'Sattelpunkapproximation:
    '-------------------------------------------------------------------------------
    
    if %Sad = "j" then
    vector(!df) taumi
    for !j = 1 to !df
        taumi(!j) = tau(!j) - mi_mi(!i)
    next

    '-------------------------------------------------------------------------------
    'Sekanten-Suchverfahren (nur für globales Moran's I):
    '-------------------------------------------------------------------------------
    
    if !i = 1 then
        !l = 1/(2*taumi(!df)) + 0.01
        !h = 1/(2*taumi(1)) - 0.01

        vector(!df) flv
        for !j = 1 to !df
            flv(!j) = taumi(!j)/(1-2*!l-taumi(!j))
        next
        !fl = @sum(flv)
        delete flv

        vector(!df) fhv
        for !j = 1 to !df
            fhv(!j) = taumi(!j)/(1-2*!h-taumi(!j))
        next
        !fh = @sum(fhv)
        delete fhv

        if !fl < 0 then
            !xl = !l
            !xh = !h
        else
        !xl = !h
        !xh = !l
        !swap = !fl
        !fl = !fh
        !fh = !swap
        endif

        !dx = !xh - !xl
        !del = 1
        !f = 1
        
        !z = 10000                               'Anzahl der möglichen Iterationen
        for !q = 1 to !z				 '(bei Bedarf erhöhen)
            !rtf = !xl + !dx*!fl/(!fl - !fh)
            vector(!df) fv
                for !j = 1 to !df
                    fv(!j) = taumi(!j)/(1-2*!rtf*taumi(!j))
                        next
                !f = @sum(fv)
                if !f < 0 then
                    !del = !xl - !rtf
                    !xl = !rtf
                    !fl = !f
                else
                    !del = !xh - !rtf
                    !xh = !rtf
                    !fh = !f
                endif
                !dx = !xh - !xl
                
					 !fehler = @abs(!f)
                statusline Sekanten-Suchverfahren für globalen Sattelpunkt: _
									 Iteration !q von maximal !z, Fehler: !fehler
                
					 if !fehler > 10000000 then    'Abbruchbedingung falls nicht konvergent
                    !rtf = 0                 '(Grund: hohe Rechenzeit)
                    exitloop
                endif

                if !fehler < 0.00001 then     'Abbruchbedingung falls konvergent
                    exitloop
                endif

                if !q = !z then               'möglicherweise zu wenig Iterationen
                    !rtf = 0
                    exitloop
                endif     
      
        next
    
        !omega = !rtf

    '-------------------------------------------------------------------------------
    'Exakte Sattelpunkt (Wurzel) für lokales MI:
    '-------------------------------------------------------------------------------

    else
        !l = tau(!df)
        !h = tau(1)
        !n = !df - 2
        !mi = mi_mi(!i)
        !aroot = !n*!mi*(!l + !h - 2*!mi) + !mi*(3*!l + 3*!h - 4*!mi) - 2*!l*!h
        !broot = (!n + 2)*!mi*(!l - !mi)*(!h - !mi)
        !c1root = !l^2*!mi^2*(!n+1)^2+!h^2*!mi^2*(!n+1)^2
        !c2root = 2*!l*!h*(2*!l*!h - 2*!l*!mi - 2*!h*!mi - 2*!n*!mi^2 - _
			         !n^2*!mi^2 + !mi^2)
        !omega = 0.25*((!aroot - @sqrt(!c1root + !c2root))/!broot)
    endif
    
    vector(!nLoop) mi_sad
    mi_sad(!i) = !omega

    '-------------------------------------------------------------------------------
    'Berechnung von r und u:
    '-------------------------------------------------------------------------------

    vector(!nLoop) mi_sad_r						'r-Parameter für Verteilung
    if !omega < 0 then
        vector(!df) r1
        for !j = 1 to !df
            r1(!j) = @log(1-2*!omega*taumi(!j))
        next
        mi_sad_r(!i) = -@sqrt(@sum(r1))
            delete r1       
    else
        vector(!df) r2
        for !j = 1 to !df
            r2(!j) = @log(1-2*!omega*taumi(!j))
        next
        mi_sad_r(!i) = @sqrt(@sum(r2))
        delete r2       
    endif
    
    vector(!nLoop) mi_sad_u						'u-Parameter für Verteilung

    vector(!df) uv
    for !j = 1 to !df
        uv(!j) = taumi(!j)^2/((1-2*!omega*taumi(!j))^2)
    next
    mi_sad_u(!i) = !omega*@sqrt(2*@sum(uv))
	 
    '-------------------------------------------------------------------------------
    'Wahrscheinlichkeit bei Sattelpunktapproximation:
    '-------------------------------------------------------------------------------


    vector(!nLoop) mi_prob_sad              'Wahrscheinlichkeitsverteilung
						   'bei Sattelpunktapprox.
    
    if mi_sad(!i) = 0 then                  'kein globaler Sattelpunkt gefunden
    else                                    'globaler Sattelpunkt gefunden
        mi_prob_sad(!i) = @cnorm(mi_sad_r(!i)-(1/mi_sad_r(!i)) _
									*log(mi_sad_r(!i)/mi_sad_u(!i)))
    endif
    
    if mi_sad(1) = 0 then					'Eintrag in Ergebnistabelle
        setcell(mi_info,30,3,"nein","l")
        if !q = !z then
            setcell(mi_info,30,4,"Iterationen (!z) erhöhen!","l")
        endif
    else
        setcell(mi_info,30,3,"ja","l")
    endif
    setcell(mi_info,31,2,"durchlaufene Iterationen:","l")
    setcell(mi_info,31,3,!q,"l")

    endif

next
endsub

'----------------------------------------------------------------------------------
'Subroutine aufrufen:
'----------------------------------------------------------------------------------

call calc_moran(mi_y, mi_x, mi_g, %code, %GlobMI, %Sad, mi_mi, mi_my, _
mi_moments_exp, mi_moments_var, mi_moments_skw, mi_moments_kur, mi_z_mi, _
mi_prob_nv, mi_sad, mi_sad_r, mi_sad_u,mi_prob_sad, mi_info)

!time = @toc                    'Beendet Timer

'----------------------------------------------------------------------------------
'Tabelle mit Ergebnissen:
'----------------------------------------------------------------------------------

setcolwidth(mi_info,1,5)
setcolwidth(mi_info,2,25)
setcolwidth(mi_info,3,25)
setcolwidth(mi_info,4,15)

setline(mi_info,2)
setcell(mi_info,3,2,"Das Programm moran.prg wurde ","l")
setcell(mi_info,3,4,@date,"r")
setcell(mi_info,4,2,"mit folgenden Argumenten aufgerufen: ","l")
setcell(mi_info,4,4,@time,"r")

setcell(mi_info,6,2,%0,"l")
setcell(mi_info,7,2,%1,"l")
setcell(mi_info,8,2,%2,"l")
setcell(mi_info,9,2,%3,"l")
setcell(mi_info,10,2,%4,"l")
setcell(mi_info,11,2,%5,"l")

setcell(mi_info,6,3,"Endogene Variable","l")
setcell(mi_info,7,3,"Exogene Variable","l")
setcell(mi_info,8,3,"Gewichtungsmatrix","l")
setcell(mi_info,9,3,"Kodierungsschema","l")
if %4 = "g" then setcell(mi_info,10,3,"globales MI","l") endif
if %4 = "gl" then setcell(mi_info,10,3,"globales und lokales MI","l") endif
if %5 = "j" then setcell(mi_info,11,3,"Sattelpunktapproximation durchführen","l") endif
if %5 = "n" then    
    setcell(mi_info,11,3,"Sattelpunktapproximation nicht durchführen","l")
    delete *sad*        'löscht nicht benötigte Variablen
endif
setline(mi_info,12)

setcell(mi_info,13,2,"RECHENZEIT","l")
setcell(mi_info,14,2,!time,"r")
setcell(mi_info,14,3," Sekunden","l")
setline(mi_info,15)

setcell(mi_info,16,2,"ERGEBNISSE NORMALVERTEILUNG","l")
setcell(mi_info,18,2,"Global","l")
setcell(mi_info,18,3,"Beschreibung","l")
setcell(mi_info,18,4,"Variable","l")
setcell(mi_info,19,2,mi_mi(1),"l")
setcell(mi_info,19,3,"Moran's I","l")
setcell(mi_info,19,4,"mi_mi","l")
setcell(mi_info,20,2,mi_moments_exp(1),"l")
setcell(mi_info,20,3,"Erwartungswert","l")
setcell(mi_info,20,4,"mi_moments_exp","l")
setcell(mi_info,21,2,mi_moments_var(1),"l")
setcell(mi_info,21,3,"Varianz","l")
setcell(mi_info,21,4,"mi_moments_var","l")
setcell(mi_info,22,2,mi_moments_skw(1),"l")
setcell(mi_info,22,3,"Schiefe","l")
setcell(mi_info,22,4,"mi_moments_skw","l")
setcell(mi_info,23,2,mi_moments_kur(1),"l")
setcell(mi_info,23,3,"Kurtosis","l")
setcell(mi_info,23,4,"mi_moments_kur","l")
setcell(mi_info,24,2,mi_z_mi(1),"l")
setcell(mi_info,24,3,"z-standardisiertes MI","l")
setcell(mi_info,24,4,"mi_z_mi","l")
setcell(mi_info,24,2,mi_prob_nv(1),"l")
setcell(mi_info,24,3,"Wahrscheinlichkeit","l")
setcell(mi_info,24,4,"mi_prob_nv","l")
setcell(mi_info,25,2,mi_z_mi(1),"l")
setcell(mi_info,25,3,"z-Wert","l")
setcell(mi_info,25,4,"mi_z_mi","l")
setline(mi_info,26)


if %5 = "j" then
    setcell(mi_info,27,2,"ERGEBNISSE SATTELPUNKTAPPROXIMATION","l")
    setcell(mi_info,29,2,"Globaler Sattelpunkt (Sekanten-Suchverfahren)","l")
    setcell(mi_info,30,2,"gefunden:","l")
    setcell(mi_info,33,2,"Global","l")
    setcell(mi_info,33,3,"Beschreibung","l")
    setcell(mi_info,33,4,"Variable","l")
    setcell(mi_info,34,2,mi_sad(1),"l")
    setcell(mi_info,34,3,"Sattelpunkt","l")
    setcell(mi_info,34,4,"mi_sad","l")
    setcell(mi_info,35,2,mi_sad_r(1),"l")
    setcell(mi_info,35,3,"r-Parameter für Verteilung","l")
    setcell(mi_info,35,4,"mi_sad_r","l")
    setcell(mi_info,36,2,mi_sad_u(1),"l")
    setcell(mi_info,36,3,"u-Parameter für Verteilung","l")
    setcell(mi_info,36,4,"mi_sad_u","l")
    setcell(mi_info,37,2,mi_prob_sad(1),"l")
    setcell(mi_info,37,3,"Wahrscheinlichkeit","l")
    setcell(mi_info,37,4,"mi_prob_sad","l")
    setline(mi_info,38)
endif

show mi_info            'Ergebnistabelle anzeigen

