'----------------------------------------------------------------------------------
'Optimales Delta für Gewichtungsmatrix bestimmen
'EViews 5
'----------------------------------------------------------------------------------

'VERWENDUNG:
'run opt_delta [%0] [%1] [%2] [%3] [%4] [%5] [%6]

'Eingabeparameter:
'%0 = Variable als Reihe
'%1 = Räumliche Distanzmatrix
'%2 = Startwert pa
'%3 = Startwert pb
'%4 =	Startwert pc
'%5 = Genauigkeit
'%6 = maximale Iterationen

'AUSGABE:

'	        opt_ergebnis = Skalar mit optimalem Delta
'	      opt_startwerte = Tabelle mit Ergebnissen der erstmaligen Schachtelung
'	opt_goldener_schnitt = Tabelle mit Ergebnissen von Methode des Goldenen Schnitts
'	     				  gi_x = Vektor mit ursprünglicher Variable
'                  gi_g = Vektor mit G_i Statistik
'  			   gi_star_g = Vektor mit G_i* Statistik
'      				 gi_eg = Vektor mit Erwartungswert der G_i Statistik
'            gi_star_eg = Vektor mit Erwartungswert der G_i* Statistik
'                 gi_sg = Vektor mit Std.abweichung der G_i Statistik
'            gi_star_sg = Vektor mit Std.abweichung der G_i* Statistik
'               gi_z_gi = z-standardisiertes Gi
'          gi_star_z_gi = z-standardisiertes Gi*
' 					   gi_x_f = Vektor mit gefilterter Variable
'               gi_z_mi =	Skalar mit z-standardisiertem globalen Moran's  für
'									gefilterte Variable (gi_x_f)
' 							  u = unstandardisierte Gewichtungsmatrix
'							  w = Zeilensummen standardisierte Gewichtungsmatrix

mode quiet ver4	'stellt Kompatibilität zu EViews 5 sicher

'----------------------------------------------------------------------------------
'Fehlerüberprüfung:
'----------------------------------------------------------------------------------

!rG = @rows({%1})
!cG = @columns({%1})
if !rG<>!cG then
	statusline Programmabbruch: Distanzmatrix nicht quadratisch
	stop
endif

'----------------------------------------------------------------------------------
'Notwendige Subroutinen laden:
'----------------------------------------------------------------------------------

include optimierung					'Subroutinen laden
include distance2weight
include getis

'----------------------------------------------------------------------------------
'Zielfunktion:
'----------------------------------------------------------------------------------

subroutine zf(scalar delta)
	call distance2weight({%1},delta)
	call getis(w,{%0})
	scalar zf = @abs(gi_z_mi)
endsub

'----------------------------------------------------------------------------------
'Anwendung:
'----------------------------------------------------------------------------------

scalar pa = {%2}					'Startwerte
scalar pb = {%3}
scalar pc = {%4}

call einschachteln(pa,pb,pc)		   'Erstmalige Schachtelung des Minimums
	
call goldener_schnitt(pa,pb,pc,{%5},{%6})	   'Methode des Goldenen Schnitts

'Hilfsvariablen löschen:
delete xa xb xc xd fa fb fc fd zf pa pb pc differenz fast_null letzte_fb _
	    praezision

