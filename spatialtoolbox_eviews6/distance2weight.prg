'----------------------------------------------------------------------------------
'Konvertiert Distanz- zu Gewichtungsmatrix:
'Eviews 5
'----------------------------------------------------------------------------------

'VERWENDUNG:

'Subroutine laden:							
'include distance2weight
'Aufruf der Subroutine mit
'call distance2weight([%0] [%1])
'%0... räumliche Distanzmatrix
'%1... Skalar mit Delta für Distancedecay Funktion

'AUSGABE:

' u = unstandardisierte Gewichtungsmatrix
' w = Zeilensummen standardisierte Gewichtungsmatrix

'----------------------------------------------------------------------------------

mode quiet ver4	'stellt Kompatibilität zu EViews 5 sicher

subroutine distance2weight(matrix di, scalar delta)

!beta = 1/delta/1000       			'Umwandlung in Kilometer
!n = @rows(di)             			'Anzahl der Beobachtungen

matrix(!n,!n) w				    
w= exp(!beta*-di)							'Distancedecay Fuktion
w = w - @identity(!n)					'Nullen in der Hauptdiagonale

matrix u = w

vector(!n) rowsums						'Zeilensummenstandardisierung
for !i = 1 to !n
	rowsums(!i) = @sum(@rowextract(w,!i))
next

w = @transpose(w*(@inverse(@makediagonal(rowsums))))

delete rowsums

endsub

