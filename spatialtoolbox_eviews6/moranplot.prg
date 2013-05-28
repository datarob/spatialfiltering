'----------------------------------------------------------------------------------
'Moran Scatterplot (mit Ausreißerdiagnose) / Moran Scatterplotmatrix
'Eviews 6
'----------------------------------------------------------------------------------

'VERWENDUNG:
'run moranplot [%0] [%1] [%2]

'Eingabeparameter:
'%0 = Variable (oder Gruppe von Variablen)
'%1 = räumliche Gewichtungsmatrix
'%2 = m (für Moran Scatterplotmatrix)

'mode quiet ver5	'stellt Kompatibilität zu EViews 5 sicher

'----------------------------------------------------------------------------------
'Fehlerüberprüfung:
'----------------------------------------------------------------------------------

!rG = @rows({%1})
!cG = @columns({%1})
if !rG<>!cG then
	statusline Programmabbruch: Gewichtungsmatrix nicht quadratisch
	stop
endif

'----------------------------------------------------------------------------------
'Subroutine: erzeugt die notwendigen Reihen
'----------------------------------------------------------------------------------

subroutine mp_calc

series mp_{%series}_y = {%series} - @mean({%series})			'Abweichungen vom Mittelwert
stom(mp_{%series}_y,mp_{%series}_y_m)					      'umwandeln in Vektor
vector mp_{%series}_wy_m = {%1}*mp_{%series}_y_m      'räumlich gelagte Variable
mtos(mp_{%series}_wy_m,mp_{%series}_wy)					   'umwandeln in Reihe	
delete mp_{%series}_y_m mp_{%series}_wy_m	
				
endsub

'----------------------------------------------------------------------------------
'Subroutine: erzeugt die Grafiken für Moran Scatterplotmatrix
'----------------------------------------------------------------------------------

subroutine mp_plot

group mp_{%0}_gr mp_{%series1}_y mp_{%series2}_wy
freeze(mp_{%0}_!i_!j)  mp_{%0}_gr.scat linefit
'graph mp_{%0}_!i_!j.scat mp_{%series1}_y mp_{%series2}_wy
mp_{%0}_!i_!j.draw(line,l) 0.0					'vertikale  0-Linie
mp_{%0}_!i_!j.draw(line,b) 0.0					'horizontale  0-Linie
%mp_string_n = "mp_" + %0 + "_" +  @str(!i) +  "_" + @str(!j)

endsub

'----------------------------------------------------------------------------------
'Subroutine: Statistiken für Ausreißerdiagnose
'----------------------------------------------------------------------------------

subroutine resdiag

'Matrix der exogenen Variablen
series const = 1
group exo const {%series}
stom(exo,x)	
			
'Hauptdiagonale der Hutmatrix
vector htt=@getmaindiagonal(x*@inverse(@inner(x))*@transpose(x))
mtos(htt, hhtt)

'Regression durchführen
equation mp_eq.ls mp_{%series}_wy mp_{%series}_y c

'intern studentisierte Resiuden
series mp_rresid=resid/(mp_eq.@se*@sqrt(1-hhtt))

'studentisierte Residuen
series mp_sresid = mp_rresid*@sqrt((@obs(mp_{%series}_wy)-2) _
                   /(@obs(mp_{%series}_y)-1-mp_rresid^2))
graph mp_sresids.line mp_sresid
mp_sresids.draw(line,b) +3.5
mp_sresids.draw(line,b) -3.5

'Cook's Distanzen
series mp_cook_d = mp_rresid^2*hhtt/(1-hhtt)/2
graph mp_{%series}_cook.line mp_cook_d

'Moran Scatterplot mit Regressionsgerade
group mp_{%series}_gr mp_{%series}_y mp_{%series}_wy

'Anstieg der Regressionsgerade
!mi = c(1)
mp_{%series}.addtext(x) MI = !mi

'Moran Scatterplot mit Kernelfit
freeze(mp_{%series}_kerfit)  mp_{%series}_gr.scat kernfit
mp_{%series}_kerfit.draw(line,l) 0.0			'vertikale 0-Linie
mp_{%series}_kerfit.draw(line,b) 0.0			'horizontale  0-Linie

'Ergebnisgrafik
graph mp_ergebnis.merge mp_{%series} mp_{%series}_kerfit mp_sresids _
	   mp_{%series}_cook

delete const exo x hhtt htt mp_eq mp_{%series}_gr mp_{%series} _
       mp_{%series}_kerfit mp_sresids mp_{%series}_cook

endsub

'----------------------------------------------------------------------------------
'Moran Scatterplotmatrix
'----------------------------------------------------------------------------------

if %2 = "m" then
	!k = {%0}.@count
	for !i = 1 to !k
		%series =  {%0}.@seriesname(!i)
		call mp_calc
	next	
	for !i = 1 to !k
		%series2 = {%0}.@seriesname(!i)
			for !j = 1 to !k	
				%series1 = {%0}.@seriesname(!j)
				call mp_plot
				%mp_string = %mp_string + " " + %mp_string_n
		next
	next
	graph mp_{%0}.merge {%mp_string}			'erzeugt Matrix aus einzelnen Grafiken
	delete %mp_string
	show mp_{%0}									'Ergebnisgrafik anzeigen

'----------------------------------------------------------------------------------
'einzelner Moran Scatterplot
'----------------------------------------------------------------------------------

else
	%series = %0
	call mp_calc								'notwendige Reihen erzeugen
	'Moran Scatterplot
                group mp_{%series}_gr mp_{%series}_y mp_{%series}_wy
	freeze(mp_{%series})  mp_{%series}_gr.scat linefit
	'graph mp_{%series}.scat mp_{%series}_y mp_{%series}_wy
	mp_{%series}.draw(line,l) 0.0		'vertikale 0-Linie
	mp_{%series}.draw(line,b) 0.0		'horizontale 0-Linie	
	call resdiag								'Ausreißerdiagnose
	show mp_ergebnis							'Ergebnisgrafik anzeigen

endif

