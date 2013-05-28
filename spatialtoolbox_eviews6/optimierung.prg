'----------------------------------------------------------------------------------
'Subroutinen zur Ermittlung eines Minimums bei einer Zielfunktion in 
'einer Variablen mit der "Methode des Goldenen Schnitts"
'EViews 5
'----------------------------------------------------------------------------------

'VERWENDUNG:

'Subroutinen laden:							
'include optimierung

'Zielfunktion definieren:		Input: ein beliebiger Skalar, 
'				        				Output: ein Skalar zf
'z.B.
'subroutine zf(scalar x)
'	scalar zf = x^2
'endsub		

'----------------------------------------------------------------------------------
'Erstmalige Schachtelung des Minimums:
'----------------------------------------------------------------------------------

'VERWENDUNG:

'call einschachteln(pa,pb,pc)
'			pa, pb, pc sind Skalare mit Startwerten

'AUSGABE:

'optimale Startwerte werden in den Skalaren pa, pb, pc gespeichert
'opt_startwerte... Tabelle mit allen durchgeführten Iterationen

'----------------------------------------------------------------------------------

mode quiet ver4	'stellt Kompatibilität zu EViews 5 sicher

subroutine einschachteln(scalar pa, scalar pb, scalar pc)
	
	scalar fa
	scalar fb
	scalar fc
	scalar xa
	scalar xb
	scalar xc

	call zf(pa)
	fa = zf	
	call zf(pb)
	fb = zf
	if fa > fb then
		xa = pa
		xb = pb
	else
		fc = fa
		fa = fb
		fb = fc
		xa = pb
		xb = pa		
	endif

	!schritt = 0
	scalar letzte_fb = fb
 	xc = xb + 1.618034*(xb - xa)
	call zf(xc)
	fc = zf
	
	table opt_startwerte
	setcell(opt_startwerte,1,1,"Schritt")
	setcell(opt_startwerte,1,2,"xa")
	setcell(opt_startwerte,1,3,"fa")
	setcell(opt_startwerte,1,4,"xb")
	setcell(opt_startwerte,1,5,"fb")
	setcell(opt_startwerte,1,6,"xc")
	setcell(opt_startwerte,1,7,"fc")

	while fc < letzte_fb
		!schritt = !schritt + 1
		statusline Erstmalige Schachtelung des Minimums: Iteration !schritt
		xc = xb + 1.618034*(xb - xa)
		call zf(xc)
		fc = zf
		letzte_fb = fb
		
		setcell(opt_startwerte,!schritt+1,1,!schritt)
		setcell(opt_startwerte,!schritt+1,2,xa)
		setcell(opt_startwerte,!schritt+1,3,fa)
		setcell(opt_startwerte,!schritt+1,4,xb)
		setcell(opt_startwerte,!schritt+1,5,fb)
		setcell(opt_startwerte,!schritt+1,6,xc)
		setcell(opt_startwerte,!schritt+1,7,fc)
		
		if fc <= fb then
			xa = xb
			fa = fb
			xb = xc
			fb = fc
		endif
	wend
	
	if xa > xc then
		pa = xc
		pb = xb
		pc = xa
	else
		pa = xa
		pb = xb
		pc = xc
	endif
	
endsub

'----------------------------------------------------------------------------------
'Methode des Goldenen Schnitts:
'----------------------------------------------------------------------------------

'VERWENDUNG:

'call goldener_schnitt(pa,pb,pc,{%5},{%6})
'		pa, pb, pc wurden vorher mit Subroutine "einschachteln" bestimmt
'		%5 = Präzision (z.B. 0.0001)
'		%6 = Maximale Iterationen (z.B. 100)	

'AUSGABE:

'opt_goldener_schnitt... Tabelle mit allen durchgeführten Iterationen
'opt_ergebnis... Skalar mit optimalem Wert, bei dem die Zielfunktion ein Minimum
'					   erreicht

'----------------------------------------------------------------------------------
subroutine goldener_schnitt(scalar a, scalar b, scalar c, scalar _
relative_praezision, scalar maximale_iterationen)
	
	scalar a
	scalar b
	scalar c
	scalar relative_praezision
	scalar maximale_iterationen
	scalar differenz
	scalar praezision
	scalar fast_null = 1E-20
	scalar xd
	scalar fd

	
	!schritt = 0
	xa = a
	xb = b
	xc = c
	call zf(xa)
	fa = zf
	call zf(xb)
	fb = zf
	call zf(xc)
	fc = zf
	differenz = xc - xa
	praezision = relative_praezision*@abs(xb)

	table opt_goldener_schnitt
	setcell(opt_goldener_schnitt,1,1,"Schritt")
	setcell(opt_goldener_schnitt,1,2,"xa")
	setcell(opt_goldener_schnitt,1,3,"fa")
	setcell(opt_goldener_schnitt,1,4,"xb")
	setcell(opt_goldener_schnitt,1,5,"fb")
	setcell(opt_goldener_schnitt,1,6,"xc")
	setcell(opt_goldener_schnitt,1,7,"fc")
	setcell(opt_goldener_schnitt,1,8,"xd")
	setcell(opt_goldener_schnitt,1,9,"fd")


	while (!schritt < maximale_iterationen) and _
			(praezision < differenz) and _
			(differenz > fast_null)
				statusline Methode des Goldenen Schnitts: Iteration !schritt
				if (xb - xa) > (xc - xb) then
				xd = xb - 0.38197*(xb - xa)
				call zf(xd)
				fd = zf
				
				setcell(opt_goldener_schnitt,!schritt+2,1,!schritt)
				setcell(opt_goldener_schnitt,!schritt+2,2,xa)
				setcell(opt_goldener_schnitt,!schritt+2,3,fa)
				setcell(opt_goldener_schnitt,!schritt+2,4,xb)
				setcell(opt_goldener_schnitt,!schritt+2,5,fb)
				setcell(opt_goldener_schnitt,!schritt+2,6,xc)
				setcell(opt_goldener_schnitt,!schritt+2,7,fc)
				setcell(opt_goldener_schnitt,!schritt+2,8,xd)
				setcell(opt_goldener_schnitt,!schritt+2,9,fd)
				
				if fd < fb then
					xc = xb
					fc = fb
					xb = xd
					fb = fd
				else
					xa = xd
					fa = fd
				endif
			else
				xd = xb + 0.38197*(xc - xb)
				call zf(xd)
				fd = zf

				setcell(opt_goldener_schnitt,!schritt+2,1,!schritt)
				setcell(opt_goldener_schnitt,!schritt+2,2,xa)
				setcell(opt_goldener_schnitt,!schritt+2,3,fa)
				setcell(opt_goldener_schnitt,!schritt+2,4,xb)
				setcell(opt_goldener_schnitt,!schritt+2,5,fb)
				setcell(opt_goldener_schnitt,!schritt+2,6,xc)
				setcell(opt_goldener_schnitt,!schritt+2,7,fc)
				setcell(opt_goldener_schnitt,!schritt+2,8,xd)
				setcell(opt_goldener_schnitt,!schritt+2,9,fd)

				if fd < fb then
					xa = xb
					fa = fb
					xb = xd
					fb = fd
				else
					xc = xd
					fc = fd
				endif
			endif
			!schritt = !schritt + 1
			differenz = xc - xa
			praezision = relative_praezision*@abs(xb)
	wend
	scalar opt_ergebnis = xb
	!result = opt_ergebnis
	statusline Goldener Schnitt liegt bei !result
endsub

