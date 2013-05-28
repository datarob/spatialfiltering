'----------------------------------------------------------------------------------
'Subroutine zur Berechnung der Getis Statistik, räumlicher Filterung und 
'globalem Moran's I für die gefilterte Variable
'Eviews 5
'----------------------------------------------------------------------------------

'VERWENDUNG:

'Subroutine laden:							
'include getis
'Aufruf der Subroutine mit
'call getis([%0] [%1])
'%0... räumliche Gewichtungsmatrix
'%1... Variable als Reihe

'AUSGABE:

'	     gi_x = Vektor mit ursprünglicher Variable
'       gi_g = Vektor mit G_i Statistik
'  gi_star_g = Vektor mit G_i* Statistik
'      gi_eg = Vektor mit Erwartungswert der G_i Statistik
' gi_star_eg = Vektor mit Erwartungswert der G_i* Statistik
'      gi_sg = Vektor mit Std.abweichung der G_i Statistik
' gi_star_sg = Vektor mit Std.abweichung der G_i* Statistik
'     gi_z_gi= z-standardisiertes Gi
'gi_star_z_gi= z-standardisiertes Gi*

' 		gi_x_f = Vektor mit gefilterter Variable
'	  gi_z_mi =	Skalar mit z-standardisiertem globalen Moran's  für
'				   gefilterte Variable (gi_x_f)

'----------------------------------------------------------------------------------

mode quiet ver4	'stellt Kompatibilität zu EViews 5 sicher

subroutine getis(matrix w, series y)

!n = @rows(w)

stom(y, gi_x)
vector(!n) gi_g
vector(!n) gi_star_g
vector(!n) gi_eg
vector(!n) gi_star_eg
vector(!n) gi_sg
vector(!n) gi_star_sg
vector(!n) gi_z_gi
vector(!n) gi_star_z_gi
vector(!n) gi_x_f

for !i = 1 to !n

	gi_g(!i) = ((@rowextract(w,!i)*gi_x)(1))/(@sum(gi_x) - gi_x(!i))	'G_i Statistik
	gi_star_g(!i) = ((@rowextract(w,!i)*gi_x)(1))/(@sum(gi_x))	'G_i* Statistik

	!wi = @sum(@rowextract(w,!i))									
	gi_eg(!i) = !wi/(!n-1)											'Erwartungswert G_i
	gi_star_eg(!i) = !wi/(!n-1)									'Erwartungswert G_i*
	
	!yi1 = (@sum(gi_x) - gi_x(!i))/(!n-1)							'Std.abweichung G_i
	!yi2 =  ((@sumsq(gi_x) - gi_x(!i)^2)/(!n-1)) - !yi1^2
	gi_sg(!i) = @sqrt((!wi*(!n - 1 - !wi))/((!n - 1)^2 * (!n - 2))*(!yi2/!yi1^2))

	!yi1_star = (@sum(gi_x))/(!n-1)								'Std.abweichung G_i*
	!yi2_star =  ((@sumsq(gi_x))/(!n-1)) - !yi1_star^2
	gi_star_sg(!i) = @sqrt((!wi*(!n - 1 - !wi))/((!n - 1)^2 * (!n - 2)) _
							*(!yi2_star/!yi1_star^2))

	gi_z_gi(!i) = (gi_g(!i)-gi_eg(!i))	/gi_sg(!i)				'z-standardisiertes Gi
	gi_star_z_gi(!i) = (gi_star_g(!i)-gi_star_eg(!i))/ _
							gi_star_sg(!i)								 'z-standardisiertes Gi*

	gi_x_f(!i) = (gi_x(!i)*gi_eg(!i))/gi_g(!i)				'gefilterte Variable

next

!i = 1

	vector(!n) mi_x
	mi_x.fill(l) 1
	!df = !n - 1
	matrix mi_m = @identity(!n) - mi_x*@inverse((@transpose(mi_x)*mi_x)) _
					*@transpose(mi_x)					    'Residualmaker
	matrix mi_my = mi_m*gi_x_f                  'Residuen
	!denom = (@transpose(mi_my)*gi_x_f)(1,1)    'Nenner von Moran's I
		
	vector(!n) mi_mi
   mi_mi(!i) = (@transpose(mi_my)*w*mi_my)(1,1)/!denom

	matrix(!n,!n) mi_k = mi_m*0.5*(w+@transpose(w))*mi_m
   vector(!n) mi_moments_exp								'Erwartungswert
   mi_moments_exp(!i) = (@trace(mi_k))/!df			
    
   vector(!n) mi_moments_var								'Varianz
   mi_moments_var(!i) = (2*(!df*@trace(mi_k*mi_k)-@trace(mi_k)^2))/(!df^2*(!df+2))
	
	'z-standardisiertes Moran's I 
	scalar gi_z_mi = (mi_mi(!i) - mi_moments_exp(!i))/@sqrt(mi_moments_var(!i))

	delete mi_*

endsub

