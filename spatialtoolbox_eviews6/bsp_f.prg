'----------------------------------------------------------------------------------
' Anwendungsbeispiel F
'----------------------------------------------------------------------------------

'Demonstriert die Verwendung des Programms distance2weight.prg und getis.prg

'----------------------------------------------------------------------------------

'Datensatz öffnen:
open Daten\eu.wf1

'Subroutinen laden:
include distance2weight
include getis

'Subroutinen aufrufen:
call distance2weight(dst,138.1611)
call getis(w,ser01)

