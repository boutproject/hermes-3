
RESTORE, "ue_bmk.idl"

nt = size(time_ue, /n_elements)

OPENW, f, "ue_bmk.csv", /GET_LUN, WIDTH=300

PRINTF, f, "Time, Ni[10], Ni[20], Ni[40], Te[10], Te[20], Te[40], Ti[10], Ti[20], Ti[40], Vi[10], Vi[20], Vi[40]"

FOR i=0, nt-1 DO PRINTF, f, time_ue[i], ", ", $
          ni_ue[i,10], ", ", ni_ue[i,20], ", ", ni_ue[i,40], ", ", $
          Te_ue[i,10], ", ", Te_ue[i,20], ", ", Te_ue[i,40], ", ", $
          Ti_ue[i,10], ", ", Ti_ue[i,20], ", ", Ti_ue[i,40], ", ", $
          Vi_ue[i,10], ", ", Vi_ue[i,20], ", ", Vi_ue[i,40]

CLOSE, f
