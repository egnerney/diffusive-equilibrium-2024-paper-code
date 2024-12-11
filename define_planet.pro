function define_planet, planet, RP, Omega, GM

;Set planet-specific variabies
;Input
;Planet		[string]	Name of planet
;Output
;RP		[m]		Radius of planet
;Omega		[1/s]		Rotation rate
;GM		[m^3/s^2]	Gravitational parameter (mass)

	status = 0
	case planet of
        	'Jupiter':      begin
                        RP   = 7.1492d7
                        Omega= 2d*!dPI/3600.d/9.9259d
                        GM   = 1.2668689059d17;126686534d9
                        end
        	'Saturn':       begin
                        RP   = 6.033d7
                        Omega= 2d*!dPI/3600.d/10.7d
                        GM   = 37940585.2d9 ;Jacobson, 2006
                        end
        	else: status=0
endcase

return,status

end
