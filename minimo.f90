module minimo


contains
function calc_minimo(a,m) result(resultado)
	
	double precision::mini,ind
	integer ::i,m
	double precision :: resultado(2)
	double precision :: a(m)
	mini=a(1)
	ind=1
	do i=1,m
		if (mini.le.a(i)) then
			mini=mini
			
			else
				mini=a(i)
				ind=i
		end if
		
	
	end do
	resultado(1)=mini
	resultado(2)=ind	
	return

end function




end module minimo