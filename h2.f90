program jastrow
use mt19937
use minimo

!gfortran -fbounds-check mt19937.f90 minimo.f90 h2.f90 -o h2.exe
! definimos las variables
double precision ::w,xtrial,p,q,el,ytrial,ztrial,x2trial,y2trial,z2trial,energia,d1,d2,d1_2
double precision :: d12,d22,inicio,final
integer :: n,i,k,m,seed,j,indice,f,l,b,caso,metodo
double precision :: res(2)
double precision :: jas_i,jas_f,jas_s,alfa_i,alfa_f,alfa_s,di_i,di_f,di_s
double precision ,allocatable :: r(:),en(:),a(:),fe(:)
double precision , allocatable :: x(:),y(:),z(:),x2(:),y2(:),z2(:)
call cpu_time(inicio)
! llamamos una vez a las funciones que nos gan a generar los numeros aleatorios

call system_clock(count=seed)
call sgrnd(seed)

!introduccion del metodo de calculo y de la funcion de onda a usar
read*, metodo
read*, caso
select case(caso)
	case(1)
		caso=1
		print*,'funcion 1'
	case(2)
		caso=2
		print*,'funcion 2'
	case(3)
		caso=3
		print*,'funcion 3'
	case(4)
		caso=4
		print*,'funcion 4'
	case(5)
		caso=5
		print*,'funcion 1 con jastrow'
	case(6)
		caso=6
		print*,'combinacion'

end select
if (metodo==1) then
	print*,'Usando diferencias centradas'
else
	print*,'Usando formula exacta'
end if


!numero de ejecuciones del montecarlo y entrada de los parametros variacionales
read(*,*) n
read(*,*) jas_i,jas_f,jas_s
read(*,*) alfa_i,alfa_f,alfa_s
read(*,*) di_i,di_f,di_s

!distancia array de distancia 
m=(di_f-di_i)/(di_s)+1
! distancia array de parametro alfa
f=(alfa_f-alfa_i)/alfa_s+1
!distancia array de parametro de jastrow
b=(jas_f-jas_i)/jas_s+1


print *,'MC run',n
print *, 'distancia,pasos',di_i,di_f,di_s,m
print *, 'alfa,pasos',alfa_i,alfa_f,alfa_s,f
print *, 'jastrow,pasos',jas_i,jas_f,jas_s,b

!condiciones que comprueban si los parametros finales e iniciales son 0
if (di_i.eq.di_f) then
	m=1
end if

if (alfa_i.eq.alfa_f) then
	f=1
end if

if (jas_i.eq.jas_f) then
	b=1
end if
	

allocate (x(n),y(n),z(n),x2(n),y2(n),z2(n))
allocate (r(m),en(m),a(f),fe(b))

!llenamos los bucles de los parametros
fe(1)=jas_i
do l=1,b-1
fe(l+1)=jas_i+l*jas_s
end do

a(1)=alfa_i
do l=1,f-1
a(l+1)=alfa_i+l*alfa_s
end do



! El fichero 1, guarda E y r
! El fichero 2, guarda las posiciones x,y,z de cada electron, 
! creando un fichero cuyo tamaño es el de las ejecuciones de MC
! El fichero 3, muestra los resultados de todos los parametros variacionales


!bucles que recorren cada combinación de parametros variacionales

do l=1,b

! open(3,file='final.txt',status='unknown')


do j=1,f
	r(1)=di_i*1.889
	do i=1,m-1		
		r(i+1)=i*di_s*1.889+di_i*1.889
	end do
		
	    
	! open(1,file="datos.txt",status='unknown')	
	
do k=1,m
   x(1)=0
   y(1)=-0.1
   z(1)=0
   x2(1)=0.2+r(k)
   y2(1)=0+r(k)
   z2(1)=0.1+r(k)
   el=0
   energia=0
      ! open(2,file="EvsR.txt",status='unknown')
!ejecucion del Monte Carlo
   do i=1,n-1

      q=grnd()
      q=-1+2*q
      xtrial=x(i)+q

      q=grnd()
      q=-1+2*q
      ytrial=y(i)+q

      q=grnd()
      q=-1+2*q
      ztrial=z(i)+q

      q=grnd()
	  q=-1+2*q     
      x2trial=x2(i)+q

      q=grnd()
	  q=-1+2*q     
      y2trial=y2(i)+q

      q=grnd()
	  q=-1+2*q
      z2trial=z2(i)+q

      w=(funcion(a(j),r(k),xtrial,ytrial,ztrial,x2trial,y2trial,z2trial,fe(l),caso)**2)/&
	  (funcion(a(j),r(k),x(i),y(i),z(i),x2(i),y2(i),z2(i),fe(l),caso)**2)
      p=grnd()
      if (w.GE.1 .or. w.GE.p)  then
         x(i+1)=xtrial
         y(i+1)=ytrial
         z(i+1)=ztrial
         x2(i+1)=x2trial
         y2(i+1)=y2trial
         z2(i+1)=z2trial
		 
      else
          x(i+1)=x(i)
          y(i+1)=y(i)
          z(i+1)=z(i)
          x2(i+1)=x2(i)
          y2(i+1)=y2(i)
          z2(i+1)=z2(i)
		  
      end if
	  ! write(2,*) x(i)/1.889,y(i)/1.889,z(i)/1.889,x2(i)/1.889,y2(i)/1.889,z2(i)/1.889
   d1=sqrt(x(i)**2+y(i)**2+z(i)**2)
   d2=sqrt(x2(i)**2+y2(i)**2+z2(i)**2)
   d12=sqrt((x(i)-r(k))**2+y(i)**2+z(i)**2)
   d22=sqrt((x2(i)-r(k))**2+y2(i)**2+z2(i)**2)
   d1_2=sqrt((x(i)-x2(i))**2+(y(i)-y2(i))**2+(z(i)-z2(i))**2)
   
    
   
   
   
select case(metodo)
		case(1)

			energia=(-0.5*centrada(a(j),r(k),x(i),y(i),z(i),x2(i),y2(i),z2(i),fe(l),caso)/&
			funcion(a(j),r(k),x(i),y(i),z(i),x2(i),y2(i),z2(i),fe(l),caso)-&
			(1/d1)-(1/d12)-(1/d2)-(1/d22)+(1/d1_2)+(1/r(k)))
		case(2)
		
			energia=(-0.5*funciones(a(j),r(k),x(i),y(i),z(i),x2(i),y2(i),z2(i),fe(l),caso)/&
			funcion(a(j),r(k),x(i),y(i),z(i),x2(i),y2(i),z2(i),fe(l),caso)-&
			(1/d1)-(1/d12)-(1/d2)-(1/d22)+(1/d1_2)+(1/r(k)))
	end select


   el=el+energia
   
   end do

   ! close(2)
   ! print*,a(j),r(k)/1.889,((el/n)+1)*27.211
   ! write(1,*) r(k)/1.889,((el/n)+1)*27.211
   r(k)=r(k)/1.889
   en(k)=((el/n)+1)*27.211

end do
close(1)

res=calc_minimo(en,m)
indice=res(2)

print *,fe(l),a(j),r(indice),res(1)
! write(3,*) fe(l),a(j),r(indice),res(1)


end do


end do
! close(3)
call cpu_time(final)

print*, final-inicial
contains



!contiene la funciones de onda
double precision function funcion(a,r,x,y,z,x2,y2,z2,fe,caso)

       double precision :: d1,d2,d12,d22,d1_2,r,fe
       double precision :: x,y,z,x2,y2,z2,a
	   integer :: caso
       
       d1=sqrt(x**2+y**2+z**2)
       d2=sqrt(x2**2+y2**2+z2**2)
       d12=sqrt((x-r)**2+y**2+z**2)
       d22=sqrt((x2-r)**2+y2**2+z2**2)
       d1_2=sqrt((x-x2)**2+(y-y2)**2+(z-z2)**2)
       
	   
	   select case(caso)
			case(1)
			   !FUNCION 1
			   funcion=exp(-a*d1)*exp(-a*d2)+&
			   exp(-a*d1)*exp(-a*d22)+&
			   exp(-a*d12)*exp(-a*d2)+&
			   exp(-a*d12)*exp(-a*d22)

			case(2)			
			   !2
			   
			   funcion=(exp(-a*d1)-exp(-a*d12))*(exp(-a*d2)-exp(-a*d22))
			case(3)
			   !3
			   
			   funcion=2*exp(-a*d1)*exp(-a*d2)-2*exp(-a*d12)*exp(-a*d22)
			case(4)
			   !4
			   
			   funcion=-2*exp(-a*d1)*exp(-a*d22)+2*exp(-a*d12)*exp(-a*d2)
	   
			case(5)
				!1 con jastrow
				
			   funcion=exp(-a*d1)*exp(-a*d2)+&
			   exp(-a*d1)*exp(-a*d22)+&
			   exp(-a*d12)*exp(-a*d2)+&
			   exp(-a*d12)*exp(-a*d22)
			   funcion=funcion*exp(-fe/(2*(1+(d1_2/fe))))
			case(6)
			   ! combinacion
			   
			   funcion=(exp(-a*d1)+exp(-a*d12))*(exp(-a*d2)+exp(-a*d22))-fe*(exp(-a*d1)-exp(-a*d12))*(exp(-a*d2)-exp(-a*d22))
       
		end select
	   
       return
       
end function
!calcula los laplacianos usando derivada numerica
double precision function centrada(a,r,x,y,z,x2,y2,z2,fe,caso)
       double precision :: h=0.001
       double precision ::x,y,z,x2,y2,z2,r,a,fe
	   integer:: caso
       centrada=&
       funcion(a,r,x+h,y,z,x2,y2,z2,fe,caso)-2*funcion(a,r,x,y,z,x2,y2,z2,fe,caso)+funcion(a,r,x-h,y,z,x2,y2,z2,fe,caso)+&
       funcion(a,r,x,y+h,z,x2,y2,z2,fe,caso)-2*funcion(a,r,x,y,z,x2,y2,z2,fe,caso)+funcion(a,r,x,y-h,z,x2,y2,z2,fe,caso)+&
       funcion(a,r,x,y,z+h,x2,y2,z2,fe,caso)-2*funcion(a,r,x,y,z,x2,y2,z2,fe,caso)+funcion(a,r,x,y,z-h,x2,y2,z2,fe,caso)+&
       funcion(a,r,x,y,z,x2+h,y2,z2,fe,caso)-2*funcion(a,r,x,y,z,x2,y2,z2,fe,caso)+funcion(a,r,x,y,z,x2-h,y2,z2,fe,caso)+&
       funcion(a,r,x,y,z,x2,y2+h,z2,fe,caso)-2*funcion(a,r,x,y,z,x2,y2,z2,fe,caso)+funcion(a,r,x,y,z,x2,y2-h,z2,fe,caso)+&
       funcion(a,r,x,y,z,x2,y2,z2+h,fe,caso)-2*funcion(a,r,x,y,z,x2,y2,z2,fe,caso)+funcion(a,r,x,y,z,x2,y2,z2-h,fe,caso)
       centrada=centrada/(h**2)

       return

end function

!calcula los laplacianos usando las formulas exactas
double precision function funciones(a,r,x,y,z,x2,y2,z2,fe,caso)
       
       double precision :: d1,d2,d12,d22,primero,segundo,tercero,cuarto,quinto
       double precision :: a,x,y,z,x2,y2,z2,r,fe
	   integer :: caso
	  
       d1=sqrt(x**2+y**2+z**2)
       d2=sqrt(x2**2+y2**2+z2**2)
       d12=sqrt((x-r)**2+y**2+z**2)
       d22=sqrt((x2-r)**2+y2**2+z2**2)
	   d1_2=sqrt((x-x2)**2+(y-y2)**2+(z-z2)**2)
	   
	   select case(caso)
			case(1)
			   !funcion 1
			   primero=exp(-a*d1)*((a**2)-(2*a/(d1)))+((a**2)-(2*a/(d12)))*exp(-a*d12)
			   primero=primero*(exp(-a*d2)+exp(-a*d22))
			   segundo=exp(-a*d2)*((a**2)-(2*a/(d2)))+((a**2)-(2*a/(d22)))*exp(-a*d22)
			   segundo=segundo*(exp(-a*d1)+exp(-a*d12))
			   funciones=primero+segundo
			case(2)
				!funcion 2
			   primero=exp(-a*d1)*((a**2)-(2*a/(d1)))-((a**2)-(2*a/(d12)))*exp(-a*d12)
			   primero=primero*(exp(-a*d2)-exp(-a*d22))
			   segundo=exp(-a*d2)*((a**2)-(2*a/(d2)))-((a**2)-(2*a/(d22)))*exp(-a*d22)
			   segundo=segundo*(exp(-a*d1)-exp(-a*d12))
			   funciones=primero+segundo
			case(3)
			   !funcion 3
			   primero=exp(-a*d1)*exp(-a*d2)*((a**2)-(2*a/(d1)))-((a**2)-(2*a/(d12)))*exp(-a*d12)*exp(-a*d22)
			   segundo=exp(-a*d2)*exp(-a*d1)*((a**2)-(2*a/(d2)))-((a**2)-(2*a/(d22)))*exp(-a*d22)*exp(-a*d12)
			   funciones=2*primero+2*segundo
			case(4)
			   !funcion 4
			   primero=-exp(-a*d1)*exp(-a*d22)*((a**2)-(2*a/(d1)))+((a**2)-(2*a/(d12)))*exp(-a*d12)*exp(-a*d2)
			   segundo=-exp(-a*d22)*exp(-a*d1)*((a**2)-(2*a/(d22)))+((a**2)-(2*a/(d2)))*exp(-a*d2)*exp(-a*d12)
			   funciones=2*primero+2*segundo
			case(5)
				!funcion 1 con jastrow
		   
				primero=((a**2-2*a/d1)*exp(-a*d1)+(a**2-2*a/d12)*exp(-a*d12))/(exp(-a*d1)+exp(-a*d12))
				segundo=(1/d1_2)*(1/(1+d1_2/fe)**3)+&
				(1/d1_2**2)*((1/(2*(1+d1_2/fe)**2))**2)*((x-x2)**2+(y-y2)**2+(z-z2)**2)
				tercero=-2*((-a*x*exp(-a*d1)/(d1)-a*(x-r)*exp(-a*d12)/d12)/(exp(-a*d1)+exp(-a*d12))*(-1/(2*d1_2*(1+d1_2/fe)**2))*(x-x2)+&
				(-a*y*exp(-a*d1)/(d1)-a*(y)*exp(-a*d12)/d12)/(exp(-a*d1)+exp(-a*d12))*(-1/(2*d1_2*(1+d1_2/fe)**2))*(y-y2)+&
				(-a*z*exp(-a*d1)/(d1)-a*(z)*exp(-a*d12)/d12)/(exp(-a*d1)+exp(-a*d12))*(-1/(2*d1_2*(1+d1_2/fe)**2))*(z-z2))				
				cuarto=((a**2-2*a/d2)*exp(-a*d2)+(a**2-2*a/d22)*exp(-a*d22))/(exp(-a*d2)+exp(-a*d22))
				quinto=-2*((-a*x2*exp(-a*d2)/(d2)-a*(x2-r)*exp(-a*d22)/d22)/(exp(-a*d2)+exp(-a*d22))*(-1/(2*d1_2*(1+d1_2/fe)**2))*(x2-x)+&
				(-a*y2*exp(-a*d2)/(d2)-a*(y2)*exp(-a*d22)/d22)/(exp(-a*d2)+exp(-a*d22))*(-1/(2*d1_2*(1+d1_2/fe)**2))*(y2-y)+&
				(-a*z2*exp(-a*d2)/(d2)-a*(z2)*exp(-a*d22)/d22)/(exp(-a*d2)+exp(-a*d22))*(-1/(2*d1_2*(1+d1_2/fe)**2))*(z2-z))				
				funciones=(primero+2*segundo+tercero+cuarto+quinto)*(exp(-a*d1)*exp(-a*d2)+&
				exp(-a*d1)*exp(-a*d22)+&
				exp(-a*d12)*exp(-a*d2)+&
				exp(-a*d12)*exp(-a*d22))*exp(-fe/(2*(1+(d1_2/fe))))
	    
			case(6)
				! combinacion
				
			   primero=exp(-a*d1)*((a**2)-(2*a/(d1)))+((a**2)-(2*a/(d12)))*exp(-a*d12)
			   primero=primero*(exp(-a*d2)+exp(-a*d22))
			   segundo=exp(-a*d2)*((a**2)-(2*a/(d2)))+((a**2)-(2*a/(d22)))*exp(-a*d22)
			   segundo=segundo*(exp(-a*d1)+exp(-a*d12))
			   
			   
			   tercero=exp(-a*d1)*((a**2)-(2*a/(d1)))-((a**2)-(2*a/(d12)))*exp(-a*d12)
			   tercero=tercero*(exp(-a*d2)-exp(-a*d22))
			   cuarto=exp(-a*d2)*((a**2)-(2*a/(d2)))-((a**2)-(2*a/(d22)))*exp(-a*d22)
			   cuarto=cuarto*(exp(-a*d1)-exp(-a*d12))
			   
			   funciones=(primero+segundo)-fe*(tercero+cuarto)
			   
		end select
	   
       return
       
end function
		
			
end program jastrow









       
       
