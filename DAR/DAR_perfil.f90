module globalDAR

	!este programa foi feito por base do modelo de exemplo do cap. 4 do livro barabasi
	!========================================================================================================
	integer(4)::vez		!vezes que o programa vai rodar. e o número de amos	
	integer(4)::L
	integer(4)::mol	!O número de moléculas que serão depositadas, as camadas serão mol/L
	integer(4),allocatable::S(:),VE(:),VD(:)
	real(8),allocatable::WA(:),W(:)
	integer(4)::J(1:3)
	integer(4)::a,n,k,i,am,ca,p,erro,i2
	real(8)::b,h,z
	real(8),parameter::m=1.3
	!---------------------------- número aleatório -------------------------------------------------------------! 
	integer, parameter ::AMAG=843314861,BMAG=453816693,im=1073741824,a1=65539
	real*8, parameter::r_1=1.0/(4.0*im)
	integer ::iseed,jseed,rcont
	!-----------------------------------------------------------------------------------------------------------!
	!========================================================================================================

	contains
	subroutine arquivos !!arquivos que aparecerão

		iseed=88925613; jseed=83991
		
		if(L==128) then
			open(1,file='perfilDARL=128mol=100000Lamo=1.dat')
			
		end if

	end subroutine
!================================================================================================!
	subroutine reg_var !!deixa o substrato circular

		do i=2, L
			VE(i)=i-1; VD(i)=i+1
		end do
		VD(L)=1; VE(1)=L; VD(1)=2

	end subroutine
!================================================================================================!	
	subroutine dados !!armazenamento dos dados
				
		do n=1, L
			write(1,*) n,S(n)
		end do		
			
	end subroutine
!================================================================================================!
	subroutine depo		!!aqui ocorre as deposições
		!-------------------------------------------
		J=0; J(1)=S(VE(a)); J(2)=S(a); J(3)=S(VD(a))
		!-------------------------------------------			
		
		if(minval(J)==J(2)) then 


			S(a)=S(a)+1


		else if(J(1)==J(3)) then
					
			iseed=(AMAG*iseed)+BMAG; rcont=rcont+1; z=r_1*iseed+0.5d0
							

			if(z>0.5) then
			 
				S(VE(a))=S(VE(a))+1			

			else 	
				S(VD(a))=S(VD(a))+1
		
			end if

		else if(minval(J)==J(1)) then 

			S(VE(a))=S(VE(a))+1
			
		else if(minval(J)==J(3)) then
 					
			S(VD(a))=S(VD(a))+1
	
		end if
	end subroutine
!================================================================================================!
	subroutine troca_var	!!aqui haverá a troca de varáveis
		 
		if(L/=0) then	!!só entra aq quando já tiver algum valor em L
			
			deallocate(S,WA,W,VE,VD,STAT=erro)
		
		end if			
			
		L=128*i2
		mol=100000*L
		allocate(S(1:L),WA(1:L),W(1:mol/L),VE(1:L),VD(1:L))
		
	end subroutine
!================================================================================================!
	subroutine programa !!programa a parte para a rugosidade
		
		vez=1	!!número máximo de camadas
		
		do i2=1, 1 !!este do auxilia na troca de L

			call troca_var	!!aqui haverá a troca do !!

			call arquivos !!configura os melhores arquivos
		
			call reg_var  !!deixa o substrato circular

			do am=1, vez		

                		WA=0; S=0; ca=0

				do n=1, mol 

					if (rcont.ge.1e7) then
    	 					rcont=0
    	 					iseed=a1*jseed
    	 					jseed=iseed
 					end if

					iseed=(AMAG*iseed)+BMAG; rcont=rcont+1; z=r_1*iseed+0.5d0
					a=int(z*L)+1

					call depo	!!regra da deposição                   
			
				end do
		
			end do

			call dados !!coleta dos dados
 
		end do		 
	
	end subroutine
!=================================================================================================!		
end module globalDAR

program DAR
	use globalDAR
	call programa	!!cálculo da rugosidade
end program DAR
