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
	real(8),parameter::m=1.19
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
			
			open(1,file='perfilDAL=128mol=100000Lamo=10000.dat')
			open(2,file='wDAL=128mol=100000Lamo=10000.dat')
			
		else if(L==256) then
				
			open(1,file='perfilDAL=256mol=100000Lamo=10000.dat')
			open(2,file='wDAL=256mol=100000Lamo=10000.dat')
			
		else if(L==384) then
				
			open(1,file='perfilDAmol=100000Lamo=10000dat')
			open(2,file='wDAL=384mol=100000Lamo=10000dat')	
			
		else if(L==512) then
				
			open(1,file='perfilDAmol=100000Lamo=10000.dat')
			open(2,file='wDAL=512mol=100000Lamo=10000.dat')
			
		else if(L==640) then
				
			open(1,file='perfilDAmol=100000Lamo=10000.dat')
			open(2,file='wDAL=1024mol=100000Lamo=10000.dat')
			
		end if
		
	end subroutine
!================================================================================================!
	subroutine calc_rug !!cáculo da rugosidade

		if(mod(n,L)==0) then !se o resto for 0 ele caracteriza a rugosidade
                	ca=ca+1
	                !=============================================================================================!
                    	!				cálculo da rugosidade de cada camada			      !
                    	!=============================================================================================!	       
			h=0; h=sum(S)/L	!média da altura da soma dos vetores
		        do i=1, L
				WA(i)=((S(i)-h)**2)	!primeira parte da esquação de rugosidade	
		        end do
			W(ca)=((sum(WA))/L)+W(ca) !aqui está contido a rugosidade de cada amostra com suas camadas
	                !==============================================================================================!
		end if

	end subroutine
!================================================================================================!	
	subroutine dados !!armazenamento dos dados
		!*****************************!
		if(vez==1) then		
			do n=1, L
				write(1,*) n,S(n)
			end do
		else		
			k=1
       			do n=1, mol/L   
                		if(n==k) then
	        		        write(2,*) n, W(n)
                		k=int(k*m)+1
                		end if
			end do
		end if
		!***********************************!
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
		
		vez=10000	!!número máximo de camadas
		
		do i2=1, 5 !!este do auxilia na troca de L

			call troca_var	!!aqui haverá a troca do !!

			call arquivos !!configura os melhores arquivos

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

					S(a)=S(a)+1
					

                	       		call calc_rug                    
			
				end do
		
			end do

			!=============================================!          
    			!		cálculo final da camadas      !       
			!=============================================!	
			do ca=1, mol/L 
        			W(ca)=sqrt(W(ca))
			end do
   			!=============================================!

			call dados !!coleta dos dados
 
		end do		 
	
	end subroutine
!=================================================================================================!
	subroutine programa_perfil !!programa a parte para o perfil

		vez=1
		
		do i2=1, 1 !!este do auxilia na troca de L	

			call troca_var	
			
			call arquivos !!configura os melhores arquivos
		
			WA=0; S=0; ca=0			
			
			do n=1, mol 
			
				if (rcont.ge.1e7) then
    	 				rcont=0
    	 				iseed=a1*jseed
    	 				jseed=iseed
 				end if
			

				iseed=(AMAG*iseed)+BMAG; rcont=rcont+1; z=r_1*iseed+0.5d0
				a=int(z*L)+1

				S(a)=S(a)+1
                   
			end do
			
			call dados !!cálculos de dados 
		
		end do 
		
	end subroutine
!=================================================================================================!		
end module globalDAR

program DAR
	use globalDAR
	call programa_perfil    !!cáculo dos perfis
	!call programa	!!cálculo da rugosidade
end program DAR
