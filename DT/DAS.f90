module globalDAS

	!este programa foi feito por base do modelo de exemplo do cap. 4 do livro barabasi
	!========================================================================================================
	integer(4)::vez		!vezes que o programa vai roDAS. e o número de amos	
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
			open(1,file='perfilDASL=128mol=100000Lamo=100000.dat')
			open(2,file='wDASL=128mol=100000Lamo=100000.dat')
		else if(L==256) then
			open(1,file='perfilDASL=256mol=100000Lamo=100000.dat')
			open(2,file='wDASL=256mol=100000Lamo=1000000.dat')
		else if(L==384) then
			open(1,file='perfilDASL=384mol=100000Lamo=1000000.dat')
			open(2,file='wDASL=384mol=100000Lamo=1000000.dat')
		else if(L==512) then
			open(1,file='perfilDASL=512mol=100000Lamo=1000000.dat')
			open(2,file='wDASL=512mol=100000Lamo=100000.dat')
		else if(L==640) then
			open(1,file='perfilDASL=640mol=100000Lamo=1000000.dat')
			open(2,file='wDASL=640mol=100000Lamo=100000.dat')
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
				
		do n=1, L
			write(1,*) n,S(n)
		end do		
			
		k=1
       		do n=1, mol/L   
                	if(n==k) then
	        		write(2,*) n, W(n)
                	k=int(k*m)+1
                	end if
		end do
		!***********************************!
	end subroutine
!================================================================================================!
	subroutine depo		!!aqui ocorre as deposições
		!**************************************************************
		lc=0; le=0; ld=0

		if(S(VE(a))>S(a)) then !quanto de ligação no meio
			lc=lc+1
		end if
		if(S(VD(a))>S(a)) then
			lc=lc+1
		end if
			
		if(S(VE(a)-1)>S(VE(a))) then !quanto de ligação na esquerda
			le=le+1
		end if
		if(S(a)>S(VE(a))) then
			le=le+1
		end if
		
		if(S(VD(a)+1)>S(VD(a))) then !quanto de ligação na direita
			ld=ld+1
		end if
		if(S(a)>S(VD(a))) then
			ld=ld+1
		end if
		!**************************************************************
		
		if((lc>=1).or.((lc+ld+le)==0)) then	!!

			S(a)=S(a)+1

		else if((le>=1).and.(ld>=1)) then !os do lado têm a mesma ligação

			iseed=(AMAG*iseed)+BMAG; rcont=rcont+1; z=r_1*iseed+0.5d0

			if(z>0.5) then !!funfa!!
		
				S(VE(a))=S(VE(a))+1

			else	!!funfa!!		

				S(VD(a))=S(VD(a))+1

			end if

		else if(le>0) then	!!funfa!! !!quando a ligação do lado esquerdo for maior

			S(VE(a))=S(VE(a))+1

		else if(ld>0) then	!!funfa!! !!quando a ligação do lado direito for maior

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
		
		vez=100000	!!número máximo de camadas
		
		do i2=1, 5 !!este do auxilia na troca de L

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
end module globalDAS

program DAS
	use globalDAS
	!call programa_perfil    !!cáculo dos perfis
	call programa	!!cálculo da rugosidade
end program DAS
