subroutine jno(cc_no, hnm, nlon, nivbas, o3t, sza, pm, tjno)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                       !
! parametrisation de la photodissociation de no         !
! d'apres minschwaner and siskind, a new calculation    !
! of nitric oxide photolysis in the stratosphere,       !
! mesosphere, and lower thermosphere, j. geophys. res., !
! 98, 20401-20412, 1993                                 !
!                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none

integer, intent(in)                       :: nlon, nivbas
integer                                   :: iniv, ilon

real, dimension(nlon,nivbas), intent(in)  :: cc_no, hnm, o3t, pm, sza
real, dimension(nlon,nivbas), intent(out) :: tjno

real, parameter                           :: factor = 3.141592/180.
real, dimension(nlon,nivbas)              :: colo2, colno, colo3
real, dimension(nlon,nivbas)              :: p, t_o3_1, t_o3_2, t_o3_3
real, dimension(nlon,nivbas)              :: bd1, bd2, bd3
real, dimension(6,3)                      :: sigma_o2, sigma_no_1, sigma_no_2
real, dimension(6,3)                      :: w_no_1, w_no_2
real, dimension(3)                        :: delta_lambda, i0, sigma_o3
real                                      :: a, d, kq, n2
real                                      :: r_o2, sec, zcmt, dp

! sections efficaces

sigma_o2(:,1) = (/1.12e-23, 2.45e-23, 7.19e-23, &
                  3.04e-22, 1.75e-21, 1.11e-20/)
sigma_o2(:,2) = (/1.35e-22, 2.99e-22, 7.33e-22, &
                  3.07e-21, 1.69e-20, 1.66e-19/)
sigma_o2(:,3) = (/2.97e-22, 5.83e-22, 2.05e-21, &
                  8.19e-21, 4.80e-20, 2.66e-19/)

sigma_no_1(:,1) = (/0.00e+00, 1.32e-18, 6.35e-19, &
                    7.09e-19, 2.18e-19, 4.67e-19/)
sigma_no_1(:,2) = (/0.00e+00, 0.00e+00, 3.05e-21, &
                    5.76e-19, 2.29e-18, 2.21e-18/)
sigma_no_1(:,3) = (/1.80e-18, 1.50e-18, 5.01e-19, &
                    7.20e-20, 6.72e-20, 1.49e-21/)

sigma_no_2(:,1) = (/0.00e+00, 4.41e-17, 4.45e-17, &
                    4.50e-17, 2.94e-17, 4.35e-17/)
sigma_no_2(:,2) = (/0.00e+00, 0.00e+00, 3.20e-21, &
                    5.71e-17, 9.09e-17, 6.00e-17/)
sigma_no_2(:,3) = (/1.40e-16, 1.52e-16, 7.00e-17, &
                    2.83e-17, 2.73e-17, 6.57e-18/)

! facteurs de ponderation

w_no_1(:,1) = (/0.00e+00, 5.12e-02, 1.36e-01, &
                1.65e-01, 1.41e-01, 4.50e-02/)
w_no_1(:,2) = (/0.00e+00, 0.00e+00, 1.93e-03, &
                9.73e-02, 9.75e-02, 3.48e-02/)
w_no_1(:,3) = (/4.50e-02, 1.80e-01, 2.25e-01, &
                2.25e-01, 1.80e-01, 4.50e-02/)

w_no_2(:,1) = (/0.00e+00, 5.68e-03, 1.52e-02, &
                1.83e-02, 1.57e-02, 5.00e-03/)
w_no_2(:,2) = (/0.00e+00, 0.00e+00, 2.14e-04, &
                1.08e-02, 1.08e-02, 3.86e-03/)
w_no_2(:,3) = (/5.00e-03, 2.00e-02, 2.50e-02, &
                2.50e-02, 2.00e-02, 5.00e-03/)

! largeurs spectrales pour les trois bandes considerees (nm)

delta_lambda = (/2.3, 1.5, 1.5/)

! flux pour les trois bandes considerees (s-1 cm-2 nm-1)

i0 = (/3.98e+11, 2.21e+11, 2.30e+11/)

! sections efficaces de l'ozone pour les trois bandes
! considerees (cm2) d'apres wmo 1985

sigma_o3 = (/4.6e-19, 6.7e-19, 7.1e-19/)
 
! zcmt: gas constant/boltzmann constant*g

zcmt = 287.0/(1.38e-23*9.81)

! d: taux de predissociation spontanee (s-1)

d = 1.65e+09

! a: taux d'emission spontanee (s-1)

a = 5.1e+07

! kq: quenching rate constant (cm3 s-1)

kq = 1.5e-09

! r_o2: rapport de melange de l'oxygene

r_o2 = 0.2095

!=============================================================
!  calcul des colonnes de o2 et no
!=============================================================

!  premier niveau du modele (mol/cm2)

colo2(:,1) = 6.8e+05*r_o2*hnm(:,1)
colno(:,1) = 6.8e+05*cc_no(:,1)
colo3(:,1) = o3t(:,1)

!  cas general

do iniv = 2,nivbas
   do ilon = 1,nlon
      dp = (pm(ilon,iniv) - pm(ilon,iniv - 1))*100.
      dp = max(dp, 0.)
      colo2(ilon,iniv) = colo2(ilon,iniv - 1)                 &
               + zcmt*(r_o2 + r_o2)*0.5*dp*1.e-4
      colno(ilon,iniv) = colno(ilon,iniv - 1) + zcmt          &
               *(cc_no(ilon,iniv - 1)/hnm(ilon,iniv - 1)      &
               + cc_no(ilon,iniv)/hnm(ilon,iniv))*0.5*dp*1.e-4
      colo3(ilon,iniv) = o3t(ilon,iniv)
   end do
end do 

!=============================================================
!  boucle sur les niveaux
!=============================================================

do iniv = 1,nivbas

!=============================================================
!  boucle sur les longitudes
!=============================================================

   do ilon = 1,nlon
      if (sza(ilon,1) <= 89.0) then

         sec = 1./cos(sza(ilon,1)*factor)

         colo2(ilon,iniv) = colo2(ilon,iniv)*sec
         colno(ilon,iniv) = colno(ilon,iniv)*sec
         colo3(ilon,iniv) = colo3(ilon,iniv)*sec

!        facteurs de transmission de l'ozone

         t_o3_1(ilon,iniv) = exp(-sigma_o3(1)*colo3(ilon,iniv))
         t_o3_2(ilon,iniv) = exp(-sigma_o3(2)*colo3(ilon,iniv))
         t_o3_3(ilon,iniv) = exp(-sigma_o3(3)*colo3(ilon,iniv))

!        calcul de la probabilite de predissociation

         n2 = hnm(ilon,iniv)*0.78
         p(ilon,iniv) = d/(a + d + kq*n2)

      end if
   end do

!  calcul proprement dit, pour les 3 bandes
!  on boucle sur chacune des bandes pour des raisons de vectorisation

   do ilon = 1,nlon
      if (sza(ilon,1) <= 89.0) then

         bd1(ilon,iniv) = delta_lambda(1)*i0(1)*t_o3_1(ilon,iniv)*p(ilon,iniv)*(  &
                exp(-sigma_o2(1,1)*colo2(ilon,iniv))   &
              *(w_no_1(1,1)*sigma_no_1(1,1)*exp(-sigma_no_1(1,1)*colno(ilon,iniv))  &
              + w_no_2(1,1)*sigma_no_2(1,1)*exp(-sigma_no_2(1,1)*colno(ilon,iniv))) &
              + exp(-sigma_o2(2,1)*colo2(ilon,iniv))   &
              *(w_no_1(2,1)*sigma_no_1(2,1)*exp(-sigma_no_1(2,1)*colno(ilon,iniv))  &
              + w_no_2(2,1)*sigma_no_2(2,1)*exp(-sigma_no_2(2,1)*colno(ilon,iniv))) &
              + exp(-sigma_o2(3,1)*colo2(ilon,iniv))   &
              *(w_no_1(3,1)*sigma_no_1(3,1)*exp(-sigma_no_1(3,1)*colno(ilon,iniv))  &
              + w_no_2(3,1)*sigma_no_2(3,1)*exp(-sigma_no_2(3,1)*colno(ilon,iniv))) &
              + exp(-sigma_o2(4,1)*colo2(ilon,iniv))   &
              *(w_no_1(4,1)*sigma_no_1(4,1)*exp(-sigma_no_1(4,1)*colno(ilon,iniv))  &
              + w_no_2(4,1)*sigma_no_2(4,1)*exp(-sigma_no_2(4,1)*colno(ilon,iniv))) &
              + exp(-sigma_o2(5,1)*colo2(ilon,iniv))   &
              *(w_no_1(5,1)*sigma_no_1(5,1)*exp(-sigma_no_1(5,1)*colno(ilon,iniv))  &
              + w_no_2(5,1)*sigma_no_2(5,1)*exp(-sigma_no_2(5,1)*colno(ilon,iniv))) &
              + exp(-sigma_o2(6,1)*colo2(ilon,iniv))  &
              *(w_no_1(6,1)*sigma_no_1(6,1)*exp(-sigma_no_1(6,1)*colno(ilon,iniv))  &
              + w_no_2(6,1)*sigma_no_2(6,1)*exp(-sigma_no_2(6,1)*colno(ilon,iniv))) &
              )
      end if
   end do

   do ilon = 1,nlon
      if (sza(ilon,1) <= 89.0) then

         bd2(ilon,iniv) = delta_lambda(2)*i0(2)*t_o3_2(ilon,iniv)*p(ilon,iniv)*(  &
                exp(-sigma_o2(1,2)*colo2(ilon,iniv))   &
              *(w_no_1(1,2)*sigma_no_1(1,2)*exp(-sigma_no_1(1,2)*colno(ilon,iniv))  &
              + w_no_2(1,2)*sigma_no_2(1,2)*exp(-sigma_no_2(1,2)*colno(ilon,iniv))) &
              + exp(-sigma_o2(2,2)*colo2(ilon,iniv))   &
              *(w_no_1(2,2)*sigma_no_1(2,2)*exp(-sigma_no_1(2,2)*colno(ilon,iniv))  &
              + w_no_2(2,2)*sigma_no_2(2,2)*exp(-sigma_no_2(2,2)*colno(ilon,iniv))) &
              + exp(-sigma_o2(3,2)*colo2(ilon,iniv))   &
              *(w_no_1(3,2)*sigma_no_1(3,2)*exp(-sigma_no_1(3,2)*colno(ilon,iniv))  &
              + w_no_2(3,2)*sigma_no_2(3,2)*exp(-sigma_no_2(3,2)*colno(ilon,iniv))) &
              + exp(-sigma_o2(4,2)*colo2(ilon,iniv))   &
              *(w_no_1(4,2)*sigma_no_1(4,2)*exp(-sigma_no_1(4,2)*colno(ilon,iniv))  &
              + w_no_2(4,2)*sigma_no_2(4,2)*exp(-sigma_no_2(4,2)*colno(ilon,iniv))) &
              + exp(-sigma_o2(5,2)*colo2(ilon,iniv))   &
              *(w_no_1(5,2)*sigma_no_1(5,2)*exp(-sigma_no_1(5,2)*colno(ilon,iniv))  &
              + w_no_2(5,2)*sigma_no_2(5,2)*exp(-sigma_no_2(5,2)*colno(ilon,iniv))) &
              + exp(-sigma_o2(6,2)*colo2(ilon,iniv))  &
              *(w_no_1(6,2)*sigma_no_1(6,2)*exp(-sigma_no_1(6,2)*colno(ilon,iniv))  &
              + w_no_2(6,2)*sigma_no_2(6,2)*exp(-sigma_no_2(6,2)*colno(ilon,iniv))) &
              )
      end if
   end do

   do ilon = 1,nlon
      if (sza(ilon,1) <= 89.0) then

         bd3(ilon,iniv) = delta_lambda(3)*i0(3)*t_o3_3(ilon,iniv)*p(ilon,iniv)*(  &
                exp(-sigma_o2(1,3)*colo2(ilon,iniv))   &
              *(w_no_1(1,3)*sigma_no_1(1,3)*exp(-sigma_no_1(1,3)*colno(ilon,iniv))  &
              + w_no_2(1,3)*sigma_no_2(1,3)*exp(-sigma_no_2(1,3)*colno(ilon,iniv))) &
              + exp(-sigma_o2(2,3)*colo2(ilon,iniv))   &
              *(w_no_1(2,3)*sigma_no_1(2,3)*exp(-sigma_no_1(2,3)*colno(ilon,iniv))  &
              + w_no_2(2,3)*sigma_no_2(2,3)*exp(-sigma_no_2(2,3)*colno(ilon,iniv))) &
              + exp(-sigma_o2(3,3)*colo2(ilon,iniv))   &
              *(w_no_1(3,3)*sigma_no_1(3,3)*exp(-sigma_no_1(3,3)*colno(ilon,iniv))  &
              + w_no_2(3,3)*sigma_no_2(3,3)*exp(-sigma_no_2(3,3)*colno(ilon,iniv))) &
              + exp(-sigma_o2(4,3)*colo2(ilon,iniv))   &
              *(w_no_1(4,3)*sigma_no_1(4,3)*exp(-sigma_no_1(4,3)*colno(ilon,iniv))  &
              + w_no_2(4,3)*sigma_no_2(4,3)*exp(-sigma_no_2(4,3)*colno(ilon,iniv))) &
              + exp(-sigma_o2(5,3)*colo2(ilon,iniv))   &
              *(w_no_1(5,3)*sigma_no_1(5,3)*exp(-sigma_no_1(5,3)*colno(ilon,iniv))  &
              + w_no_2(5,3)*sigma_no_2(5,3)*exp(-sigma_no_2(5,3)*colno(ilon,iniv))) &
              + exp(-sigma_o2(6,3)*colo2(ilon,iniv))  &
              *(w_no_1(6,3)*sigma_no_1(6,3)*exp(-sigma_no_1(6,3)*colno(ilon,iniv))  &
              + w_no_2(6,3)*sigma_no_2(6,3)*exp(-sigma_no_2(6,3)*colno(ilon,iniv))) &
              )

         tjno(ilon,iniv) = max(bd1(ilon,iniv) &
                             + bd2(ilon,iniv) &
                             + bd3(ilon,iniv), 1.e-30)

      end if
   end do

end do

end subroutine jno
