subroutine altitude(nlon, nlat, niv, relief, psol, aaa, bbb, t, h2o, &
                    trajlon, trajlat, trajpre, alt)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                              !
!     calcul de l'altitude geometrique en km                   !
!                                                              !
! d'apres http://mtp.jpl.nasa.gov/notes/altitude/altitude.html !
!                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none
 
integer, intent(in)                         :: nlon, nlat, niv
integer                                     :: ilon, ilat, iniv
integer                                     :: ilonreg, ilatreg, inivreg
integer, dimension(nlon,nlat,niv)           :: indlon, indlat, indpre

real, dimension(nlon,nlat,niv), intent(in)  :: h2o, t
real, dimension(nlon,nlat,niv), intent(in)  :: trajlon, trajlat, trajpre
real, dimension(nlon,nlat)    , intent(in)  :: psol, relief
real, dimension(niv)          , intent(in)  :: aaa, bbb
real, dimension(nlon,nlat,niv), intent(out) :: alt

real, dimension(nlon,nlat,niv)              :: p, georeg, geo, tv
real, dimension(nlon,nlat,niv)              :: pirr
real                                        :: ctegas, dz, pi, xlat

! constante de l'air sec en m2 s-2 k-1, divisee par
! g a 45 degres

ctegas = 287.04/9.80665
pi = 2.*asin(1.)

! pressions de la grille reguliere
!

do iniv = 1,niv
   do ilat = 1,nlat
      do ilon = 1,nlon
         p(ilon,ilat,iniv) = aaa(iniv) + bbb(iniv)*psol(ilon,ilat)
      end do
   end do
end do

! calcul de la temperature virtuelle (avec conversion rapport
! de melange -> humidite specifique) sur la grille reguliere

tv = t*(1. + 0.608*h2o*18./28.97)

! geopotentiel sur la grille reguliere

! dernier niveau

georeg(:,:,niv) = relief(:,:) + ctegas*tv(:,:,niv) &
                 *log(psol(:,:)/p(:,:,niv))*1.e-3

! cas general

do iniv = niv-1, 1, -1
   georeg(:,:,iniv) =  georeg(:,:,iniv+1)                          &
                    + ctegas*((tv(:,:,iniv+1) + tv(:,:,iniv))/2.)  &
                    *log(p(:,:,iniv+1)/p(:,:,iniv))*1.e-3
end do

! indlon, indlat: indices du point le plus proche sur la grille reguliere

do iniv = 1,niv
   do ilat = 1,nlat
      do ilon = 1,nlon

         indlon(ilon,ilat,iniv) = nint(trajlon(ilon,ilat,iniv)/2.) + 1
         if (indlon(ilon,ilat,iniv) > nlon) then
            indlon(ilon,ilat,iniv) = 1
         end if
         indlat(ilon,ilat,iniv) = nint((90. - trajlat(ilon,ilat,iniv))/2.) + 1

      end do
   end do
end do

do iniv = 1,niv
   do ilat = 1,nlat
      do ilon = 1,nlon

         ilonreg = indlon(ilon,ilat,iniv)
         ilatreg = indlat(ilon,ilat,iniv)

         pirr(ilon,ilat,iniv) = max(trajpre(ilon,ilat,iniv), p(ilonreg,ilatreg,1))
         pirr(ilon,ilat,iniv) = min(trajpre(ilon,ilat,iniv), p(ilonreg,ilatreg,niv))

      end do
   end do
end do


! on attribue le geopotentiel jusqu'au niveau inivreg immediatement en-dessous,
! puis on ajoute dz

do iniv = 1,niv
   do inivreg = niv,2,-1
      do ilat = 1,nlat
         do ilon = 1,nlon

            ilonreg = indlon(ilon,ilat,iniv)
            ilatreg = indlat(ilon,ilat,iniv)

            if (p(ilonreg,ilatreg,inivreg) >= pirr(ilon,ilat,iniv)) then
               geo(ilon,ilat,iniv) = georeg(ilonreg,ilatreg,inivreg)
               indpre(ilon,ilat,iniv) = inivreg 
            end if

         end do
      end do
   end do
end do 

do iniv = 1,niv
   do ilat = 1,nlat
      do ilon = 1,nlon

         ilonreg = indlon(ilon,ilat,iniv)
         ilatreg = indlat(ilon,ilat,iniv)
         inivreg = indpre(ilon,ilat,iniv)

         dz = ctegas*(tv(ilonreg,ilatreg,inivreg)                     &
                    + tv(ilonreg,ilatreg,inivreg-1))/2.               &
              *log(p(ilonreg,ilatreg,inivreg)/pirr(ilon,ilat,iniv))   &
              *1.e-3
         geo(ilon,ilat,iniv) = geo(ilon,ilat,iniv) + dz
      end do
   end do
end do

! conversion geopotentiel -> altitude geometrique

do iniv = 1,niv
   do ilat = 1,nlat
      do ilon = 1,nlon
         xlat = trajlat(ilon,ilat,iniv)*pi/180.
         alt(ilon,ilat,iniv) = (1. + 0.002644*cos(2.*xlat))*geo(ilon,ilat,iniv)   &
                             + (1. + 0.0089*cos(2.*xlat))*geo(ilon,ilat,iniv)**2  &
                                                         /6245.
      end do
   end do
end do

end subroutine altitude
