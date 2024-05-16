      program reprobus
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     REPROBUS - Three-dimensional chemistry-transport model of     c
c                the stratosphere.                                  c
c                                                                   c
c     PURPOSE                                                       c
c     -------                                                       c
c           The REPROBUS model describes the distribution and       c
c     evolution of the main stratospheric species, using a compre-  c
c     hensive photochemical scheme coupled to a semi-lagrangian     c
c     transport model. Temperatures and winds driving the photo-    c
c     chemistry and transport are taken from the ECMWF analysis.    c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      parameter(nlon = 180, nlat = 91, niv = 137, nivbas = 94)
      parameter(nivh2oecmwf = 59)
      parameter(lolani = nlon*nlat*niv, lola = nlon*nlat)
      parameter(nbcon = 43, ncm = 15)
      parameter(nphot = 41, nz = 101, nsza = 27, nozo = 7)
c
      parameter(nrappmax = 12)
c
      real dt, adapdt
      real qj1(nlon,nlat,niv,nbcon), qj1m(nlon,nlat,niv,nbcon)
      real hc(nlon,nlat,niv,ncm), hcm(nlon,nlat,niv,ncm)
      real alt(nlon,nlat,niv), altm(nlon,nlat,niv)
      real relief(nlon,nlat)
      real vsed3d(nlon,nlat,nivbas), sza3d(nlon,nlat,nivbas)
      real parthno3sed(nlon,nlat,nivbas), parth2osed(nlon,nlat,nivbas)
      real surfarea(nlon,nlat,niv), surfaream(nlon,nlat,niv)
      real trajlon(nlon,nlat,niv,nrappmax)
      real trajlat(nlon,nlat,niv,nrappmax)
      real trajpre(nlon,nlat,niv,nrappmax)
      real trajtem(nlon,nlat,niv,nrappmax)
      real pj1m(nlon,nlat)
      real uj1m(nlon,nlat,niv), vj1m(nlon,nlat,niv)
      real tj1m(nlon,nlat,niv)
      real aa(niv+1), bb(niv+1)
c
      integer idate(18)
c
      character*6 namexp
c
      common /constant/
     $     a0hox, a1hox, a0nox1,
     $     a1nox1, a0nox2, a1nox2, a0nox3, a1nox3,
     $     a0clx1, a1clx1, a0clx2, a1clx2, a0clx3,
     $     a1clx3, a0brx, a1brx, a0c1, a1c1
c
      common /photodis/ ajl(nphot-1,nz,nsza,nozo), o3up(0:79)
c
      common/grille/dlon(0:nlon+2),ddlon,dlat(0:nlat+1),ddlat
     +             ,slon(nlon),clon(nlon),slat(nlat),clat(nlat)
     +             ,slonr(nlon,nlat),clonr(nlon,nlat)
     +             ,slatr(nlon,nlat),clatr(nlon,nlat)
     +             ,xcoef(-1:2),ycoef(-1:2),ip(nlon)
c
      common/verniv/ aaa(niv), bbb(niv)
c
      common/forcj1/uj1(nlon,nlat,niv),vj1(nlon,nlat,niv)
     +             ,wj1(nlon,nlat,niv),tj1(nlon,nlat,niv)
     +             ,hj1(nlon,nlat,niv),pj1(nlon,nlat)
c
      common/forcm1/um1(nlon,nlat,niv),vm1(nlon,nlat,niv)
     +             ,wm1(nlon,nlat,niv),tm1(nlon,nlat,niv)
     +             ,hm1(nlon,nlat,niv),pm1(nlon,nlat),daym1
c
      common/forcp1/up1(nlon,nlat,niv),vp1(nlon,nlat,niv)
     +             ,wp1(nlon,nlat,niv),tp1(nlon,nlat,niv)
     +             ,hp1(nlon,nlat,niv),pp1(nlon,nlat),dayp1
c
      common/cstes/ gg, ra, cp, omega, rter, xp00, rter2, deuomg
     +            , xpi, conv, xpi2, xpih, dlatlim
c
c     coefficients aaa et bbb des niveaux hybrides
c
      open(unit = 20, file = 'ecmwf_137_levels.txt', form = 'formatted')
c
      read(20,*)
      read(20,*)
      do iniv = 1,niv+1
         read(20,*) ibid, aa(iniv), bb(iniv)
      end do
      close(20)
c
      do iniv = 1,niv
         aaa(iniv) = (aa(iniv) + aa(iniv+1))*0.5*0.01
         bbb(iniv) = (bb(iniv) + bb(iniv+1))*0.5
      end do
c
c     call profile
c
ccc   ouverture du fichier historique
c
      open(unit = 78, file = 'history', form = 'unformatted')
      rewind(78)
c
ccc   ouverture du fichier stations
c
      open(unit = 75, file = 'stations.txt', form = 'formatted')
      rewind(75)
c
cccccccccccccccccccccccc tableau de bord cccccccccccccccccccccccccccc
c                                                                   c
c     maxdays    : nombre de jours d'integration                    c
c                                                                   c
c     nstart     : conditions initiales                             c
c                =0 lecture de moyennes zonales                     c
c                =1 lecture d'un restart                            c
c                                                                   c
c     init_o3    : initialisation de o3 et du traceur               c
c                =0 d'apres le restart                              c
c                =1 d'apres l'analyse ecmwf                         c
c                                                                   c
c     init_h2so4 : forcage de h2so4                                 c
c                =0 d'apres le restart                              c
c                =1 reinitialisation chaque mois d'apres modele 2d  c
c                                                                   c
c     intmean    : sauvegarde des moyennes sur l'integration        c
c                =0 non                                             c
c                =1 oui                                             c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      maxdays    = 1
      nstart     = 1
      init_o3    = 0
      init_h2so4 = 1
      intmean    = 0
c                                                                   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     ntime    : nombre de pas de temps par heure                   c
c     lengdy   : nombre de pas de temps par jour                    c
c     nend     : nombre total de pas de temps                       c
c     nout     : frequence d'ecriture des resultats (pas de temps)  c 
c     nstep    : compteur de pas de temps                           c
c     dtdyn    : pas de temps pour l'advection                      c 
c     dt       : pas de temps pour la chimie                        c
c     nrapp    : rapport entre les deux pas de temps (dtdyn/dt)     c
c     nintecmwf: intervalle (en pas de temps) entre deux echeances  c
c                de forcage                                         c
c     dtsedi   : pas de temps pour le calcul de la sedimentation    c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ntime      = 4
      lengdy     = ntime*24
      nend       = maxdays*lengdy
      nout       = maxdays*lengdy
      nstep      = 0
      dt         = 3600./real(ntime)
      nrapp      = min(nrappmax,nend)
      dtdyn      = dt*real(nrapp)
      nintecmwf  = 12
      nrappsedi  = 12
      dtsedi     = dt*real(nrappsedi)
      adapdt     = 0
      nech       = 0 
c
      print *
      print *,'le pas de temps chimique est de ',dt/60.,' minutes.'
      print *,'le pas de temps dynamique est de ',dtdyn/60.,' minutes.'
      print *,'le pas de temps de la sedimentation est de ',
     $         dtsedi/60.,' minutes.'
      print *
c
      if( nstart .eq. 0 ) then
c
ccc      date initiale de la simulation
c
         ian    = 2022
         imois  = 12
         ijour  = 01
         iheure = 12
         imin   = 0 
c
ccc      calcul du jour julien
c
         call gregtojul(ian,imois,ijour,iheure,imin,day0)
c
c        initialisation: lecture du modele 2d
c
         namexp = '001442'
         print*, 'numero experience : ', namexp
c
         print*
         print*,'initialisation : ouverture fichier 2d'
         open(unit = 13, file = 'qinit2d.txt',form = 'formatted')
         read(13,1301)
         read(13,1301)(((qj1(1,ilat,iniv,nb),ilat = 1,nlat),
     $                   iniv = 1,niv),nb = 1,nbcon)
         print*,'lecture constituants longue duree de vie: ok'
         read(13,1301)(((hc(1,ilat,iniv,nb),ilat = 1,nlat),
     $                   iniv = 1,niv),nb = 1,ncm)
         print*,'lecture constituants courte duree de vie: ok'
         print *
 1301    format(1x,6e13.5)
         close(13)
c
         do iniv = 1,niv
            do ilat = 1,nlat
               do nb = 1,nbcon - 1
                  do ilon = 1,nlon
                     qj1(ilon,ilat,iniv,nb) = max(
     $                                        qj1(1,ilat,iniv,nb),
     $                                        1.e-30)
                  end do
               end do
               do nb = 1,ncm
                  do ilon = 1,nlon
                     hc(ilon,ilat,iniv,nb) = max(hc(1,ilat,iniv,nb),
     $                                           1.e-30)
                  end do
               end do
            end do
         end do
c
ccc      initialisation: h2so4
c
         if (init_h2so4 .eq. 1) then
            call h2so4_2d(ian,imois,ijour,iheure,imin,qj1)
         end if
c
      else
c
c        initialisation: lecture d'un restart 3d
c
         open(unit = 90, form = 'unformatted')
         read(90) namexp, ian, imois, ijour, iheure, imin, 
     $            pj1, uj1, vj1, alt, tj1, qj1, surfarea, hc
         print*,'lecture restart experience ',namexp,' ok'
         close(90)
         namexp = '001442'
         call gregtojul(ian,imois,ijour,iheure,imin,day0)
c
      end if
c
      if (init_o3 .eq. 1) then
c
c        lecture du champ d'ozone ecmwf
c
         open(unit = 13, file = 'ecmwf_o3_20230501',
     $        form = 'unformatted')
c
         read(13)
         read(13)
         print*
         print*,'lecture du champ ozone ecmwf:'
         print*
         do iniv = 1,niv
            read(13) idate
            write(6,*) (idate(i),i = 5,16)
            read(13) ((qj1(ilon,ilat,iniv,8),
     $                ilon = 1,nlon), ilat = 1,nlat)
         end do
         close(13)
c
c        conversion en rapport de melange volumique
c        et initialisation de o3 et de ox passif
c
         do iniv = 1,niv
            do ilat = 1,nlat
               do ilon = 1,nlon
                  qj1(ilon,ilat,iniv,8)  = qj1(ilon,ilat,iniv,8)*29./48.
                  hc(ilon,ilat,iniv,5)   = qj1(ilon,ilat,iniv,8)
                  qj1(ilon,ilat,iniv,11) = qj1(ilon,ilat,iniv,8)
               end do
            end do
         end do
      end if
c
      ian0    = ian
      imois0  = imois
      ijour0  = ijour
      iheure0 = iheure
      imin0   = imin
c
      write(6,*)
      write(6,101)' date initiale  : ',
     $            ian,imois,ijour,iheure,'h',imin,'mn'
      write(6,102)' jour julien    : ',day0
      write(6,*)
c
 101  format(a17,i4,i3,i3,i4,a1,i3,a2)
 102  format(a17,f10.6)
c
      daynum = day0
c
ccc   initialisation des constantes chimiques
c
      call const(relief)
c
ccc   initialisation du transport
c
      call initransp
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                              c
c     debut de la boucle temporelle                            c
c                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  100 continue
c
ccc   forcage vent, temperature et psol
c
      if (nstep .eq. 0 .or. mod(nstep,nintecmwf) .eq. 0) then
         call ecmwf(ian, daynum, ntime, nintecmwf)
      end if
c
ccc   forcage h2so4
c
      if (init_h2so4 .eq. 1) then
         call h2so4_2d(ian,imois,ijour,iheure,imin,qj1)
      end if
c
ccc   avance temporelle
c
      nstep = nstep + 1
      daynum = daynum + dt/86400.
c
ccc   modification eventuelle du pas de temps dynamique
c     si la duree de la simulation n'est pas un multiple
c     de dtdyn
cnamexp
      if (mod(nstep,nrapp).ne. 0 .and. nstep .eq. nend) then
         nrapp  = mod(nend,nrapp)
         dtdyn  = dt*real(nrapp)
         dtsedi = dtdyn
         adapdt = 1 
         print*,'adaptation du pas de temps dynamique'
         print*,'nouveau dtdyn = ',dtdyn/60.,' minutes'
      end if
c
      if (mod(nstep,nrapp).eq.0 .or. adapdt .eq. 1) then
c
ccc      calcul de la date
c       
         call jultogreg(ian,imois,ijour,iheure,imin,daynum)
c
ccc      transport semi-lagrangien
c
         print*
         write(6,103)nstep,' calcul du transport a la date : ',
     $                  ian,imois,ijour,iheure,'h',imin,'mn'
 103     format(i6,a33,i4,i3,i3,i4,a1,i3,a2)
c
         call semilag(daynum,dtdyn,nrapp,qj1,
     $                trajlon,trajlat,trajpre,trajtem,
     $                ntime,nintecmwf)
c
ccc      nrapp iterations chimiques le long des trajectoires
c
         daynum = daynum - dtdyn/86400.
c
         do irapp = 1,nrapp
c
            daynum = daynum + dt/86400.
c
ccc         calcul de la date
c    
            call jultogreg(ian,imois,ijour,iheure,imin,daynum)
c
            write(6,103)nstep,
     $                  ' calcul de la chimie a la date : ',
     $                  ian,imois,ijour,iheure,'h',imin,'mn'
c
ccc         relaxation des especes dans la troposphere
c
            call tropo(ian,imois,dt,trajpre,irapp,qj1,hc)
c
ccc         forcage de co dans la troposphere
c
            call co_tropo(imois,dt,qj1)
c
ccc         traitement de la vapeur d'eau tropospherique:
c           analyses ecmwf a partir du niveau nivh2oecmwf
c
            call ecmwfh2o (daynum, qj1, nivh2oecmwf,
     $                     ntime,nintecmwf)
c
ccc         calcul de l'altitude de chaque point de grille
c
            call instant1(daynum,tj1,tm1,daym1,tp1,dayp1,lolani,
     $                    ntime,nintecmwf)
            call instant2(daynum,pj1,pm1,daym1,pp1,dayp1,lola,
     $                    ntime,nintecmwf)
            call altitude(nlon, nlat, niv, relief, pj1, aaa, bbb,
     $                    tj1, qj1(:,:,:,3), trajlon(:,:,:,irapp), 
     $                    trajlat(:,:,:,irapp), trajpre(:,:,:,irapp),
     $                    alt)
c
ccc         appel de la chimie
c
            do ilat = 1,nlat
               call chem3d(ilat,
     $                     daynum,
     $                     dt,
     $                     alt,
     $                     trajlon,
     $                     trajlat,
     $                     trajpre,
     $                     trajtem,
     $                     irapp,
     $                     qj1,
     $                     hc,
     $                     vsed3d,
     $                     parth2osed,
     $                     parthno3sed,
     $                     surfarea,
     $                     sza3d)
            end do
c
c           sedimentation
c
            if (irapp .eq. nrapp) then
               call sedimentation(trajpre, trajtem, irapp,
     $                            qj1, vsed3d, nivh2oecmwf,
     $                            parth2osed, parthno3sed, dtsedi)
            end if
c
c           correction pour conservation
c
            call conserv(qj1, hc)
c
c           valeur minimale des traceurs a 1.e-30
c
            do n = 1,nbcon-1
               do iniv = 1,niv
                  do ilat = 1,nlat
                     do ilon = 1,nlon
                        qj1(ilon,ilat,iniv,n) 
     $                  = max(qj1(ilon,ilat,iniv,n), 1.e-30)
                     end do
                  end do
               end do
            end do
c
c           hno3 gas dans la tropo
c
            do iniv = nivbas+1,niv
               do ilat = 1,nlat
                  do ilon = 1,nlon
                     qj1(ilon,ilat,iniv,nbcon) = qj1(ilon,ilat,iniv,5)
                  end do
               end do
            end do
c
c           call odin(namexp, daynum, ian, imois, ijour,
c    $                iheure, imin, alt, sza3d, trajlon, trajlat,
c    $                trajpre, trajtem, irapp, qj1, hc)
c
         end do
      end if
c
ccc   calcul des colonnes au-dessus des stations choisies
c
      if(mod(nstep,nrapp).eq.0) then
         call stations(daynum, qj1, hc, trajpre, trajtem, nrapp)
      end if
c
ccc   extraction des profils au-dessus des stations choisies
c
      if(mod(nstep,nrapp).eq.0) then
         call profils(namexp, daynum, ian, imois, ijour,
     $                iheure, imin, alt, trajpre,
     $                trajtem, trajlon, trajlat, 
     $                qj1, hc, nrapp)
      end if
c
ccc   calcul des moyennes sur la duree de l'integration
c
      if (intmean .eq. 1) then
         if (mod(nstep,nrapp).eq.0) then
            nech = nech + 1
            do ilat = 1,nlat
               do ilon = 1,nlon
                  pj1m(ilon,ilat) = pj1m(ilon,ilat) + pj1(ilon,ilat) 
               end do
            end do
            do iniv = 1,niv
               do ilat = 1,nlat
                  do ilon = 1,nlon
                     uj1m(ilon,ilat,iniv) = uj1m(ilon,ilat,iniv)
     $                                    + uj1(ilon,ilat,iniv) 
                     vj1m(ilon,ilat,iniv) = vj1m(ilon,ilat,iniv)
     $                                    + vj1(ilon,ilat,iniv) 
                     tj1m(ilon,ilat,iniv) = tj1m(ilon,ilat,iniv)
     $                                    + tj1(ilon,ilat,iniv) 
                     altm(ilon,ilat,iniv) = altm(ilon,ilat,iniv)
     $                                    + alt(ilon,ilat,iniv) 
                     surfaream(ilon,ilat,iniv) = 
     $                                     surfaream(ilon,ilat,iniv)
     $                                   + surfarea(ilon,ilat,iniv) 
                  end do
               end do
            end do
            do is = 1,nbcon
               do iniv = 1,niv
                  do ilat = 1,nlat
                     do ilon = 1,nlon
                        qj1m(ilon,ilat,iniv,is) = 
     $                                qj1m(ilon,ilat,iniv,is)
     $                              + qj1(ilon,ilat,iniv,is)
                     end do
                  end do
               end do
            end do 
            do is = 1,ncm
               do iniv = 1,niv
                  do ilat = 1,nlat
                     do ilon = 1,nlon
                        hcm(ilon,ilat,iniv,is) = 
     $                               hcm(ilon,ilat,iniv,is)
     $                             + hc(ilon,ilat,iniv,is)
                     end do
                  end do
               end do
            end do 
         end if
         if (mod(nstep,nout) .eq. 0) then
            do ilat = 1,nlat
               do ilon = 1,nlon
                  pj1m(ilon,ilat) = pj1m(ilon,ilat)/real(nech)
               end do
            end do
            do iniv = 1,niv
               do ilat = 1,nlat
                  do ilon = 1,nlon
                     uj1m(ilon,ilat,iniv) = uj1m(ilon,ilat,iniv)
     $                                     /real(nech)
                     vj1m(ilon,ilat,iniv) = vj1m(ilon,ilat,iniv)
     $                                     /real(nech)
                     tj1m(ilon,ilat,iniv) = tj1m(ilon,ilat,iniv)
     $                                     /real(nech)
                     altm(ilon,ilat,iniv) = altm(ilon,ilat,iniv)
     $                                     /real(nech)
                     surfaream(ilon,ilat,iniv) = 
     $                                      surfaream(ilon,ilat,iniv)
     $                                     /real(nech)
                  end do
               end do
            end do
            do is = 1,nbcon
               do iniv = 1,niv
                  do ilat = 1,nlat
                     do ilon = 1,nlon
                        qj1m(ilon,ilat,iniv,is) = 
     $                                      qj1m(ilon,ilat,iniv,is)
     $                                     /real(nech)
                     end do
                  end do
               end do
            end do
            do is = 1,ncm
               do iniv = 1,niv
                  do ilat = 1,nlat
                     do ilon = 1,nlon
                        hcm(ilon,ilat,iniv,is) = hcm(ilon,ilat,iniv,is)
     $                                          /real(nech)
                     end do
                  end do
               end do
            end do
            write(76) namexp,
     $                ian0, imois0, ijour0, iheure0, imin0,
     $                ian, imois, ijour, iheure, imin,
     $                pj1m, uj1m, vj1m, altm, tj1m, 
     $                qj1m, surfaream, hcm
     
         end if
      end if 
c
ccc   ecriture du fichier historique
c
      if( mod(nstep,nout) .eq. 0 ) then
         call instant1(daynum,uj1,um1,daym1,up1,dayp1,lolani,
     $                 ntime,nintecmwf)
         call instant1(daynum,tj1,tm1,daym1,tp1,dayp1,lolani,
     $                 ntime,nintecmwf)
         call instant2(daynum,pj1,pm1,daym1,pp1,dayp1,lola,
     $                 ntime,nintecmwf)
c
         write(78) namexp, ian, imois, ijour, iheure, imin,
     $             pj1, uj1, vj1, alt, tj1, qj1, surfarea, hc
      end if
c
      if( nstep .lt. nend ) then
         go to 100
      end if
c
      stop
c
      end program reprobus
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine chem3d( ilat,
     $                   daynum,
     $                   dt,
     $                   alt,
     $                   trajlon,
     $                   trajlat,
     $                   trajpre,
     $                   trajtem,
     $                   irapp,
     $                   qj1,
     $                   hc,
     $                   vsed3d,
     $                   parth2osed,
     $                   parthno3sed,
     $                   surfarea,
     $                   sza3d)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this routine calculates the chemical rate coefficients,         c
c     the concentration of short-lived species (assumed to            c
c     be in photochemical equilibrium) and the production             c
c     and destruction rates of the long-lived species.                c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      parameter(nlon = 180, nlat = 91, niv = 137, nivbas = 94)
c
      parameter(nbcon = 43, ncm = 15)
      parameter(niter = 5)
      parameter(nphot = 41, nz = 101, nsza = 27, nozo = 7)
      parameter(nrappmax = 12)
c
      real qj1(nlon,nlat,niv,nbcon)
      real ck0(nlon,nivbas,nbcon)
      real hc(nlon,nlat,niv,ncm)
      real alt(nlon,nlat,niv)
      real daynum, dt
      real trajlon(nlon,nlat,niv,nrappmax)
      real trajlat(nlon,nlat,niv,nrappmax)
      real trajpre(nlon,nlat,niv,nrappmax)
      real trajtem(nlon,nlat,niv,nrappmax)
      real het1(nlon,nivbas), het2(nlon,nivbas), het3(nlon,nivbas)
      real het4(nlon,nivbas), het5(nlon,nivbas), het6(nlon,nivbas)
      real het7(nlon,nivbas), het8(nlon,nivbas), het9(nlon,nivbas)
      real parthno3(nlon,nivbas), parthno3sl(nlon,nivbas)
      real parth2o(nlon,nivbas), parthcl(nlon,nivbas)
      real parthbr(nlon,nivbas)
      real lntemp, xpo
c
      real parthno3sed(nlon,nlat,nivbas), parth2osed(nlon,nlat,nivbas)
      real vsed(nlon,nivbas), vsed3d(nlon,nlat,nivbas)
      real surfarea(nlon,nlat,niv)
c
      real pr(nlon,nivbas,nbcon), pe(nlon,nivbas,nbcon)
      real hnm(nlon,nivbas), n2(nlon,nivbas), o2(nlon,nivbas),
     $     co2(nlon,nivbas)
      real ck(nlon,nivbas,nbcon), t(nlon,nivbas), tj(nlon,nivbas,nphot)
      real sza(nlon,nivbas), sza3d(nlon,nlat,nivbas)
      real h2(nlon,nivbas), cc(nlon,nivbas,ncm)
      real coea, coeb, coec, argum
      real a1(nlon,nivbas), a1et(nlon,nivbas), a2(nlon,nivbas), 
     $     a5(nlon,nivbas),
     $     a6(nlon,nivbas), a6b(nlon,nivbas), a7(nlon,nivbas),
     $     a17(nlon,nivbas), a19(nlon,nivbas),
     $     a26(nlon,nivbas), a27(nlon,nivbas), a30(nlon,nivbas),
     $     a31(nlon,nivbas), a36(nlon,nivbas)
      real b3(nlon,nivbas), b4(nlon,nivbas), b6(nlon,nivbas),
     $     b7(nlon,nivbas),
     $     b9(nlon,nivbas), b12(nlon,nivbas), b22(nlon,nivbas),
     $     b23(nlon,nivbas), b24(nlon,nivbas), b27(nlon,nivbas),
     $     b28(nlon,nivbas), b32(nlon,nivbas), b38(nlon,nivbas),
     $     b39(nlon,nivbas)
      real c1a(nlon,nivbas),c1b(nlon,nivbas),c1c(nlon,nivbas),
     $     c2(nlon,nivbas), c3(nlon,nivbas), c4(nlon,nivbas),
     $     c5(nlon,nivbas), c7(nlon,nivbas), c8(nlon,nivbas),
     $     c9(nlon,nivbas), c10(nlon,nivbas), c12(nlon,nivbas),
     $     c14(nlon,nivbas), c15(nlon,nivbas), c17a(nlon,nivbas),
     $     c17b(nlon,nivbas), c18(nlon,nivbas), c19(nlon,nivbas)
      real d2(nlon,nivbas), d3(nlon,nivbas),
     $     d4(nlon,nivbas), d5(nlon,nivbas), d6(nlon,nivbas),
     $     d7(nlon,nivbas), d8(nlon,nivbas), d9(nlon,nivbas),
     $     d10(nlon,nivbas), d11(nlon,nivbas), d31(nlon,nivbas),
     $     d32(nlon,nivbas), d33(nlon,nivbas), d34(nlon,nivbas),
     $     d35(nlon,nivbas), d36(nlon,nivbas), d37(nlon,nivbas),
     $     d60(nlon,nivbas), d61(nlon,nivbas), d62(nlon,nivbas),
     $     d90(nlon,nivbas), d91(nlon,nivbas), d92(nlon,nivbas),
     $     d93(nlon,nivbas), d94(nlon,nivbas), d95(nlon,nivbas),
     $     d96(nlon,nivbas), d97(nlon,nivbas)
      real e2(nlon,nivbas), e3(nlon,nivbas), e4(nlon,nivbas),
     $     e5a(nlon,nivbas),e5b(nlon,nivbas),e5c(nlon,nivbas),
     $     e6(nlon,nivbas), e7(nlon,nivbas), e8(nlon,nivbas),
     $     e9(nlon,nivbas), e10(nlon,nivbas), e11(nlon,nivbas),
     $     e13(nlon,nivbas), e15(nlon,nivbas), e16(nlon,nivbas),
     $     e50(nlon,nivbas), e51(nlon,nivbas),e52(nlon,nivbas),
     $     e53(nlon,nivbas), e54(nlon,nivbas)
      real hk1(nlon,nivbas), hk2(nlon,nivbas), hk3(nlon,nivbas),
     $     hk4(nlon,nivbas), hk5(nlon,nivbas)
      real rain_h2o2(nlon,nivbas), rain_hcl(nlon,nivbas)
      real rain_hbr(nlon,nivbas), rain_hno3(nlon,nivbas)
      real rain_hno4(nlon,nivbas), rain_ch2o(nlon,nivbas)
      real rain_ch3o2h(nlon,nivbas), rain_h2o, rain_rate
      real o3t(nlon,nivbas), day(nlon,nivbas)
      real pm(nlon,nivbas)
      real ch3o(nlon,nivbas), hco(nlon,nivbas)
c
      integer ilat, irapp
c
      common /constant/
     $     a0hox, a1hox, a0nox1,
     $     a1nox1, a0nox2, a1nox2, a0nox3, a1nox3,
     $     a0clx1, a1clx1, a0clx2, a1clx2, a0clx3,
     $     a1clx3, a0brx, a1brx, a0c1, a1c1
      common /photodis/ ajl(nphot-1,nz,nsza,nozo), o3up(0:79)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c...  especes transportees:
c
c     ck( 1): n2o       ck( 2): ch4      ck( 3): h2o       ck( 4): noy
c     ck( 5): hno3      ck( 6): n2o5     ck( 7): cly       ck( 8): ox
c     ck( 9): co        ck(10): oclo     ck(11): ox passif ck(12): h2so4
c     ck(13): hcl       ck(14): clono2   ck(15): hocl      ck(16): cl2
c     ck(17): h2o2      ck(18): clno2    ck(19): hbr       ck(20): brono2
c     ck(21): nox       ck(22): ho2no2   ck(23): clox      ck(24): brox
c     ck(25): cl2o2     ck(26): hobr     ck(27): brcl      ck(28): ch2o
c     ck(29): ch3o2     ck(30): ch3o2h   ck(31): cfc-11    ck(32): cfc-12*
c     ck(33): cfc-113   ck(34): ccl4     ck(35): ch3ccl3*  ck(36): ch3cl
c     ck(37): hcfc-22*  ck(38): ch3br    ck(39): h-1211    ck(40): h-1301*
c     ck(41): bry       ck(42): ch2br2*  ck(43): hno3 gas
c
c     cfc-12*  inclut une contribution de cfc-114 et cfc-115
c     ch3ccl3* inclut une contribution de hcfc-141b
c     hcfc-22* inclut une contribution de hcfc-142b
c     h-1301*  inclut une contribution de h-2402
c     ch2br2*  inclut les contributions de toutes les especes bromees
c              a courte duree de vie
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c...  especes a l'equilibre:
c
c     cc( 1): o(1d)     cc( 2): oh        cc( 3): cl       cc( 4): o(3p)
c     cc( 5): o3        cc( 6): ho2       cc( 7): no2      cc( 8): no
c     cc( 9): br        cc(10): n         cc(11): clo      cc(12): bro
c     cc(13): no3       cc(14): h         cc(15): ch3
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c...  numerotation des coefficients de photodissociation:
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      j_o2         =  1
      j_no2        =  2
      j_o3_o1d     =  3
      j_o3_o3p     =  4
      j_no3_no     =  5
      j_no3_no2    =  6
      j_n2o5       =  7
      j_hno3       =  8
      j_hno4       =  9
      j_h2o2       = 10
      j_h2o        = 11
      j_n2o        = 12
      j_hocl       = 13
      j_clono2_cl  = 14
      j_clono2_clo = 15
      j_cl2o2      = 16
      j_hcl        = 17
      j_cl2        = 18
      j_clno2      = 19
      j_bro        = 20
      j_brono2     = 21
      j_hobr       = 22
      j_oclo       = 23
      j_brcl       = 24
      j_ch3o2h     = 25
      j_ch2o_hco   = 26
      j_ch2o_co    = 27
      j_cfc11      = 28
      j_cfc12      = 29
      j_cfc113     = 30
      j_ccl4       = 31
      j_ch3ccl3    = 32
      j_ch3cl      = 33
      j_hcfc22     = 34
      j_ch3br      = 35
      j_h1211      = 36
      j_h1301      = 37
      j_co2        = 38
      j_ch2br2     = 39
      j_hbr        = 40
      j_no         = 41
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     definition of total air density and concentration of h2, o2 and n2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do iniv = 1,nivbas
         do ilon = 1,nlon
            t(ilon,iniv)  = trajtem(ilon,ilat,iniv,irapp)
            pm(ilon,iniv) = trajpre(ilon,ilat,iniv,irapp)
            hnm(ilon,iniv) = trajpre(ilon,ilat,iniv,irapp)
     $                     /(trajtem(ilon,ilat,iniv,irapp)
     $                        *1.38e-19)
            n2(ilon,iniv)  = 0.78*hnm(ilon,iniv)
            o2(ilon,iniv)  = 0.21*hnm(ilon,iniv)
            co2(ilon,iniv) = 360.e-6*hnm(ilon,iniv)
            h2(ilon,iniv)  = 5.e-7*hnm(ilon,iniv)
         end do
      end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     definition of (initial) number densities of species
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do n = 1,nbcon
         do iniv = 1,nivbas
            do ilon = 1,nlon
               ck0(ilon,iniv,n) = qj1(ilon,ilat,iniv,n)
     $                            *hnm(ilon,iniv)
               ck(ilon,iniv,n) = ck0(ilon,iniv,n)
            end do
         end do
      end do
      do n = 1,ncm
         do iniv = 1,nivbas
            do ilon = 1,nlon
               cc(ilon,iniv,n) = hc(ilon,ilat,iniv,n)*hnm(ilon,iniv)
            end do
         end do
      end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calcul de l'ozone total                                        c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call total( trajpre, irapp, ilat, hc, hnm, o3t )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calcul de l'angle zenithal solaire                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call zenith( trajlon, trajlat, irapp, ilat, daynum, sza )
c
      do iniv = 1,nivbas
         do ilon = 1,nlon
            sza3d(ilon,ilat,iniv) = sza(ilon,iniv)
         end do
      end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calcul des coefficients de photodissociation                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call phot( pm, ilat, alt, sza, o3t, tj,
     $           cc(:,:,8), hnm )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialisation of heterogeneous reactions                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do iniv = 1,nivbas
         do ilon = 1,nlon
            het1(ilon,iniv) = 0.
            het2(ilon,iniv) = 0.
            het3(ilon,iniv) = 0.
            het4(ilon,iniv) = 0.
            het5(ilon,iniv) = 0.
            het6(ilon,iniv) = 0.
            het7(ilon,iniv) = 0.
            het8(ilon,iniv) = 0.
            het9(ilon,iniv) = 0.
            parth2o(ilon,iniv)  = 1.
            parthno3(ilon,iniv) = 1.
            parthcl(ilon,iniv)  = 1.
            parthbr(ilon,iniv)  = 1.
         end do
      end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     heterogeneous chemistry module                                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                    c
c  1)    clono2 + h2o -> hocl + hno3         het1                    c
c  2)    clono2 + hcl -> cl2  + hno3         het2                    c
c  3)    n2o5   + h2o ->       2hno3         het3                    c
c  4)    n2o5   + hcl -> clno2 +hno3         het4                    c
c  5)    hocl   + hcl -> cl2 + h2o           het5                    c
c  6)    brono2 + h2o -> hobr + hno3         het6                    c
c  7)    hobr   + hcl -> brcl + h2o          het7                    c
c  8)    hobr   + hbr -> br2 + h2o           het8                    c
c  9)    hocl   + hbr -> brcl + h2o          het9                    c
c                                                                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call analytic(t,pm,hnm,ck,vsed,
     $              het1,het2,het3,het4,het5,het6,het7,het8,het9,
     $              parthno3,parthno3sl,parth2o,parthcl,parthbr,
     $              surfarea,ilat)
c
c     stockage des vitesses de sedimentation 
c     et partitions non sedimentees
c
      do iniv = 1,nivbas
         do ilon = 1,nlon
            vsed3d(ilon,ilat,iniv)      = vsed(ilon,iniv)
            parth2osed(ilon,ilat,iniv)  = parth2o(ilon,iniv)
            parthno3sed(ilon,ilat,iniv) = parthno3sl(ilon,iniv)
         end do
      end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     1ere Partie: calcul des vitesses de reactions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      xln6 = log(0.6)
c
      do iniv = 1,nivbas
      do ilon = 1,nlon
c
         lntemp = log(t(ilon,iniv))
c
c...reactions des especes hydrogenees
c
c        a1: h + o2 + m -> ho2 + m
         ak0 = a0hox*exp(lntemp*(-1.3))                     !jpl 2006
         ak1 = a1hox*exp(lntemp*( 0.2))                     !jpl 2011
         rate = (ak0*hnm(ilon,iniv))
     $          /(1. + ak0*hnm(ilon,iniv)/ak1)
         xpo = 1./(1. + log10((ak0*hnm(ilon,iniv))/ak1)**2)
         a1(ilon,iniv) = rate*exp(xln6*xpo)
c
c        a1et : o(1d) + h2o   -> oh + oh
         a1et(ilon,iniv) = 1.63e-10*exp(60./t(ilon,iniv))   ! jpl 2006
c
c        a2 : h + o3      -> oh + o2
         a2(ilon,iniv) = 1.4e-10*exp(-470./t(ilon,iniv))
c
c        a3et : o(1d) + h2    -> oh + h
         a3et = 1.2e-10                                     ! jpl 2011
c
c        a5 : o + oh      -> o2 + h
         a5(ilon,iniv) = 1.8e-11*exp(180./t(ilon,iniv))     ! jpl 2011
c
c        a6 : oh + o3     -> ho2 + o2                       ! jpl 2003
         a6(ilon,iniv) = 1.7e-12*exp(-940./t(ilon,iniv))
c
c        a6b: ho2 + o3    -> oh + o2 + o2 (jpl 2003)
         a6b(ilon,iniv) = 1.0e-14*exp(-490./t(ilon,iniv))
c
c        a7 : o + ho2     -> oh + o2
         a7(ilon,iniv) = 3.0e-11*exp(200./t(ilon,iniv))
c
c        a17: oh + ho2    -> h2o + o2
         a17(ilon,iniv) = 4.8e-11*exp(250./t(ilon,iniv))
c
c        a19: oh + h2     -> h2o + h                        ! jpl 2006
         a19(ilon,iniv) = 2.8e-12*exp(-1800./t(ilon,iniv))
c
c        a23a : h + ho2     -> oh + oh
         a23a = 7.2e-11                 ! jpl 2006
c
c        a23b : h + ho2     -> h2 + o2
         a23b = 6.9e-12                 ! jpl 2006
c
c        a23c : h + ho2     -> h2o + o
         a23c = 1.6e-12                 ! jpl 2006
c
c        a26: no + ho2    -> no2 + oh
         a26(ilon,iniv) = 3.3e-12*exp(270./t(ilon,iniv))     ! jpl 2011
c
c        a27: ho2 + ho2   -> h2o2 + o2
         a27(ilon,iniv) = 1.5e-12*exp(19./t(ilon,iniv))  ! christensen et al., 2002
     $                  + 2.1e-33*exp(920./t(ilon,iniv))*hnm(ilon,iniv) ! jpl 2011
c
c        a30: oh + h2o2   -> h2o + ho2
         a30(ilon,iniv) = 1.8e-12          ! jpl 2006
c
c        a31: oh + oh -> h2o + o
         a31(ilon,iniv) = 1.8e-12          ! jpl 2006
c
c        a36: oh + co -> co2 + h
         a36(ilon,iniv) = 1.5e-13*(1. + .6*pm(ilon,iniv)/1013.)
c
c...reactions des especes azotees
c
c        b3 : o + no2    -> no + o2
         b3(ilon,iniv) = 5.1e-12*exp(210./t(ilon,iniv))  ! jpl 2006
c
c        b4 : o3 + no    -> no2 + o2
         b4(ilon,iniv) = 3.0e-12*exp((-1500./t(ilon,iniv)))
c
c        b6 : n + no    -> n2 + o
         b6(ilon,iniv) = 2.1e-11*exp((100./t(ilon,iniv)))
c
c        b7 : n + o2     -> no + o
         b7(ilon,iniv) = 1.5e-11*exp((-3600./t(ilon,iniv)))
c
c        b9 : o3 + no2   -> no3 + o2
         b9(ilon,iniv) = 1.2e-13*exp((-2450./t(ilon,iniv)))
c
c        b12: no2 + no3 + m -> n2o5 + m
         ak0 = a0nox1*exp(lntemp*(-4.4))
         ak1 = a1nox1*exp(lntemp*(-0.7))
         rate = (ak0*hnm(ilon,iniv))
     $          /(1. + ak0*hnm(ilon,iniv)/ak1)
         xpo = 1./(1. + log10((ak0*hnm(ilon,iniv))/ak1)**2)
         b12(ilon,iniv) = rate*exp(xln6*xpo)
c
c        b22: oh + no2 + m  -> hno3 + m
         ak0 = a0nox2*exp(lntemp*(-3.0))    ! jpl 2006
         ak1 = a1nox2*exp(lntemp*( 0.0))    ! jpl 2006
         rate = (ak0*hnm(ilon,iniv))
     $       /(1. + ak0*hnm(ilon,iniv)/ak1)
         xpo = 1./(1. + log10((ak0*hnm(ilon,iniv))/ak1)**2)
         b22(ilon,iniv) = rate*exp(xln6*xpo)
c
c        b23: ho2 + no2 + m -> ho2no2 + m
         ak0 = a0nox3*exp(lntemp*(-3.4))    ! jpl 2006
         ak1 = a1nox3*exp(lntemp*(-1.1))    ! jpl 2006
         rate = (ak0*hnm(ilon,iniv))
     $             /(1. + ak0*hnm(ilon,iniv)/ak1)
         xpo = 1./(1. + log10((ak0*hnm(ilon,iniv))/ak1)**2)
         b23(ilon,iniv) = rate*exp(xln6*xpo)
c
c        b24: ho2no2 + m -> ho2 + no2 + m
         deq = 2.1e-27*exp(10900./t(ilon,iniv))
         b24(ilon,iniv) = b23(ilon,iniv)/(deq*hnm(ilon,iniv))
c
c        b27: hno3 + oh  -> h2o + no3
         aux1 = 2.7e-17*exp(2199./t(ilon,iniv))
         aux2 = 1. 
     $        + aux1/(6.5e-34*exp(1335./t(ilon,iniv))*hnm(ilon,iniv))
         b27(ilon,iniv) = 2.4e-14*exp(460./t(ilon,iniv))+ aux1/aux2
c
c        b28: oh + ho2no2 -> h2o + no2 + o2
         b28(ilon,iniv) = 1.3e-12*exp(380./t(ilon,iniv))
c
c        b32: n2o5 + m   -> no2 + no3 + m
         deq = 2.7e-27*exp(11000./t(ilon,iniv))   ! jpl 2006
         b32(ilon,iniv) = b12(ilon,iniv)/(deq*hnm(ilon,iniv))
c
c        b38: o(1d) + n2o -> n2 + o2
         b38(ilon,iniv) = 4.63e-11*exp(20./t(ilon,iniv)) ! jpl 2009
c
c        b39: o(1d) + n2o -> no + no
         b39(ilon,iniv) = 7.25e-11*exp(20./t(ilon,iniv)) ! jpl 2009
c
c...reactions des especes carbonees
c
c        c1a : o(1d) + ch4-> oh + ch3
         c1a(ilon,iniv) = 1.31e-10  ! jpl 2009
c
c        c1b : o(1d) + ch4-> ch3o + h
         c1b(ilon,iniv) = 0.35e-10  ! jpl 2009
c
c        c1c : o(1d) + ch4-> ch2o + h2
         c1c(ilon,iniv) = 0.09e-10  ! jpl 2009
c
c        c2  : ch4 + oh -> ch3 + h2o
         c2(ilon,iniv) = 2.45e-12*exp((-1775./t(ilon,iniv)))
c
c        c3  : ch3 + o -> ch2o + h
         c3(ilon,iniv) = 1.1e-10
c
c        c4  : ch3 + o2 + m -> ch3o2 + m
         ak0 = a0c1*exp(lntemp*(-3.6))            ! jpl 2006
         ak1 = a1c1*exp(lntemp*( 1.1))            ! jpl 2006
         rate = (ak0*hnm(ilon,iniv))
     $       /(1. + ak0*hnm(ilon,iniv)/ak1)
         xpo = 1./(1. + log10((ak0*hnm(ilon,iniv))/ak1)**2)
         c4(ilon,iniv) = rate*exp(xln6*xpo)
c
c        c5  : ch3o2 + no -> ch3o + no2 (jpl 2003)
         c5(ilon,iniv) = 2.8e-12*exp((300./t(ilon,iniv)))
c
c        c7  : ch3o2 + ho2 -> ch3o2h + o2
         c7(ilon,iniv) = 3.8e-13*exp((800./t(ilon,iniv)))
c
c        c8  : ch2o + oh -> hco + h2o
         c8(ilon,iniv) = 5.5e-12*exp(125./t(ilon,iniv)) ! jpl 2006
c
c        c9  : ch2o + o -> hco + oh
         c9(ilon,iniv) = 3.4e-11*exp((-1600./t(ilon,iniv)))
c
c        c10 : ch2o + cl -> hco + hcl
         c10(ilon,iniv) = 8.1e-11*exp((-30./t(ilon,iniv)))
c
c        c12 : hco + o2 -> co + ho2 (jpl 2003)
         c12(ilon,iniv) = 5.2e-12
c
c        c14 : ch3o2 + ch3o2 -> 2 ch3o + o2 (jpl 2003)
         c14(ilon,iniv) = 9.5e-14*exp((390./t(ilon,iniv)))
c
c        c15 : ch3o + o2 -> ch2o + ho2
         c15(ilon,iniv) = 3.9e-14*exp((-900./t(ilon,iniv)))
c
c        c17a: ch3o2h + oh -> ch3o2 + h2o
         c17a(ilon,iniv) = 3.8e-12*exp((200./t(ilon,iniv)))*0.7
c
c        c17b: ch3o2h + oh -> ch2o + oh + h2o
         c17b(ilon,iniv) = 3.8e-12*exp((200./t(ilon,iniv)))*0.3
c
c        c18  : ch3o2 + clo -> ch3o + cl + o2
         c18(ilon,iniv) = 3.3e-12*exp((-115./t(ilon,iniv)))
c
c        c19 : ch2o + br -> hco + hbr
         c19(ilon,iniv) = 1.7e-11*exp((-800./t(ilon,iniv)))
c
c...reactions des especes chlorees
c
c        d2 : cl + o3      -> clo + o2
         d2(ilon,iniv) = 2.3e-11*exp((-200./t(ilon,iniv)))
c
c        d3 : clo + o      -> cl + o2
         d3(ilon,iniv) = 2.8e-11*exp(85./t(ilon,iniv))       ! jpl 2006
c
c        d4 : clo + no     -> no2 + cl
         d4(ilon,iniv) = 6.4e-12*exp(290./t(ilon,iniv))
c
c        d5 : cl + ch4     -> hcl + ch3
         d5(ilon,iniv) = 7.3e-12*exp((-1280./t(ilon,iniv)))  ! jpl 2006
c
c        d6 : cl + h2      -> hcl + h
         d6(ilon,iniv) = 3.05e-11*exp((-2270./t(ilon,iniv))) ! jpl 2006
c
c        d7 : cl + ho2     -> hcl + o2
         d7(ilon,iniv) = 1.4e-11*exp(269./t(ilon,iniv)) ! jpl 2009
c
c        d8 : clo + oh     -> cl + ho2
         d8(ilon,iniv) = 7.4e-12*exp(270./t(ilon,iniv))
c
c        d9 : clo + oh     -> hcl + o2 (jpl 2003)
         d9(ilon,iniv) = 6.0e-13*exp(230./t(ilon,iniv))
c
c        d10: cl + ho2     -> oh + clo
         d10(ilon,iniv) = 3.6e-11*exp(-375./t(ilon,iniv)) ! jpl 2009
c
c        d11: oh + hcl     -> h2o + cl
         d11(ilon,iniv) = 1.8e-12*exp((-250./t(ilon,iniv))) ! jpl 2009
c
c        d31: clo + no2 + m -> clono2 + m
         ak0 = a0clx1*exp(lntemp*(-3.4))
         ak1 = a1clx1*exp(lntemp*(-1.9))
         rate = (ak0*hnm(ilon,iniv))
     $          /(1. + ak0*hnm(ilon,iniv)/ak1)
         xpo = 1./(1. + log10((ak0*hnm(ilon,iniv))/ak1)**2)
         d31(ilon,iniv) = rate*exp(xln6*xpo)
c
c        d32: o + clono2   -> products
         d32(ilon,iniv) = 3.6e-12*exp((-840./t(ilon,iniv))) ! jpl 2011
c
c        d33: clo + ho2    -> hocl + o2
         d33(ilon,iniv) = 2.6e-12*exp(290./t(ilon,iniv)) ! jpl 2009
c
c        d34: oh + hocl    -> h2o + clo
         d34(ilon,iniv) = 3.0e-12*exp((-500./t(ilon,iniv)))
c
c        d35: o + hocl     -> oh + clo
         d35(ilon,iniv) = 1.7e-13
c
c        d36: cl + no2 + m -> clno2 + m
         ak0 = a0clx2*exp(lntemp*(-2))
         ak1 = a1clx2*exp(lntemp*(-1))
         rate = (ak0*hnm(ilon,iniv))
     $             /(1. + ak0*hnm(ilon,iniv)/ak1)
         xpo = 1./(1. + log10((ak0*hnm(ilon,iniv))/ak1)**2)
         d36(ilon,iniv) = rate*exp(xln6*xpo)
c
c        d37: cl + hocl    -> oh + cl2
         d37(ilon,iniv) = 3.4e-12*exp((-130./t(ilon,iniv))) ! jpl 2011
c
c        d60: clo + clo + m -> cl2o2 + m
         ak0 = a0clx3*exp(lntemp*(-4.5))
         ak1 = a1clx3*exp(lntemp*(-2.0))         ! jpl 2009
         rate = (ak0*hnm(ilon,iniv))
     $          /(1. + ak0*hnm(ilon,iniv)/ak1)
         xpo = 1./(1. + log10((ak0*hnm(ilon,iniv))/ak1)**2)
         d60(ilon,iniv) = rate*exp(xln6*xpo)
c
c        d61: cl2o2 + m     -> clo + clo + m
         deq = 1.72e-27*exp(8649./t(ilon,iniv)) ! jpl 2009
c
         d61(ilon,iniv) = d60(ilon,iniv)/(deq*hnm(ilon,iniv))
c
c        d62: clo + clo -> cl2 + o2
         d62(ilon,iniv) = 1.0e-12*exp(-1590./t(ilon,iniv))
     $                  + 3.0e-11*exp(-2450./t(ilon,iniv))
c
c        d90: o(1d) + cfc-11  -> products
         d90(ilon,iniv) = 2.3e-10
c
c        d91: o(1d) + cfc-12 -> products
         d91(ilon,iniv) = 1.4e-10
c
c        d92: o(1d) + cfc-113 -> products
         d92(ilon,iniv) = 2.0e-10
c
c        d93: o(1d) + ccl4 -> products
         d93(ilon,iniv) = 3.3e-10
c
c        d94: oh + ch3ccl3 -> ch2ccl3 + h2o
         d94(ilon,iniv) = 1.64e-12*exp((-1520./t(ilon,iniv)))  ! jpl 2006
c
c        d95: oh + ch3cl -> ch2cl + h2o (jpl 2003)
         d95(ilon,iniv) = 2.4e-12*exp((-1250./t(ilon,iniv)))
c
c        d96: o(1d) + hcfc-22 -> products
         d96(ilon,iniv) = 1.0e-10
c
c        d97: oh + hcfc-22 -> cf2cl + h2o
         d97(ilon,iniv) = 1.05e-12*exp((-1600./t(ilon,iniv)))
c
c...reactions des especes bromees
c
c        e2 : br + o3    -> bro + o2
         e2(ilon,iniv) = 1.6e-11*exp((-780./t(ilon,iniv))) ! jpl 2009
c
c        e3 : bro + o    -> br + o2
         e3(ilon,iniv) = 1.9e-11*exp(230./t(ilon,iniv))
c
c        e4 : bro + no   -> no2 + br
         e4(ilon,iniv) = 8.8e-12*exp(260./t(ilon,iniv))
c
c        e5a: bro + clo  -> oclo + br
         e5a(ilon,iniv) = 9.5e-13*exp(550./t(ilon,iniv))
c
c        e5b: bro + clo  -> br + cl + o2
         e5b(ilon,iniv) = 2.3e-12*exp(260./t(ilon,iniv))
c
c        e5c: bro + clo  -> brcl + o2
         e5c(ilon,iniv) = 4.1e-13*exp(290./t(ilon,iniv))
c
c        e6 : bro + bro  -> 2br + o2
         e6(ilon,iniv) = 1.5e-12*exp(230./t(ilon,iniv))
c
c        e7 : br + ho2   -> hbr + o2
         e7(ilon,iniv) = 4.8e-12*exp(-310./t(ilon,iniv))   ! jpl 2006
c
c        e8 : brono2 + o -> no3 + bro            
         e8(ilon,iniv) = 1.9e-11*exp(215./t(ilon,iniv))    ! jpl 2006
c
c        e9: o(1d) + hbr -> oh + br
         e9(ilon,iniv) =  1.5e-10
c
c        e10: bro + oh  -> br + ho2
         e10(ilon,iniv) = 1.7e-11*exp(250./t(ilon,iniv))   ! jpl 2006
c
c        e11: oh + hbr   -> h2o + br
         e11(ilon,iniv) = 5.5e-12*exp(200./t(ilon,iniv))   ! jpl 2006
c
c        e13: bro + no2 + m -> brono2 + m
         ak0 = a0brx*exp(lntemp*(-3.2))
         ak1 = a1brx*exp(lntemp*(-2.9))
         rate = (ak0*hnm(ilon,iniv))
     $          /(1. + ak0*hnm(ilon,iniv)/ak1)
         xpo = 1./(1. + log10((ak0*hnm(ilon,iniv))/ak1)**2)
         e13(ilon,iniv) = rate*exp(xln6*xpo)
c
c        e15: bro + ho2  -> hobr + o2
         e15(ilon,iniv) = 4.5e-12*exp(460./t(ilon,iniv))   ! jpl 2006
c
c        e16: hobr + o   -> oh + bro
         e16(ilon,iniv) = 1.2e-10*exp((-430./t(ilon,iniv)))
c
c        e50: o(1d) + ch3br -> products
         e50(ilon,iniv) = 1.8e-10
c
c        e51: oh + ch3br  -> ch2br + h2o (jpl 2003)
         e51(ilon,iniv) = 2.35e-12*exp((-1300./t(ilon,iniv)))
c
c        e52: o(1d) + h-1211 -> products
         e52(ilon,iniv) = 1.5e-10
c
c        e53: o(1d) + h-1301 -> products
         e53(ilon,iniv) = 1.0e-10
c       
c        e54: ch2br2 + oh -> products
         e54(ilon,iniv) = 2.0e-12*exp(-840./t(ilon,iniv))
c 
c
c...reactions de l'oxygene
c
c        hk1: o + o + m  -> o2 + m
         hk1(ilon,iniv) = 4.23e-28*hnm(ilon,iniv)
     $                           /(t(ilon,iniv)*t(ilon,iniv))
c
c        hk2: o + o2 + m -> o3 + m
         hk2(ilon,iniv) = 6.0e-34*(300./t(ilon,iniv))**2.4
     $                   *hnm(ilon,iniv)
c
c        hk3: o + o3     -> o2 + o2
         hk3(ilon,iniv) = 8.e-12*exp((-2060./t(ilon,iniv)))
c
c        hk4: o(1d) + n2   -> o + n2
         hk4(ilon,iniv) = 2.15e-11*exp(110./t(ilon,iniv)) ! jpl 2006
c
c        hk5: o(1d) + o2   -> o + o2
         hk5(ilon,iniv) = 3.3e-11*exp(55./t(ilon,iniv))   ! jpl 2006
c
      end do
      end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     lessivage des especes solubles dans la troposphere
c
c     actif a partir du seuil de rapport de melange de h2o: rain_h2o
c     taux de lessivage : rain_rate
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      rain_h2o  = 100.e-6
      rain_rate = 1.e-6    ! 10 jours
c
      do iniv = 1,nivbas
         do ilon = 1,nlon
            if (ck(ilon,iniv,3)/hnm(ilon,iniv) .ge. rain_h2o) then
               rain_hno3(ilon,iniv)   = rain_rate
               rain_hno4(ilon,iniv)   = 0.
               rain_hcl(ilon,iniv)    = rain_rate
               rain_hbr(ilon,iniv)    = rain_rate
               rain_h2o2(ilon,iniv)   = rain_rate
               rain_ch2o(ilon,iniv)   = 0.
               rain_ch3o2h(ilon,iniv) = 0.
            else
               rain_hno3(ilon,iniv)   = 0.
               rain_hno4(ilon,iniv)   = 0.
               rain_hcl(ilon,iniv)    = 0.
               rain_hbr(ilon,iniv)    = 0.
               rain_h2o2(ilon,iniv)   = 0.
               rain_ch2o(ilon,iniv)   = 0.
               rain_ch3o2h(ilon,iniv) = 0.
            end if
         end do
      end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     2eme Partie: chimie des composes a courte duree de vie
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccc   boucle d'iteration
c
      do 10050 iter = 1,niter
c
         do iniv = 1,nivbas
         do ilon = 1,nlon
            if( sza(ilon,iniv) .le. 95. ) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c              jour
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
               day(ilon,iniv) = 1.
c
c              ch3
c
               cc(ilon,iniv,15) =
     $               ((c2(ilon,iniv)*cc(ilon,iniv,2)
     $               + c1a(ilon,iniv)*cc(ilon,iniv,1)
     $               + d5(ilon,iniv)*cc(ilon,iniv,3))*ck(ilon,iniv,2))
     $               /(c4(ilon,iniv)*o2(ilon,iniv)
     $               + c3(ilon,iniv)*cc(ilon,iniv,4))
c
c              ch3o
c
               ch3o(ilon,iniv) =
     $                (c1b(ilon,iniv)*cc(ilon,iniv,1)*ck(ilon,iniv,2)
     $              + (c5(ilon,iniv)*cc(ilon,iniv,8)
     $              +  2.*c14(ilon,iniv)*ck(ilon,iniv,29) 
     $              +  c18(ilon,iniv)*cc(ilon,iniv,11))*ck(ilon,iniv,29)
     $              +  tj(ilon,iniv,j_ch3o2h)*ck(ilon,iniv,30))
     $              / (c15(ilon,iniv)*o2(ilon,iniv))
c
c              hco
c
               hco(ilon,iniv) = 
     $                ( c8(ilon,iniv)*cc(ilon,iniv,2)
     $                + c9(ilon,iniv)*cc(ilon,iniv,4)
     $                +c10(ilon,iniv)*cc(ilon,iniv,3)
     $                +c19(ilon,iniv)*cc(ilon,iniv,9)
     $                + tj(ilon,iniv,j_ch2o_hco))*ck(ilon,iniv,28)
     $                / (c12(ilon,iniv)*o2(ilon,iniv))
c
c              o/o3
c
               r1 = (tj(ilon,iniv,j_o3_o1d) + tj(ilon,iniv,j_o3_o3p)
     $             + a2(ilon,iniv)*cc(ilon,iniv,14)
     $             + a6(ilon,iniv)*cc(ilon,iniv,2)
     $             +a6b(ilon,iniv)*cc(ilon,iniv,6)
     $             + b4(ilon,iniv)*cc(ilon,iniv,8)
     $             + d2(ilon,iniv)*cc(ilon,iniv,3))
     $            /(hk2(ilon,iniv)*o2(ilon,iniv))
c
c              o3/ox
c
               r2 = 1./(1. + r1)
c
c              ozone
c
               cc(ilon,iniv,5) = r2*ck(ilon,iniv,8)
c
c              atomic oxygen o(3p)
c
               cc(ilon,iniv,4) = r1*cc(ilon,iniv,5)
c
c              o(1d)/o3
c
               r3 = tj(ilon,iniv,j_o3_o1d)
     $             /(hk4(ilon,iniv)*n2(ilon,iniv)
     $             + hk5(ilon,iniv)*o2(ilon,iniv))
c
c              excited atomic oxygen o(1d)
c
               cc(ilon,iniv,1) = r3*cc(ilon,iniv,5)
c
c              h/oh
c
               r4 = (a5(ilon,iniv)*cc(ilon,iniv,4)
     $             +a36(ilon,iniv)*ck(ilon,iniv,9)
     $             +a19(ilon,iniv)*h2(ilon,iniv))
     $             /(a1(ilon,iniv)*o2(ilon,iniv)
     $             + a2(ilon,iniv)*cc(ilon,iniv,5)
     $             + (a23a + a23b + a23c)*cc(ilon,iniv,6))
c
c              hydrogen atom
c
               cc(ilon,iniv,14) = r4*cc(ilon,iniv,2)
c
c              n/no
c
               r7 = tj(ilon,iniv,j_no)
     $             /(b7(ilon,iniv)*o2(ilon,iniv)
     $             + b6(ilon,iniv)*cc(ilon,iniv,8))
c
c              no3/no2
c
               r8 = (b9(ilon,iniv)*cc(ilon,iniv,5))
     $              /(tj(ilon,iniv,j_no3_no) 
     $               +tj(ilon,iniv,j_no3_no2)
     $               +b12(ilon,iniv)*cc(ilon,iniv,7))
c
c              no/no2
c
               r5 = (tj(ilon,iniv,j_no2)
     $             + b3(ilon,iniv)*cc(ilon,iniv,4)
     $             + tj(ilon,iniv,j_no3_no)*r8)
     $             /(b4(ilon,iniv)*cc(ilon,iniv,5)
     $             + a26(ilon,iniv)*cc(ilon,iniv,6)
     $             + d4(ilon,iniv)*cc(ilon,iniv,11)
     $             + e4(ilon,iniv)*cc(ilon,iniv,12)
     $             + c5(ilon,iniv)*ck(ilon,iniv,29))
c
c              no2
c
               cc(ilon,iniv,7) = ck(ilon,iniv,21)/(1.+ r5 + r7*r5 + r8)
c
c              no
c
               cc(ilon,iniv,8) = r5*cc(ilon,iniv,7)
c
c              n
c
               cc(ilon,iniv,10) = r7*cc(ilon,iniv,8)
c
c              no3
c
               cc(ilon,iniv,13) = r8*cc(ilon,iniv,7)
c
c              ho2/oh
c
               r6 = (a6(ilon,iniv)*cc(ilon,iniv,5)
     $             +a30(ilon,iniv)*ck(ilon,iniv,17)
     $             +c8(ilon,iniv)*ck(ilon,iniv,28)
     $             +d8(ilon,iniv)*cc(ilon,iniv,11)
     $             +e10(ilon,iniv)*cc(ilon,iniv,12)
     $             +a1(ilon,iniv)*o2(ilon,iniv)*r4)
     $            /((a23a + a23b + a23c)*cc(ilon,iniv,14)
     $             +a6b(ilon,iniv)*cc(ilon,iniv,5)
     $             +a7(ilon,iniv)*cc(ilon,iniv,4)
     $             +a26(ilon,iniv)*cc(ilon,iniv,8)
     $             +2.*a27(ilon,iniv)*cc(ilon,iniv,6)
     $             +a17(ilon,iniv)*cc(ilon,iniv,2)
     $             +d7(ilon,iniv)*cc(ilon,iniv,3)
     $             +d10(ilon,iniv)*cc(ilon,iniv,3)
     $             +d33(ilon,iniv)*cc(ilon,iniv,11)
     $             +e7(ilon,iniv)*cc(ilon,iniv,9)
     $             +e15(ilon,iniv)*cc(ilon,iniv,12))
c
c              ho2
c
               a = .5
               b = 1. - a
               cc(ilon,iniv,6) = a*cc(ilon,iniv,6) 
     $                         + b*cc(ilon,iniv,2)*r6
c
c              clono2/clo
c
               r13 = d31(ilon,iniv)*cc(ilon,iniv,7)
     $               /(tj(ilon,iniv,j_clono2_cl)
     $               + tj(ilon,iniv,j_clono2_clo)
     $               + d32(ilon,iniv)*cc(ilon,iniv,4)
     $               + het1(ilon,iniv)*ck(ilon,iniv,3)
     $               + het2(ilon,iniv)*ck(ilon,iniv,13))
c
c              hocl/clo
c
               r9 = (d33(ilon,iniv)*cc(ilon,iniv,6)
     $              + het1(ilon,iniv)*ck(ilon,iniv,3)*r13)
     $              /(het5(ilon,iniv)*ck(ilon,iniv,13)
     $               + d35(ilon,iniv)*cc(ilon,iniv,4)
     $               + tj(ilon,iniv,j_hocl))
c
c              cl2o2/clo
c
               r10 = (d60(ilon,iniv)*cc(ilon,iniv,11))
     $              /(d61(ilon,iniv)*hnm(ilon,iniv)
     $               + tj(ilon,iniv,j_cl2o2))
c
c              brcl/clo
c
               r12 =  e5c(ilon,iniv)*cc(ilon,iniv,12)
     $                /tj(ilon,iniv,j_brcl)
c
c              cl/clo
c
               r11 = (d3(ilon,iniv)*cc(ilon,iniv,4)
     $              + d4(ilon,iniv)*cc(ilon,iniv,8)
     $              + d8(ilon,iniv)*cc(ilon,iniv,2)
     $              + c18(ilon,iniv)*ck(ilon,iniv,29)
     $              + e5b(ilon,iniv)*cc(ilon,iniv,12)
     $              + tj(ilon,iniv,j_clono2_cl)*r13
     $              + tj(ilon,iniv,j_hocl)*r9
     $              + 2.*tj(ilon,iniv,j_cl2o2)*r10
     $              + tj(ilon,iniv,j_brcl)*r12)
     $              /(d2(ilon,iniv)*cc(ilon,iniv,5)
     $               +d5(ilon,iniv)*ck(ilon,iniv,2)
     $               +d7(ilon,iniv)*cc(ilon,iniv,6)
     $               +d10(ilon,iniv)*cc(ilon,iniv,6)
     $               +d37(ilon,iniv)*ck(ilon,iniv,15))
c
c              clo
c
               a = .5
               b = 1. - a
               cc(ilon,iniv,11) = a*cc(ilon,iniv,11)
     $                          + b*ck(ilon,iniv,23)/(1. + r11)
c
c              cl
c
               cc(ilon,iniv,3) = r11*cc(ilon,iniv,11)
c
c              brono2/bro
c
               r27 = e13(ilon,iniv)*cc(ilon,iniv,7)
     $               /(tj(ilon,iniv,j_brono2)
     $               + e8(ilon,iniv)*cc(ilon,iniv,4))
c
c              hobr/bro
c
               r28 = e15(ilon,iniv)*cc(ilon,iniv,6)
     $               /(tj(ilon,iniv,j_hobr) 
     $               + e16(ilon,iniv)*cc(ilon,iniv,4))
c
c              brcl/bro
c
               r29 =  e5c(ilon,iniv)*cc(ilon,iniv,11)
     $                /tj(ilon,iniv,j_brcl)
c
c              br/bro
c
               r30 = (e3(ilon,iniv)*cc(ilon,iniv,4)
     $              + e4(ilon,iniv)*cc(ilon,iniv,8)
     $              + e5a(ilon,iniv)*cc(ilon,iniv,11)
     $              + e5b(ilon,iniv)*cc(ilon,iniv,11)
     $              + 2.*e6(ilon,iniv)*cc(ilon,iniv,12)
     $              + e10(ilon,iniv)*cc(ilon,iniv,2)
     $              + tj(ilon,iniv,j_bro)
     $              + tj(ilon,iniv,j_brono2)*0.85*r27
     $              + tj(ilon,iniv,j_hobr)*r28
     $              + tj(ilon,iniv,j_brcl)*r29)
     $               /(e2(ilon,iniv)*cc(ilon,iniv,5)
     $                +e7(ilon,iniv)*cc(ilon,iniv,6)
     $                +c19(ilon,iniv)*ck(ilon,iniv,28))
c
c              bro
c
               a = .5
               b = 1. - a
               cc(ilon,iniv,12) = a*cc(ilon,iniv,12)
     $                          + b*ck(ilon,iniv,24)/(1. + r30)
c
c              br
c
               cc(ilon,iniv,9) = r30*cc(ilon,iniv,12)
c
c              hydroxyl radical oh
c
c              quadratic loss:
c
               coea =
     $         2.*a31(ilon,iniv)
     $         + 2.*(a17(ilon,iniv)
     $         + a27(ilon,iniv)*r6 + (a23b + a23c)*r4)*r6
c
c              linear loss:
c
               coeb =
     $           b27(ilon,iniv)*ck(ilon,iniv,5)*parthno3(ilon,iniv)
     $         + b22(ilon,iniv)*cc(ilon,iniv,7)
     $         + b28(ilon,iniv)*ck(ilon,iniv,22)
     $         + b23(ilon,iniv)*cc(ilon,iniv,7)*r6
     $         + c2(ilon,iniv)*ck(ilon,iniv,2)
     $         + c7(ilon,iniv)*ck(ilon,iniv,29)*r6
     $         + c8(ilon,iniv)*ck(ilon,iniv,28)
     $         + c17a(ilon,iniv)*ck(ilon,iniv,30)
     $         + d7(ilon,iniv)*cc(ilon,iniv,3)*r6
     $         + d9(ilon,iniv)*cc(ilon,iniv,11)
     $         + d11(ilon,iniv)*ck(ilon,iniv,13)*parthcl(ilon,iniv)
     $         + d33(ilon,iniv)*cc(ilon,iniv,11)*r6
     $         + d34(ilon,iniv)*ck(ilon,iniv,15)
     $         + e11(ilon,iniv)*ck(ilon,iniv,19)*parthbr(ilon,iniv)
     $         + e7(ilon,iniv)*cc(ilon,iniv,9)*r6
     $         + e15(ilon,iniv)*cc(ilon,iniv,12)*r6
c
c              production:
c
               t01 = 2.*(a1et(ilon,iniv)*cc(ilon,iniv,1) 
     $                  + tj(ilon,iniv,j_h2o))*ck(ilon,iniv,3)
     $                  *parth2o(ilon,iniv)
               t02 = ((c1a(ilon,iniv) + c1b(ilon,iniv))*ck(ilon,iniv,2)
     $               + 2.*a3et*h2(ilon,iniv))*cc(ilon,iniv,1)
               t03 = c3(ilon,iniv)*cc(ilon,iniv,15)*cc(ilon,iniv,4)
               t04 = tj(ilon,iniv,j_ch3o2h)*ck(ilon,iniv,30)
               t05 = c15(ilon,iniv)*ch3o(ilon,iniv)*o2(ilon,iniv)
               t06 = (tj(ilon,iniv,j_ch2o_hco)
     $               + c9(ilon,iniv)*cc(ilon,iniv,4))*ck(ilon,iniv,28)
               t07 = c12(ilon,iniv)*hco(ilon,iniv)*o2(ilon,iniv)
               t08 = tj(ilon,iniv,j_hno3)*ck(ilon,iniv,5)
     $               *parthno3(ilon,iniv)
               t09 = (b24(ilon,iniv)*hnm(ilon,iniv)
     $               + tj(ilon,iniv,j_hno4))*ck(ilon,iniv,22)
               t10 = (d35(ilon,iniv)*cc(ilon,iniv,4)
     $              + d37(ilon,iniv)*cc(ilon,iniv,3)
     $              + tj(ilon,iniv,j_hocl))*ck(ilon,iniv,15)
               t11 = 2.*tj(ilon,iniv,j_h2o2)*ck(ilon,iniv,17)
               t12 = d6(ilon,iniv)*h2(ilon,iniv)*cc(ilon,iniv,3)
               t13 = tj(ilon,iniv,j_hobr)*ck(ilon,iniv,26)
               t14 = e16(ilon,iniv)*cc(ilon,iniv,4)*ck(ilon,iniv,26)
               t15 = tj(ilon,iniv,j_hbr)*ck(ilon,iniv,19)
c
               coec = t01 + t02 + t03 + t04 + t05
     $              + t06 + t07 + t08 + t09 + t10 
     $              + t11 + t12 + t13 + t14 + t15
c
               argum = coeb*coeb + 4.*coea*coec
               argum = max(argum,0.)
c
               a = .5
               b = 1. - a
               cc(ilon,iniv,2) = a*cc(ilon,iniv,2)
     $                          +b*(sqrt(argum) - coeb)/(2.*coea)
c
            else
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c              nuit
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
               day(ilon,iniv) = 0.
c
c              h
c
               cc(ilon,iniv,14) = 0.
c
c              ch3
c
               cc(ilon,iniv,15) = 0.
c
c              ch3o
c
               ch3o(ilon,iniv) = 0.
c
c              hco
c
               hco(ilon,iniv)  = 0.
c
c              o(1d)
c
               cc(ilon,iniv,1) = 0.
c
c              o3p
c
               cc(ilon,iniv,4) = 0.
c
c              o3
c
               cc(ilon,iniv,5) = ck(ilon,iniv,8)
c
c              oh
c
               cc(ilon,iniv,2) = 0.
c
c              cl
c
               cc(ilon,iniv,3) = 0.
c
c              ho2
c
               cc(ilon,iniv,6) = 0.
c
c              no
c
               cc(ilon,iniv,8) = 0.
c
c              no3/no2
c
               r8 = b9(ilon,iniv)*cc(ilon,iniv,5)
     $             /(b12(ilon,iniv)*max(cc(ilon,iniv,7),1.e-30))
c
c              no2
c
               cc(ilon,iniv,7) = ck(ilon,iniv,21)/(1. + r8)
c
c              no3
c
               cc(ilon,iniv,13) = r8*cc(ilon,iniv,7)
c
c              n
c
               cc(ilon,iniv,10) = 0.
c
c              clo
c
               cc(ilon,iniv,11) = max(ck(ilon,iniv,23), 0.)
c
c              br
c
               cc(ilon,iniv,9) = 0.
c
c              bro
c
               cc(ilon,iniv,12) = max(ck(ilon,iniv,24), 0.)
c
            end if
         end do
         end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     3eme Partie: chimie des composes a longue duree de vie
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         do iniv = 1,nivbas
         do ilon = 1,nlon
c
c           n2o
c
            pr(ilon,iniv,1) = 0.e0
c
            pe(ilon,iniv,1) = tj(ilon,iniv,j_n2o)
     $                      + (b38(ilon,iniv) 
     $                      +  b39(ilon,iniv))*cc(ilon,iniv,1)
c
c           ch4
c
            pr(ilon,iniv,2) = 0.e0
c
            pe(ilon,iniv,2) =
     $      c2(ilon,iniv)*cc(ilon,iniv,2)
     $    + c1a(ilon,iniv)*cc(ilon,iniv,1)
     $    + c1b(ilon,iniv)*cc(ilon,iniv,1)
     $    + c1c(ilon,iniv)*cc(ilon,iniv,1)
     $    + d5(ilon,iniv)*cc(ilon,iniv,3)
c
c           h2o
c
            pr(ilon,iniv,3) =
     $      c2(ilon,iniv)*cc(ilon,iniv,2)*ck(ilon,iniv,2)
     $    + c8(ilon,iniv)*cc(ilon,iniv,2)*ck(ilon,iniv,28)
     $    + c17a(ilon,iniv)*ck(ilon,iniv,30)*cc(ilon,iniv,2)
     $    + c17b(ilon,iniv)*ck(ilon,iniv,30)*cc(ilon,iniv,2)
     $    + (a17(ilon,iniv)*cc(ilon,iniv,2)
     $    +        a23c*cc(ilon,iniv,14))*cc(ilon,iniv,6)
     $    + a19(ilon,iniv)*cc(ilon,iniv,2)*h2(ilon,iniv)
     $    + a30(ilon,iniv)*cc(ilon,iniv,2)*ck(ilon,iniv,17)
     $    + a31(ilon,iniv)*cc(ilon,iniv,2)*cc(ilon,iniv,2)
     $    + b27(ilon,iniv)*cc(ilon,iniv,2)*ck(ilon,iniv,5)
     $      *parthno3(ilon,iniv)
     $    + d11(ilon,iniv)*cc(ilon,iniv,2)*ck(ilon,iniv,13)
     $      *parthcl(ilon,iniv)
c
            pe(ilon,iniv,3) =
     $      (tj(ilon,iniv,j_h2o) + a1et(ilon,iniv)*cc(ilon,iniv,1))
     $      *parth2o(ilon,iniv)
c
c           noy
c
            pr(ilon,iniv,4) = 
     $      2.*b39(ilon,iniv)*cc(ilon,iniv,1)*ck(ilon,iniv,1)
c
            pe(ilon,iniv,4) = 
     $      2.*b6(ilon,iniv)*cc(ilon,iniv,10)*cc(ilon,iniv,8)
     $    + rain_hno3(ilon,iniv)*ck(ilon,iniv,5)
     $    + rain_hno4(ilon,iniv)*ck(ilon,iniv,22)
            pe(ilon,iniv,4) = pe(ilon,iniv,4)/ck(ilon,iniv,4)
c
c           hno3
c
            pr(ilon,iniv,5) =
     $      b22(ilon,iniv)*cc(ilon,iniv,2)*cc(ilon,iniv,7)
     $    + het1(ilon,iniv)*ck(ilon,iniv,14)*ck(ilon,iniv,3)
     $    + het2(ilon,iniv)*ck(ilon,iniv,14)*ck(ilon,iniv,13)
     $    + 2.*het3(ilon,iniv)*ck(ilon,iniv,6)*ck(ilon,iniv,3)
     $    + het4(ilon,iniv)*ck(ilon,iniv,6)*ck(ilon,iniv,13)
     $    + het6(ilon,iniv)*ck(ilon,iniv,20)*ck(ilon,iniv,3)
c
            pe(ilon,iniv,5) =
     $      (tj(ilon,iniv,j_hno3) + b27(ilon,iniv)*cc(ilon,iniv,2))
     $      *parthno3(ilon,iniv)
     $    + rain_hno3(ilon,iniv)
c
c           n2o5
c
            pr(ilon,iniv,6) =
     $      b12(ilon,iniv)*cc(ilon,iniv,13)*cc(ilon,iniv,7)
c
            pe(ilon,iniv,6) =
     $      tj(ilon,iniv,j_n2o5)
     $    + b32(ilon,iniv)*hnm(ilon,iniv)
     $    + het3(ilon,iniv)*ck(ilon,iniv,3)
     $    + het4(ilon,iniv)*ck(ilon,iniv,13)
c
c           cly
c
            pr(ilon,iniv,7) =
     $      3.*(tj(ilon,iniv,j_cfc11)
     $         +d90(ilon,iniv)*cc(ilon,iniv,1))*ck(ilon,iniv,31)
     $    + 2.*(tj(ilon,iniv,j_cfc12)
     $         +d91(ilon,iniv)*cc(ilon,iniv,1))*ck(ilon,iniv,32)
     $    + 3.*(tj(ilon,iniv,j_cfc113)
     $         +d92(ilon,iniv)*cc(ilon,iniv,1))*ck(ilon,iniv,33)
     $    + 4.*(tj(ilon,iniv,j_ccl4)
     $         +d93(ilon,iniv)*cc(ilon,iniv,1))*ck(ilon,iniv,34)
     $    + 3.*(tj(ilon,iniv,j_ch3ccl3)
     $         +d94(ilon,iniv)*cc(ilon,iniv,2))*ck(ilon,iniv,35)
     $    + 1.*(tj(ilon,iniv,j_ch3cl)
     $         +d95(ilon,iniv)*cc(ilon,iniv,2))*ck(ilon,iniv,36)
     $    + 1.*(tj(ilon,iniv,j_hcfc22)
     $         +d96(ilon,iniv)*cc(ilon,iniv,1)
     $         +d97(ilon,iniv)*cc(ilon,iniv,2))*ck(ilon,iniv,37)
     $    + 1.*(tj(ilon,iniv,j_h1211)
     $         +e52(ilon,iniv)*cc(ilon,iniv,1))*ck(ilon,iniv,39)
c
            pe(ilon,iniv,7) = rain_hcl(ilon,iniv)*ck(ilon,iniv,13)
            pe(ilon,iniv,7) = pe(ilon,iniv,7)/ck(ilon,iniv,7)
c
c           ox
c
            ratio = tj(ilon,iniv,j_cl2o2)
     $            /(tj(ilon,iniv,j_cl2o2)
     $              + d61(ilon,iniv)*hnm(ilon,iniv))
c
c           production:
c
            pr(ilon,iniv,8) =
     $      2.*tj(ilon,iniv,j_o2)*o2(ilon,iniv)
     $    + a26(ilon,iniv)*cc(ilon,iniv,6)*cc(ilon,iniv,8)
     $    + c5(ilon,iniv)*cc(ilon,iniv,8)*ck(ilon,iniv,29)
c
c           perte:
c
c           recombinaison o-o:
c
            t01 = 2.*hk1(ilon,iniv)*cc(ilon,iniv,4)*cc(ilon,iniv,4)
c
c           recombinaison o-o3:
c
            t02 = 2.*hk3(ilon,iniv)*cc(ilon,iniv,5)*cc(ilon,iniv,4)
c
c           cycle no-no2:
c
            t03 = 2.*b3(ilon,iniv)*cc(ilon,iniv,7)*cc(ilon,iniv,4)
c
c           cycle clo-o:
c
            t04 = 2.*d3(ilon,iniv)*cc(ilon,iniv,11)*cc(ilon,iniv,4)
c
c           cycle clo-clo:
c
            t05 = 2.*d60(ilon,iniv)*ratio
     $            *cc(ilon,iniv,11)*cc(ilon,iniv,11)
c
c           cycle clo-ho2:
c
            t06 = min(tj(ilon,iniv,j_hocl)*ck(ilon,iniv,15),
     $                d33(ilon,iniv)*cc(ilon,iniv,6)*cc(ilon,iniv,11))
c
c           cycle bro-clo (via brcl):
c
            t07 = 2.*min(e5c(ilon,iniv)*cc(ilon,iniv,12)
     $                  *cc(ilon,iniv,11),
     $                   tj(ilon,iniv,j_brcl)*ck(ilon,iniv,27))
c
c           cycle bro-clo (via clo2):
c
            t08 = 2.*e5b(ilon,iniv)
     $            *cc(ilon,iniv,12)*cc(ilon,iniv,11)*day(ilon,iniv)
c
c           cycle bro-ho2:
c
            t09 = min(e15(ilon,iniv)*cc(ilon,iniv,12)*cc(ilon,iniv,6),
     $                tj(ilon,iniv,j_hobr)*ck(ilon,iniv,26))
c
c           cycle bro-bro:
c
            t10 = 2.*e6(ilon,iniv)
     $            *cc(ilon,iniv,12)*cc(ilon,iniv,12)*day(ilon,iniv)
c
c           cycle clono2-hv:
c
            t11 = 2.*(tj(ilon,iniv,j_no3_no)
     $            /max(tj(ilon,iniv,j_no3_no)
     $               + tj(ilon,iniv,j_no3_no2), 1.e-30))
     $             *tj(ilon,iniv,j_clono2_cl)*ck(ilon,iniv,14)
c
c           destruction par hox:
c
            t12 = (a2(ilon,iniv)*cc(ilon,iniv,14)
     $           + a6(ilon,iniv)*cc(ilon,iniv,2)
     $           + a6b(ilon,iniv)*cc(ilon,iniv,6))*cc(ilon,iniv,5)
     $           + (a5(ilon,iniv)*cc(ilon,iniv,2)
     $           + a7(ilon,iniv)*cc(ilon,iniv,6))*cc(ilon,iniv,4)
c
c           cycle no2-no3:
c
            t13 = 2.*b9(ilon,iniv)*cc(ilon,iniv,7)*cc(ilon,iniv,5)
     $           *tj(ilon,iniv,j_no3_no)/max(tj(ilon,iniv,j_no3_no)
     $           + tj(ilon,iniv,j_no3_no2), 1.e-30)
c
            pe(ilon,iniv,8) = t01 + t02 + t03 + t04
     $                      + t05 + t06 + t07 + t08
     $                      + t09 + t10 + t11 + t12  
     $                      + t13
c
            pe(ilon,iniv,8) = pe(ilon,iniv,8)/ck(ilon,iniv,8)
c
c           co
c
            pr(ilon,iniv,9) =
     $        tj(ilon,iniv,j_co2)*co2(ilon,iniv)
     $      + tj(ilon,iniv,j_ch2o_co)*ck(ilon,iniv,28)
     $      + c12(ilon,iniv)*hco(ilon,iniv)*o2(ilon,iniv)
c
            pe(ilon,iniv,9) = a36(ilon,iniv)*cc(ilon,iniv,2)
c
c           oclo
c
            pr(ilon,iniv,10) = 
     $      e5a(ilon,iniv)*cc(ilon,iniv,12)*cc(ilon,iniv,11)
c
            pe(ilon,iniv,10) = tj(ilon,iniv,j_oclo)
c
c           ox passif
c
            pr(ilon,iniv,11) = 0.
            pe(ilon,iniv,11) = 0.
c
c           h2so4
c
            pr(ilon,iniv,12) = 0.
            pe(ilon,iniv,12) = 0.
c
c           hcl
c
            pr(ilon,iniv,13) =
     $      (c10(ilon,iniv)*ck(ilon,iniv,28)
     $     + d5(ilon,iniv)*ck(ilon,iniv,2)
     $     + d7(ilon,iniv)*cc(ilon,iniv,6)
     $     + d6(ilon,iniv)*h2(ilon,iniv))*cc(ilon,iniv,3)
     $     + d9(ilon,iniv)*cc(ilon,iniv,2)*cc(ilon,iniv,11)
c
            pe(ilon,iniv,13) =
     $      (tj(ilon,iniv,j_hcl) + d11(ilon,iniv)*cc(ilon,iniv,2))
     $    * parthcl(ilon,iniv)
     $    + het2(ilon,iniv)*ck(ilon,iniv,14)
     $    + het4(ilon,iniv)*ck(ilon,iniv,6)
     $    + het5(ilon,iniv)*ck(ilon,iniv,15)
     $    + het7(ilon,iniv)*ck(ilon,iniv,26)
     $    + rain_hcl(ilon,iniv)
c
c           clono2
c
            pr(ilon,iniv,14) =
     $      d31(ilon,iniv)*cc(ilon,iniv,7)*cc(ilon,iniv,11)
c
            pe(ilon,iniv,14) =
     $      tj(ilon,iniv,j_clono2_cl)
     $    + tj(ilon,iniv,j_clono2_clo)
     $    + d32(ilon,iniv)*cc(ilon,iniv,4)
     $    + het1(ilon,iniv)*ck(ilon,iniv,3)
     $    + het2(ilon,iniv)*ck(ilon,iniv,13)
c
c           hocl
c
            pr(ilon,iniv,15) =
     $      d33(ilon,iniv)*cc(ilon,iniv,6)*cc(ilon,iniv,11)
     $    + het1(ilon,iniv)*ck(ilon,iniv,14)*ck(ilon,iniv,3)
c
            pe(ilon,iniv,15) =
     $      tj(ilon,iniv,j_hocl)
     $    + d34(ilon,iniv)*cc(ilon,iniv,2)
     $    + d35(ilon,iniv)*cc(ilon,iniv,4)
     $    + d37(ilon,iniv)*cc(ilon,iniv,3)
     $    + het5(ilon,iniv)*ck(ilon,iniv,13)
     $    + het9(ilon,iniv)*ck(ilon,iniv,19)
c
c           cl2
c
            pr(ilon,iniv,16) =
     $      d37(ilon,iniv)*ck(ilon,iniv,15)*cc(ilon,iniv,3)
     $    + d62(ilon,iniv)*cc(ilon,iniv,11)*cc(ilon,iniv,11)
     $    + het2(ilon,iniv)*ck(ilon,iniv,14)*ck(ilon,iniv,13)
     $    + het5(ilon,iniv)*ck(ilon,iniv,15)*ck(ilon,iniv,13)
c
            pe(ilon,iniv,16) = tj(ilon,iniv,j_cl2)
c
c           h2o2
c
            pr(ilon,iniv,17) =
     $      a27(ilon,iniv)*cc(ilon,iniv,6)*cc(ilon,iniv,6)
c
            pe(ilon,iniv,17) = tj(ilon,iniv,j_h2o2)
     $                    + a30(ilon,iniv)*cc(ilon,iniv,2)
     $                    + rain_h2o2(ilon,iniv)
c
c           clno2
c
            pr(ilon,iniv,18) =
     $      d36(ilon,iniv)*cc(ilon,iniv,3)*cc(ilon,iniv,7)
     $    + het4(ilon,iniv)*ck(ilon,iniv,6)*ck(ilon,iniv,13)
c
            pe(ilon,iniv,18) = tj(ilon,iniv,j_clno2)
c
c           hbr
c
            pr(ilon,iniv,19) =
     $      e7(ilon,iniv)*cc(ilon,iniv,6)*cc(ilon,iniv,9)
     $    + c19(ilon,iniv)*ck(ilon,iniv,28)*cc(ilon,iniv,9)
c
            pe(ilon,iniv,19) = (e11(ilon,iniv)*cc(ilon,iniv,2)
     $                     +e9(ilon,iniv)*cc(ilon,iniv,1)
     $                     +tj(ilon,iniv,j_hbr))
     $                     *parthbr(ilon,iniv)
     $                     +het8(ilon,iniv)*ck(ilon,iniv,26)
     $                     +het9(ilon,iniv)*ck(ilon,iniv,15)
     $                     +rain_hbr(ilon,iniv)
c
c           brono2
c
            pr(ilon,iniv,20) =
     $      e13(ilon,iniv)*cc(ilon,iniv,12)*cc(ilon,iniv,7)
c
            pe(ilon,iniv,20) = 
     $      tj(ilon,iniv,j_brono2)
     $    + e8(ilon,iniv)*cc(ilon,iniv,4)
     $    + het6(ilon,iniv)*ck(ilon,iniv,3)
c
c           nox
c
            pr(ilon,iniv,21) =
     $     (tj(ilon,iniv,j_hno4)
     $    + b28(ilon,iniv)*cc(ilon,iniv,2)
     $    + b24(ilon,iniv)*hnm(ilon,iniv))*ck(ilon,iniv,22)
     $    + tj(ilon,iniv,j_hno3)*ck(ilon,iniv,5)*parthno3(ilon,iniv)
     $    + 2.*(tj(ilon,iniv,j_n2o5)
     $    + b32(ilon,iniv)*hnm(ilon,iniv))*ck(ilon,iniv,6)
     $    + 2.*b39(ilon,iniv)*cc(ilon,iniv,1)*ck(ilon,iniv,1)
     $    + tj(ilon,iniv,j_clno2)*ck(ilon,iniv,18)
     $    + b27(ilon,iniv)*ck(ilon,iniv,5)
     $      *parthno3(ilon,iniv)*cc(ilon,iniv,2)
     $    + d32(ilon,iniv)*ck(ilon,iniv,14)*cc(ilon,iniv,4)
     $    + (tj(ilon,iniv,j_clono2_cl)
     $    + tj(ilon,iniv,j_clono2_clo))*ck(ilon,iniv,14)
     $    + tj(ilon,iniv,j_brono2)*ck(ilon,iniv,20)
     $    + e8(ilon,iniv)*ck(ilon,iniv,20)*cc(ilon,iniv,4)
c
            pe(ilon,iniv,21) =
     $     (d31(ilon,iniv)*cc(ilon,iniv,11)
     $    + b22(ilon,iniv)*cc(ilon,iniv,2)
     $    + b23(ilon,iniv)*cc(ilon,iniv,6)
     $    + 2.*b12(ilon,iniv)*cc(ilon,iniv,13)
     $    + d36(ilon,iniv)*cc(ilon,iniv,3)
     $    + e13(ilon,iniv)*cc(ilon,iniv,12))*cc(ilon,iniv,7)
     $    + 2.*b6(ilon,iniv)*cc(ilon,iniv,10)*cc(ilon,iniv,8)
            pe(ilon,iniv,21) = pe(ilon,iniv,21)/ck(ilon,iniv,21)
c
c           ho2no2
c
            pr(ilon,iniv,22) =
     $      b23(ilon,iniv)*cc(ilon,iniv,6)*cc(ilon,iniv,7)
c
            pe(ilon,iniv,22) =
     $      tj(ilon,iniv,j_hno4)
     $    + b28(ilon,iniv)*cc(ilon,iniv,2)
     $    + b24(ilon,iniv)*hnm(ilon,iniv)
     $    + rain_hno4(ilon,iniv)
c
c           clox
c
            t01 = 2.* (d61(ilon,iniv)*hnm(ilon,iniv)
     $          + tj(ilon,iniv,j_cl2o2))*ck(ilon,iniv,25)
            t02 = (tj(ilon,iniv,j_clono2_cl) 
     $          + tj(ilon,iniv,j_clono2_clo)
     $          + d32(ilon,iniv)*cc(ilon,iniv,4))*ck(ilon,iniv,14)
            t03 = (tj(ilon,iniv,j_hocl) + d34(ilon,iniv)*cc(ilon,iniv,2)
     $          + d35(ilon,iniv)*cc(ilon,iniv,4))*ck(ilon,iniv,15)
            t04 = tj(ilon,iniv,j_oclo)*ck(ilon,iniv,10)
     $          + (tj(ilon,iniv,j_hcl) + d11(ilon,iniv)*cc(ilon,iniv,2))
     $          *ck(ilon,iniv,13)*parthcl(ilon,iniv)
            t05 = 2.*tj(ilon,iniv,j_cl2)*ck(ilon,iniv,16)
            t06 = tj(ilon,iniv,j_clno2)*ck(ilon,iniv,18)
            t07 = tj(ilon,iniv,j_brcl)*ck(ilon,iniv,27)
c
c         source gases
c
            t08 = 3.*(tj(ilon,iniv,j_cfc11)
     $          + d90(ilon,iniv)*cc(ilon,iniv,1))*ck(ilon,iniv,31)
            t09 = 2.*(tj(ilon,iniv,j_cfc12)
     $          + d91(ilon,iniv)*cc(ilon,iniv,1))*ck(ilon,iniv,32)
            t10 = 3.*(tj(ilon,iniv,j_cfc113)
     $          + d92(ilon,iniv)*cc(ilon,iniv,1))*ck(ilon,iniv,33)
            t11 = 4.*(tj(ilon,iniv,j_ccl4)
     $          + d93(ilon,iniv)*cc(ilon,iniv,1))*ck(ilon,iniv,34)
            t12 = 3.*(tj(ilon,iniv,j_ch3ccl3)
     $          + d94(ilon,iniv)*cc(ilon,iniv,2))*ck(ilon,iniv,35)
            t13 = 1.*(tj(ilon,iniv,j_ch3cl)
     $          + d95(ilon,iniv)*cc(ilon,iniv,2))*ck(ilon,iniv,36)
            t14 = 1.*(tj(ilon,iniv,j_hcfc22)
     $          + d96(ilon,iniv)*cc(ilon,iniv,1)
     $          + d97(ilon,iniv)*cc(ilon,iniv,2))*ck(ilon,iniv,37)
c
            pr(ilon,iniv,23) = t01 + t02 + t03 + t04
     $                       + t05 + t06 + t07 + t08
     $                       + t09 + t10 + t11 + t12
     $                       + t13 + t14
c
            pe(ilon,iniv,23) =
     $      (d31(ilon,iniv)*cc(ilon,iniv,7)
     $    + 2.*d60(ilon,iniv)*cc(ilon,iniv,11)
     $    + 2.*d62(ilon,iniv)*cc(ilon,iniv,11)
     $    + e5a(ilon,iniv)*cc(ilon,iniv,12)
     $    + e5c(ilon,iniv)*cc(ilon,iniv,12)
     $    + d9(ilon,iniv)*cc(ilon,iniv,2)
     $    + d33(ilon,iniv)*cc(ilon,iniv,6))*cc(ilon,iniv,11)
     $    + (c10(ilon,iniv)*ck(ilon,iniv,28)
     $    + d5(ilon,iniv)*ck(ilon,iniv,2)
     $    + d7(ilon,iniv)*cc(ilon,iniv,6)
     $    + d6(ilon,iniv)*h2(ilon,iniv)
     $    + d36(ilon,iniv)*cc(ilon,iniv,7)
     $    + d37(ilon,iniv)*ck(ilon,iniv,15))*cc(ilon,iniv,3)
            pe(ilon,iniv,23) = pe(ilon,iniv,23)/ck(ilon,iniv,23)
c
c           brox
c
            pr(ilon,iniv,24) =
     $      tj(ilon,iniv,j_hobr)*ck(ilon,iniv,26)
     $    + tj(ilon,iniv,j_brono2)*ck(ilon,iniv,20)
     $    + e8(ilon,iniv)*ck(ilon,iniv,20)*cc(ilon,iniv,4)
     $    + e16(ilon,iniv)*ck(ilon,iniv,26)*cc(ilon,iniv,4)
     $    + tj(ilon,iniv,j_brcl)*ck(ilon,iniv,27)
     $    + (e11(ilon,iniv)*cc(ilon,iniv,2) 
     $    + e9(ilon,iniv)*cc(ilon,iniv,1)
     $    + tj(ilon,iniv,j_hbr))*ck(ilon,iniv,19)*parthbr(ilon,iniv)
     $    + 2.*het8(ilon,iniv)*ck(ilon,iniv,19)*ck(ilon,iniv,26)
c
            pe(ilon,iniv,24) =
     $     (e13(ilon,iniv)*cc(ilon,iniv,7)
     $    + e5c(ilon,iniv)*cc(ilon,iniv,11)
     $    + e15(ilon,iniv)*cc(ilon,iniv,6))*cc(ilon,iniv,12)
     $    + c19(ilon,iniv)*cc(ilon,iniv,9)*ck(ilon,iniv,28)
     $    + e7(ilon,iniv)*cc(ilon,iniv,9)*cc(ilon,iniv,6)
            pe(ilon,iniv,24) = pe(ilon,iniv,24)/ck(ilon,iniv,24)
c
c           cl2o2
c
            pr(ilon,iniv,25) =
     $      d60(ilon,iniv)*cc(ilon,iniv,11)*cc(ilon,iniv,11)
c
            pe(ilon,iniv,25) = tj(ilon,iniv,j_cl2o2)
     $                       + d61(ilon,iniv)*hnm(ilon,iniv)
c
c           hobr
c
            pr(ilon,iniv,26) =
     $      e15(ilon,iniv)*cc(ilon,iniv,12)*cc(ilon,iniv,6)
     $    + het6(ilon,iniv)*ck(ilon,iniv,20)*ck(ilon,iniv,3)
c
            pe(ilon,iniv,26) = 
     $      tj(ilon,iniv,j_hobr)
     $    + e16(ilon,iniv)*cc(ilon,iniv,4)
     $    + het7(ilon,iniv)*ck(ilon,iniv,13)
     $    + het8(ilon,iniv)*ck(ilon,iniv,19)
c
c           brcl
c
            pr(ilon,iniv,27) =
     $      e5c(ilon,iniv)*cc(ilon,iniv,11)*cc(ilon,iniv,12)
     $    + het7(ilon,iniv)*ck(ilon,iniv,26)*ck(ilon,iniv,13)
     $    + het9(ilon,iniv)*ck(ilon,iniv,19)*ck(ilon,iniv,15)
c
            pe(ilon,iniv,27) = tj(ilon,iniv,j_brcl)
c
c           ch2o
c
            pr(ilon,iniv,28) =
     $      (c1c(ilon,iniv)*cc(ilon,iniv,1)*ck(ilon,iniv,2)
     $    + c3(ilon,iniv)*cc(ilon,iniv,4)*cc(ilon,iniv,15)
     $    + c15(ilon,iniv)*ch3o(ilon,iniv)*o2(ilon,iniv)
     $    + c17b(ilon,iniv)*cc(ilon,iniv,2)*ck(ilon,iniv,30))
c
            pe(ilon,iniv,28) =
     $      c8(ilon,iniv)*cc(ilon,iniv,2)
     $    + c9(ilon,iniv)*cc(ilon,iniv,4)
     $    + c10(ilon,iniv)*cc(ilon,iniv,3)
     $    + c19(ilon,iniv)*cc(ilon,iniv,9)
     $    + tj(ilon,iniv,j_ch2o_hco)
     $    + tj(ilon,iniv,j_ch2o_co)
     $    + rain_ch2o(ilon,iniv)
c
c           ch3o2
c
            pr(ilon,iniv,29) =
     $      c4(ilon,iniv)*cc(ilon,iniv,15)*o2(ilon,iniv)
     $    + c17a(ilon,iniv)*cc(ilon,iniv,2)*ck(ilon,iniv,30)
c
            pe(ilon,iniv,29) = 
     $      2.*c14(ilon,iniv)*ck(ilon,iniv,29)
     $    + c7(ilon,iniv)*cc(ilon,iniv,6)
     $    + c5(ilon,iniv)*cc(ilon,iniv,8)
     $    + c18(ilon,iniv)*cc(ilon,iniv,11)
c
c           ch3o2h
c
            pr(ilon,iniv,30) = 
     $      c7(ilon,iniv)*ck(ilon,iniv,29)*cc(ilon,iniv,6)
c
            pe(ilon,iniv,30) = 
     $      c17a(ilon,iniv)*cc(ilon,iniv,2)
     $    + c17b(ilon,iniv)*cc(ilon,iniv,2)
     $    + tj(ilon,iniv,j_ch3o2h)
     $    + rain_ch3o2h(ilon,iniv)
c
c           cfc-11
c
            pr(ilon,iniv,31) = 0.
c
            pe(ilon,iniv,31) = tj(ilon,iniv,j_cfc11)
     $                       + d90(ilon,iniv)*cc(ilon,iniv,1)
c
c           cfc-12*
c
            pr(ilon,iniv,32) = 0.
c
            pe(ilon,iniv,32) = tj(ilon,iniv,j_cfc12)
     $                       + d91(ilon,iniv)*cc(ilon,iniv,1)
c
c           cfc-113
c
            pr(ilon,iniv,33) = 0.
c
            pe(ilon,iniv,33) = tj(ilon,iniv,j_cfc113)
     $                       + d92(ilon,iniv)*cc(ilon,iniv,1)
c
c           ccl4
c
            pr(ilon,iniv,34) = 0.
c
            pe(ilon,iniv,34) = tj(ilon,iniv,j_ccl4)
     $                       + d93(ilon,iniv)*cc(ilon,iniv,1)
c
c           ch3ccl3*
c
            pr(ilon,iniv,35) = 0.
            pe(ilon,iniv,35) = tj(ilon,iniv,j_ch3ccl3)
     $                       + d94(ilon,iniv)*cc(ilon,iniv,2)
c
c           ch3cl
c
            pr(ilon,iniv,36) = 0.
            pe(ilon,iniv,36) = tj(ilon,iniv,j_ch3cl)
     $                       + d95(ilon,iniv)*cc(ilon,iniv,2)
c
c           hcfc-22*
c
            pr(ilon,iniv,37) = 0.
c
            pe(ilon,iniv,37) = tj(ilon,iniv,j_hcfc22)
     $                       + d96(ilon,iniv)*cc(ilon,iniv,1)
     $                       + d97(ilon,iniv)*cc(ilon,iniv,2)
c
c           ch3br
c
            pr(ilon,iniv,38) = 0.
c
            pe(ilon,iniv,38) = tj(ilon,iniv,j_ch3br)
     $                       + e50(ilon,iniv)*cc(ilon,iniv,1)
     $                       + e51(ilon,iniv)*cc(ilon,iniv,2)
c
c           h-1211
c
            pr(ilon,iniv,39) = 0.
c
            pe(ilon,iniv,39) = tj(ilon,iniv,j_h1211)
     $                       + e52(ilon,iniv)*cc(ilon,iniv,1)
c
c           h-1301*
c
            pr(ilon,iniv,40) = 0.
c
            pe(ilon,iniv,40) = tj(ilon,iniv,j_h1301)
     $                       + e53(ilon,iniv)*cc(ilon,iniv,1)
c
c           bry
c
            pr(ilon,iniv,41) = 
     $      (tj(ilon,iniv,j_ch3br)
     $    + e50(ilon,iniv)*cc(ilon,iniv,1) 
     $    + e51(ilon,iniv)*cc(ilon,iniv,2))*ck(ilon,iniv,38)
     $    +(tj(ilon,iniv,j_h1211)
     $    + e52(ilon,iniv)*cc(ilon,iniv,1))*ck(ilon,iniv,39)
     $    +(tj(ilon,iniv,j_h1301)
     $    + e53(ilon,iniv)*cc(ilon,iniv,1))*ck(ilon,iniv,40)
     $    + 2.*(e54(ilon,iniv)*cc(ilon,iniv,2)
     $    + tj(ilon,iniv,j_ch2br2))*ck(ilon,iniv,42)
c
            pe(ilon,iniv,41) = rain_hbr(ilon,iniv)*ck(ilon,iniv,19)
            pe(ilon,iniv,41) = pe(ilon,iniv,41)/ck(ilon,iniv,41)
c
c           ch2br2*
c
            pr(ilon,iniv,42) = 0.
c
            pe(ilon,iniv,42) = tj(ilon,iniv,j_ch2br2)
     $                       + e54(ilon,iniv)*cc(ilon,iniv,2)
c
c           hno3 phase gazeuse
c
            pr(ilon,iniv,nbcon) = 0.
c
            pe(ilon,iniv,nbcon) = 0.
c
         end do
         end do
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     tendance chimique pour les constituents a longue duree de vie
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         do n = 1,nbcon - 1
            do iniv = 1,nivbas
               do ilon = 1,nlon
                  ck(ilon,iniv,n) = 
     $            (ck0(ilon,iniv,n) + dt*pr(ilon,iniv,n))
     $            /(1. + dt*pe(ilon,iniv,n))
               end do
            end do
         end do
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     fin de la boucle sur les iterations
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
10050 continue
c
      do n = 1,nbcon - 1
         do iniv = 1,nivbas
            do ilon = 1,nlon
               qj1(ilon,ilat,iniv,n) = ck(ilon,iniv,n)/hnm(ilon,iniv)
            end do
         end do
      end do
c
c     hno3 phase gazeuse
c
      do iniv = 1,nivbas
         do ilon = 1,nlon
            qj1(ilon,ilat,iniv,nbcon) = 
     $                               ck(ilon,iniv,5)*parthno3(ilon,iniv)
     $                              /hnm(ilon,iniv)
         end do
      end do
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     stockage de la concentration des especes
c     a courte duree de vie
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do n = 1,ncm
         do iniv = 1,nivbas
            do ilon = 1,nlon
               hc(ilon,ilat,iniv,n) =
     $         cc(ilon,iniv,n)/hnm(ilon,iniv)
            end do
         end do
      end do
c
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine sedimentation( trajpre, trajtem, irapp, qj1, 
     $                          vsed3d, nivh2oecmwf, parth2osed,
     $                          parthno3sed, dt )
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calcule les quantites de h2o et hno3 sedimentees
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      parameter(nlon = 180, nlat = 91, niv = 137, nivbas = 94,
     $          nivhau = 14)
      parameter(nbcon = 43)
      parameter(nrappmax = 12)
c
      real trajpre(nlon,nlat,niv,nrappmax)
      real trajtem(nlon,nlat,niv,nrappmax)
      real qj1(nlon,nlat,niv,nbcon)
      real dt
      real hnm(nlon,nivbas)
      real vsed3d(nlon,nlat,nivbas)
      real parth2osed(nlon,nlat,nivbas), parthno3sed(nlon,nlat,nivbas)
      real dz, fracup, fracdown
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     boucle sur les latitudes 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do ilat = 1,nlat
c
ccc      calcul de la concentration
c
         do iniv = 1,nivbas
            do ilon = 1,nlon
               hnm(ilon,iniv) = trajpre(ilon,ilat,iniv,irapp)
     $                        /(trajtem(ilon,ilat,iniv,irapp)*
     $                          1.38e-19)
            end do
         end do
c
ccc      calcul des fractions sedimentees
c
         do iniv = nivh2oecmwf-1,nivhau,-1
            do ilon = 1,nlon
c
               dz = 7000.*log(trajpre(ilon,ilat,iniv  ,irapp)
     $                       /trajpre(ilon,ilat,iniv-1,irapp))
c
               fracup   = vsed3d(ilon,ilat,iniv-1)*dt/dz
               fracdown = vsed3d(ilon,ilat,iniv  )*dt/dz
c
               fracup   = min(fracup,1.0)
               fracdown = min(fracdown,1.0)
c
ccc            sedimentation de h2o
c
               qj1(ilon,ilat,iniv,3) = 
     $                      qj1(ilon,ilat,iniv,3)
     $                    + qj1(ilon,ilat,iniv-1,3)
     $                    *(1. - parth2osed(ilon,ilat,iniv-1))
     $                    *hnm(ilon,iniv-1)/hnm(ilon,iniv)*fracup
     $                    - qj1(ilon,ilat,iniv,3)
     $                    *(1. - parth2osed(ilon,ilat,iniv))*fracdown
c
            end do
         end do
c
         do iniv = nivbas,nivhau,-1
            do ilon = 1,nlon
c
               dz = 7000.*log(trajpre(ilon,ilat,iniv  ,irapp)
     $                       /trajpre(ilon,ilat,iniv-1,irapp))
c
               fracup   = vsed3d(ilon,ilat,iniv-1)*dt/dz
               fracdown = vsed3d(ilon,ilat,iniv  )*dt/dz
c
               fracup   = min(fracup,1.0)
               fracdown = min(fracdown,1.0)
c
C BUG corrige le 01/12/2012 inversion noy et hno3
ccc            sedimentation de noy
c
               qj1(ilon,ilat,iniv,4) = 
     $                      qj1(ilon,ilat,iniv,4)
     $                    + qj1(ilon,ilat,iniv-1,5)
     $                    *(1. - parthno3sed(ilon,ilat,iniv-1))
     $                    *hnm(ilon,iniv-1)/hnm(ilon,iniv)*fracup
     $                    - qj1(ilon,ilat,iniv,5)
     $                    *(1. - parthno3sed(ilon,ilat,iniv))*fracdown
c
ccc            sedimentation de hno3
c
               qj1(ilon,ilat,iniv,5) = 
     $                      qj1(ilon,ilat,iniv,5)
     $                    + qj1(ilon,ilat,iniv-1,5)
     $                    *(1. - parthno3sed(ilon,ilat,iniv-1))
     $                    *hnm(ilon,iniv-1)/hnm(ilon,iniv)*fracup
     $                    - qj1(ilon,ilat,iniv,5)
     $                    *(1. - parthno3sed(ilon,ilat,iniv))*fracdown
c
            end do
         end do
c
cc       fin de la boucle sur les latitudes
c
      end do
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine conserv( qj1, hc )
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     correction for mass conservation                                c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      parameter(nlon = 180, nlat = 91, niv = 137, nivbas = 94)
      parameter(nbcon = 43, ncm = 15)
c
      real qj1(nlon,nlat,niv,nbcon)
      real hc(nlon,nlat,niv,ncm)
      real rapcly, rapbry, rapnoy, totmas
c
      do iniv = 1,nivbas
      do ilat = 1,nlat
      do ilon = 1,nlon
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        chlorine                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         totmas =   qj1(ilon,ilat,iniv,10)
     $            + qj1(ilon,ilat,iniv,13) 
     $            + qj1(ilon,ilat,iniv,14)
     $            + qj1(ilon,ilat,iniv,15) 
     $            + 2.*qj1(ilon,ilat,iniv,16)
     $            + qj1(ilon,ilat,iniv,18)
     $            + qj1(ilon,ilat,iniv,23) 
     $            + 2.*qj1(ilon,ilat,iniv,25)
c
         rapcly = (qj1(ilon,ilat,iniv,7) - qj1(ilon,ilat,iniv,27))
     $            /totmas
c
         hc(ilon,ilat,iniv,  3) =  hc(ilon,ilat,iniv, 3)*rapcly
         hc(ilon,ilat,iniv, 11) =  hc(ilon,ilat,iniv,11)*rapcly
         qj1(ilon,ilat,iniv,10) = qj1(ilon,ilat,iniv,10)*rapcly
         qj1(ilon,ilat,iniv,13) = qj1(ilon,ilat,iniv,13)*rapcly
         qj1(ilon,ilat,iniv,14) = qj1(ilon,ilat,iniv,14)*rapcly
         qj1(ilon,ilat,iniv,15) = qj1(ilon,ilat,iniv,15)*rapcly
         qj1(ilon,ilat,iniv,16) = qj1(ilon,ilat,iniv,16)*rapcly
         qj1(ilon,ilat,iniv,18) = qj1(ilon,ilat,iniv,18)*rapcly
         qj1(ilon,ilat,iniv,23) = qj1(ilon,ilat,iniv,23)*rapcly
         qj1(ilon,ilat,iniv,25) = qj1(ilon,ilat,iniv,25)*rapcly
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        nitrogen                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         totmas = qj1(ilon,ilat,iniv,5) 
     $          + 2.*qj1(ilon,ilat,iniv,6) 
     $          + qj1(ilon,ilat,iniv,21)
     $          + qj1(ilon,ilat,iniv,22)
c
         rapnoy = max((qj1(ilon,ilat,iniv,4) - qj1(ilon,ilat,iniv,14)
     $           - qj1(ilon,ilat,iniv,18) - qj1(ilon,ilat,iniv,20)), 0.)
     $           /totmas
c
         hc(ilon,ilat,iniv,  7) =  hc(ilon,ilat,iniv, 7)*rapnoy
         hc(ilon,ilat,iniv,  8) =  hc(ilon,ilat,iniv, 8)*rapnoy
         hc(ilon,ilat,iniv, 10) =  hc(ilon,ilat,iniv,10)*rapnoy
         hc(ilon,ilat,iniv, 13) =  hc(ilon,ilat,iniv,13)*rapnoy
         qj1(ilon,ilat,iniv ,5) = qj1(ilon,ilat,iniv, 5)*rapnoy
         qj1(ilon,ilat,iniv ,6) = qj1(ilon,ilat,iniv, 6)*rapnoy
         qj1(ilon,ilat,iniv,21) = qj1(ilon,ilat,iniv,21)*rapnoy
         qj1(ilon,ilat,iniv,22) = qj1(ilon,ilat,iniv,22)*rapnoy
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        bromine                                                      c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         totmas = qj1(ilon,ilat,iniv,19)
     $          + qj1(ilon,ilat,iniv,20)
     $          + qj1(ilon,ilat,iniv,24)
     $          + qj1(ilon,ilat,iniv,26)
     $          + qj1(ilon,ilat,iniv,27)
c
         rapbry = qj1(ilon,ilat,iniv,41)/totmas
c
         hc(ilon,ilat,iniv,  9) =  hc(ilon,ilat,iniv, 9)*rapbry
         hc(ilon,ilat,iniv, 12) =  hc(ilon,ilat,iniv,12)*rapbry
         qj1(ilon,ilat,iniv,19) = qj1(ilon,ilat,iniv,19)*rapbry
         qj1(ilon,ilat,iniv,20) = qj1(ilon,ilat,iniv,20)*rapbry
         qj1(ilon,ilat,iniv,24) = qj1(ilon,ilat,iniv,24)*rapbry
         qj1(ilon,ilat,iniv,26) = qj1(ilon,ilat,iniv,26)*rapbry
         qj1(ilon,ilat,iniv,27) = qj1(ilon,ilat,iniv,27)*rapbry
c
      end do
      end do
      end do
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine analytic(temp,ptot,hnm,ck,vsed,
     $                    k1,k2,k3,k4,k5,k6,k7,k8,k9,
     $                    parthno3,parthno3sl,parth2o,parthcl,
     $                    parthbr,
     $                    surfarea,ilat)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      parameter(nlon = 180, nlat = 91, niv = 137, nivbas = 94,
     $          nivhau = 14)
      parameter (nbcon = 43)
      parameter (r = 8.205e-5)
      parameter (sigma = 1.8)
c
      real surfarea(nlon,nlat,niv)
      real temp(nlon,nivbas), ptot(nlon,nivbas), hnm(nlon,nivbas)
      real ck(nlon,nivbas,nbcon)
      real k1(nlon,nivbas), k2(nlon,nivbas), k3(nlon,nivbas),
     $     k4(nlon,nivbas), k5(nlon,nivbas), k6(nlon,nivbas),
     $     k7(nlon,nivbas), k8(nlon,nivbas), k9(nlon,nivbas)
c
      real parthno3l(nlon,nivbas), parthno3s(nlon,nivbas)
      real parthno3sl(nlon,nivbas)
      real parthno3(nlon,nivbas)
      real parth2o(nlon,nivbas), parthcl(nlon,nivbas)
c
      real a, b, det, msb,mnb,mnf,msf
      real qn(10),qs(10),kn(7),ks(7)
      save qn,qs,kn,ks
c
      real ndrop
      real dens(nlon,nivbas)
      real ns(nlon,nivbas)
      real pn0(nlon,nivbas), phcl0(nlon,nivbas) 
      real pw(nlon,nivbas), ms(nlon,nivbas), mn(nlon,nivbas),
     $     ws(nlon,nivbas), wn(nlon,nivbas)
      real hhcl(nlon,nivbas), hhocl(nlon,nivbas),
     $     hhobr(nlon,nivbas), hclono2(nlon,nivbas),
     $     hhbr(nlon,nivbas)
      real tice(nlon,nivbas)
      real mh2so4c, bh2so4, ph2so4, muh2so4
c
      integer lnat(nlon,nivbas), lice(nlon,nivbas)
      integer condliq(nlon,nivbas)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ajouts jpl 2000
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real aw(nlon,nivbas)
      real wt(nlon,nivbas), mh2so4(nlon,nivbas)
      real m, y1, y2, bigx
      real a1(3), b1(3), c1(3), d1(3)
      real a2(3), b2(3), c2(3), d2(3)
c
      data a1/12.37208932, 11.820654354, -180.06541028/
      data b1/-0.16125516114, -0.20786404244, -0.38601102592/
      data c1/-30.490657554, -4.807306373, -93.317846778/
      data d1/-2.1133114241, -5.1727540348, 273.88132245/
      data a2/13.455394705, 12.891938068, -176.95814097/
      data b2/-0.1921312255,-0.232338847708,-0.36257048154/
      data c2/-34.285174607, -6.4261237757, -90.469744201/
      data d2/-1.7620073078, -4.9005471319, 267.45509988/
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real rmean(nlon,nivbas)
      real aliq(nlon,nivbas)
c
      real nnat_small_max, vnat_small_max
      real nnats(nlon,nivbas)
      real nnatl(nlon,nivbas)
      real rnats, rnatl
      real anat(nlon,nivbas) 
      real vnat, vnatl
c
      real nice(nlon,nivbas)
      real rice
      real aice(nlon,nivbas)
      real vice
c
c     sedimentation
c
      real vsed(nlon,nivbas)
      real rhoair, rhoice, rhonat
      real netha, lambda, lambda0, nre, us, cdnre2
      real bnre(0:6), tcel
c
      logical pscliq, pscnat, pscice
c
c     pruppacher and klett, 1988
c
      data bnre /-3.18657,0.992696,-0.153193e-2,
     $           -0.987059e-3,-0.578878e-3,0.855176e-4,
     $           -0.327815e-5/
c
c     carslaw et al., 1995
c
      data qn/14.5734,0.0615994,-1.14895,0.691693, -0.098863,
     $     0.0051579,0.123472,-0.115574,0.0110113,0.0097914/
      data qs/14.4700,0.0638795,-3.29597,1.778224,-0.223244,
     $     0.0086486,0.536695,-0.335164,0.0265153,0.0157550/
      data kn/-39.136,6358.4,83.29,-17650.0,198.53,-11948,-28.469/
      data ks/-21.661,2724.2,51.81,-15732.0,47.004,-6969.0,-4.6183/
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     choix du schema de psc:
c   
c     pscliq = true : autoriser les aerosols liquides
c     pscnat = true : autoriser le nat
c     pscice = true : autoriser la glace
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      pscliq = .true.
      pscnat = .true.
      pscice = .true.
c
      ndrop = 10.                 ! densite d'aerosols liquides (cm-3)
c
      nnat_small_max  = 1.        ! densite maximale de nat (cm-3)
      rnats           = 0.5e-4    ! rayon du small mode de nat (cm)
      rnatl           = 6.5e-4    ! rayon du large mode de nat (cm)
c
      rice = 10.e-4               ! rayon de la glace
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialisations
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      pi = 2.*asin(1.0)
      amu  = 1.66e-24
c
      do iniv = 1,nivbas
         do ilon = 1,nlon
            rmean(ilon,iniv)      = 0.
            aliq(ilon,iniv)       = 0.
            lnat(ilon,iniv)       = 0
            nnats(ilon,iniv)      = 0.
            nnatl(ilon,iniv)      = 0.
            vnat                  = 0.
            vnatl                 = 0.
            anat(ilon,iniv)       = 0.
            lice(ilon,iniv)       = 0
            nice(ilon,iniv)       = 0.
            vice                  = 0.
            aice(ilon,iniv)       = 0.
            vsed(ilon,iniv)       = 0.
            condliq(ilon,iniv)    = 0
            parthno3(ilon,iniv)   = 1.
            parthno3l(ilon,iniv)  = 1.
            parthno3s(ilon,iniv)  = 1.
            parthno3sl(ilon,iniv) = 1.
            parth2o(ilon,iniv)    = 1.
            parthcl(ilon,iniv)    = 1.
         end do
      end do
c
      do iniv = nivhau,nivbas
         do ilon = 1,nlon
            qh2o = ck(ilon,iniv,3)/hnm(ilon,iniv)
            pw(ilon,iniv) = qh2o*ptot(ilon,iniv)/1013.0
            pw(ilon,iniv) = min(pw(ilon,iniv),2.e-3/1013.)
c
            qhno3 = ck(ilon,iniv,5)/hnm(ilon,iniv)
            pn0(ilon,iniv)= qhno3*ptot(ilon,iniv)/1013.0
         end do
      end do
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     traitement des aerosols solides
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cc    traitement du nat
c
      if (pscnat) then
         do iniv = nivhau,nivbas
         do ilon = 1,nlon
            t = min(temp(ilon,iniv),273.)
            pn0t = pn0(ilon,iniv)*1013.*0.75
            zmt = -2.7836 - 0.00088*t
            zbt = 38.9855 - 11397.0/t + 0.009179*t
            zhno3eq = 10.0**(zmt*log10(pw(ilon,iniv)*1013.*0.75)
     $                       + zbt)
            if (pn0t/zhno3eq .ge. 1.) then
               lnat(ilon,iniv) = 1
               pn = zhno3eq/(0.75*1013.)
               hno3s = (pn0(ilon,iniv) - pn)*1013.
     $                 *hnm(ilon,iniv)/ptot(ilon,iniv)
               hno3s = max(hno3s,0.)
c
c              volume de nat condense
c
               vnat = hno3s*amu*117./1.35
c
c              volume maximal du small mode
c
               vnat_small_max = nnat_small_max
     $                         *(4./3.*pi*rnats**3.)
c
c              densite de nat (cm-3) pour les deux modes
c
               if (vnat .gt. vnat_small_max) then
                  nnats(ilon,iniv) = nnat_small_max
                  nnatl(ilon,iniv) = (vnat - vnat_small_max)
     $                           /(4./3.*pi*rnatl**3.)
                  vnatl = vnat - vnat_small_max
               else
                  nnats(ilon,iniv) = vnat
     $                           /(4./3.*pi*rnats**3.)
                  nnatl(ilon,iniv) = 0.
                  vnatl = 0.
               end if
c
c              surface du nat
c
               anat(ilon,iniv) = nnats(ilon,iniv)*4*pi*rnats**2
     $                         + nnatl(ilon,iniv)*4*pi*rnatl**2
c
c              parthno3s: part en phase gazeuse
c
               parthno3s(ilon,iniv) = 1.0 - (pn0(ilon,iniv)-pn)
     $                                      /pn0(ilon,iniv)
               parthno3s(ilon,iniv) = max(parthno3s(ilon,iniv),0.)
               parthno3s(ilon,iniv) = min(parthno3s(ilon,iniv),1.)
c
c              parthno3sl: part non sedimentee 
c                          (gaz + petit mode de nat)
c
               parthno3sl(ilon,iniv) = 1.0
     $                            - (pn0(ilon,iniv)-pn)*vnatl/vnat
     $                              /pn0(ilon,iniv)
               parthno3sl(ilon,iniv) = max(parthno3sl(ilon,iniv),0.)
               parthno3sl(ilon,iniv) = min(parthno3sl(ilon,iniv),1.)
            end if
         end do
         end do
      end if
c
cc    traitement de la glace
c
      if (pscice) then
         do iniv = nivhau,nivbas
         do ilon = 1,nlon
            t = min(temp(ilon,iniv),273.)
            tice(ilon,iniv) =
     $      2668.70/(10.4310-(log(pw(ilon,iniv))+log(760.0))/log(10.0))
            if (t .le. tice(ilon,iniv)) then
               lice(ilon,iniv) = 1
               lnat(ilon,iniv) = 0
               anat(ilon,iniv) = 0.
               zh2oeq = 10.0**(-2668.7/t + 10.4310)
               ph2o = zh2oeq/760.0
               h2os = max(pw(ilon,iniv) - ph2o,0.)
     $                *1013.*hnm(ilon,iniv)/ptot(ilon,iniv)
c
c              volume de glace condensee
c
               vice = h2os*amu*18./0.928
c
c              densite de glace (cm-3)
c
               nice(ilon,iniv) = vice/(4./3.*pi*rice**3.)
c
c              surface de la glace
c
               aice(ilon,iniv) = nice(ilon,iniv)*4*pi*rice**2
c
c              parth2o: part en phase gazeuse
c
               parth2o(ilon,iniv) = 1.0 - max(pw(ilon,iniv)-ph2o,0.)
     $                                       /pw(ilon,iniv)
               parth2o(ilon,iniv) = max(parth2o(ilon,iniv),0.)
               parth2o(ilon,iniv) = min(parth2o(ilon,iniv),1.)
c
c              on retire hno3 de la phase gazeuse sous forme de nat
c
               zmt = -2.7836 - 0.00088*t
               zbt = 38.9855 - 11397.0/t + 0.009179*t
               zhno3eq = 10.0**(zmt*log10(pw(ilon,iniv)
     $                          *parth2o(ilon,iniv)*1013.*0.75)+ zbt)
               pn = zhno3eq/(0.75*1013.)
c
c              parthno3s: part en phase gazeuse
c
               parthno3s(ilon,iniv) = 1.0 - (pn0(ilon,iniv)-pn)
     $                                      /pn0(ilon,iniv)
               parthno3s(ilon,iniv) = max(parthno3s(ilon,iniv),0.)
               parthno3s(ilon,iniv) = min(parthno3s(ilon,iniv),1.)
               parthno3sl(ilon,iniv) = parthno3s(ilon,iniv)
            end if
         end do
         end do
      end if
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     traitement des aerosols liquides
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (pscliq) then
c
      do iniv = nivhau,nivbas
      do ilon = 1,nlon
c
         condliq(ilon,iniv) = 1
         t = max(temp(ilon,iniv),max(tice(ilon,iniv)-3.0,185.))
         qh2o = ck(ilon,iniv,3)*parth2o(ilon,iniv)/hnm(ilon,iniv)
         pw(ilon,iniv) = qh2o*ptot(ilon,iniv)/1013.0
         pw(ilon,iniv) = min(pw(ilon,iniv),2.e-3/1013.)
         pw(ilon,iniv) = max(pw(ilon,iniv),2.e-5/1013.)
c
         qh2so4 = max(ck(ilon,iniv,12)/hnm(ilon,iniv),1.e-18)
         ns(ilon,iniv) = ptot(ilon,iniv)*100.0*qh2so4/8.314/t
c
         qhno3 = ck(ilon,iniv,5)*parthno3s(ilon,iniv)/hnm(ilon,iniv)
         pn0(ilon,iniv)= ptot(ilon,iniv)/1013.0*qhno3
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc         aerosols binaires h2so4/h2o 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         a = ks(3) + ks(4)/t
         b = ks(1) + ks(2)/t
         c = ks(5) + ks(6)/t + ks(7)*log(t) - log(pw(ilon,iniv))
c
         det = b**2 - 4.*a*c
c
         if (det .gt. 0.) then
            xsb = (-b - sqrt(det))/(2.*a)
            msb = 55.51*xsb/(1.0 - xsb)
         else
c           print*,'determinant negatif'
c           print*,'ptot = ', ptot(ilon,iniv), 't = ',t
c           print*,'qj1(h2o) = ', qh2o
            msb = 0.
            condliq(ilon,iniv) = 0
         end if
c
         ms(ilon,iniv) = msb
         ws(ilon,iniv) = msb*0.098076/(1.0 + msb*0.098076)
c
c        xsb = 1.0/(2.0*(ks(3)+ks(4)/t))*( -ks(1)-ks(2)/t-
c    $        ((ks(1)+ks(2)/t)**2-4.0*(ks(3)+ks(4)/t)
c    $         *(ks(5)+ks(6)/t+ks(7)*log(t)-log(pw(ilon,iniv))))**0.5)
c
         mn(ilon,iniv) = 0.0
         wn(ilon,iniv) = 0.0
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c           pression saturante de h2so4
c          
c           d'apres ayers et al., grl, 1980
c           et giauque et al., j. am. chem. soc., 1960
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         wts = max(41., ws(ilon,iniv)*100.)
         muh2so4 = 4.184*(1.514e4 - 286.*(wts - 40.)
     $           + 1.08*(wts - 40.)**2
     $           - 3941./(wts - 40.)**0.1)
c
c        pression saturante en atmosphere
c
         ph2so4 = exp(-10156./t + 16.2590 - muh2so4/(8.314*t)) 
c
         if (qh2so4*ptot(ilon,iniv)/1013. .le. ph2so4) then
            ms(ilon,iniv) = 0.
            ws(ilon,iniv) = 0.
            condliq(ilon,iniv) = 0
         end if

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc         aerosols ternaires hno3/h2so4/h2o
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         if (condliq(ilon,iniv) .eq. 1 .and. t .lt. 215.) then
c
            tf = t
            tt = r*tf*ns(ilon,iniv)
            tr = 1.0e4/tf-43.4782608
            pr = log(pw(ilon,iniv))+18.4
c
            hsb = qs(1)+qs(2)*tr**2+(qs(3)+qs(4)*tr+qs(5)*tr**2
     $           +qs(6)*tr**3)*pr + (qs(7)+qs(8)*tr+qs(9)*tr**2)
     $           *pr**2+qs(10)*tr*pr**3
            hsb = exp(hsb)
            xnb = 1.0/(2.0*(kn(3)+kn(4)/tf))*(-kn(1)-kn(2)/tf-
     $           ((kn(1)+kn(2)/tf)**2-4.0*(kn(3)+kn(4)/tf)
     $           *(kn(5)+kn(6)/tf +kn(7)*log(tf)
     $           -log(pw(ilon,iniv))))**0.5)
            mnb = 55.51*xnb/(1.0-xnb)
            hnb = qn(1)+qn(2)*tr**2+(qn(3)+qn(4)*tr+qn(5)*tr**2
     $           +qn(6)*tr**3)*pr + (qn(7)+qn(8)*tr+qn(9)*tr**2)
     $           *pr**2+qn(10)*tr*pr**3
            hnb = exp(hnb)
c
            a = (tt*hnb*mnb**2 - tt*hsb*mnb*msb - 2.0*mnb**2*msb
     $           + mnb*msb**2 + hnb*mnb*msb*pn0(ilon,iniv)
     $          - hsb*msb**2*pn0(ilon,iniv))/(mnb**2 - mnb*msb)
            b = msb*(-2.0*tt*hnb*mnb+tt*hsb*msb+mnb*msb
     $         -hnb*msb*pn0(ilon,iniv))/(mnb-msb)
            c = (tt*hnb*mnb*msb**2)/(mnb - msb)
c
            phi = atan(sqrt(
     $            max(4.0*(a**2-3.0*b)**3-(-2.0*a**3+9.0*a*b
     $           -27.0*c)**2,0.))/(-2.0*a**3+9.0*a*b-27.0*c) )
            if (phi .lt. 0.) then
               phi = phi + pi
            end if
c
            msf = -1.0/3.0*(a+2.0*sqrt(a**2-3.0*b)*cos((pi+phi)/3.0))
            msf = max(msf,0.)
c
            mnf = mnb*(1.0-msf/msb)
            mnf = max(mnf,0.)
            pn = mnf/(hnb*mnf/(mnf+msf)+hsb*msf/(mnf+msf))
            wsf = msf*0.098076/(1.0+msf*0.098076+mnf*0.063012)
            wnf = mnf*0.063012/(1.0+msf*0.098076+mnf*0.063012)
c 
            ms(ilon,iniv) = msf
            mn(ilon,iniv) = mnf
            ws(ilon,iniv) = wsf
            wn(ilon,iniv) = wnf
c
c           parthno3l: part en phase gazeuse
c
            parthno3l(ilon,iniv) = 1.0 - (pn0(ilon,iniv)-pn)
     $                                   /pn0(ilon,iniv)
            parthno3l(ilon,iniv) = min(parthno3l(ilon,iniv),1.)
c
         end if
c
c        dilution minimale de h2so4: 0.5%
c
         ws(ilon,iniv) = max(ws(ilon,iniv), 0.005)
c
      end do
      end do
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   densite des aerosols ternaires (g/cm3)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call density(ws,wn,temp,dens)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc      total specific liquid aerosol volume (dimensionless)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do iniv = nivhau,nivbas
      do ilon = 1,nlon
c
         if (condliq(ilon,iniv) .eq. 1) then
            vliq = ns(ilon,iniv)*98.076e-6/ws(ilon,iniv)/dens(ilon,iniv)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc         liquid aerosol radius and surface area
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c           rmean=mean radius of liquid droplets (units = cm)
c           this varies with the total liquid volume per unit volume of air
c           (called vliq in main program).  rmean can be calculated
c           from vliq using
c
c           rmean=rmode*exp(0.5*(log(sigma))**2), where
c
c           rmode=(3.0*vliq/(4.*pi*ndrop)*
c                 exp(-9.0/2.0*(log(sigma)**2)))**(1./3.)
c
c           ndrop=number of droplets per cm3 air
c           sigma=width of lognormal (use sigma=1.8).
c
            rmode = (3.0*vliq/(4.*pi*ndrop)*
     $              exp(-9.0/2.0*(log(sigma)**2)))**(1./3.)
            rmean(ilon,iniv) = rmode*exp(0.5*(log(sigma))**2)
c
            aliq(ilon,iniv) = 4.0*pi*rmode**2*ndrop
     $                        *exp(2.0*(log(sigma))**2)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc         the solubility of hcl
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c           hhcl est calcule directement a partir de jpl 2000,
c           lui-meme utilisant carslaw et al., j. phys. chem., 1995.
c           l'ensemble des calculs suppose un aerosol binaire:
c           la molalite de h2so4 (wt) est donc recalculee dans
c           ces conditions.
c
            t = max(temp(ilon,iniv),max(tice(ilon,iniv)-3.0,185.))
c
c           ph2o0 : saturation water vapor pressure, hpa
c
            ph2o0 = exp(18.452406985 - 3505.1578807/t
     $                  - 330918.55082/t**2 + 12725068.262/t**3)
c
c           aw : water activity
c
            aw(ilon,iniv) = pw(ilon,iniv)*1013./ph2o0
            aw(ilon,iniv) = min(aw(ilon,iniv), 1.)
c
            if (aw(ilon,iniv) .le. 0.05) then
               y1 = a1(1)*aw(ilon,iniv)**b1(1) 
     $            + c1(1)*aw(ilon,iniv) + d1(1)
               y2 = a2(1)*aw(ilon,iniv)**b2(1)
     $            + c2(1)*aw(ilon,iniv) + d2(1)
            else if (aw(ilon,iniv) .lt. 0.85) then
               y1 = a1(2)*aw(ilon,iniv)**b1(2)
     $            + c1(2)*aw(ilon,iniv) + d1(2)
               y2 = a2(2)*aw(ilon,iniv)**b2(2)
     $            + c2(2)*aw(ilon,iniv) + d2(2)
            else
               y1 = a1(3)*aw(ilon,iniv)**b1(3)
     $            + c1(3)*aw(ilon,iniv) + d1(3)
               y2 = a2(3)*aw(ilon,iniv)**b2(3) 
     $            + c2(3)*aw(ilon,iniv) + d2(3)
            end if
c
c           m: h2so4 molality, mol/kg
c
            m = y1 + (t - 190.)*(y2 - y1)/70.
c
c           wt: h2so4 weight percentage, %
c
            wt(ilon,iniv) = 9800.*m/(98.*m + 1000.)
            wt(ilon,iniv) = max(wt(ilon,iniv), 0.5)
c
c           mh2so4: h2so4 molarity, mol l-1
c
            mh2so4(ilon,iniv) = dens(ilon,iniv)*wt(ilon,iniv)/9.8
c
            bigx = wt(ilon,iniv)/(wt(ilon,iniv)
     $           + (100. - wt(ilon,iniv))*98./18.)
c 
            hhcl(ilon,iniv) = (0.094 - 0.61*bigx + 1.2*bigx**2)
     $                     *exp(-8.68 + (8515. - 10718.*(bigx**0.7))/t)
c
c           phcl1 : pression de hcl restant en phase gazeuse
c           ns/ms : kilogrammes (litres) d'eau par metre cube 
c
            qhcl = ck(ilon,iniv,13)/hnm(ilon,iniv)
            phcl0(ilon,iniv) = qhcl*ptot(ilon,iniv)/1013.0
c
            phcl1 = (phcl0(ilon,iniv)/(r*t))
     $          /(1./(r*t) + hhcl(ilon,iniv)*ns(ilon,iniv)
     $                      /ms(ilon,iniv)) 
c
            parthcl(ilon,iniv) = 1.0 - (phcl0(ilon,iniv)-phcl1)
     $                                 /phcl0(ilon,iniv)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc         the solubility of clono2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c           on peut negliger la diminution de la concentration
c           en phase gazeuse.
c
c           sclono2 : setchenow coefficient m-1
c
            sclono2 = 0.306 + 24./t
c
c           hclono2 : shi et al., m atm-1
c
            hclono2(ilon,iniv) = 1.6e-6*exp(4710./t)
     $                           *exp(-sclono2*mh2so4(ilon,iniv))
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc         the solubility of hocl
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c           hhocl d'apres jpl 2000 (shi et al.)
c           on peut negliger la diminution de la concentration
c           en phase gazeuse.
c
c           shocl : setchenow coefficient m-1
c
            shocl = 0.0776 + 59.18/t
c
c           hhocl : shi et al.
c
            hhocl(ilon,iniv) = 1.91e-6*exp(5862.4/t)
     $                   *exp(-shocl*mh2so4(ilon,iniv))
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc         the solubility of hobr
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c           hhobr d'apres waschewsky and abbatt, j. phys. chem. a,
c           5312-5320, 1999.
c           on peut negliger la diminution de la concentration
c           en phase gazeuse.
c
            hhobr(ilon,iniv) = 4.6e-4*exp(4.5e3/t) 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc         the solubility of hbr
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c           hhbr d'apres kleffmann et al., j. phys. chem. a,
c           8489-8495, 2000. attention a l'erreur typographique
c           contenue dans la version imprimee de l'article
c           (m3 = - 4.445 imprime au lieu de m3 = + 4.445)
c
            mh2so4c = -1.977e-4*wt(ilon,iniv)**2.
     $                -2.096e-2*wt(ilon,iniv) + 4.445
c
            bh2so4 = -8.979e-5*wt(ilon,iniv)**2.
     $               +2.141e-2*wt(ilon,iniv) - 6.067
c
            hhbr(ilon,iniv) = 10.**(mh2so4c*1000./t + bh2so4)
c
         end if
      end do
      end do
c
      end if        ! pscliq
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     surface totale des psc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do iniv = nivhau,nivbas
         do ilon = 1,nlon
            surfarea(ilon,ilat,iniv) = aliq(ilon,iniv)
     $                       + anat(ilon,iniv)*lnat(ilon,iniv)      
     $                       + aice(ilon,iniv)*lice(ilon,iniv)      
         end do
      end do
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     partition totale de hno3
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do iniv = nivhau,nivbas
         do ilon = 1,nlon
            parthno3(ilon,iniv) = parthno3l(ilon,iniv)
     $                           *parthno3s(ilon,iniv)
            parthno3(ilon,iniv) = min(parthno3(ilon,iniv), 1.)
         end do
      end do
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     sedimentation
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      rhoice   = 0.928*1000.0
      rhonat   = 1.350*1000.0
c
      do iniv = nivhau,nivbas
      do ilon = 1,nlon
c
         if (lice(ilon,iniv) .eq. 1 .or. lnat(ilon,iniv) .eq. 1) then
c
            t = temp(ilon,iniv)
c
            rpart = (lice(ilon,iniv)*rice
     $             + lnat(ilon,iniv)*rnatl)*1.e-2
            rpart = min(rpart,535.e-6)
c
            if (rpart .ge. 0.5e-6) then
c
               rhopart = lice(ilon,iniv)*rhoice
     $                 + lnat(ilon,iniv)*rhonat
               rhoair=(hnm(ilon,iniv)*28.97/6.022e23)*1000.0
c
               tcel = t - 273.15
               netha = (1.718 + 0.0049*tcel
     $                 -1.2e-5*tcel*tcel)*1.e-5
c
               cdnre2 = 32.0*(rpart)**3*(rhopart-rhoair)
     $                  *rhoair*9.80665
     $                  /(3.0*(netha)**2)
c
c              pruppacher and klett 1998
c
               ynre = bnre(0)+bnre(1)*log(cdnre2)
     $               +bnre(2)*(log(cdnre2))**2
     $               +bnre(3)*(log(cdnre2))**3
     $               +bnre(4)*(log(cdnre2))**4
     $               +bnre(5)*(log(cdnre2))**5
     $               +bnre(6)*(log(cdnre2))**6
               nre = exp(ynre)
c
               if (nre .ge. 1.e-2) then
c
c                 big particles
c                (1.e-2 < reynolds number < 300)
c                (10.0  < radius  < 535 microns)
c
                  vsed(ilon,iniv) = netha*nre
     $                         /(2*rhoair*rpart)
c
               else
c
c                 small particules
c                 (1.e-6 < reynolds number < 1.e-2)
c                 (0.5   < radius     < 10 microns)
c
c                 lambda0 : libre parcours moyen (m) a 1013.25 hpa
c                           et 293.15k. voir pruppacher and klett
c                           page 416.
c      
                  lambda0 = 6.6e-8
                  lambda  = lambda0
     $                     *(1013.25/ptot(ilon,iniv))
     $                     *(t/293.15)
 
                  us = 2.0*rpart**2*9.80665
     $                 *(rhopart-rhoair)/(9.0*netha)
c
ccc               solid particles:
c
c                 niels larsen - 1991
c                 columnar or elipsoid shape
c                 with a shape factor of 1.2 (shapek)
c
                  shapek = 1.2
c
                  vsedsol = us/shapek
     $                      *(1.0+1.246*lambda/rpart
     $                      +0.42*lambda/rpart
     $                      *exp(-0.87*rpart/lambda))
c
                  vsed(ilon,iniv) = vsedsol
c
               end if        ! nre .ge. 1.e-2
            end if           ! rpart .ge. 0.5e-6
         end if              ! lice .eq. 1 .or. lnat .eq. 1
      end do
      end do
c
c     only the parameters that are needed in 'hetero' are 
c     included in the variable list.  output variables should
c     be included in the list as needed!
c
      call hetero(aw,temp,ptot,hnm,ck,ws,
     $       hhcl,hhocl,hhobr,hclono2,hhbr,
     $       rice,nice,rnats,rnatl,nnats,nnatl,
     $       aliq,rmean,condliq,lnat,lice,
     $       k1,k2,k3,k4,k5,k6,k7,k8,k9,wt)
c
 100  return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      subroutine hetero(aw,temp,ptot,hnm,ck,ws,
     $                  hhcl,hhocl,hhobr,hclono2,hhbr,
     $                  rice,nice,rnats,rnatl,nnats,nnatl,
     $                  aliq,rmean,condliq,lnat,lice,
     $                  k1,k2,k3,k4,k5,k6,k7,k8,k9,wt)
C
C ROUTINE TO CALCULATE UPTAKE COEFFICIENTS (GAMMA VALUES).
C GAMMA VALUES ARE INDICATED BY VARIABLES WITH PREFIX 'G', FOR
C EXAMPLE GHOCLHCL IS THE GAMMA VALUE OF HOCL DUE TO REACTION WITH
C HCL IN THE DROPLETS.
C FROM THE GAMMA VALUES, SECOND ORDER RATE CONSTANTS ARE CALCULATED.
C THESE HAVE THE PREFIX 'R' AND HAVE UNITS CM3 MOLECULE-1 S-1. FOR
C EXAMPLE, THE LOSS OF CLNO3 AND HCL DUE TO THE HETEROGENEOUS REACTION
C CLNO3+HCL -> CL2+HNO3 IS D(CLNO3)/DT (UNITS MOLECULE CM-3 S-1) =
C -RCLNO3HCL.[CLNO3].[HCL], WHERE [CLNO3] AND [HCL] ARE THE
C ****TOTAL**** AMOUNTS OF THESE SPECIES IN UNITS MOLECULE CM-3.
C
      parameter(nlon = 180, nlat = 91, niv = 137, nivbas = 94,
     $          nivhau = 14)
      PARAMETER (nbcon = 43)
C
      REAL TEMP(NLON,NIVBAS), PTOT(NLON,NIVBAS), HNM(NLON,NIVBAS)
      REAL CK(NLON,NIVBAS,NBCON)
C
      REAL K1(NLON,NIVBAS), K2(NLON,NIVBAS), K3(NLON,NIVBAS),
     $     K4(NLON,NIVBAS), K5(NLON,NIVBAS), K6(NLON,NIVBAS),
     $     K7(NLON,NIVBAS), K8(NLON,NIVBAS), K9(NLON,NIVBAS)
C
      REAL K1L(NLON,NIVBAS), K2L(NLON,NIVBAS), K3L(NLON,NIVBAS),
     $     K4L(NLON,NIVBAS), K5L(NLON,NIVBAS), K6L(NLON,NIVBAS),
     $     K7L(NLON,NIVBAS), K8L(NLON,NIVBAS), K9L(NLON,NIVBAS)
C
      REAL K1S(NLON,NIVBAS), K2S(NLON,NIVBAS), K3S(NLON,NIVBAS),
     $     K4S(NLON,NIVBAS), K5S(NLON,NIVBAS), K6S(NLON,NIVBAS),
     $     K7S(NLON,NIVBAS), K8S(NLON,NIVBAS), K9S(NLON,NIVBAS)
C
      REAL WS(NLON,NIVBAS)
      REAL HHCL(NLON,NIVBAS), HHOCL(NLON,NIVBAS),
     $     HHOBR(NLON,NIVBAS), HCLONO2(NLON,NIVBAS), 
     $     HHBR(NLON,NIVBAS)
      REAL RICE, ALIQ(NLON,NIVBAS), RMEAN(NLON,NIVBAS)
      REAL NNATS(NLON,NIVBAS), NNATL(NLON,NIVBAS)
      REAL RNATS, RNATL
      REAL NICE(NLON,NIVBAS)
C
      INTEGER LNAT(NLON,NIVBAS), LICE(NLON,NIVBAS), CONDLIQ(NLON,NIVBAS)
c
c     additions pour jpl 2000
c
      real aw(nlon,nivbas)
      real wt(nlon,nivbas)
      real kh, kh2o, khydr
      real phcl, mhcl, khcl
      real pclono2, lclono2, fclono2
      real phobr, mhobr, k1hobrhcl, k2hobrhcl
      real phbr, mhbr, lhbr, k1hoclhbr, k2hoclhbr
      real khoclhcl, lhocl, phocl, mhocl
      real lhobr, lhcl
C
C ******************************************************
C THE FOLLOWING PARAMETERS ARE NEEDED IN THIS ROUTINE!!
C ******************************************************
C
C RNAT=MEAN RADIUS OF NAT PARTICLES (CM)
C NNAT=NUMBER OF NAT PARTICLES PER CM3 AIR
C RICE=MEAN RADIUS OF ICE PARTICLES (CM)
C NICE=NUMBER OF ICE PARTICLES PER CM3 AIR
C RMEAN=MEAN RADIUS OF LIQUID AEROSOL DISTRIBUTION (CM)
C ALIQ=TOTAL AREA OF LIQUID AEROSOLS (CM2 CM-3 AIR)
C
C
      do iniv = nivhau,nivbas
      do ilon = 1,nlon
c
      if (condliq(ilon,iniv) .eq. 1) then
C
      T = TEMP(ILON,INIV)
C
C Z VARIABLES ARE SET TO AVOID DIVISION BY ZERO
C
      ZH2O    = MAX(CK(ILON,INIV, 3),1.E-08*HNM(ILON,INIV))
      ZHCL    = MAX(CK(ILON,INIV,13),1.E-13*HNM(ILON,INIV))
      ZCLONO2 = MAX(CK(ILON,INIV,14),1.E-13*HNM(ILON,INIV))
      ZHOCL   = MAX(CK(ILON,INIV,15),1.E-13*HNM(ILON,INIV))
      ZHOBR   = MAX(CK(ILON,INIV,26),1.E-15*HNM(ILON,INIV))
      ZHBR    = MAX(CK(ILON,INIV,19),1.E-15*HNM(ILON,INIV))
c
c     pressions de hcl, clono2, hobr, et hbr (atm)
c
      phcl    = (zhcl/hnm(ilon,iniv))*ptot(ilon,iniv)/1013.
      pclono2 = (zclono2/hnm(ilon,iniv))*ptot(ilon,iniv)/1013.
      phocl   = (zhocl/hnm(ilon,iniv))*ptot(ilon,iniv)/1013.
      phobr   = (zhobr/hnm(ilon,iniv))*ptot(ilon,iniv)/1013.
      phbr    = (zhbr/hnm(ilon,iniv))*ptot(ilon,iniv)/1013.
c
C =====================================================================
C ====================== CLONO2+HCL -> CL2+HNO3 =======================
C ====================== CLONO2+H2O -> HOCL+HNO3 ======================
C ========================= ON LIQUID AEROSOL =========================
C =====================================================================
c
c     d'apres jpl 2000. On ne tient pas compte de la deposition
c     de hno3 (solutions ternaires) sur les probabilites de reaction.
c
      t0   = 144.11 + 0.166*wt(ilon,iniv) - 0.015*wt(ilon,iniv)**2
     $     + 2.18e-4*wt(ilon,iniv)**3
      biga = 169.5 + 5.18*wt(ilon,iniv) - 0.0825*wt(ilon,iniv)**2
     $     + 3.27e-3*wt(ilon,iniv)**3
c
c     eta : viscosity of h2so4 solution, cp
c
      eta = biga*t**(-1.43)*exp(448./(t - t0))
c
c     ah : acid activity in molarity
c
      ah = exp(60.51 - 0.095*wt(ilon,iniv) + 0.0077*wt(ilon,iniv)**2
     $         -1.61e-5*wt(ilon,iniv)**3 - (1.76 +
     $          2.52e-4*wt(ilon,iniv)**2)*sqrt(t)
     $         + (-805.89 + 253.05*wt(ilon,iniv)**0.076)/sqrt(t))
c
c     dclono2 : clono2 diffusivity, klassen et al., cm2 cm-1
c
      dclono2 = 5.e-8*t/eta
c
c     cbar
c
      cclono2 = 1474.*sqrt(t)
c
c     kh2o : shi et al., s-1
c
      kh2o = 1.95e10*exp(-2800./t)
c
c     kh : shi et al., s-1
c
      kh = 1.22e12*exp(-6200./t)
c
c     khydr : shi et al., s-1
c
      khydr = kh2o*aw(ilon,iniv) + kh*ah*aw(ilon,iniv)
c
      gbh2o = 4*hclono2(ilon,iniv)*0.082*t*sqrt(dclono2*khydr)
     $      / cclono2
c
      mhcl = hhcl(ilon,iniv)*phcl
c
c     khcl : shi et al., s-1
c
      khcl = 7.9e11*ah*dclono2*mhcl
c
c     lclono2 : reacto-diffusive length, cm
c
      lclono2 = sqrt(dclono2/(khydr + khcl))
c
c     fclono2
c
      fclono2 = 1./tanh(rmean(ilon,iniv)/lclono2)
     $        - lclono2/rmean(ilon,iniv)
c
      grxn  = fclono2*gbh2o*sqrt(1. + khcl/khydr)
      gbhcl = grxn*khcl/(khcl + khydr)
      gs    = 66.12*hclono2(ilon,iniv)*mhcl*exp(-1374./t)
c
      fhcl  = 1.
     $      /(1 + 0.612*(gs + gbhcl)*pclono2/phcl)
c
      gps    = fhcl*gs
      gbhclp = fhcl*gbhcl
      gb     = gbhclp + grxn*khydr/(khcl + khydr)
c
      gclono2 = 1./(1. + 1./(gps + gb))
c
      gclono2hcl = gclono2*(gps + gbhclp)/(gps + gb)
      gclono2h2o = gclono2 - gclono2hcl
c
      rclono2h2o = 0.25*gclono2h2o*cclono2*aliq(ilon,iniv)/zh2o
      k1l(ilon,iniv) = rclono2h2o
c
      rclono2hcl = 0.25*gclono2hcl*cclono2*aliq(ilon,iniv)/zhcl
      k2l(ilon,iniv) = rclono2hcl
c
C =====================================================================
C ========================== HOCL + HCL ===============================
C ======================= ON LIQUID AEROSOL ===========================
C =====================================================================
c
c     d'apres jpl 2000
c
c     cbar
c
      chocl = 2009.*sqrt(t)
c
c     dhocl : hocl diffusivity, klassen et al., cm2 cm-1
c
      dhocl = 6.4e-8*t/eta
c
c     khoclhcl : shi et al., s-1
c
      khoclhcl = 1.25e9*ah*dhocl*mhcl
c
c     lhocl : reacto-diffusive length, cm
c
      lhocl = sqrt(dhocl/khoclhcl)
c
c     fhocl
c
      fhocl = 1./tanh(rmean(ilon,iniv)/lhocl)
     $        - lhocl/rmean(ilon,iniv)
c
      grxn = 4.*hhocl(ilon,iniv)*0.082*t*sqrt(dhocl*khoclhcl)/chocl
c
      ghoclhcl = 1./(1. + 1./(fhocl*grxn*fhcl))
c
      rhoclhcl = 0.25*ghoclhcl*chocl*aliq(ilon,iniv)/zhcl
      k5l(ilon,iniv) = rhoclhcl
C
C ====================================================================
C =========================== HOBR + HCL =============================
C ======================= ON LIQUID AEROSOL ==========================
C ====================================================================
c
c     d'apres waschewsky and abbatt, j. phys. chem. a,
c     5312-5320, 1999.
c
c     mhobr : concentration de hobr en phase liquide
c
      mhobr = hhobr(ilon,iniv)*phobr
c
      k2hobrhcl = exp(0.542*wt(ilon,iniv) - 6.44e3/t + 10.3)
c
      ratio = mhobr/mhcl
c
      if (ratio .lt. 1.) then
c
c        hcl est en exces dans l'aerosol.
c        hobr est le reactif limitant.
c
         k1hobrhcl = k2hobrhcl*mhcl
c
c        dhobr : hobr diffusivity, klassen et al., cm2 cm-1
c
         dhobr = 6.2e-8*t/eta
c
c        lhobr : reacto-diffusive length, cm
c
         lhobr = sqrt(dhobr/k1hobrhcl)
c
         fhobr = 1./tanh(rmean(ilon,iniv)/lhobr)
     $         - lhobr/rmean(ilon,iniv)
c
         chobr = 1477.*sqrt(t)
c
         ghobrhcl = 4.*hhobr(ilon,iniv)*0.082*t
     $              *sqrt(dhobr*k1hobrhcl)*fhobr/chobr
c
         rhobrhcl = 0.25*ghobrhcl*chobr*aliq(ilon,iniv)/zhcl
      else
c
c        hobr est en exces dans l'aerosol.
c        hcl est le reactif limitant.
c
         k1hobrhcl = k2hobrhcl*mhobr
c
c        dhcl : hcl diffusivity, klassen et al., cm2 cm-1
c
         dhcl = 7.8e-8*t/eta
c
c        lhcl : reacto-diffusive length, cm
c
         lhcl = sqrt(dhcl/k1hobrhcl)
c
         fhcl = 1./tanh(rmean(ilon,iniv)/lhcl)
     $        - lhcl/rmean(ilon,iniv)
c
         chcl  = 2408.*sqrt(t)
c
         ghobrhcl = 4.*hhcl(ilon,iniv)*0.082*t
     $              *sqrt(dhcl*k1hobrhcl)*fhcl/chcl
c
         rhobrhcl = 0.25*ghobrhcl*chcl*aliq(ilon,iniv)/zhobr
      end if
c
      k7l(ilon,iniv) = rhobrhcl
C
C ====================================================================
C =========================== HOBR + HBR =============================
C ======================= ON LIQUID AEROSOL ==========================
C ====================================================================
c
      k8l(ilon,iniv) = 0.
C
C ====================================================================
C =========================== HOCL + HBR =============================
C ======================= ON LIQUID AEROSOL ==========================
C ====================================================================
c
c     mhbr  : concentration de hbr  en phase liquide
c     mhocl : concentration de hocl en phase liquide
c
c     mhbr  = hhbr(ilon,iniv)*phbr
c     mhocl = hhocl(ilon,iniv)*phocl
c
c     k2hoclhbr: abbatt and nowak, j. phys. chem. a, 2131, 1997.
c     valeur a 228k et 69.3% h2so4
c
c     k2hoclhbr = 2.e6
c
c     ratio = mhbr/mhocl 
c
c     if (ratio .lt. 1.) then
c
c        hocl est en exces dans l'aerosol.
c        hbr est le reactif limitant.
c
c        k1hoclhbr = k2hoclhbr*mhocl
c
c        dhbr : hbr diffusivity, klassen et al., cm2 cm-1
c
c        dhbr = 7.9e-8*t/eta
c        chbr = 1616.*sqrt(t)
c
c        lhbr : reacto-diffusive length, cm
c
c        lhbr = sqrt(dhbr/k1hoclhbr)
c
c        fhbr = 1./tanh(rmean(ilon,iniv)/lhbr)
c    $         - lhbr/rmean(ilon,iniv)
c
c        ghoclhbr = 4.*hhbr(ilon,iniv)*0.082*t
c    $              *sqrt(dhbr*k1hoclhbr)*fhbr/chbr
c
c     else
c
c        hbr est en exces dans l'aerosol.
c        hocl est le reactif limitant.
c
c        k1hoclhbr = k2hoclhbr*mhbr
c
c        lhocl : reacto-diffusive length, cm
c
c        lhocl = sqrt(dhocl/k1hoclhbr)
c
c        fhocl = 1./tanh(rmean(ilon,iniv)/lhocl)
c    $         - lhocl/rmean(ilon,iniv)
c
c        ghoclhbr = 4.*hhocl(ilon,iniv)*0.082*t
c    $              *sqrt(dhocl*k1hoclhbr)*fhocl/chocl
c     end if
c
      k9l(ilon,iniv) = 0.
C
C =====================================================================
C ======================== N2O5 + H2O =================================
C ====================== ON LIQUID AEROSOL ============================
C =====================================================================
C
      GN2O5H2O=0.1
      CBAR=1400.1*SQRT(T)
C
      RN2O5H2O=0.25*GN2O5H2O*CBAR*ALIQ(ILON,INIV)/ZH2O
      K3L(ILON,INIV) = RN2O5H2O
c
C =====================================================================
C ======================== N2O5 + HCl =================================
C ====================== ON LIQUID AEROSOL ============================
C =====================================================================
C
      K4L(ILON,INIV) = 0.
C
C =====================================================================
C ======================== BrONO2 + H2O ===============================
C ====================== ON LIQUID AEROSOL ============================
C =====================================================================
C
c     parametrisation jpl 2000 (hanson, private comm.)
c
      gxrn = exp(29.24 - 0.396*100.*ws(ilon,iniv)) + 0.114
      gbrno3h2o = 1./(1./0.805 + 1./gxrn)
c
      CBAR=1221.4*SQRT(T)
C
      RBRNO3H2O=0.25*GBRNO3H2O*CBAR*ALIQ(ILON,INIV)/ZH2O
      K6L(ILON,INIV) = RBRNO3H2O
C
      else
         k1l(ilon,iniv) = 0.
         k2l(ilon,iniv) = 0.
         k3l(ilon,iniv) = 0.
         k4l(ilon,iniv) = 0.
         k5l(ilon,iniv) = 0.
         k6l(ilon,iniv) = 0.
         k7l(ilon,iniv) = 0.
         k8l(ilon,iniv) = 0.
         k9l(ilon,iniv) = 0.
      end if
c
      end do
      end do
C
      do iniv = nivhau,nivbas
      do ilon = 1,nlon
C
      IF (LNAT(ILON,INIV) .EQ. 1) THEN
C
      T = TEMP(ILON,INIV)
C
C Z VARIABLES ARE SET TO AVOID DIVISION BY ZERO
C
      ZH2O  = MAX(CK(ILON,INIV, 3),1.E-08*HNM(ILON,INIV))
      ZHCL  = MAX(CK(ILON,INIV,13),1.E-13*HNM(ILON,INIV))
      ZHOCL = MAX(CK(ILON,INIV,15),1.E-13*HNM(ILON,INIV))
      ZHOBR = MAX(CK(ILON,INIV,26),1.E-13*HNM(ILON,INIV))
C
C =====================================================================
C ====================== CLNO3 + HCL ==================================
C ====================== CLNO3 + H2O ==================================
C ====================== N2O5  + H2O ==================================
C ====================== N2O5  + HCL ==================================
C ====================== HOCL  + HCL ==================================
C ====================== BRNO3 + H2O  =================================
C ====================== HOBR  + HCL  =================================
C ====================== HOBR  + HBR  =================================
C ====================== HOCL  + HBR  =================================
C ======================= ON NAT ======================================
C =====================================================================
C
C GAMMA VALUES
C
      GCLNO3HCLNAT=0.2
      GCLNO3H2ONAT=0.004
      GN2O5H2ONAT=4.0E-4
      GN2O5HCLNAT=3.0E-3
      GHOCLHCLNAT=0.1
      GBRNO3H2ONAT=0.
      GHOBRHCLNAT=0.
      GHOBRHBRNAT=0.
      GHOCLHBRNAT=0.
C
C BECAUSE NAT PARTICLES CAN BE VERY LARGE, GAS DIFFUSION LIMITATION
C IS TAKEN INTO ACCOUNT. THE SECOND ORDER RATE CONSTANTS ARE
C CALCULATED FROM THE EQUATION R=G*PI*R**2*CBAR*N/(1+3G.R/(4.L)), WHERE
C G=GAMMA, R=PARTICLE RADIUS, CBAR=MEAN MOLECULAR SPEED, L=MEAN FREE
C PATH (TURCO ET AL, JGR, 94, 16493, 1989).
C
      RCLNO3H2ONAT =
     >     GCLNO3H2ONAT*4.56E4*SQRT(T/97.45)*RNATS**2
     >     *NNATS(ILON,INIV)/
     >     (1.0+3.3E4*GCLNO3H2ONAT*RNATS*PTOT(ILON,INIV)/T)
     >     /ZH2O
     >   + GCLNO3H2ONAT*4.56E4*SQRT(T/97.45)*RNATL**2
     >     *NNATL(ILON,INIV)/
     >     (1.0+3.3E4*GCLNO3H2ONAT*RNATL*PTOT(ILON,INIV)/T)
     >     /ZH2O
      K1S(ILON,INIV) = RCLNO3H2ONAT
C
      RCLNO3HCLNAT=
     >     GCLNO3HCLNAT*4.56E4*SQRT(T/97.45)*RNATS**2
     >     *NNATS(ILON,INIV)/
     >     (1.0+3.3E4*GCLNO3HCLNAT*RNATS*PTOT(ILON,INIV)/T)
     >     /ZHCL
     >   + GCLNO3HCLNAT*4.56E4*SQRT(T/97.45)*RNATL**2
     >     *NNATL(ILON,INIV)/
     >     (1.0+3.3E4*GCLNO3HCLNAT*RNATL*PTOT(ILON,INIV)/T)
     >     /ZHCL
      K2S(ILON,INIV) = RCLNO3HCLNAT
C
      RN2O5H2ONAT  =
     >     GN2O5H2ONAT*4.56E4*SQRT(T/108.0)*RNATS**2
     >     *NNATS(ILON,INIV)/
     >     (1.0+3.3E4*GN2O5H2ONAT*RNATS*PTOT(ILON,INIV)/T)
     >     /ZH2O
     >   + GN2O5H2ONAT*4.56E4*SQRT(T/108.0)*RNATL**2
     >     *NNATL(ILON,INIV)/
     >     (1.0+3.3E4*GN2O5H2ONAT*RNATL*PTOT(ILON,INIV)/T)
     >     /ZH2O
      K3S(ILON,INIV) = RN2O5H2ONAT
C
      RBRNO3H2ONAT =
     >     GBRNO3H2ONAT*4.56E4*SQRT(T/142.0)*RNATS**2
     >     *NNATS(ILON,INIV)/
     >     (1.0+3.3E4*GBRNO3H2ONAT*RNATS*PTOT(ILON,INIV)/T)
     >     /ZH2O
     >   + GBRNO3H2ONAT*4.56E4*SQRT(T/142.0)*RNATL**2
     >     *NNATL(ILON,INIV)/
     >     (1.0+3.3E4*GBRNO3H2ONAT*RNATL*PTOT(ILON,INIV)/T)
     >     /ZH2O
      K6S(ILON,INIV) = RBRNO3H2ONAT
C
      RN2O5HCLNAT=
     >     GN2O5HCLNAT*4.56E4*SQRT(T/108.0)*RNATS**2
     >     *NNATS(ILON,INIV)/
     >     (1.0+3.3E4*GN2O5HCLNAT*RNATS*PTOT(ILON,INIV)/T)
     >     /ZHCL
     >   + GN2O5HCLNAT*4.56E4*SQRT(T/108.0)*RNATL**2
     >     *NNATL(ILON,INIV)/
     >     (1.0+3.3E4*GN2O5HCLNAT*RNATL*PTOT(ILON,INIV)/T)
     >     /ZHCL
      K4S(ILON,INIV) = RN2O5HCLNAT
C
      RHOCLHCLNAT=
     >     GHOCLHCLNAT*4.56E4*SQRT(T/52.5)*RNATS**2
     >     *NNATS(ILON,INIV)/
     >     (1.0+3.3E4*GHOCLHCLNAT*RNATS*PTOT(ILON,INIV)/T)
     >     /ZHCL
     >   + GHOCLHCLNAT*4.56E4*SQRT(T/52.5)*RNATL**2
     >     *NNATL(ILON,INIV)/
     >     (1.0+3.3E4*GHOCLHCLNAT*RNATL*PTOT(ILON,INIV)/T)
     >     /ZHCL
      K5S(ILON,INIV) = RHOCLHCLNAT
C
      RHOBRHCLNAT=
     >     GHOBRHCLNAT*4.56E4*SQRT(T/96.9)*RNATS**2
     >     *NNATS(ILON,INIV)/
     >     (1.0+3.3E4*GHOBRHCLNAT*RNATS*PTOT(ILON,INIV)/T)
     >     /ZHCL
     >   + GHOBRHCLNAT*4.56E4*SQRT(T/96.9)*RNATL**2
     >     *NNATL(ILON,INIV)/
     >     (1.0+3.3E4*GHOBRHCLNAT*RNATL*PTOT(ILON,INIV)/T)
     >     /ZHCL
      K7S(ILON,INIV) = RHOBRHCLNAT
C
      K8S(ILON,INIV) = 0.
      K9S(ILON,INIV) = 0.
C
      ENDIF
C
      end do
      end do
C
      do iniv = nivhau,nivbas
      do ilon = 1,nlon
C
      IF (LICE(ILON,INIV) .EQ. 1) THEN
C
      T = TEMP(ILON,INIV)
C
C Z VARIABLES ARE SET TO AVOID DIVISION BY ZERO
C
      ZH2O  = MAX(CK(ILON,INIV, 3),1.E-08*HNM(ILON,INIV))
      ZHCL  = MAX(CK(ILON,INIV,13),1.E-13*HNM(ILON,INIV))
      ZHOCL = MAX(CK(ILON,INIV,15),1.E-13*HNM(ILON,INIV))
      ZHOBR = MAX(CK(ILON,INIV,26),1.E-13*HNM(ILON,INIV))
C
C =====================================================================
C ====================== CLNO3 + HCL ==================================
C ====================== CLNO3 + H2O ==================================
C ====================== N2O5  + H2O ==================================
C ====================== N2O5  + HCL ==================================
C ====================== HOCL  + HCL ==================================
C ====================== BRNO3 + H2O  =================================
C ====================== HOBR  + HCL  =================================
C ====================== HOBR  + HBR  =================================
C ====================== HOCL  + HBR  =================================
C ======================= ON ICE ======================================
C =====================================================================
C GAMMA VALUES
C
      GCLNO3HCLICE=0.3
      GCLNO3H2OICE=0.3
      GN2O5H2OICE=0.02
      GN2O5HCLICE=0.03
      GHOCLHCLICE=0.2
      GBRNO3H2OICE=0.3
      GHOBRHCLICE=0.3
      GHOBRHBRICE=0.1
      GHOCLHBRICE=0.
C
C BECAUSE ICE PARTICLES CAN BE VERY LARGE, GAS DIFFUSION LIMITATION
C IS TAKEN INTO ACCOUNT. THE SECOND ORDER RATE CONSTANTS ARE
C CALCULATED FROM THE EQUATION R=G*PI*R**2*CBAR*N/(1+3G.R/(4.L)), WHERE
C G=GAMMA, R=PARTICLE RADIUS, CBAR=MEAN MOLECULAR SPEED, L=MEAN FREE
C PATH (TURCO ET AL, JGR, 94, 16493, 1989).
C
      RCLNO3HCLICE=
     >     GCLNO3HCLICE*4.56E4*SQRT(T/97.45)*RICE**2
     >     *NICE(ILON,INIV)/
     >     (1.0+3.3E4*GCLNO3HCLICE*RICE*PTOT(ILON,INIV)/T)
     >     /ZHCL
      K2S(ILON,INIV) = RCLNO3HCLICE
C
      RCLNO3H2OICE=
     >     GCLNO3H2OICE*4.56E4*SQRT(T/97.45)*RICE**2
     >     *NICE(ILON,INIV)/
     >     (1.0+3.3E4*GCLNO3H2OICE*RICE*PTOT(ILON,INIV)/T)
     >     /ZH2O
      K1S(ILON,INIV) = RCLNO3H2OICE
C
      RN2O5H2OICE=
     >     GN2O5H2OICE*4.56E4*SQRT(T/108.0)*RICE**2
     >     *NICE(ILON,INIV)/
     >     (1.0+3.3E4*GN2O5H2OICE*RICE*PTOT(ILON,INIV)/T)
     >     /ZH2O
      K3S(ILON,INIV) = RN2O5H2OICE
C
      RN2O5HCLICE=
     >     GN2O5HCLICE*4.56E4*SQRT(T/108.0)*RICE**2
     >     *NICE(ILON,INIV)/
     >     (1.0+3.3E4*GN2O5HCLICE*RICE*PTOT(ILON,INIV)/T)
     >     /ZHCL
      K4S(ILON,INIV) = RN2O5HCLICE
C
      RHOCLHCLICE=
     >     GHOCLHCLICE*4.56E4*SQRT(T/52.5)*RICE**2
     >     *NICE(ILON,INIV)/
     >     (1.0+3.3E4*GHOCLHCLICE*RICE*PTOT(ILON,INIV)/T)
     >     /ZHCL
      K5S(ILON,INIV) = RHOCLHCLICE
C
      RBRNO3H2OICE=
     >     GBRNO3H2OICE*4.56E4*SQRT(T/142.0)*RICE**2
     >     *NICE(ILON,INIV)/
     >     (1.0+3.3E4*GBRNO3H2OICE*RICE*PTOT(ILON,INIV)/T)
     >     /ZH2O
      K6S(ILON,INIV) = RBRNO3H2OICE
C
      RHOBRHCLICE=
     >     GHOBRHCLICE*4.56E4*SQRT(T/96.9)*RICE**2
     >     *NICE(ILON,INIV)/
     >     (1.0+3.3E4*GHOBRHCLICE*RICE*PTOT(ILON,INIV)/T)
     >     /ZHCL
      K7S(ILON,INIV) = RHOBRHCLICE
C
      RHOBRHBRICE=
     >     GHOBRHBRICE*4.56E4*SQRT(T/96.9)*RICE**2
     >     *NICE(ILON,INIV)/
     >     (1.0+3.3E4*GHOBRHBRICE*RICE*PTOT(ILON,INIV)/T)
     >     /ZHOBR
      K8S(ILON,INIV) = RHOBRHBRICE
C
      END IF
C
      end do
      end do
C
      do iniv = nivhau,nivbas
         do ilon = 1,nlon
            K1(ILON,INIV) = K1L(ILON,INIV)
     $                   +(LNAT(ILON,INIV) 
     $                   + LICE(ILON,INIV))*K1S(ILON,INIV)
            K2(ILON,INIV) = K2L(ILON,INIV)
     $                   +(LNAT(ILON,INIV)
     $                   + LICE(ILON,INIV))*K2S(ILON,INIV)
            K3(ILON,INIV) = K3L(ILON,INIV)
     $                   +(LNAT(ILON,INIV)
     $                   + LICE(ILON,INIV))*K3S(ILON,INIV)
            K4(ILON,INIV) = K4L(ILON,INIV)
     $                   +(LNAT(ILON,INIV) 
     $                   + LICE(ILON,INIV))*K4S(ILON,INIV)
            K5(ILON,INIV) = K5L(ILON,INIV)
     $                   +(LNAT(ILON,INIV) 
     $                   + LICE(ILON,INIV))*K5S(ILON,INIV)
            K6(ILON,INIV) = K6L(ILON,INIV)
     $                   +(LNAT(ILON,INIV) 
     $                   + LICE(ILON,INIV))*K6S(ILON,INIV)
            K7(ILON,INIV) = K7L(ILON,INIV)
     $                   +(LNAT(ILON,INIV) 
     $                   + LICE(ILON,INIV))*K7S(ILON,INIV)
            K8(ILON,INIV) = K8L(ILON,INIV)
     $                   +(LNAT(ILON,INIV) 
     $                   + LICE(ILON,INIV))*K8S(ILON,INIV)
            K9(ILON,INIV) = 0.
         end do
      end do
C
      return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE DENSITY(WS,WN,TEMP,DENS)
C
C     DENSITY OF TERNARY SOLUTION IN G/CM3
C     WS ,WN ARE WT FRACTION,
C     FITTED TO 0.05<WS+WN<0.70 WT FRACTION, BUT EXTRAPOLATES WELL
C     185 < T (K)
C
      parameter(nlon = 180, nlat = 91, niv = 137, nivbas = 94,
     $          nivhau = 14)
C
      REAL WS(NLON,NIVBAS),WN(NLON,NIVBAS),TEMP(NLON,NIVBAS),
     $     DENS(NLON,NIVBAS)
      REAL X(22)
      SAVE X
C
      DATA X/
     $     2.393284E-02,-4.359335E-05,7.961181E-08,0.0,-0.198716351,
     $     1.39564574E-03,-2.020633E-06,0.51684706,-3.0539E-03,
     $     4.505475E-06,-0.30119511,1.840408E-03,-2.7221253742E-06,
     $     -0.11331674116,8.47763E-04,-1.22336185E-06,0.3455282,
     $     -2.2111E-03,3.503768245E-06,-0.2315332,1.60074E-03,
     $     -2.5827835E-06/
C
      do iniv = nivhau,nivbas
         do ilon = 1,nlon
            T = TEMP(ILON,INIV)
            W = WS(ILON,INIV) + WN(ILON,INIV)
            WH= 1.0 - W
            V1=X(1)+X(2)*T+X(3)*T**2+X(4)*T**3
            VS=X(5)+X(6)*T+X(7)*T**2+(X(8)+X(9)*T+X(10)*T**2)*W
     $         +(X(11)+X(12)*T+X(13)*T**2)*W*W
            VN=X(14)+X(15)*T+X(16)*T**2+(X(17)+X(18)*T+X(19)*T**2)*W
     $         +(X(20)+X(21)*T+X(22)*T**2)*W*W
            VMCAL=WH/18.0160*V1 + VS*WS(ILON,INIV)/98.080
     $            + VN*WN(ILON,INIV)/63.0160
            DENS(ILON,INIV) = 1.0/VMCAL*0.001
         end do
      end do
C
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine tropo (ian,imois,dt,trajpre,irapp,qj1,hc)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     forcage des especes au dessous du niveau le plus bas de la chimie (nivbas)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      parameter(nlon = 180, nlat = 91, niv = 137, nivbas = 94)
c
      parameter(nbcon = 43, ncm = 15)
      parameter(nrappmax = 12)
c
c     nccmval : nombre d'annees prises en compte pour le forcage
c
      parameter(nccmval = 37)
c
      real    qj1(nlon,nlat,niv,nbcon)
      real    hc(nlon,nlat,niv,ncm)
      real    trajpre(nlon,nlat,niv,nrappmax)
      real    dt, tau, ilim
c
      real      s_noy(nccmval),      s_cly(nccmval)
      real      s_bry(nccmval)
      real      s_n2o(nccmval),      s_ch4(nccmval), s_cfc11(nccmval)
      real    s_cfc12(nccmval),   s_cfc113(nccmval),  s_ccl4(nccmval)
      real  s_ch3ccl3(nccmval),   s_hcfc22(nccmval), s_h1211(nccmval)
      real    s_h1301(nccmval),    s_ch3br(nccmval), s_ch3cl(nccmval)
      real s_hcfc141b(nccmval), s_hcfc142b(nccmval)
      real   s_cfc114(nccmval),   s_cfc115(nccmval)
      real    s_h2402(nccmval),   s_ch2br2(nccmval)

c     valeurs de forcage au niveau nivbas pour noy, cly, bry
c     les valeurs sont donnees pour chaque annee (mois de juin) entre
c     1988 et 2023
c
c                           1988
c                           1989        1990        1991        1992
c                           1993        1994        1995        1996 
c                           1997        1998        1999        2000
c                           2001        2002        2003        2004
c                           2005        2006        2007        2008
c                           2009        2010        2011        2012
c                           2013        2014        2015        2016
c                           2017        2018        2019        2020
c                           2021        2022        2023        2024

      data s_noy /        0.3e-9,
     $                    0.3e-9,     0.3e-9,     0.3e-9,     0.3e-9,
     $                    0.3e-9,     0.3e-9,     0.3e-9,     0.3e-9,
     $                    0.3e-9,     0.3e-9,     0.3e-9,     0.3e-9,
     $                    0.3e-9,     0.3e-9,     0.3e-9,     0.3e-9,
     $                    0.3e-9,     0.3e-9,     0.3e-9,     0.3e-9,
     $                    0.3e-9,     0.3e-9,     0.3e-9,     0.3e-9,
     $                    0.3e-9,     0.3e-9,     0.3e-9,     0.3e-9,
     $                    0.3e-9,     0.3e-9,     0.3e-9,     0.3e-9,
     $                    0.3e-9,     0.3e-9,     0.3e-9,     0.3e-9 /
c
      data s_cly /       40.e-12,
     $                   40.e-12,    40.e-12,    40.e-12,    40.e-12,
     $                   40.e-12,    40.e-12,    40.e-12,    40.e-12,
     $                   40.e-12,    40.e-12,    40.e-12,    40.e-12,
     $                   40.e-12,    40.e-12,    40.e-12,    40.e-12,
     $                   40.e-12,    40.e-12,    40.e-12,    40.e-12,
     $                   40.e-12,    40.e-12,    40.e-12,    40.e-12,
     $                   40.e-12,    40.e-12,    40.e-12,    40.e-12,
     $                   40.e-12,    40.e-12,    40.e-12,    40.e-12,
     $                   40.e-12,    40.e-12,    40.e-12,    40.e-12 /
c
      data s_bry /       0.5e-12,
     $                   0.5e-12,    0.5e-12,    0.5e-12,    0.5e-12,
     $                   0.5e-12,    0.5e-12,    0.5e-12,    0.5e-12,
     $                   0.5e-12,    0.5e-12,    0.5e-12,    0.5e-12,
     $                   0.5e-12,    0.5e-12,    0.5e-12,    0.5e-12,
     $                   0.5e-12,    0.5e-12,    0.5e-12,    0.5e-12,
     $                   0.5e-12,    0.5e-12,    0.5e-12,    0.5e-12,
     $                   0.5e-12,    0.5e-12,    0.5e-12,    0.5e-12,
     $                   0.5e-12,    0.5e-12,    0.5e-12,    0.5e-12,
     $                   0.5e-12,    0.5e-12,    0.5e-12,    0.5e-12 /
c
c     valeurs de forcage au niveau nivbas pour les gaz sources d'apres ccmval
c     les valeurs sont donnees pour chaque annee (mois de juin) entre
c     1988 et 2022
c
c                           1988
c                           1989        1990        1991        1992
c                           1993        1994        1995        1996 
c                           1997        1998        1999        2000
c                           2001        2002        2003        2004
c                           2005        2006        2007        2008
c                           2009        2010        2011        2012
c                           2013        2014        2015        2016
c                           2017        2018        2019        2020
c                           2021        2022        2023        2024
c
!     n2o : combined hats n2o data file june data
!           noaa (https://gml.noaa.gov/hats/combined/N2O.html)

      data s_n2o /       306.2e-9,
     $                   307.6e-9,   308.1e-9,   309.0e-9,   309.7e-9,
     $                   309.9e-9,   310.7e-9,   310.9e-9,   311.8e-9,
     $                   312.5e-9,   313.3e-9,   314.5e-9,   315.6e-9,
     $                   316.1e-9,   316.9e-9,   317.6e-9,   318.1e-9,
     $                   318.8e-9,   319.9e-9,   320.4e-9,   321.5e-9,
     $                   322.0e-9,   323.0e-9,   323.9e-9,   324.8e-9,
     $                   325.8e-9,   326.8e-9,   327.9e-9,   328.8e-9,
     $                   329.6e-9,   330.9e-9,   332.2e-9,   332.8e-9,
     $                   334.3e-9,   334.6e-9,   335.5e-9,   336.7e-9 /

!     ch4 : gml.noaa.gov from 2022 onwards (annual mean)

      data s_ch4 /       1.693e-6,
     $                   1.704e-6,   1.714e-6,   1.725e-6,   1.735e-6,
     $                   1.736e-6,   1.742e-6,   1.749e-6,   1.751e-6,
     $                   1.754e-6,   1.765e-6,   1.772e-6,   1.773e-6,
     $                   1.771e-6,   1.773e-6,   1.777e-6,   1.777e-6,
     $                   1.774e-6,   1.775e-6,   1.781e-6,   1.787e-6,
     $                   1.793e-6,   1.799e-6,   1.803e-6,   1.808e-6,
     $                   1.813e-6,   1.822e-6,   1.833e-6,   1.842e-6,
     $                   1.849e-6,   1.856e-6,   1.863e-6,   1.879e-6,
     $                   1.892e-6,   1.895e-6,   1.911e-6,   1.956e-6/

!     ch3cl : gml.noaa.gov (flask data) from 2022 onwards

      data s_ch3cl /   549.72e-12,
     $                 549.78e-12, 549.83e-12, 549.87e-12, 549.90e-12,
     $                 549.93e-12, 549.95e-12, 549.97e-12, 549.99e-12,
     $                 550.00e-12, 550.00e-12, 550.01e-12, 550.01e-12,
     $                 550.01e-12, 550.01e-12, 550.01e-12, 550.01e-12,
     $                 550.01e-12, 550.01e-12, 550.01e-12, 550.01e-12,
     $                 550.01e-12, 550.01e-12, 550.01e-12, 550.01e-12,
     $                 550.01e-12, 550.01e-12, 550.01e-12, 530.e-12,
     $                 520.e-12,   510.e-12,   510.e-12,   510.e-12,
     $                 510.e-12,   560.e-12,   560.e-12,   560.e-12/

!     cfc-11 : gml.noaa.gov (flask data) from 2022 onwards

      data s_cfc11 /   245.21e-12,
     $                 253.06e-12, 260.45e-12, 265.07e-12, 267.88e-12,
     $                 269.50e-12, 269.50e-12, 268.94e-12, 268.16e-12,
     $                 266.84e-12, 265.08e-12, 263.36e-12, 259.0e-12,
     $                 257.3e-12, 255.7e-12, 253.4e-12, 252.0e-12,
     $                 249.4e-12, 247.4e-12, 245.0e-12, 243.2e-12,
     $                 241.3e-12, 239.4e-12, 237.2e-12, 235.1e-12,
     $                 233.1e-12, 231.8e-12, 230.7e-12, 229.3e-12,
     $                 228.6e-12, 227.9e-12, 226.4e-12, 224.0e-12,
     $                 221.9e-12, 219.7e-12, 217.5e-12, 215.5e-12 /

!     cfc-12 : gml.noaa.gov (flask data) from 2022 onwards

      data s_cfc12 /   445.e-12,
     $                 466.78e-12, 481.96e-12, 494.85e-12, 505.25e-12,
     $                 513.37e-12, 519.82e-12, 525.26e-12, 529.34e-12,
     $                 532.53e-12, 535.06e-12, 537.17e-12, 538.96e-12,
     $                 540.00e-12, 540.03e-12, 539.87e-12, 539.34e-12,
     $                 537.38e-12, 534.56e-12, 531.64e-12, 528.59e-12,
     $                 525.41e-12, 522.10e-12, 518.67e-12, 515.42e-12,
     $                 511.78e-12, 508.04e-12, 504.23e-12, 512.e-12,
     $                 509.e-12,   507.e-12,   501.e-12,   497.e-12,
     $                 493.e-12,   490.e-12,   487.e-12,   484.e-12/

!     cfc-114

      data s_cfc114 /   14.57e-12,
     $                  15.10e-12,  15.57e-12,  16.00e-12, 16.36e-12,
     $                  16.64e-12,  16.83e-12,  16.97e-12, 17.05e-12,
     $                  17.08e-12,  17.07e-12,  17.05e-12, 17.01e-12,
     $                  16.98e-12,  16.94e-12,  16.89e-12, 16.85e-12,
     $                  16.81e-12,  16.76e-12,  16.70e-12, 16.65e-12,
     $                  16.60e-12,  16.54e-12,  16.49e-12, 16.44e-12,
     $                  16.38e-12,  16.33e-12,  16.28e-12, 16.22e-12,
     $                  16.e-12,    16.e-12,    16.e-12,   16.e-12  ,
     $                  16.e-12,    16.e-12,    16.e-12,   16.e-12 /

!     cfc-115

      data s_cfc115 /    4.62e-12,
     $                   5.08e-12,   5.56e-12,   6.04e-12,  6.51e-12,
     $                   6.98e-12,   7.42e-12,   7.80e-12,  8.11e-12,
     $                   8.35e-12,   8.5e-12,    8.5e-12,   8.5e-12,
     $                   8.5e-12,    8.5e-12,    8.5e-12,   8.5e-12,
     $                   8.5e-12,    8.5e-12,    8.5e-12,   8.5e-12,
     $                   8.5e-12,    8.5e-12,    8.5e-12,   8.5e-12,
     $                   8.5e-12,    8.5e-12,    8.5e-12,   8.5e-12,
     $                   8.5e-12,    8.5e-12,    8.5e-12,   8.5e-12,
     $                   8.5e-12,    8.5e-12,    8.5e-12,   8.5e-12 /

!     cfc-113: : gml.noaa.gov (flask data) from 2022 onwards

      data s_cfc113 /   58.23e-12,
     $                  64.74e-12,  70.54e-12,  76.20e-12,  80.36e-12,
     $                  82.40e-12,  83.50e-12,  83.83e-12,  83.68e-12,
     $                  83.44e-12,  83.00e-12,  82.49e-12,  82.05e-12,
     $                  81.57e-12,  80.89e-12,  80.08e-12,  79.28e-12,
     $                  78.44e-12,  77.56e-12,  76.68e-12,  75.81e-12,
     $                  74.95e-12,  74.09e-12,  73.23e-12,  72.45e-12,
     $                  71.61e-12,  70.77e-12,  69.94e-12,  71.3e-12,
     $                  70.8e-12,   70.3e-12,   71.0e-12,   70.8e-12,
     $                  68.4e-12,   67.8e-12,   67.2e-12,   66.6e-12 /
 
!     ccl4: combined hats ccl4 data file (2016 onwards, june data)
!           noaa (ftp://ftp.cmdl.noaa.gov/hats/solvents/CCl4/combined/HATS_global_CCl4.txt

      data s_ccl4 /    103.52e-12,
     $                 104.50e-12, 105.66e-12, 105.69e-12, 104.87e-12,
     $                 104.32e-12, 103.53e-12, 102.78e-12, 101.89e-12,
     $                 100.86e-12, 100.15e-12,  99.32e-12,  98.36e-12,
     $                  97.33e-12,  96.28e-12,  95.34e-12,  94.39e-12,
     $                  93.39e-12,  92.25e-12,  90.89e-12,  89.34e-12,
     $                  87.60e-12,  85.67e-12,  83.56e-12,  81.47e-12,
     $                  79.04e-12,  76.44e-12,  73.69e-12,  81.3e-12,
     $                  80.3e-12,   79.3e-12,   78.3e-12,   77.0e-12,
     $                  76.3e-12,   75.1e-12,   73.9e-12,   72.7e-12/

!     ch3ccl3 : in situ global file (2016 onwards, june data)
!               noaa (ftp://ftp.cmdl.noaa.gov/hats/solvents/CH3CCl3/insituGCs/CATS/global/insitu_global_MC.txt)

      data s_ch3ccl3 / 122.04e-12,
     $                 124.51e-12, 128.30e-12, 130.89e-12, 130.42e-12,
     $                 125.42e-12, 115.07e-12, 103.30e-12,  90.27e-12,
     $                  76.97e-12,  64.82e-12,  54.18e-12,  45.23e-12,
     $                  37.68e-12,  31.40e-12,  26.15e-12,  21.81e-12,
     $                  18.31e-12,  15.44e-12,  13.09e-12,  11.17e-12,
     $                   9.59e-12,   8.22e-12,   7.03e-12,   6.13e-12,
     $                   5.32e-12,   4.65e-12,   3.96e-12,   3.0e-12,
     $                   2.6e-12,    2.3e-12,    2.0e-12,    1.6e-12,
     $                   1.3e-12,    1.0e-12,    1.0e-12,    1.0e-12 /

!     hcfc-22 : niwot ridge flask data (2016 onwards, june data)
!               noaa (ftp://ftp.cmdl.noaa.gov/hats/hcfcs/hcfc22/flasks/HCFC22_GCMS_flask.txt)

      data s_hcfc22 /   78.53e-12,
     $                  84.62e-12,  90.76e-12,  97.04e-12, 102.39e-12,
     $                 106.82e-12, 110.63e-12, 116.25e-12, 122.16e-12,
     $                 126.54e-12, 132.02e-12, 137.61e-12, 142.13e-12,
     $                 147.77e-12, 153.17e-12, 157.68e-12, 162.50e-12,
     $                 168.01e-12, 174.21e-12, 180.86e-12, 187.88e-12,
     $                 195.20e-12, 202.67e-12, 210.17e-12, 217.07e-12,
     $                 224.44e-12, 231.44e-12, 237.78e-12, 249.e-12,
     $                 251.e-12,   253.e-12,   253.e-12,   253.e-12,
     $                 255.e-12,   255.e-12,   255.e-12,   255.e-12 /

!     hcfc-141b : niwot ridge flask data (2016 onwards, june data)
!               noaa (ftp://ftp.cmdl.noaa.gov/hats/hcfcs/hcfc141b/HCFC141B_GCMS_flask.txt)

      data s_hcfc141b/         0.,
     $                         0.,         0.,         0.,         0.,
     $                   0.63e-12,   2.00e-12,   3.69e-12,   5.64e-12,
     $                   7.47e-12,   9.32e-12,  11.17e-12,  12.78e-12,
     $                  14.26e-12,  15.64e-12,  16.65e-12,  17.22e-12,
     $                  17.71e-12,  18.25e-12,  18.84e-12,  19.46e-12,
     $                  20.08e-12,  20.70e-12,  21.31e-12,  21.84e-12,
     $                  22.39e-12,  22.90e-12,  23.37e-12,  26.0e-12,
     $                  25.7e-12,   25.4e-12,   25.4e-12,   25.4e-12,
     $                  26.e-12,    25.e-12,    25.e-12,    25.e-12/

!     hcfc-142b : northern hemisphere monthly mean (2016 onwards, june data)
!               noaa (ftp://ftp.cmdl.noaa.gov/hats/hcfcs/hcfc142b/insituGCs/CATS/global/insitu_global_HCFC142b.txt)

      data s_hcfc142b/   1.00e-12,
     $                   1.26e-12,   1.72e-12,   2.47e-12,   3.45e-12,
     $                   4.54e-12,   5.72e-12,   6.87e-12,   7.99e-12,
     $                   8.99e-12,   9.99e-12,  11.09e-12,  12.08e-12,
     $                  13.02e-12,  13.78e-12,  14.38e-12,  14.97e-12,
     $                  15.38e-12,  15.64e-12,  15.88e-12,  16.08e-12,
     $                  16.25e-12,  16.38e-12,  16.45e-12,  16.48e-12,
     $                  16.47e-12,  16.40e-12,  16.27e-12,  22.5e-12,
     $                  22.6e-12,   22.7e-12,   22.8e-12,   22.9e-12,
     $                  22.e-12,    21.e-12,    21.e-12,    21.e-12 /
 
!     ch3br : gml.noaa.gov (flask data) from 2022 onwards

      data s_ch3br /     8.92e-12,
     $                   9.01e-12,   9.09e-12,   9.18e-12,   9.28e-12,
     $                   9.37e-12,   9.47e-12,   9.22e-12,   9.07e-12,
     $                   9.19e-12,   9.26e-12,   9.01e-12,   8.60e-12,
     $                   8.29e-12,   8.18e-12,   8.05e-12,   7.90e-12,
     $                   7.70e-12,   7.58e-12,   7.54e-12,   7.53e-12,
     $                   7.52e-12,   7.52e-12,   7.52e-12,   7.52e-12,
     $                   7.52e-12,   7.52e-12,   7.27e-12,   7.8e-12,
     $                   7.6e-12,    7.4e-12,    7.2e-12,    7.4e-12,
     $                   7.4e-12,    7.0e-12,    7.0e-12,    7.0e-12 /

!     h-1211 : gml.noaa.gov (flask data) from 2022 onwards

      data s_h1211 /     1.88e-12,
     $                   2.18e-12,   2.46e-12,   2.72e-12,   2.93e-12,
     $                   3.10e-12,   3.25e-12,   3.40e-12,   3.56e-12,
     $                   3.71e-12,   3.84e-12,   3.96e-12,   4.05e-12,
     $                   4.11e-12,   4.17e-12,   4.21e-12,   4.23e-12,
     $                   4.24e-12,   4.24e-12,   4.22e-12,   4.19e-12,
     $                   4.14e-12,   4.09e-12,   4.02e-12,   3.95e-12,
     $                   3.86e-12,   3.77e-12,   3.68e-12,   3.5e-12,
     $                   3.4e-12,    3.3e-12,    3.2e-12,    3.1e-12,
     $                   3.05e-12,   2.95e-12,   2.85e-12,   2.75e-12 /

!     h-1301

      data s_h1301 /     1.37e-12,
     $                   1.57e-12,   1.76e-12,   1.92e-12,   2.02e-12,
     $                   2.12e-12,   2.25e-12,   2.39e-12,   2.49e-12,
     $                   2.56e-12,   2.62e-12,   2.68e-12,   2.73e-12,
     $                   2.78e-12,   2.83e-12,   2.87e-12,   2.91e-12,
     $                   2.94e-12,   2.97e-12,   3.00e-12,   3.02e-12,
     $                   3.04e-12,   3.06e-12,   3.07e-12,   3.08e-12,
     $                   3.09e-12,   3.10e-12,   3.10e-12,   3.30e-12,
     $                   3.30e-12,   3.30e-12,   3.30e-12,   3.30e-12,
     $                   3.30e-12,   3.3e-12,    3.3e-12,    3.3e-12/

!     h-2402 : gml.noaa.gov (flask data) from 2022 onwards

      data s_h2402 /     0.28e-12,
     $                   0.30e-12,   0.33e-12,   0.35e-12,   0.37e-12,
     $                   0.39e-12,   0.41e-12,   0.43e-12,   0.43e-12,
     $                   0.43e-12,   0.42e-12,   0.42e-12,   0.41e-12,
     $                   0.40e-12,   0.39e-12,   0.38e-12,   0.37e-12,
     $                   0.36e-12,   0.35e-12,   0.34e-12,   0.33e-12,
     $                   0.32e-12,   0.31e-12,   0.30e-12,   0.28e-12,
     $                   0.27e-12,   0.26e-12,   0.25e-12,   0.42e-12,
     $                   0.41e-12,   0.41e-12,   0.41e-12,   0.41e-12,
     $                   0.40e-12,   0.39e-12,   0.39e-12,   0.39e-12 /

!     ch2br2

      data s_ch2br2/      2.5e-12,
     $                    2.5e-12,    2.5e-12,    2.5e-12,     2.5e-12,
     $                    2.5e-12,    2.5e-12,    2.5e-12,     2.5e-12,
     $                    2.5e-12,    2.5e-12,    2.5e-12,     2.5e-12,
     $                    2.5e-12,    2.5e-12,    2.5e-12,     2.5e-12,
     $                    2.5e-12,    2.5e-12,    2.5e-12,     2.5e-12,
     $                    2.5e-12,    2.5e-12,    2.5e-12,     2.5e-12,
     $                    2.5e-12,    2.5e-12,    2.5e-12,     2.5e-12,
     $                    2.5e-12,    2.5e-12,    2.5e-12,     2.5e-12,
     $                    2.5e-12,    2.5e-12,    2.5e-12,     2.5e-12 /

c
c     calcul des valeurs de forcage pour la date courante
c
      yearccmval = real(ian) + real(imois - 1)/12.
c     print*,'yearccmval = ',yearccmval
c
      icount = 0
      do iyear = 1988,2024
         icount = icount + 1
         if (real(iyear) + 0.5 .gt. yearccmval) then
            cinf = (yearccmval - real(iyear - 1) - 0.5)/1.
            csup = 1. - cinf
            v_noy       = cinf*s_noy(icount)
     $                  + csup*s_noy(icount - 1)     
            v_cly       = cinf*s_cly(icount)
     $                  + csup*s_cly(icount - 1)     
            v_bry       = cinf*s_bry(icount)
     $                  + csup*s_bry(icount - 1)     
            v_n2o       = cinf*s_n2o(icount)
     $                  + csup*s_n2o(icount - 1)     
            v_ch4       = cinf*s_ch4(icount)
     $                  + csup*s_ch4(icount - 1)     
            v_ch3cl     = cinf*s_ch3cl(icount)
     $                  + csup*s_ch3cl(icount - 1)     
            v_cfc11     = cinf*s_cfc11(icount)
     $                  + csup*s_cfc11(icount - 1)     
            v_cfc12     = cinf*s_cfc12(icount)
     $                  + csup*s_cfc12(icount - 1)     
            v_cfc113    = cinf*s_cfc113(icount)
     $                  + csup*s_cfc113(icount - 1)     
            v_cfc114    = cinf*s_cfc114(icount)
     $                  + csup*s_cfc114(icount - 1)     
            v_cfc115    = cinf*s_cfc115(icount)
     $                  + csup*s_cfc115(icount - 1)     
            v_ccl4      = cinf*s_ccl4(icount)
     $                  + csup*s_ccl4(icount - 1)     
            v_ch3ccl3   = cinf*s_ch3ccl3(icount)
     $                  + csup*s_ch3ccl3(icount - 1)     
            v_hcfc22    = cinf*s_hcfc22(icount)
     $                  + csup*s_hcfc22(icount - 1)     
            v_hcfc141b  = cinf*s_hcfc141b(icount)
     $                  + csup*s_hcfc141b(icount - 1)     
            v_hcfc142b  = cinf*s_hcfc142b(icount)
     $                  + csup*s_hcfc142b(icount - 1)     
            v_ch3br     = cinf*s_ch3br(icount)
     $                  + csup*s_ch3br(icount - 1)     
            v_h1211     = cinf*s_h1211(icount)
     $                  + csup*s_h1211(icount - 1)     
            v_h1301     = cinf*s_h1301(icount)
     $                  + csup*s_h1301(icount - 1)     
            v_h2402     = cinf*s_h2402(icount)
     $                  + csup*s_h2402(icount - 1)     
            v_ch2br2    = cinf*s_ch2br2(icount)
     $                  + csup*s_ch2br2(icount - 1)     
            go to 300
         end if
      end do
 300  continue
c
c     on ajoute cfc-114 (2 atomes de chlore) et 0.5*cfc-115
c     (1 atome de chlore) a cfc-12
c
      v_cfc12 = v_cfc12 + v_cfc114 + 0.5*v_cfc115
c
c     on ajoute hcfc-141b a ch3ccl3 car de durees de vie proche
c     avec un facteur 2/3 car ch3ccl3 a trois atomes de chlore
c     et hcfc-141b en a deux
c
      v_ch3ccl3 = v_ch3ccl3 + v_hcfc141b*2./3.
c
c     on ajoute hcfc-142b a hcfc22 car de durees de vie proche
c
      v_hcfc22 = v_hcfc22 + v_hcfc142b
c
c     on ajoute h-2402 a h-1301
c     avec un facteur 2 car h-2402 a deux atomes de brome
c     et h-1301 en a un
c
      v_h1301 = v_h1301 + v_h2402*2.
c
c     print*, 'n2o        : ', v_n2o     
c     print*, 'ch4        : ', v_ch4    
c     print*, 'noy        : ', v_noy    
c     print*, 'cly        : ', v_cly    
c     print*, 'bry        : ', v_bry    
c     print*, 'ch3cl      : ', v_ch3cl  
c     print*, 'cfc-11     : ', v_cfc11
c     print*, 'cfc-12*    : ', v_cfc12
c     print*, 'cfc-113    : ', v_cfc113
c     print*, 'ccl4       : ', v_ccl4
c     print*, 'ch3ccl3*   : ', v_ch3ccl3
c     print*, 'hcfc-22*   : ', v_hcfc22
c     print*, 'ch3br      : ', v_ch3br
c     print*, 'h-1211     : ', v_h1211
c     print*, 'h-1301*    : ', v_h1301
c     print*, 'ch2br2     : ', v_ch2br2
c     print*
c
c     forcage a partir de nivbas + 1
c
      do iniv = nivbas + 1, niv
         do ilat = 1,nlat
            do ilon = 1,nlon
               qj1(ilon,ilat,iniv, 1) = v_n2o
               qj1(ilon,ilat,iniv, 2) = v_ch4
               qj1(ilon,ilat,iniv, 4) = v_noy
               qj1(ilon,ilat,iniv, 7) = v_cly
               qj1(ilon,ilat,iniv,31) = v_cfc11
               qj1(ilon,ilat,iniv,32) = v_cfc12
               qj1(ilon,ilat,iniv,33) = v_cfc113
               qj1(ilon,ilat,iniv,34) = v_ccl4
               qj1(ilon,ilat,iniv,35) = v_ch3ccl3
               qj1(ilon,ilat,iniv,36) = v_ch3cl
               qj1(ilon,ilat,iniv,37) = v_hcfc22
               qj1(ilon,ilat,iniv,38) = v_ch3br
               qj1(ilon,ilat,iniv,39) = v_h1211
               qj1(ilon,ilat,iniv,40) = v_h1301
               qj1(ilon,ilat,iniv,41) = v_bry
               qj1(ilon,ilat,iniv,42) = v_ch2br2
            end do
         end do 
      end do
c
c     cas particulier de ox
c
c     tau : temps de relaxation en jours
c     dtt : pas de temps en jours
c     ilim: valeur de forcage a t infini
c
      tau = 1.
      dtt = dt/(24.*60.*60.)
c
      do iniv = nivbas + 1,niv
         do ilat = 1,nlat
              do ilon = 1,nlon
                 if (trajpre(ilon,ilat,iniv,irapp) .gt. 370.) then
                   ilim = -2.3825e-8
     $                  *log(trajpre(ilon,ilat,iniv,irapp)) + 1.9457e-7
                 else
                    ilim = qj1(ilon,ilat,iniv,8)
                 end if
                 qj1(ilon,ilat,iniv,8)= 
     $          max(ilim+(qj1(ilon,ilat,iniv,8)-ilim)*exp(-dtt/tau), 0.)
             end do
          end do
      end do
c
c     en dessous de nivbas:
c
c     o3        = ox
c     ox passif = ox
c
      do iniv = nivbas,niv
         do ilat = 1,nlat
            do ilon = 1,nlon
               hc(ilon,ilat,iniv,5)   = qj1(ilon,ilat,iniv,8)
               qj1(ilon,ilat,iniv,11) = qj1(ilon,ilat,iniv,8)
            end do
         end do
      end do
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine co_tropo (imois,dt,qj1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     forcage de co au dessous du niveau le plus bas de la chimie (nivbas)
c     par les mesures de mopitt
c     (jean-luc attie, laboratoire d'aerologie)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      parameter(nlon = 180, nlat = 91, niv = 137, nivbas = 94)
      parameter(nbcon = 43)
c
      parameter(nmois = 12)
c
      real    qj1(nlon,nlat,niv,nbcon)
      real    co(nlon,nlat,nmois)
      real    tau, dt, dtt, ilim
c
      character*35 head
      logical firstcall
c
      data    firstcall /.true./
      save    firstcall
c
c     lecture du fichier co, 12 mois
c
      if (firstcall) then
         firstcall = .false. 
         print*
         print*,'ouverture du fichier co mopitt'
         open(33, form = 'formatted', file = 'mopitt_corrected.txt')
         do i = 1,nmois
            read(33,'(a35)') head
c           print*, head
            do ilat = 1,nlat
               read(33,'(9e13.5)') (co(ilon,ilat,i), ilon = 1,nlon)
            end do
         end do
         print*,'lecture du fichier co mopitt: ok'
         print*
         close(33)
      end if
c
c     tau : temps de relaxation en heures
c     dtt : pas de temps en heures
c     ilim: valeur de forcage a t infini
c
      tau = 24.
      dtt = dt/(60.*60.)
c
      do iniv = nivbas + 1,niv
         do ilat = 1,nlat
            do ilon = 1,nlon
               ilim = co(ilon,ilat,imois)
               qj1(ilon,ilat,iniv,9) =
     $         max(ilim+(qj1(ilon,ilat,iniv,9)-ilim)*exp(-dtt/tau),
     $             1.e-30)
            end do
         end do
      end do
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine h2so4_2d(ian,imois,ijour,iheure,imin,qj1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     lecture de la distribution de h2so4 total a partir
c     d'une simulation 2d avec chimie du soufre
c     et microphysique (Slimane Bekki).
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      parameter(nlon = 180, nlat = 91, niv = 137)
c
      parameter(nbcon = 43)
      parameter(nyear_f = 1989, nyear_l = 2005, nmonth = 12)
c
      real    h2so4(nlat,niv,nyear_f:nyear_l,nmonth)
      real    qj1(nlon,nlat,niv,nbcon)
c
      integer yearinit, monthinit, ian_h2so4, imois_h2so4
c
      logical firstcall
c
      data firstcall /.true./
      save firstcall
c
      if (firstcall) then
c    
c     lecture du fichier 2d
c
         open(52, form='formatted', file='h2so4.txt')
         print*,'lecture du fichier h2so4...'
         do iyear = nyear_f,nyear_l
            do imonth = 1,nmonth
               read(52,'(i4,1x,i4)') yearinit, monthinit
c              write(6,'(i4,1x,i4)') yearinit, monthinit
               read(52,'(1x,6e13.5)')((h2so4(ilat,iniv,iyear,imonth),
     $                                                 ilat = 1,nlat),
     $                                                 iniv = 1,niv)
            end do
         end do
         print*,'lecture du fichier h2so4: ok'
         close(52)
         firstcall = .false.
      end if
c
c     forcage de h2so4 au 1er de chaque mois
c     pour l'annee courante
c
      ian_h2so4   = min(ian,nyear_l) 
      imois_h2so4 = imois 
c
      if (ijour .eq. 1 .and. iheure .eq. 12 
     $    .and. imin .eq. 0) then
         print*
         print*,'initialisation de h2so4:'
         print*,'ian_h2so4   = ',ian_h2so4
         print*,'imois_h2so4 = ',imois_h2so4
c
         do iniv = 1, niv
            do ilat = 1,nlat
               do ilon = 1,nlon
                  qj1(ilon,ilat,iniv,12) =
     $            h2so4(ilat,iniv,ian_h2so4,imois_h2so4)
               end do
            end do 
         end do
      end if
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ecmwf(ian, day, ntime, nintecmwf)
c
      parameter(nlon = 180, nlat = 91, niv = 137)
c
      real day
c
      integer idate(18)
c
      logical entered
c
      character filename*14
c
      common/forcm1/um1(nlon,nlat,niv),vm1(nlon,nlat,niv)
     +             ,wm1(nlon,nlat,niv),tm1(nlon,nlat,niv)
     +             ,hm1(nlon,nlat,niv),pm1(nlon,nlat),daym1
c
      common/forcp1/up1(nlon,nlat,niv),vp1(nlon,nlat,niv)
     +             ,wp1(nlon,nlat,niv),tp1(nlon,nlat,niv)
     +             ,hp1(nlon,nlat,niv),pp1(nlon,nlat),dayp1
c
      save entered
      data entered / .false. /
c
 103  format(1x,a13,i4,i3,i3,i4,a1,i3,a2)
 104  format('ecmwf_',i8)
c
      iunite = 14
c
      call jultogreg(ian,imois,ijour,iheure,imin,day)
      print*
      write(6,103)'date modele: ',
     $             ian,imois,ijour,iheure,'h',imin,'mn'
c
      if( .not. entered ) then
c
c        cas particulier du premier pas de temps: lecture de
c        deux echeances consecutives
c 
         entered = .true.
         print*,
     $   '1er pas de temps, lecture de deux echeances consecutives'
c
c        ouverture du fichier ecmwf
c
         ifilename = ian*10000 + imois*100 + ijour
         write(filename,104) ifilename
         print*,'ouverture du fichier: ',filename
         open(iunite, form = 'unformatted', file = filename)
c
c        saut des echeances anterieures a la date modele
c
         do iday = 1,iheure/(nintecmwf/ntime)
            read(iunite) idate
            read(iunite)
            do iniv = 1,niv
               do i = 1,10
                  read(iunite)
               end do
            end do
         end do
c
c        lecture du forcage ecmwf au temps t0
c        avec verification de la date
c
         read(iunite)  idate
         ianecmwf    = idate(10)
         imoisecmwf  = idate(11)
         ijourecmwf  = idate(12)
         iheureecmwf = idate(13) + idate(16)
         iminecmwf   = 0
c
         call gregtojul(ianecmwf,imoisecmwf,ijourecmwf,iheureecmwf,
     $                  iminecmwf,dayecmwf)
c
         if (dayecmwf .eq. day) then
            daym1 = dayecmwf
            print*,'lecture du forcage ecmwf au temps t0'
            write(6,103)'date ecmwf : ',ianecmwf,imoisecmwf,
     $                  ijourecmwf,iheureecmwf,'h',iminecmwf,'mn'
            print*
            read(iunite) pm1
c
            write(6,*)(idate(i),i = 5,16)
            do iniv = 1,niv
               read(iunite) idate
               write(6,*)(idate(i),i = 5,16)
               read(iunite)((tm1(ilon,ilat,iniv),
     $                           ilon = 1,nlon),ilat = 1,nlat) 
               read(iunite) idate
               write(6,*)(idate(i),i = 5,16)
               read(iunite)((wm1(ilon,ilat,iniv),
     $                           ilon = 1,nlon),ilat = 1,nlat) 
               read(iunite) idate
               write(6,*)(idate(i),i = 5,16)
               read(iunite)((hm1(ilon,ilat,iniv),
     $                           ilon = 1,nlon),ilat = 1,nlat) 
               do ilat = 1,nlat
                  do ilon = 1,nlon
                     hm1(ilon,ilat,iniv) = 
     $               max(hm1(ilon,ilat,iniv),1.e-30*18./28.97)
                  end do
               end do
               read(iunite) idate
               write(6,*)(idate(i),i = 5,16)
               read(iunite)((um1(ilon,ilat,iniv),
     $                           ilon = 1,nlon),ilat = 1,nlat) 
               read(iunite) idate
               write(6,*)(idate(i),i = 5,16)
               read(iunite)((vm1(ilon,ilat,iniv),
     $                           ilon = 1,nlon),ilat = 1,nlat) 
            end do
         else
            print*,'la date ecmwf ne correspond pas, arret du programme'
            stop
         end if
c
c        lecture du forcage ecmwf au temps t0 + delta t
c        avec verification de la date
c
         read(iunite)  idate
         ianecmwf    = idate(10)
         imoisecmwf  = idate(11)
         ijourecmwf  = idate(12)
         iheureecmwf = idate(13) + idate(16)
         iminecmwf   = 0
c
         call gregtojul(ianecmwf,imoisecmwf,ijourecmwf,iheureecmwf,
     $                  iminecmwf,dayecmwf)
c
         if (dayecmwf - day .eq. real(nintecmwf)
     $                          /real(24*ntime)) then
            dayp1 = dayecmwf
            print*
            print*,'lecture du forcage ecmwf au temps t0 + delta t'
            write(6,103)'date ecmwf : ',ianecmwf,imoisecmwf,
     $                  ijourecmwf,iheureecmwf,'h',iminecmwf,'mn'
            print*
            read(iunite) pp1
            write(6,*)(idate(i),i = 5,16)
            do iniv = 1,niv
               read(iunite) idate
               write(6,*)(idate(i),i = 5,16)
               read(iunite)((tp1(ilon,ilat,iniv),
     $                           ilon = 1,nlon),ilat = 1,nlat)
               read(iunite) idate
               write(6,*)(idate(i),i = 5,16)
               read(iunite)((wp1(ilon,ilat,iniv),
     $                           ilon = 1,nlon),ilat = 1,nlat)
               read(iunite) idate
               write(6,*)(idate(i),i = 5,16)
               read(iunite)((hp1(ilon,ilat,iniv),
     $                           ilon = 1,nlon),ilat = 1,nlat)
               do ilat = 1,nlat
                  do ilon = 1,nlon
                     hp1(ilon,ilat,iniv) = 
     $               max(hp1(ilon,ilat,iniv),1.e-30*18./28.97)
                  end do
               end do
               read(iunite) idate
               write(6,*)(idate(i),i = 5,16)
               read(iunite)((up1(ilon,ilat,iniv),
     $                           ilon = 1,nlon),ilat = 1,nlat)
               read(iunite) idate
               write(6,*)(idate(i),i = 5,16)
               read(iunite)((vp1(ilon,ilat,iniv),
     $                           ilon = 1,nlon),ilat = 1,nlat)
            end do
            goto 100 
         else
            print*,'la date ecmwf ne correspond pas, arret du programme'
            stop
         end if
c
      else 
c
c        cas general: lecture du forcage ecmwf au temps t + delta t
c
         do iniv = 1,niv
            do ilat = 1,nlat
               do ilon = 1,nlon
                  tm1(ilon,ilat,iniv) = tp1(ilon,ilat,iniv)
                  um1(ilon,ilat,iniv) = up1(ilon,ilat,iniv)
                  vm1(ilon,ilat,iniv) = vp1(ilon,ilat,iniv)
                  wm1(ilon,ilat,iniv) = wp1(ilon,ilat,iniv)
                  hm1(ilon,ilat,iniv) = hp1(ilon,ilat,iniv)
               end do
            end do
         end do
         do ilat = 1,nlat
            do ilon = 1,nlon
               pm1(ilon,ilat) = pp1(ilon,ilat)
            end do
         end do
c
         day   = dayp1
         daym1 = day
         dayp1 = day + real(nintecmwf)/real(24*ntime)
         ianp1 = ian
         call jultogreg(ianp1,imoisp1,ijourp1,iheurep1,
     $                  iminp1,dayp1)
         ifilename = ianp1*10000 + imoisp1*100 + ijourp1
         write(filename,104) ifilename
c
c        ouverture du fichier ecmwf
c
         print*,'ouverture du fichier: ',filename
         open(iunite, form = 'unformatted', file = filename)
c
c        saut des echeances anterieures a la date modele
c
         do iday = 1,iheurep1/(nintecmwf/ntime)
            read(iunite) idate
            read(iunite)
            do iniv = 1,niv
               do i = 1,10
                  read(iunite)
               end do
            end do
         end do
c
c        lecture du forcage ecmwf au temps t + delta t
c        avec verification de la date
c
         read(iunite)  idate
         ianecmwf    = idate(10)
         imoisecmwf  = idate(11)
         ijourecmwf  = idate(12)
         iheureecmwf = idate(13) + idate(16)
         iminecmwf   = 0
c
         call gregtojul(ianecmwf,imoisecmwf,ijourecmwf,iheureecmwf,
     $                  iminecmwf,dayecmwf)
c 
         if (dayecmwf .eq. dayp1) then
            print*,'lecture du forcage ecmwf au temps t + delta t'
            write(6,103)'date ecmwf : ',ianecmwf,imoisecmwf,
     $                  ijourecmwf,iheureecmwf,'h',iminecmwf,'mn'
            print*
            read(iunite) pp1
            write(6,*)(idate(i),i = 5,16)
            do iniv = 1,niv
               read(iunite) idate
               write(6,*)(idate(i),i = 5,16)
               read(iunite)((tp1(ilon,ilat,iniv),
     $                       ilon = 1,nlon),ilat = 1,nlat) 
               read(iunite) idate
               write(6,*)(idate(i),i = 5,16)
               read(iunite)((wp1(ilon,ilat,iniv),
     $                       ilon = 1,nlon),ilat = 1,nlat) 
               read(iunite) idate
               write(6,*)(idate(i),i = 5,16)
               read(iunite)((hp1(ilon,ilat,iniv),
     $                       ilon = 1,nlon),ilat = 1,nlat) 
               do ilat = 1,nlat
                  do ilon = 1,nlon
                     hp1(ilon,ilat,iniv) = 
     $               max(hp1(ilon,ilat,iniv),1.e-30*18./28.97)
                  end do
               end do
               read(iunite) idate
               write(6,*)(idate(i),i = 5,16)
               read(iunite)((up1(ilon,ilat,iniv),
     $                       ilon = 1,nlon),ilat = 1,nlat) 
               read(iunite) idate
               write(6,*)(idate(i),i = 5,16)
               read(iunite)((vp1(ilon,ilat,iniv),
     $                       ilon = 1,nlon),ilat = 1,nlat) 
            end do
         else
            print*,'la date ecmwf ne correspond pas, arret du programme'
            print*,'dayecmwf = ',dayecmwf,' dayp1 = ',dayp1
            stop
         end if
c
      end if
c
 100  close(iunite)
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine phot( pm,
     $                 ilat,
     $                 alt,
     $                 sza,
     $                 o3t,
     $                 tj,
     $                 cc_no,
     $                 hnm ) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculates photodissociation frequencies (J coefficients) [sec-1]
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     input : ilat= latitude index
c             sza = solar zenith angle
c             o3t = ozone column
c             day = day/night switch
c     output: tj  = photodiss. freq.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      parameter(nlon = 180, nlat = 91, niv = 137, nivbas = 94)
      parameter(nphot = 41, nz = 101, nsza = 27, nozo = 7)
      parameter(dz = 1.)
c
      real    cc_no(nlon,nivbas)
      real    pm(nlon,nivbas)
      real    alt(nlon,nlat,niv)
      real    tj(nlon,nivbas,nphot)
      real    sza(nlon,nivbas)
      real    o3t(nlon,nivbas), hnm(nlon,nivbas)
      real    xv3(nozo), stepo3, dsza(nsza)
      real    zin(nlon), v3rat(nlon)
      real    weight(nlon,2,2,2)
      real    cizin, cisza, cio3t
c
      integer is(nlon,nivbas), iv(nlon), iz(nlon)
c
      common /photodis/ ajl(nphot-1,nz,nsza,nozo), o3up(0:79)
c
      data dsza/0.,  5., 10., 15., 20., 25., 30., 35., 
     $         40., 45., 50., 55., 60., 65., 70., 75., 
     $         80., 82., 84., 86., 88., 90., 91., 92.,
     $         93., 94., 95./
      data xv3 /0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6/
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialisation des coefficients de photodissociation            c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do n = 1,nphot
         do iniv = 1,nivbas
            do ilon = 1,nlon
               tj(ilon,iniv,n) = 0.
            end do
         end do
      end do
c
      do ilon = 1,nlon
         iv(ilon) = nozo - 1
      end do
      do iniv = 1,nivbas
         do ilon = 1,nlon
            is(ilon,iniv) = nsza - 1
         end do
      end do
c
      do iniv = 1,nivbas
         do i = 1,nsza 
            do ilon = 1,nlon
               if (dsza(i) .ge. sza(ilon,iniv)) then
                  is(ilon,iniv) = min(is(ilon,iniv),i - 1)
               end if
            end do
         end do
      end do
c
      stepo3 = (xv3(nozo) - xv3(1))/(nozo - 1)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     boucle sur l'altitude                                           c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do iniv = nivbas,1,-1
c
         do ilon = 1,nlon
            zin(ilon) = alt(ilon,ilat,iniv)
            zin(ilon) = min(zin(ilon),78.9)
c
            iz(ilon) = int(zin(ilon)/dz) + 1
c
            i = int(zin(ilon))
            v3std = o3up(i)-
     $             (zin(ilon)-real(i))*(o3up(i)-o3up(i+1))
            v3rat(ilon) = o3t(ilon,iniv)/v3std
            v3rat(ilon) = max(v3rat(ilon),xv3(1))
            v3rat(ilon) = min(v3rat(ilon),xv3(nozo))
c
            iv(ilon) = aint((v3rat(ilon) - xv3(1))
     $                       /stepo3) + 1
            iv(ilon) = min(iv(ilon),nozo - 1)
         end do
c
         do ilon = 1,nlon
            if (sza(ilon,iniv) .le. 95.) then
               j = max(iz(ilon),1)
               k = max(is(ilon,iniv),1)
               l = max(iv(ilon),1)
c
               cizin = (zin(ilon)-real(j-1)*dz)/dz
               cisza = (sza(ilon,iniv)-dsza(k))/(dsza(k+1)-dsza(k))
               cio3t = (v3rat(ilon)-xv3(l))/(xv3(l+1)-xv3(l))
c
               weight(ilon,1,1,1) = (1.-cizin)*(1.-cisza)*(1.-cio3t)
               weight(ilon,1,1,2) = (1.-cizin)*(1.-cisza)*cio3t
               weight(ilon,1,2,1) = (1.-cizin)*cisza*(1.-cio3t)
               weight(ilon,1,2,2) = (1.-cizin)*cisza*cio3t
               weight(ilon,2,1,1) = cizin*(1.-cisza)*(1.-cio3t)
               weight(ilon,2,1,2) = cizin*(1.-cisza)*cio3t
               weight(ilon,2,2,1) = cizin*cisza*(1.-cio3t)
               weight(ilon,2,2,2) = cizin*cisza*cio3t
            end if
         end do
c
         do nn = 1,nphot-1
            do ilon = 1,nlon
               if (sza(ilon,iniv) .le. 95.) then
                  j = max(iz(ilon),1)
                  k = max(is(ilon,iniv),1)
                  l = max(iv(ilon),1)
                  tj(ilon,iniv,nn) = 
     $                   weight(ilon,1,1,1)*ajl(nn,j  ,k  ,l  )
     $                 + weight(ilon,1,1,2)*ajl(nn,j  ,k  ,l+1)
     $                 + weight(ilon,1,2,1)*ajl(nn,j  ,k+1,l  )
     $                 + weight(ilon,1,2,2)*ajl(nn,j  ,k+1,l+1)
     $                 + weight(ilon,2,1,1)*ajl(nn,j+1,k  ,l  )
     $                 + weight(ilon,2,1,2)*ajl(nn,j+1,k  ,l+1)
     $                 + weight(ilon,2,2,1)*ajl(nn,j+1,k+1,l  )
     $                 + weight(ilon,2,2,2)*ajl(nn,j+1,k+1,l+1)
               end if
            end do
         end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        fin de la boucle sur les niveaux
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     photodissociation de no
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call jno(cc_no, hnm, nlon, nivbas, o3t, sza, pm, tj(:,:,nphot))
c
      return
      end
c
c************************************************************************
c
      subroutine total( trajpre, irapp, ilat, hc, hnm, o3t )
c
c------------------------------------------------------------------------
c     calcul de la colonne d'ozone [cm-2]
c------------------------------------------------------------------------
c
      parameter(nlon = 180, nlat = 91, niv = 137, nivbas = 94)
      parameter(ncm = 15)
      parameter(nrappmax = 12)
c
      real trajpre(nlon,nlat,niv,nrappmax)
      real hc(nlon,nlat,niv,ncm)
      real hnm(nlon,nivbas)
      real o3t(nlon,nivbas)
      integer ilat
c
c     zcmt = Gas constant/(Boltzmann Cst.*g)
      zcmt = 2.120e20
c
      do ilon = 1,nlon
         o3t(ilon,1) = hc(ilon,ilat,1,5)*hnm(ilon,1)*6.8e5
      end do
      do iniv = 2,nivbas
         do ilon = 1,nlon
            dp = (trajpre(ilon,ilat,iniv  ,irapp)
     $          - trajpre(ilon,ilat,iniv-1,irapp))*100.
            dp = max(dp, 0.)
            o3t(ilon,iniv) = o3t(ilon,iniv-1)
     $                       + zcmt*(hc(ilon,ilat,iniv,5)
     $                              +hc(ilon,ilat,iniv-1,5))*.5*dp
         end do
      end do
c
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine gregtojul(ian,imois,ijour,iheure,imin,jul)
c
c     calcul du jour julien decimal en fonction de la date
c     et de l'heure.
c     on utilise la convention suivante: 
c     1 janvier    0h tu =   0.0
c     31 decembre 24h tu = 365.0 (ou 366.0 si bissextile)
c
      real jul
c
      integer imn(12)
c
      data imn/31,28,31,30,31,30,31,31,30,31,30,31/
c
c     correction du format ecmwf/grib
c 
      if (ian .le. 100) then
         if (ian .le. 40) then
            ian = 2000 + ian
         else
            ian = 1900 + ian
         end if
      end if
c
c     annees bissextiles
c
      if (mod(ian,4) .eq. 0) then
         imn(2) = 29
      else
         imn(2) = 28
      end if
c
      jul = 0.
      do im = 1,imois - 1
         jul = jul + real(imn(im))
      end do  
      jul = jul + real(ijour - 1) + real(iheure)/24. 
     $    + real(imin)/1440.
c   
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine jultogreg(ian,imois,ijour,iheure,imin,jul)
c
c     calcul de la date et de l'heure en fonction du
c     jour julien decimal.
c     on utilise la convention suivante: 
c     1 janvier    0h tu =   0.0
c     31 decembre 24h tu = 365.0 (ou 366.0 si bissextile)
c
      real jul, heure
c
      integer imn(12)
c
      data imn/31,28,31,30,31,30,31,31,30,31,30,31/
c
      if (jul .lt. 0.) then
         ian = ian - 1
         if (mod(ian,4) .eq. 0) then
            jul = jul + 366.
         else
            jul = jul + 365.
         end if
      end if
c         
      ijul   = aint(jul)
      heure  = (jul - real(ijul))*24.
      iheure = aint(heure)
      imin   = nint((heure - real(iheure))*60.)
      if (imin .eq. 60) then
         imin   = 0
         iheure = iheure + 1
         if (iheure .eq. 24) then
            iheure = 0
            ijul   = ijul + 1
         end if
      end if
c
      if (mod(ian,4) .eq. 0) then
         imn(2) = 29
         if (ijul .ge. 366) then
            ijul = ijul - 366
            ian = ian + 1
         end if
      else
         imn(2) = 28
         if (ijul .ge. 365) then
            ijul = ijul - 365
            ian = ian + 1
         end if
      end if
c
      jcumul = 0
      do im = 1,12
         jcumul = jcumul + imn(im)
         if (jcumul - 1 .ge. ijul) then
            imois = im
            goto 10
         end if 
      end do  
 10   continue
c
      jcumul = 0
      do im = 1,imois - 1
         jcumul = jcumul + imn(im)
      end do
      ijour = ijul - (jcumul - 1)
c
      jul = real(jcumul) + real(ijour - 1) + real(iheure)/24.
     $    + real(imin)/1440.
c
      return
      end
c
c*************************************************************************
c
      subroutine zenith(trajlon,
     $                  trajlat,
     $                  irapp,
     $                  ilat,
     $                  day,
     $                  sza )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calcul de l'angle zenithal solaire                                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      parameter(nlon = 180, nlat = 91, niv = 137, nivbas = 94)
      parameter(nrappmax = 12)
c
c input
      real trajlon(nlon,nlat,niv,nrappmax)
      real trajlat(nlon,nlat,niv,nrappmax)
      real day
      integer ilat, irapp
c output
      real sza(nlon,nivbas)
c local
      real ut
      real rlt
c
      real lbut,lzut
      real tz, rdecl, eqr, eqh, zpt
      real csz, zr
      real sintz, costz, sin2tz, cos2tz, sin3tz, cos3tz
c
      pi = 2.*asin(1.0)
c
      ut = (day - aint(day))*24.
c
* Equation 3.8 for "day-angle"

c     tz = 2.*pi*(day - 1.)/365.   bug ?
c
      tz = 2.*pi*day/365.
c
* Calculate sine and cosine from addition theoremes for
* better performance;  the computation of sin2tz,
* sin3tz, cos2tz and cos3tz is about 5-6 times faster
* than the evaluation of the intrinsic functions
*
* It is SIN(x+y) = SIN(x)*COS(y)+COS(x)*SIN(y)
* and   COS(x+y) = COS(x)*COS(y)-SIN(x)*SIN(y)
*
* sintz  = SIN(tz)      costz  = COS(tz)
* sin2tz = SIN(2.*tz)   cos2tz = SIN(2.*tz)
* sin3tz = SIN(3.*tz)   cos3tz = COS(3.*tz)
*
      sintz = SIN(tz)
      costz = COS(tz)
      sin2tz = 2.*sintz*costz
      cos2tz = costz*costz-sintz*sintz
      sin3tz = sintz*cos2tz + costz*sin2tz
      cos3tz = costz*cos2tz - sintz*sin2tz

* Equation 3.7 for declination in radians

      rdecl = 0.006918 - 0.399912*costz  + 0.070257*sintz
     $                 - 0.006758*cos2tz + 0.000907*sin2tz
     $                 - 0.002697*cos3tz + 0.001480*sin3tz

* Equation 3.11 for Equation of time  in radians

      eqr   = 0.000075 + 0.001868*costz  - 0.032077*sintz
     $                 - 0.014615*cos2tz - 0.040849*sin2tz

* convert equation of time to hours:

      eqh = eqr*24./(2.*pi)
c
      do iniv = 1,nivbas
         do ilon = 1,nlon
            rlt = trajlat(ilon,ilat,iniv,irapp)*pi/180.
c
c           calculate local hour angle (hours):
c
            lbut = 12. - eqh - trajlon(ilon,ilat,iniv,irapp)
     $                         *24./360.
c
c           convert to angle from UT
c
            lzut = 15.*(ut - lbut)
            zpt = lzut*pi/180.
c
c           Equation 2.4 for cosine of zenith angle
c
            csz = SIN(rlt)*SIN(rdecl) + COS(rlt)*COS(rdecl)*COS(zpt)
            zr = ACOS(csz)
            sza(ilon,iniv) = zr*180./pi
            sza(ilon,iniv) = max(sza(ilon,iniv),0.)
         end do
      end do
c
      return
      end
c
c***********************************************************************
c
      subroutine stations( day, qj1, hc, trajpre, trajtem, irapp )
c
c------------------------------------------------------------------------
c     calcul des colonnes de o3, o3 passif, no2, oclo
c     au dessus de stations de mesure
c------------------------------------------------------------------------
c
      parameter(nlon = 180, nlat = 91, niv = 137, nivbas = 94)
      parameter(nbcon = 43, ncm = 15)
      parameter(nstation = 40)
      parameter(nrappmax = 12)
c
      real trajpre(nlon,nlat,niv,nrappmax)
      real trajtem(nlon,nlat,niv,nrappmax)
c
      real qj1(nlon,nlat,niv,nbcon)
      real hc(nlon,nlat,niv,ncm)
      real ozone(niv)
      real xla(nstation), xlon(nstation)
      real day
c
      real colo3, colno2, coloclo
c
      character*12 name(nstation)
c
c      1  Scoresbysund     ( 70.50N,  -22.00W) 
c      2  Thule            ( 76.50N,  -68.80W) 
c      3  Ny Aalesund      ( 78.93N,   11.95E) 
c      4  Sodankyla        ( 67.40N,   26.60E) 
c      5  Zhigansk         ( 67.40N,   26.60E) 
c      6  Haute-Provence   ( 43.93N,   05.75E) 
c      7  Aberystwyth      ( 52.00N,  -04.00W)
c      8  Salekhard        ( 66.60N,   66.70E)
c      9  Harestua         ( 60.20N,   10.80E)
c     10  Jungfraujoch     ( 46.32N,    7.60E)
c     11  Verrieres        ( 48.75N,    2.23E)
c     12  Dumont D'Urville (-66.70S,  140.01E)
c     13  Reunion          (-20.88S,   56.00E)
c     14  Bauru            (-22.00S,  -50.00W)
c     15  Kerguelen        (-50.00S,   70.00E)
c     16  Eureka           ( 79.59N,  -85.48W)
c     17  Campbell Island  (-52.60S,  169.20E)
c     18  Macquarie Island (-54.50S,  158.90E)
c     19  Marambio         (-64.27S,  303.28E)
c     20  Faraday          (-65.30S,  295.70E)
c     21  Rothera          (-67.60S,  291.90E)
c     22  Syowa            (-69.01S,   39.60E)
c     23  Neumayer         (-70.70S,  351.80E)
c     24  Terra Nova Bay   (-74.00S,  164.00E)
c     25  Halley Bay       (-75.60S,  333.20E)
c     26  Arrival Heights  (-77.80S,  166.70E)
c     27  Sanae            (-71.67S,  357.16E)
c     28  Concordia        (-75.10S,  123.35E)
c     29  South pole       (-90.00S,  258.00E) 
c     30  Rio Gallegos     (-51.92S,  -69.23W) 
c     31  Belgrano         (-77.85S,  -34.55W) 
c     32  Anadyr           (64.7N,    177.5E)
c     33  Davis            (-68.55,    78.48E)
c     34  Maitri           (-70.5,     11.4E)
c     35  Lauder           (-45.0,    169.7E)
c     36  McMurdo          (-77.83,   166.6E)
c     37  Ushuaia          (-54.9,     68.3W)
c Ajout au 01/01/2018
c     38 So Stromfjord     (66.99N, 50.62W)
c     39 St Petersbourg    (59.88N, 29.82E)
c     40 Villeneuve dAscq  (50.61N, 3.14E)
c


      data name /'scoresbysund','thule       ','nyaalesund  ',
     $           'sodankyla   ','zhigansk    ','ohp         ',
     $           'aberystwyth ','salekhard   ','harestua    ',
     $           'jungfraujoch',
     $           'verrieres   ','ddu         ','reunion     ',
     $           'bauru       ','kerguelen   ','eureka      ',
     $           'campbell    ','macquarie   ','marambio    ',
     $           'faraday     ',
     $           'rothera     ','syowa       ','neumayer    ',
     $           'terra nova  ','halley bay  ','arrival hgt ',
     $           'sanae       ','concordia   ','south pole  ',
     $           'rio gallegos','belgrano    ','anadyr      ',
     $           'davis       ','maitri      ','lauder      ',
     $           'mc murdo    ','ushuaia     ','sstromfjord ', 
     $           'st petersbg ','villeneuveda'/
c
      data xla  /  70.50, 76.50, 78.93,
     $             67.40, 66.72, 43.93,
     $             52.00, 66.60, 60.20,
     $             46.32, 
     $             48.75,-66.70,-20.88,
     $            -22.00,-50.00, 79.59,
     $            -52.60,-54.50,-64.27,
     $            -65.30,
     $            -67.60,-69.01,-70.70,
     $            -74.00,-75.60,-77.80,
     $            -71.67,-75.10,-90.00,
     $            -51.92,-77.85,64.70,
     $            -68.55,-70.5,-45.0,
     $            -77.83,-54.9,66.99,
     $            59.88,50.61/
c
      data xlon / -22.00, -68.80,  11.95,
     $             26.60, 123.40,  05.75,
     $            -04.00,  66.70,  10.80,
     $              7.60,
     $             02.23, 140.01,  56.00,
     $            -50.00,  70.00, -85.48,
     $            169.20, 158.90, 303.28,
     $            295.70,
     $            291.90,  39.60, 351.80,
     $            164.00, 333.20, 166.70,
     $            357.16, 123.35, 258.00,
     $            -69.23, -34.55, 177.5,
     $             78.48,  11.40, 169.7,
     $            166.6,  -68.3, -50.62,
     $            29.82, 3.14 /
c
c     zcmt = gas constant/(boltzmann cst.*g)
c
      zcmt = 2.120e20
      dlon = 360./real(nlon)
      dlat = 180./real(nlat - 1)
c
      do istation = 1,nstation
         if(xlon(istation) .lt. 0.) then
            xlon(istation) = xlon(istation) + 360.
         end if
c
         ilon = nint(xlon(istation)/dlon) + 1
         if (ilon .gt. nlon) then
            ilon = 1
         end if
         ilat = nint((90. - xla(istation))/dlat) + 1
c
c        estimation des colonnes au dessus du modele
c
         hnm     = trajpre(ilon,ilat,1,irapp)
     $           /(trajtem(ilon,ilat,1,irapp)*1.38e-19)
         colo3   = 6.8e5*hc(ilon,ilat,1,5)*hnm
         coloclo = 6.8e5*qj1(ilon,ilat,1,10)*hnm
         colno2  = 6.8e5*hc(ilon,ilat,1,7)*hnm
         colpass = colo3
c
c        equivalence o3 <-> ox en-dessous de nivbas
c
         do iniv = 1,niv
            if (iniv .gt. nivbas) then
               ozone(iniv) = qj1(ilon,ilat,iniv,8)
            else
               ozone(iniv) = hc(ilon,ilat,iniv,5)
            end if
         end do
c
c        colonnes d'ozone et d'ozone passif
c
         do iniv = 2,niv
            dp = (trajpre(ilon,ilat,  iniv,irapp) 
     $          - trajpre(ilon,ilat,iniv-1,irapp))*100.
            colo3   = colo3
     $              + zcmt*(ozone(iniv) + ozone(iniv-1))*.5*dp
            colpass = colpass
     $              + zcmt*(qj1(ilon,ilat,iniv,11)
     $                    + qj1(ilon,ilat,iniv-1,11))*.5*dp
         end do
c
c        colonnes de no2 et oclo
c
         do iniv = 2,nivbas
            dp = (trajpre(ilon,ilat,  iniv,irapp) 
     $          - trajpre(ilon,ilat,iniv-1,irapp))*100.
            colno2  = colno2
     $              + zcmt*(hc(ilon,ilat,iniv,7)
     $                    + hc(ilon,ilat,iniv-1,7))*.5*dp
            coloclo = coloclo
     $              + zcmt*(qj1(ilon,ilat,iniv,10)
     $                    + qj1(ilon,ilat,iniv-1,10))*.5*dp
         end do
c
c        passage en unites dobson
c
         colo3   = colo3*3.72e-17
         colpass = colpass*3.72e-17
c
c        ecriture dans le fichier stations
c
         write(75,1000) name(istation),day,
     $                  colo3,' DU',colpass,' DU',
     $                  colno2,' mol cm-2',
     $                  coloclo,' mol cm-2'
      end do
c
 1000 format(a12,f7.2,f7.2,a3,f7.2,a3,e12.4,a9,e12.4,a9)
c
      return
      end
c
c*************************************************************************
c
      subroutine mir(namexp, daynum, dt, tj1, qj1, hc)
c
      parameter(nlon = 180, nlat = 91, niv = 137)
      parameter(nbcon = 43, ncm = 15)
c
      real qj1(nlon,nlat,niv,nbcon)
      real hc(nlon,nlat,niv,ncm)
      real tj1(nlon,nlat,niv)
      real pression(niv), tempe(niv), theta(niv)
      real xlon(nlon)
      real long(niv,nbcon), short(niv,ncm)
c
      character*6 namexp
      character*12 lab3d(nbcon)
      character*12 labsh(ncm)
c
      data lab3d/'N2O         ','CH4         ','H2O         ',
     +           'NOy         ','HNO3        ','N2O5        ',
     +           'Cly         ','Ox          ','CO          ',
     +           'OClO        ','Passive Ox  ','H2SO4       ',
     +           'HCl         ','ClONO2      ','HOCl        ',
     +           'Cl2         ','H2O2        ','ClNO2       ',
     +           'HBr         ','BrONO2      ','NOx         ',
     +           'HNO4        ','ClOx        ','BrOx        ',
     +           'Cl2O2       ','HOBr        ','BrCl        ',
     +           'CH2O        ','CH3O2       ','CH3O2H      ',
     +           'CFC-11      ','CFC-12      ','CFC-113     ',
     +           'CCl4        ','CH3CCl3*    ','CH3Cl       ',
     +           'HCFC-22*    ','CH3Br       ','H-1211*     ',
     +           'H-1301      ','Bry         ','CH2Br2*     ',
     +           'HNO3 GAS    '/
c
      data labsh/'O(1D)       ','OH          ','Cl          ',
     +           'O(3P)       ','O3          ','HO2         ',
     +           'NO2         ','NO          ','Br          ',
     +           'N           ','ClO         ','BrO         ',
     +           'NO3         ','H           ','CH3         '/
c
      common /t21l30/ pmb(nlon,nlat,niv), xlat(nlat)
c
ccc   analyse mir
c
      daynumi = daynum/86400.
      if (daynumi .ge. 367.) then
         daynumi = daynumi - 366.
      end if
c
      rewind(91)
c
      do i = 1,10000
         read(91,*,end=459) daymir, xlatmir, xlonmir 
         if (xlonmir .lt. 0.) then
            xlonmir = xlonmir + 360.
         end if
c
         delta = dt/(2.*86400.)
         if (abs(daynumi - daymir) .le. delta) then
c
c        calcul des poids
c
            do ilat = 1,nlat
               if (xlat(ilat) .lt. xlatmir) then
                  cinfy = (xlatmir - xlat(ilat))/2.
                  csupy = 1. - cinfy
                  indlat = ilat
                  goto 457
               end if
            end do
 457        continue
            do ilon = 1,nlon
               xlon(ilon) = real(ilon - 1)*2.
            end do
            do ilon = 2,nlon
               if (xlon(ilon) .gt. xlonmir) then
                  cinfx = (xlonmir - xlon(ilon - 1))/2.
                  csupx = 1. - cinfx
                  indlon = ilon
                  goto 458
               end if
            end do
 458        continue
            if (xlonmir .ge. 358.) then
               cinfx = (xlonmir - xlon(nlon))/2.
               csupx = 1. - cinfx
               indlon = 1
            end if
c
            if (xlonmir .lt. 358.) then
c
c     cas general: longitudes inferieures a 358 degres
c
               do iniv = 1,niv
                  temp1 = cinfy*pmb(indlon-1,indlat-1,iniv)
     $                  + csupy*pmb(indlon-1,indlat  ,iniv)
                  temp2 = cinfy*pmb(indlon  ,indlat-1,iniv)
     $                  + csupy*pmb(indlon  ,indlat  ,iniv)
                  pression(iniv) = cinfx*temp2 + csupx*temp1
               end do
c
               do iniv = 1,niv
                  temp1 = cinfy*tj1(indlon-1,indlat-1,iniv)
     $                  + csupy*tj1(indlon-1,indlat  ,iniv)
                  temp2 = cinfy*tj1(indlon  ,indlat-1,iniv)
     $                  + csupy*tj1(indlon  ,indlat  ,iniv)
                  tempe(iniv) = cinfx*temp2 + csupx*temp1
                  theta(iniv) = tempe(iniv)
     $                          *(1000./pression(iniv))**(2./7.)
               end do
c
               do ic = 1,nbcon
                  do iniv = 1,niv
                     temp1 = cinfy*qj1(indlon-1,indlat-1,iniv,ic)
     $                     + csupy*qj1(indlon-1,indlat  ,iniv,ic)
                     temp2 = cinfy*qj1(indlon  ,indlat-1,iniv,ic)
     $                     + csupy*qj1(indlon  ,indlat  ,iniv,ic)
                     long(iniv,ic) = cinfx*temp2 + csupx*temp1
                     long(iniv,ic) = max(long(iniv,ic),1.e-30)
                  end do
               end do
c
               do ic = 1,ncm
                  do iniv = 1,niv
                     temp1 = cinfy*hc(indlon-1,indlat-1,iniv,ic)
     $                     + csupy*hc(indlon-1,indlat  ,iniv,ic)
                     temp2 = cinfy*hc(indlon  ,indlat-1,iniv,ic)
     $                     + csupy*hc(indlon  ,indlat  ,iniv,ic)
                     short(iniv,ic) = cinfx*temp2 + csupx*temp1
                     short(iniv,ic) = max(short(iniv,ic),1.e-30)
                  end do
               end do
            else
c
c     longitudes superieures a 358 degres
c
               do iniv = 1,niv
                  temp1 = cinfy*pmb(nlon,indlat-1,iniv)
     $                  + csupy*pmb(nlon,indlat  ,iniv)
                  temp2 = cinfy*pmb(1   ,indlat-1,iniv)
     $                  + csupy*pmb(1   ,indlat  ,iniv)
                  pression(iniv) = cinfx*temp2 + csupx*temp1
c
                  temp1 = cinfy*tj1(nlon,indlat-1,iniv)
     $                  + csupy*tj1(nlon,indlat  ,iniv)
                  temp2 = cinfy*tj1(1   ,indlat-1,iniv)
     $                  + csupy*tj1(1   ,indlat  ,iniv)
                  tempe(iniv) = cinfx*temp2 + csupx*temp1
                  theta(iniv) = tempe(iniv)
     $                          *(1000./pression(iniv))**(2./7.)
               end do
c
               do ic = 1,nbcon
                  do iniv = 1,niv
                     temp1 = cinfy*qj1(nlon,indlat-1,iniv,ic)
     $                     + csupy*qj1(nlon,indlat  ,iniv,ic)
                     temp2 = cinfy*qj1(1   ,indlat-1,iniv,ic)
     $                     + csupy*qj1(1   ,indlat  ,iniv,ic)
                     long(iniv,ic) = cinfx*temp2 + csupx*temp1
                     long(iniv,ic) = max(long(iniv,ic),1.e-30)
                  end do
               end do
c
               do ic = 1,ncm
                  do iniv = 1,niv
                     temp1 = cinfy*hc(nlon,indlat-1,iniv,ic)
     $                     + csupy*hc(nlon,indlat  ,iniv,ic)
                     temp2 = cinfy*hc(1   ,indlat-1,iniv,ic)
     $                     + csupy*hc(1   ,indlat  ,iniv,ic)
                     short(iniv,ic) = cinfx*temp2 + csupx*temp1
                     short(iniv,ic) = max(short(iniv,ic),
     $                                     1.e-30)
                  end do
               end do
            end if
            write(6,460) 'EXTRACTION PROFILS MIR '
            write(6,461) ' JOUR = ',daymir,
     $                   ' LAT  = ',xlatmir,
     $                   ' LON  = ',xlonmir
            write(92,332) 'RUN ',namexp
            write(92,460) 'EXTRACTION PROFILS MIR'
            write(92,461) ' JOUR = ',daymir,
     $                    ' LAT  = ',xlatmir,
     $                    ' LON  = ',xlonmir
            write(92,*)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c           ecriture proprement dite, en colonnes
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
            write(92,1000)     'PRESSURE    ',
     $                         'TEMPERATURE ',
     $                         'THETA       ',
     $                         (lab3d(is),is = 1,3)
            write(92,*)
            do iniv = 1,niv
               write(92,1001) pression(iniv),
     $                        tempe(iniv),
     $                        theta(iniv),
     $                        (long(iniv,is),is=1,3)
            end do
            write(92,*)
            write(92,1000)'PRESSION    ',(lab3d(is),is = 4,8)
            write(92,*)
            do iniv = 1,niv
               write(92,1001) pression(iniv),
     $                       (long(iniv,is),is=4,8)
            end do
            write(92,*)
            write(92,1000)'PRESSION    ',(lab3d(is),is = 9,13)
            write(92,*)
            do iniv = 1,niv
               write(92,1001) pression(iniv),
     $                       (long(iniv,is),is=9,13)
            end do
            write(92,*)
            write(92,1000)'PRESSION    ',(lab3d(is),is = 14,18)
            write(92,*)
            do iniv = 1,niv
               write(92,1001) pression(iniv),
     $                       (long(iniv,is),is=14,18)
            end do
            write(92,*)
            write(92,1000)'PRESSION    ',(lab3d(is),is = 19,23)
            write(92,*)
            do iniv = 1,niv
               write(92,1001) pression(iniv),
     $                       (long(iniv,is),is=19,23)
            end do
            write(92,*)
            write(92,1000)'PRESSION    ',(lab3d(is),is = 24,28)
            write(92,*)
            do iniv = 1,niv
               write(92,1001) pression(iniv),
     $                       (long(iniv,is),is=24,28)
            end do
            write(92,*)
            write(92,1000)'PRESSION    ',(lab3d(is),is = 29,33)
            write(92,*)
            do iniv = 1,niv
               write(92,1001) pression(iniv),
     $                       (long(iniv,is),is=29,33)
            end do
            write(92,*)
            write(92,1000)'PRESSION    ',(lab3d(is),is = 34,38)
            write(92,*)
            do iniv = 1,niv
               write(92,1001) pression(iniv),
     $                       (long(iniv,is),is=34,38)
            end do
            write(92,*)
            write(92,1000)'PRESSION    ',(lab3d(is),is = 39,40)
            write(92,*)
            do iniv = 1,niv
               write(92,1001) pression(iniv),
     $                       (long(iniv,is),is=39,40)
            end do
            write(92,*)
            write(92,1000)'PRESSION    ',(lab3d(is),is = 41,43)
            write(92,*)
            do iniv = 1,niv
               write(92,1001) pression(iniv),
     $                       (long(iniv,is),is=41,43)
            end do
            write(92,*)
            write(92,1000)'PRESSION    ',(labsh(is),is = 1,5)
            write(92,*)
            do iniv = 1,niv
               write(92,1001) pression(iniv),
     $                       (short(iniv,is),is=1,5)
            end do
            write(92,*)
            write(92,1000)'PRESSION    ',(labsh(is),is = 6,10)
            write(92,*)
            do iniv = 1,niv
               write(92,1001) pression(iniv),
     $                       (short(iniv,is),is=6,10)
            end do
            write(92,*)
            write(92,1000)'PRESSION    ',(labsh(is),is = 11,15)
            write(92,*)
            do iniv = 1,niv
               write(92,1001) pression(iniv),
     $                       (short(iniv,is),is=11,15)
            end do
            write(92,*)
c
         end if
      end do
 459  continue
c
 332  format(1x,a4,a6)
 460  format(a23)
 461  format(3(a8,f10.5))
 462  format(8e12.5)
 1000 format(2x,6a12)
 1001 format(6e12.4)
c
      return
      end
c
c*************************************************************************
c
      subroutine const(relief)
c
      parameter(nlon = 180, nlat = 91, niv = 137)
      parameter(nphot = 41, nz = 101, nsza = 27, nozo = 7)
c
      real relief(nlon,nlat)
c
      common /constant/
     $     a0hox, a1hox, a0nox1,
     $     a1nox1, a0nox2, a1nox2, a0nox3, a1nox3,
     $     a0clx1, a1clx1, a0clx2, a1clx2, a0clx3,
     $     a1clx3, a0brx, a1brx, a0c1, a1c1
c
      common /photodis/ ajl(nphot-1,nz,nsza,nozo), o3up(0:79)
c
ccc   assign variables that have constant values
c
c...hox
         a0hox = 4.4e-32*300.**1.3      ! jpl 2006
         a1hox = 7.5e-11*300.**(-0.2)   ! jpl 2011
c...nox
         a0nox1 = 2.0e-30*300.**4.4
         a1nox1 = 1.4e-12*300.**0.7
         a0nox2 = 1.8e-30*300.**3.0  ! jpl 2006
         a1nox2 = 2.8e-11*300.**0.0  ! jpl 2006
         a0nox3 = 2.0e-31*300.**3.4  ! jpl 2006
         a1nox3 = 2.9e-12*300.**1.1  ! jpl 2006
c...hydrocarbon
         a0c1 = 4.0e-31*300.**3.6    ! jpl 2006
         a1c1 = 1.2e-12*300.**(-1.1) ! jpl 2006
c...clx
         a0clx1 = 1.8e-31*300.**3.4
         a1clx1 = 1.5e-11*300.**1.9
         a0clx2 = 1.8e-31*300.**2.0
         a1clx2 = 1.0e-10*300.**1.0
         a0clx3 = 1.6e-32*300.**4.5  ! jpl 2003
         a1clx3 = 3.0e-12*300.**2.0  ! jpl 2009
c...brx
         a0brx = 5.2e-31*300.**3.2
         a1brx = 6.9e-12*300.**2.9
c
c...lecture des parametres lies au calcul des
c   coefficients de photodissociation
c
         open( unit = 30,
     $        file = 'jstrato.txt',
     $        form = 'formatted')
         do l = 1,nozo
            do k = 1,nsza
               do j = nz,1,-1
                  read(30,*) iz
                  read(30,100) (ajl(i,j,k,l),i=1,nphot-1)
               end do
            end do
         end do
c
         print*,'lecture coefficients de photodissociation: ok'
c
         read(30,*)
         do i = 0,79
            read(30,*) z, o3up(i)
c           write(6,200) z,o3up(i)
         end do
         print*,'lecture de la colonne d ozone standard: ok'
         close(30)
c
c        lecture du relief en km
c
         open(40,form='formatted',file='relief.txt')
         read(40,*) relief
         close(40)
         print*,'lecture du relief: ok'
         do ilat = 1,nlat
            do ilon = 1,nlon
               relief(ilon,ilat) = max(relief(ilon,ilat),0.)
            end do
         end do
c
  100    format(7e11.4)
  200    format(f5.1,e10.3)
c
         return
         end
c
ccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initransp
ccccccccccccccccccccccccccccccccccccccccccccccc
c                                             c
c  initialisations du modele de transport     c
c                                             c
ccccccccccccccccccccccccccccccccccccccccccccccc
c
c  parametres du modele de transport
c
      parameter(nlon = 180, nlat = 91, niv = 137)
c
      real pmod(niv)
c
c  definition de la grille du modele de transport
c
      common/grille/dlon(0:nlon+2),ddlon,dlat(0:nlat+1),ddlat
     +             ,slon(nlon),clon(nlon),slat(nlat),clat(nlat)
     +             ,slonr(nlon,nlat),clonr(nlon,nlat)
     +             ,slatr(nlon,nlat),clatr(nlon,nlat)
     +             ,xcoef(-1:2),ycoef(-1:2),ip(nlon)
c
      common/verniv/aaa(niv),bbb(niv)
c
c  common des constantes universelles
c
      common/cstes/ gg, ra, cp, omega, rter, xp00, rter2, deuomg
     +            , xpi, conv, xpi2, xpih, dlatlim
c
      gg = 9.80665
      ra = 287.05
      cp = 1005.46
      omega = 7.292e-05
      rter = 6371229.
      xp00 = 1000.
c
      rter2 = rter * rter
      deuomg = 2. * omega
      xpi = 2. * asin(1.)
      xpi2 = 2. * xpi
      xpih = xpi / 2.
      conv = 180. / xpi
c
c  initialisation de la grille (globale avec points aux poles)
c
c  dlat : latitudes des points de grille (du n vers le s)
c  dlon : longitudes des points de grille (vers l'est, a partir de 0)
c
      ddlat=xpi/real(nlat-1)
      do 21 k=0,nlat+1
         dlat(k)=xpih-(real(k)-1.)*ddlat
 21   continue
c
      do 22 k=1,nlat
         clat(k)=cos(dlat(k))
         slat(k)=sin(dlat(k))
 22   continue
c
      ddlon=xpi2/real(nlon)
      do 23 i=0,nlon+2
         dlon(i)=ddlon*(real(i)-1.)
 23   continue
c
      do 24 i=1,nlon
         slon(i)=sin(dlon(i))
         clon(i)=cos(dlon(i))
 24   continue
c
      do 25 k=1,nlat
         if (clat(k) .eq. 1.) goto 25
         do 26 i=1,nlon
            slatr(i,k)=-clon(i)*clat(k)
            clatr(i,k)=sqrt(1.-slatr(i,k)*slatr(i,k))
            clonr(i,k)=max(-1.,min(slat(k)/clatr(i,k),1.))
            sl=sqrt(1.-clonr(i,k)*clonr(i,k))
            slonr(i,k)=sign(sl,slon(i))
 26      continue
 25   continue
c
c  coefficients pour le calcul des poids des interpolations cubiques
c     (interpolateur de lagrange)
c
      coefx=1./(ddlon*ddlon*ddlon)
      xcoef(-1)=-coefx/6.
      xcoef( 0)= coefx/2.
      xcoef( 1)=-coefx/2.
      xcoef( 2)= coefx/6.
c
      coefy=1./(ddlat*ddlat*ddlat)
      ycoef(-1)=-coefy/6.
      ycoef( 0)= coefy/2.
      ycoef( 1)=-coefy/2.
      ycoef( 2)= coefy/6.
c
c  indices des points symetriques par rapport aux poles
c
      do 27 i=1,nlon
         ip(i)=mod(i+nlon/2-1,nlon)+1
 27   continue
c
c  niveaux du modele (pour psol=1000 hPa)
c
      do l = 1,niv
         pmod(l) = aaa(l) + bbb(l)*xp00
      end do
c
c  limite pour changer de repere
c
      dlatlim=45./conv
c
c     print *
c     print *,'   latitudes'
c     print *
c     do 3 k=0,nlat+1
c        print *,k,dlat(k)*conv
c3    continue
c
      print *
      print *,'   niveaux (hPa)'
      print *
      do l = 1,niv
         print *,l,pmod(l)
      end do
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine semilag(daynum,dtdyn,nrapp,qj1,
     $                   trajlon,trajlat,trajpre,trajtem,
     $                   ntime,nintecmwf)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c  transport semi-lagrangien                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      parameter(nlon = 180, nlat = 91, niv = 137)
      parameter(nbcon = 43)
      parameter(nrappmax = 12)
c
      real qj1(nlon,nlat,niv,nbcon)
      real trajlon(nlon,nlat,niv,nrappmax)
      real trajlat(nlon,nlat,niv,nrappmax)
      real trajpre(nlon,nlat,niv,nrappmax)
      real trajtem(nlon,nlat,niv,nrappmax)
c
c  calcul des trajectoires arrieres pour les noeuds de la grille
c
      call backtra(daynum,dtdyn,nrapp,
     $             trajlon,trajlat,trajpre,trajtem,
     $             ntime,nintecmwf)
c
c  interpolation des champs de traceurs
c
      do jc = 1,nbcon
         call internew(qj1,jc)
      end do
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine backtra(daynum,dtdyn,nrapp,
     $                   trajlon,trajlat,trajpre,trajtem,
     $                   ntime,nintecmwf)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c  calcul des trajectoires arrieres pour les noeuds de la grille  c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  parametres du modele de transport
c
      parameter(nlon = 180, nlat = 91, niv = 137)
      parameter(nrappmax = 12)
c
      parameter(lola=nlon*nlat)
      parameter(lpro=nlon*nlat*niv)
c
      real trajlon(nlon,nlat,niv,nrappmax)
      real trajlat(nlon,nlat,niv,nrappmax)
      real trajpre(nlon,nlat,niv,nrappmax)
      real trajtem(nlon,nlat,niv,nrappmax)
c
      real t0(nlon,nlat,niv)
c
      real x0r(lpro),y0r(lpro),xmr(lpro),ymr(lpro)
      real dx(lpro),dy(lpro)
c
c  definition de la grille du modele de transport
c
      common/grille/dlon(0:nlon+2),ddlon,dlat(0:nlat+1),ddlat
     +             ,slon(nlon),clon(nlon),slat(nlat),clat(nlat)
     +             ,slonr(nlon,nlat),clonr(nlon,nlat)
     +             ,slatr(nlon,nlat),clatr(nlon,nlat)
     +             ,xcoef(-1:2),ycoef(-1:2),ip(nlon)
c
      common/verniv/aaa(niv),bbb(niv)
c
      common/forcj1/uj1(nlon,nlat,niv),vj1(nlon,nlat,niv)
     +             ,wj1(nlon,nlat,niv),tj1(nlon,nlat,niv)
     +             ,hj1(nlon,nlat,niv),pj1(nlon,nlat)
c
      common/forcm1/um1(nlon,nlat,niv),vm1(nlon,nlat,niv)
     +             ,wm1(nlon,nlat,niv),tm1(nlon,nlat,niv)
     +             ,hm1(nlon,nlat,niv),pm1(nlon,nlat),daym1
c
      common/forcp1/up1(nlon,nlat,niv),vp1(nlon,nlat,niv)
     +             ,wp1(nlon,nlat,niv),tp1(nlon,nlat,niv)
     +             ,hp1(nlon,nlat,niv),pp1(nlon,nlat),dayp1
c
c  vitesse
c
      common/vit/u1s(0:nlon+2,0:nlat+1,0:niv+1)
     +          ,v1s(0:nlon+2,0:nlat+1,0:niv+1)
     +          ,w1s(0:nlon+2,0:nlat+1,0:niv+1)
     +          ,u1m(0:nlon+2,0:nlat+1,0:niv+1)
     +          ,v1m(0:nlon+2,0:nlat+1,0:niv+1)
     +          ,u2(lpro),v2(lpro),w2(lpro)
c
c  position des points origine et medians
c
      common/depo/x0(nlon,nlat,niv),y0(nlon,nlat,niv),z0(nlon,nlat,niv)
     +           ,p0(nlon,nlat,niv)
     +           ,xm(lpro),ym(lpro),zm(lpro),pm(lpro),ind(lpro)
c
c  common des constantes universelles
c
      common/cstes/ gg, ra, cp, omega, rter, xp00, rter2, deuomg
     +            , xpi, conv, xpi2, xpih, dlatlim
c
      common/poids/alfxm1(lpro),alfxp0(lpro),alfxp1(lpro),
     +              alfxp2(lpro),alfym1(lpro),alfyp0(lpro),
     +              alfyp1(lpro),alfyp2(lpro),alfzm1(lpro),
     +              alfzp0(lpro),alfzp1(lpro),alfzp2(lpro)
c
      common/indices/ig0(lpro),kg0(lpro),lg0(lpro)
c
      xpi = 2.*asin(1.)
c
c  dtn est le pas de temps effectif pour le calcul des trajectoires
c
      dtn = dtdyn/real(nrapp)
c     print*,'le pas de temps de calcul des trajectoires est de ',
c    $        dtn/60.,' minutes'
c
c  nrapp iterations de dtn chacune
c
      do 2 iter = 1,nrapp
c
      daypre = daynum - dtn*real(iter-1)/86400.
c
c  interpolation temporelle de la vitesse et de psol
c
      call instant1(daypre,uj1,um1,daym1,up1,dayp1,lpro,
     $              ntime,nintecmwf)
      call instant1(daypre,vj1,vm1,daym1,vp1,dayp1,lpro,
     $              ntime,nintecmwf)
      call instant1(daypre,tj1,tm1,daym1,tp1,dayp1,lpro,
     $              ntime,nintecmwf)
      call instant1(daypre,wj1,wm1,daym1,wp1,dayp1,lpro,
     $              ntime,nintecmwf)
      call instant2(daypre,pj1,pm1,daym1,pp1,dayp1,lola,
     $              ntime,nintecmwf)
c
c  positions initiales: points de la grille de transport
c
      if (iter .eq. 1) then
         do l = 1,niv
            do k = 1,nlat
               do i = 1,nlon
                  x0(i,k,l)=dlon(i)
                  y0(i,k,l)=dlat(k)
                  z0(i,k,l)=(aaa(l)+bbb(l)*pj1(i,k))*100.
                  p0(i,k,l)=pj1(i,k)
                  trajlon(i,k,l,nrapp) = x0(i,k,l)*180./xpi
                  trajlat(i,k,l,nrapp) = y0(i,k,l)*180./xpi
                  trajpre(i,k,l,nrapp) = z0(i,k,l)/100.
                  trajtem(i,k,l,nrapp) = tj1(i,k,l)
               end do
            end do
         end do
      end if
c
c  initialisation du champ de vitesse sur la grille etendue
c
      call inivit( uj1, vj1, wj1, dtn )
c
c  first-guess du point median
c
         do l = 1,niv
            do k = 1,nlat
               do i = 1,nlon
                  ikl = i+(k-1)*nlon+(l-1)*lola
                  xm(ikl) = x0(i,k,l)
                  ym(ikl) = y0(i,k,l)
                  zm(ikl) = z0(i,k,l)
                  pm(ikl) = p0(i,k,l)
                  ind(ikl) = 0
                  if(abs(ym(ikl)).ge.dlatlim) ind(ikl) = 1
               end do
            end do
         end do
c
         do l = 1,niv
            do k = 1,nlat
               do i = 1,nlon
                  ikl = i+(k-1)*nlon+(l-1)*lola
                  if (ind(ikl) .ne. 0) then
                     sy0r = -cos(x0(i,k,l))*cos(y0(i,k,l))
                     y0r(ikl) = asin(sy0r)
                     cy0r = sqrt(1. - sy0r*sy0r)
                     if (cy0r .eq. 0.) then
                        cx0r = 0.
                     else
                        cx0r = max(-1.,min(1.,sin(y0(i,k,l))/cy0r))
                     end if
                     x0r(ikl) = sign(acos(cx0r),sin(x0(i,k,l)))
                  end if
               end do
            end do
         end do
c
c  nitap approximations successives du point median (xm,ym)
c
         nitap = 2
         do 23 it = 1,nitap-1
c
c  initialisation du champ de vitesse:
c
            isum=it+iter
c
c  a) cas particulier de la premiere iteration
c
            if (isum.eq.2) then
c
               do l = 1,niv
                  do k = 1,nlat
                     do i = 1,nlon
                        ikl=i+(k-1)*nlon+(l-1)*lola
                        if (ind(ikl) .eq. 0) then
                           u2(ikl) = u1s(i,k,l)
                           v2(ikl) = v1s(i,k,l)
                        else
                           u2(ikl) = u1m(i,k,l)
                           v2(ikl) = v1m(i,k,l)
                        end if
                        w2(ikl) = w1s(i,k,l)
                     end do
                  end do
               end do
c
            else
c
c  b) cas general: interpolation
c
               call interv
c
            end if
c
c  calcul du point median
c
            do l = 1,niv
               do k = 1,nlat
                  do i = 1,nlon
                     ikl = i+(k-1)*nlon+(l-1)*lola
                     if (ind(ikl) .eq. 0) then
                        xm(ikl) = x0(i,k,l) - 0.5*u2(ikl)/cos(ym(ikl))
                        ym(ikl) = y0(i,k,l) - 0.5*v2(ikl)
                     else
                        symr = -cos(xm(ikl))*cos(ym(ikl))
                        cymr = sqrt(1. - symr*symr)
                        xmr(ikl) = x0r(ikl) - 0.5*u2(ikl)/cymr
                        ymr(ikl) = y0r(ikl) - 0.5*v2(ikl)
c
                        sym = cos(xmr(ikl))*cos(ymr(ikl))
                        ym(ikl) = asin(sym)
                        cym = sqrt(1. - sym*sym)
                        if (cym .eq. 0.) then
                           cxm = 0.
                        else
                           cxm = max(-1.,min(1.,-sin(ymr(ikl))/cym))
                        end if
                        xm(ikl) = sign(acos(cxm),sin(xmr(ikl)))
                     end if
                     zm(ikl) = z0(i,k,l) - 0.5*w2(ikl)
                  end do
               end do
            end do
c
            call corlat(xm,ym,lpro)
c
c  interpolation de la pression sol
c
            call interp(pj1,xm,ym,pm)
c
            call corniv(zm,pm,lpro)
c
 23      continue
c
c  point de depart de la trajectoire arriere pour l'iteration suivante
c
c  interpolation du champ de vitesse
c
         call interv
c
         do l = 1,niv
            do k = 1,nlat
               do i = 1,nlon
                  ikl = i+(k-1)*nlon+(l-1)*lola
                  if (ind(ikl) .eq. 0) then
                     x0(i,k,l) = x0(i,k,l) - u2(ikl)/cos(ym(ikl))
                     y0(i,k,l) = y0(i,k,l) - v2(ikl)
                  else
                     symr = -cos(xm(ikl))*cos(ym(ikl))
                     cymr = sqrt(1. -symr*symr)
                     x0r(ikl) = x0r(ikl) - u2(ikl)/cymr
                     y0r(ikl) = y0r(ikl) - v2(ikl)
c
                     sy0 = cos(x0r(ikl))*cos(y0r(ikl))
                     y0(i,k,l) = asin(sy0)
                     cy0 = sqrt(1. - sy0*sy0)
                     if (cy0 .eq. 0.) then
                        cx0 = 0.
                     else
                        cx0 = max(-1.,min(1.,-sin(y0r(ikl))/cy0))
                     end if
                     x0(i,k,l) = sign(acos(cx0),sin(x0r(ikl)))
                  end if
                  z0(i,k,l) = z0(i,k,l) - w2(ikl)
               end do
            end do
         end do
c
         call corlat(x0,y0,lpro)
c
c  interpolation de la pression sol
c
         call interp(pj1,x0,y0,p0)
c
         call corniv(z0,p0,lpro)
c
c  determination du coin inferieur de la maille contenant (x0,y0,z0)
c
         call posx(x0,ig0,dx,lpro)
         call posy(y0,kg0,dy,lpro)
         call posz(z0,p0,lg0,lpro)
c
c  calcul des poids pour l'interpolateur de lagrange
c
         call poidx(dx,alfxm1,alfxp0,alfxp1,alfxp2,lpro)
         call poidy(dy,alfym1,alfyp0,alfyp1,alfyp2,lpro)
         call poidz(z0,p0,lg0,alfzm1,alfzp0,alfzp1,alfzp2,lpro)
c
c  interpolation de la temperature
c
         call intert(tj1,t0)
c
c  stockage des positions et temperatures le long de la trajectoire
c
         if (iter .lt. nrapp) then
            do l = 1,niv
               do k = 1,nlat
                  do i = 1,nlon
                     trajlon(i,k,l,nrapp-iter) = x0(i,k,l)*180./xpi
                     trajlat(i,k,l,nrapp-iter) = y0(i,k,l)*180./xpi
                     trajpre(i,k,l,nrapp-iter) = z0(i,k,l)/100.
                     trajtem(i,k,l,nrapp-iter) = t0(i,k,l)
                  end do
               end do
            end do
         end if
c
 2    continue
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine instant1(daynum,xj1,xm1,dm1,xp1,dp1,ndim,
     $                    ntime,nintecmwf)
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                  c
c  sous-programme d'interpolation temporelle pour  c
c  avoir les champs du forcage a l'instant t       c
c  entre deux intervalles de forcage               c
c                                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  variable a l'instant t
c
      real xj1(ndim)
      real day
c
c  variable dans les forcages n et n+1 (encadrant l'instant t)
c
      real xm1(ndim),xp1(ndim)
      real dm1, dp1
c
      day = daynum
c
c  correction de l'eventuelle valeur negative 
c  de day (au changement d'annee)
c
      if (day .lt. 1.e-5) then
         day = day + aint(dm1) + 1.
      end if
c
      cp1 = (day - dm1)/(real(nintecmwf)/real(24*ntime))
      cp1 = min(cp1, 1.)
      cp1 = max(cp1, 0.)
      cm1 = 1. - cp1
c
      do 1 n=1,ndim
         xj1(n) = cp1 * xp1(n) + cm1 * xm1(n)
 1    continue
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine instant2(daynum,xj1,xm1,dm1,xp1,dp1,ndim,
     $                    ntime,nintecmwf)
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                  c
c  sous-programme d'interpolation temporelle pour  c
c  avoir la pression sol (en hPa) a l'instant t    c
c  entre deux intervalles de forcage               c
c                                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  variable a l'instant t
c
      real xj1(ndim)
      real day
c
c  variable dans les forcages n et n+1 (encadrant l'instant t)
c
      real xm1(ndim),xp1(ndim)
      real dm1, dp1
c
      day = daynum
c
c  correction de l'eventuelle valeur negative 
c  de day (au changement d'annee)
c
      if (day .lt. 1.e-5) then
         day = day + aint(dm1) + 1.
      end if
c
      cp1 = (day - dm1)/(real(nintecmwf)/real(24*ntime))
      cp1 = min(cp1, 1.)
      cp1 = max(cp1, 0.)
      cm1 = 1. - cp1
c
      do 1 n=1,ndim
         xj1(n) = .01*( cp1*exp(xp1(n)) + cm1*exp(xm1(n)) )
 1    continue
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine inivit(u3d,v3d,w3d,dtn)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                        c
c  sous-programme d'initialisation de la vitesse sur la  c
c  grille etendue pour les interpolations horizontales   c
c                                                        c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  parametres du modele de transport
c
      parameter(nlon = 180, nlat = 91, niv = 137)
c
      parameter(lpro=nlon*nlat*niv)
c
      real u3d(nlon,nlat,niv),v3d(nlon,nlat,niv),w3d(nlon,nlat,niv)
      real dtn
c
c  definition de la grille du modele de transport
c
      common/grille/dlon(0:nlon+2),ddlon,dlat(0:nlat+1),ddlat
     +             ,slon(nlon),clon(nlon),slat(nlat),clat(nlat)
     +             ,slonr(nlon,nlat),clonr(nlon,nlat)
     +             ,slatr(nlon,nlat),clatr(nlon,nlat)
     +             ,xcoef(-1:2),ycoef(-1:2),ip(nlon)
c
c  vitesse
c
      common/vit/u1s(0:nlon+2,0:nlat+1,0:niv+1)
     +          ,v1s(0:nlon+2,0:nlat+1,0:niv+1)
     +          ,w1s(0:nlon+2,0:nlat+1,0:niv+1)
     +          ,u1m(0:nlon+2,0:nlat+1,0:niv+1)
     +          ,v1m(0:nlon+2,0:nlat+1,0:niv+1)
     +          ,u2(lpro),v2(lpro),w2(lpro)
c
c  common des constantes universelles
c
      common/cstes/ gg, ra, cp, omega, rter, xp00, rter2, deuomg
     +            , xpi, conv, xpi2, xpih, dlatlim
c
      coefh=dtn/rter
      coefv=dtn
c
c  noeuds de la grille du modele
c
      do l = 1,niv
         do k = 1,nlat
            do i = 1,nlon
               u1s(i,k,l)=u3d(i,k,l)*coefh
               v1s(i,k,l)=v3d(i,k,l)*coefh
               w1s(i,k,l)=w3d(i,k,l)*coefv
            end do
         end do
      end do 
c
      do l=1,niv
         do k=1,(nlat-1)/2
            do i=1,nlon
               v1m(i,k,l)= (u1s(i,k,l)*slon(i)
     +                     +v1s(i,k,l)*clon(i)*slat(k))
     +                     / clatr(i,k)
               u1m(i,k,l)= (u1s(i,k,l)*clon(i)
     +                     -v1s(i,k,l)*slon(i)*slat(k)
     +                     +v1m(i,k,l)*slonr(i,k)*slatr(i,k))
     +                     / clonr(i,k)
            end do
         end do
         do k=(nlat-1)/2+2,nlat
            do i=1,nlon
               v1m(i,k,l)= (u1s(i,k,l)*slon(i)
     +                     +v1s(i,k,l)*clon(i)*slat(k))
     +                     / clatr(i,k)
               u1m(i,k,l)= (u1s(i,k,l)*clon(i)
     +                     -v1s(i,k,l)*slon(i)*slat(k)
     +                     +v1m(i,k,l)*slonr(i,k)*slatr(i,k))
     +                     / clonr(i,k)
            end do
         end do
      end do
c
c  direction meridienne: symetrie par rapport au pole nord (k=0)
c                      + symetrie par rapport au pole sud  (k=nlat+1)
c
      do l = 1,niv
         do i = 1,nlon
            u1s(i,    0,l)=-u1s(ip(i),  2,l)
            v1s(i,    0,l)=-v1s(ip(i),  2,l)
            w1s(i,    0,l)= w1s(ip(i),  2,l)
            u1s(i,nlat+1,l)=-u1s(ip(i),nlat-1,l)
            v1s(i,nlat+1,l)=-v1s(ip(i),nlat-1,l)
            w1s(i,nlat+1,l)= w1s(ip(i),nlat-1,l)
         end do
      end do
c
      do l = 1,niv
         do i = 1,nlon
            u1m(i,    0,l)=u1m(ip(i),  2,l)
            v1m(i,    0,l)=v1m(ip(i),  2,l)
            u1m(i,nlat+1,l)=u1m(ip(i),nlat-1,l)
            v1m(i,nlat+1,l)=v1m(ip(i),nlat-1,l)
         end do
      end do
c
c  direction zonale (periodicite)
c
      do l = 1,niv
         do k = 0,nlat+1
            u1s(0    ,k,l)=u1s(nlon,k,l)
            v1s(0    ,k,l)=v1s(nlon,k,l)
            w1s(0    ,k,l)=w1s(nlon,k,l)
            u1s(nlon+1,k,l)=u1s(1  ,k,l)
            v1s(nlon+1,k,l)=v1s(1  ,k,l)
            w1s(nlon+1,k,l)=w1s(1  ,k,l)
            u1s(nlon+2,k,l)=u1s(2  ,k,l)
            v1s(nlon+2,k,l)=v1s(2  ,k,l)
            w1s(nlon+2,k,l)=w1s(2  ,k,l)
         end do
      end do
c
      do l=1,niv
         do k=0,(nlat-1)/2
            u1m(0     ,k,l)=u1m(nlon,k,l)
            v1m(0     ,k,l)=v1m(nlon,k,l)
            u1m(nlon+1,k,l)=u1m(1   ,k,l)
            v1m(nlon+1,k,l)=v1m(1   ,k,l)
            u1m(nlon+2,k,l)=u1m(2   ,k,l)
            v1m(nlon+2,k,l)=v1m(2   ,k,l)
         end do
         do k=(nlat-1)/2+2,nlat+1
            u1m(0     ,k,l)=u1m(nlon,k,l)
            v1m(0     ,k,l)=v1m(nlon,k,l)
            u1m(nlon+1,k,l)=u1m(1   ,k,l)
            v1m(nlon+1,k,l)=v1m(1   ,k,l)
            u1m(nlon+2,k,l)=u1m(2   ,k,l)
            v1m(nlon+2,k,l)=v1m(2   ,k,l)
         end do
      end do
c
c  direction verticale
c
      do i = 0,nlon+2
         do k = 0,nlat+1
            u1s(i,k,    0)=u1s(i,k,1)
            v1s(i,k,    0)=v1s(i,k,1)
            w1s(i,k,    0)=w1s(i,k,1)
            u1s(i,k,niv+1)=u1s(i,k,niv)
            v1s(i,k,niv+1)=v1s(i,k,niv)
            w1s(i,k,niv+1)=w1s(i,k,niv)
         end do
      end do
c
      do i=0,nlon+2
         do k=0,(nlat-1)/2
            u1m(i,k,    0)=u1m(i,k,1)
            v1m(i,k,    0)=v1m(i,k,1)
            u1m(i,k,niv+1)=u1m(i,k,niv)
            v1m(i,k,niv+1)=v1m(i,k,niv)
         end do
         do k=(nlat-1)/2+2,nlat+1
            u1m(i,k,    0)=u1m(i,k,1)
            v1m(i,k,    0)=v1m(i,k,1)
            u1m(i,k,niv+1)=u1m(i,k,niv)
            v1m(i,k,niv+1)=v1m(i,k,niv)
         end do
      end do
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine interv
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                              c
c  sous-programme d'interpolation pour le champ de vitesse 3d  c
c  interpolateur de Lagrange: bicubique sur l'horizontale      c
c                             lineaire sur la verticale        c
c  (ecriture valable pour une grille horizontale reguliere )   c
c                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  parametres du modele de transport
c
      parameter(nlon = 180, nlat = 91, niv = 137)
c
      parameter(lpro=nlon*nlat*niv)
c
c  definition de la grille du modele de transport
c
      common/grille/dlon(0:nlon+2),ddlon,dlat(0:nlat+1),ddlat
     +             ,slon(nlon),clon(nlon),slat(nlat),clat(nlat)
     +             ,slonr(nlon,nlat),clonr(nlon,nlat)
     +             ,slatr(nlon,nlat),clatr(nlon,nlat)
     +             ,xcoef(-1:2),ycoef(-1:2),ip(nlon)
c
c  position des points origine et medians
c
      common/depo/x0(nlon,nlat,niv),y0(nlon,nlat,niv),z0(nlon,nlat,niv)
     +           ,p0(nlon,nlat,niv)
     +           ,xm(lpro),ym(lpro),zm(lpro),pm(lpro),ind(lpro)
c
c  vitesse
c
      common/vit/u1s(0:nlon+2,0:nlat+1,0:niv+1)
     +          ,v1s(0:nlon+2,0:nlat+1,0:niv+1)
     +          ,w1s(0:nlon+2,0:nlat+1,0:niv+1)
     +          ,u1m(0:nlon+2,0:nlat+1,0:niv+1)
     +          ,v1m(0:nlon+2,0:nlat+1,0:niv+1)
     +          ,u2(lpro),v2(lpro),w2(lpro)
c
c  common des constantes universelles
c
      common/cstes/ gg, ra, cp, omega, rter, xp00, rter2, deuomg
     +            , xpi, conv, xpi2, xpih, dlatlim
c
      common/poids/alfxm1(lpro),alfxp0(lpro),alfxp1(lpro),
     +              alfxp2(lpro),alfym1(lpro),alfyp0(lpro),
     +              alfyp1(lpro),alfyp2(lpro),alfzm1(lpro),
     +              alfzp0(lpro),alfzp1(lpro),alfzp2(lpro)
c
      common/indices/ig0(lpro),kg0(lpro),lg0(lpro)
c
      dimension dx(lpro),dy(lpro)
      dimension yu(lpro,-1:2)
      dimension zu(lpro,-1:2)
c
c  determination du coin inferieur de la maille contenant (xm,ym,zm)
c
      call posx(xm,ig0,dx,lpro)
      call posy(ym,kg0,dy,lpro)
      call posz(zm,pm,lg0,lpro)
c
c  calcul des poids pour l'interpolateur de lagrange
c
      call poidx(dx,alfxm1,alfxp0,alfxp1,alfxp2,lpro)
      call poidy(dy,alfym1,alfyp0,alfyp1,alfyp2,lpro)
      call poidz(zm,pm,lg0,alfzm1,alfzp0,alfzp1,alfzp2,lpro)
c
c  boucle sur les niveaux
c
      do ls=-1,2
c
c  boucle sur les latitudes
c
         do ks=-1,2
c
c  interpolation selon x du vent zonal
c
            do n=1,lpro
               i=ig0(n)
               k=kg0(n)+ks
               l=lg0(n)+ls
c
               if (ind(n) .eq. 0) then
                  xum1 = u1s(i-1,k,l)
                  xup0 = u1s(i  ,k,l)
                  xup1 = u1s(i+1,k,l)
                  xup2 = u1s(i+2,k,l)
               else
                  xum1 = u1m(i-1,k,l)
                  xup0 = u1m(i  ,k,l)
                  xup1 = u1m(i+1,k,l)
                  xup2 = u1m(i+2,k,l)
               end if
c
               yu(n,ks)=alfxm1(n)*xum1+alfxp0(n)*xup0
     +                 +alfxp1(n)*xup1+alfxp2(n)*xup2
c
            end do
c
         end do
c
c  interpolation selon y du vent zonal
c
         do n=1,lpro
c
            yum1=yu(n,-1)
            yup0=yu(n, 0)
            yup1=yu(n, 1)
            yup2=yu(n, 2)
c
            zu(n,ls)=alfym1(n)*yum1+alfyp0(n)*yup0
     +              +alfyp1(n)*yup1+alfyp2(n)*yup2
c
         end do
c
      end do
c
c  interpolation selon z du vent zonal
c
      do n=1,lpro
c
         zum1=zu(n,-1)
         zup0=zu(n, 0)
         zup1=zu(n, 1)
         zup2=zu(n, 2)
c
         u2(n)=alfzm1(n)*zum1+alfzp0(n)*zup0
     +        +alfzp1(n)*zup1+alfzp2(n)*zup2
c
      end do
c
c  boucle sur les niveaux
c
      do ls=-1,2
c
c  boucle sur les latitudes
c
         do ks=-1,2
c
c  interpolation selon x du vent meridien
c
            do n=1,lpro
               i=ig0(n)
               k=kg0(n)+ks
               l=lg0(n)+ls
c
               if (ind(n) .eq. 0) then
                  xvm1 = v1s(i-1,k,l)
                  xvp0 = v1s(i  ,k,l)
                  xvp1 = v1s(i+1,k,l)
                  xvp2 = v1s(i+2,k,l)
               else
                  xvm1 = v1m(i-1,k,l)
                  xvp0 = v1m(i  ,k,l)
                  xvp1 = v1m(i+1,k,l)
                  xvp2 = v1m(i+2,k,l)
               end if
c
               yu(n,ks)=alfxm1(n)*xvm1+alfxp0(n)*xvp0
     +                 +alfxp1(n)*xvp1+alfxp2(n)*xvp2
c
            end do
c
         end do
c
c  interpolation selon y du vent meridien
c
         do n=1,lpro
c
            yvm1=yu(n,-1)
            yvp0=yu(n, 0)
            yvp1=yu(n, 1)
            yvp2=yu(n, 2)
c
            zu(n,ls)=alfym1(n)*yvm1+alfyp0(n)*yvp0
     +              +alfyp1(n)*yvp1+alfyp2(n)*yvp2
c
         end do
c
      end do
c
c  interpolation selon z du vent meridien
c
      do n=1,lpro
c
         zvm1=zu(n,-1)
         zvp0=zu(n, 0)
         zvp1=zu(n, 1)
         zvp2=zu(n, 2)
c
         v2(n)=alfzm1(n)*zvm1+alfzp0(n)*zvp0
     +        +alfzp1(n)*zvp1+alfzp2(n)*zvp2
c
      end do
c
c  boucle sur les niveaux
c
      do ls=-1,2
c
c  boucle sur les latitudes
c
         do ks=-1,2
c
c  interpolation selon x du vent vertical
c
            do n=1,lpro
               i=ig0(n)
               k=kg0(n)+ks
               l=lg0(n)+ls
c
               xwm1=w1s(i-1,k,l)
               xwp0=w1s(i  ,k,l)
               xwp1=w1s(i+1,k,l)
               xwp2=w1s(i+2,k,l)
c
               yu(n,ks)=alfxm1(n)*xwm1+alfxp0(n)*xwp0
     +                 +alfxp1(n)*xwp1+alfxp2(n)*xwp2
c
            end do
c
         end do
c
c  interpolation selon y du vent vertical
c
         do n=1,lpro
c
            ywm1=yu(n,-1)
            ywp0=yu(n, 0)
            ywp1=yu(n, 1)
            ywp2=yu(n, 2)
c
            zu(n,ls)=alfym1(n)*ywm1+alfyp0(n)*ywp0
     +              +alfyp1(n)*ywp1+alfyp2(n)*ywp2
c
         end do
c
      end do
c
c  interpolation selon z du vent vertical
c
      do n=1,lpro
c
         zwm1=zu(n,-1)
         zwp0=zu(n, 0)
         zwp1=zu(n, 1)
         zwp2=zu(n, 2)
c
         w2(n)=alfzm1(n)*zwm1+alfzp0(n)*zwp0
     +        +alfzp1(n)*zwp1+alfzp2(n)*zwp2
c
      end do
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine interp(psol,xx,yy,pres)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                              c
c  sous-programme d'interpolation pour la pression de surface  c
c  interpolateur de Lagrange: bicubique sur l'horizontale      c
c  (ecriture valable pour une grille horizontale reguliere )   c
c                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  parametres du modele de transport
c
      parameter(nlon = 180, nlat = 91, niv = 137)
c
      parameter(lpro=nlon*nlat*niv)
c
c  definition de la grille du modele de transport
c
      common/grille/dlon(0:nlon+2),ddlon,dlat(0:nlat+1),ddlat
     +             ,slon(nlon),clon(nlon),slat(nlat),clat(nlat)
     +             ,slonr(nlon,nlat),clonr(nlon,nlat)
     +             ,slatr(nlon,nlat),clatr(nlon,nlat)
     +             ,xcoef(-1:2),ycoef(-1:2),ip(nlon)
c
      common/poids/alfxm1(lpro),alfxp0(lpro),alfxp1(lpro),
     +              alfxp2(lpro),alfym1(lpro),alfyp0(lpro),
     +              alfyp1(lpro),alfyp2(lpro),alfzm1(lpro),
     +              alfzp0(lpro),alfzp1(lpro),alfzp2(lpro)
c
      common/indices/ig0(lpro),kg0(lpro),lg0(lpro)
c
      dimension dx(lpro),dy(lpro)
      dimension y(lpro,-1:2)
      dimension x(0:nlon+2,0:nlat+1)
c
      real xx(lpro),yy(lpro),psol(nlon,nlat)
      real pres(lpro)
c
c  determination du coin inferieur de la maille contenant (x0,y0)
c
      call posx(xx,ig0,dx,lpro)
      call posy(yy,kg0,dy,lpro)
c
c  calcul des poids pour l'interpolateur de lagrange
c
      call poidx(dx,alfxm1,alfxp0,alfxp1,alfxp2,lpro)
      call poidy(dy,alfym1,alfyp0,alfyp1,alfyp2,lpro)
c
c  extension de la grille pour les interpolations cubiques
c
      do 11 k=1,nlat
         do 110 i=1,nlon
            x(i,k)=psol(i,k)
 110     continue
 11   continue
c
c  direction meridienne: symetrie par rapport au pole nord (k=0)
c                      + symetrie par rapport au pole sud  (k=nlat+1)
c
      do 12 i=1,nlon
         x(i,     0)=x(ip(i),  2)
         x(i,nlat+1)=x(ip(i),nlat-1)
 12   continue
c
c  direction zonale (periodicite)
c
      do 13 k=0,nlat+1
         x(0    ,k)=x(nlon,k)
         x(nlon+1,k)=x(1  ,k)
         x(nlon+2,k)=x(2  ,k)
 13   continue
c
c  interpolations proprement dites
c
c  boucle sur les latitudes
c
      do 14 ks=-1,2
c
c  interpolation selon x
c
         do 140 n=1,lpro
c
            i=ig0(n)
            k=kg0(n)+ks
c
            xm1=x(i-1,k)
            xp0=x(i  ,k)
            xp1=x(i+1,k)
            xp2=x(i+2,k)
c
            ppy=alfxm1(n)*xm1+alfxp0(n)*xp0
     +         +alfxp1(n)*xp1+alfxp2(n)*xp2
c
c  limitation pour eviter l'apparition de nouveaux extremas
c
            rmx=max(xp0,xp1)
            rmn=min(xp0,xp1)
            ppy=max(min(ppy,rmx),rmn)
c
            y(n,ks)=ppy
c
 140     continue
c
 14   continue
c
c  interpolation selon y
c
      do 15 n=1,lpro
c
         ym1=y(n,-1)
         yp0=y(n, 0)
         yp1=y(n, 1)
         yp2=y(n, 2)
c
         ppz=alfym1(n)*ym1+alfyp0(n)*yp0
     +      +alfyp1(n)*yp1+alfyp2(n)*yp2
c
c  limitation pour eviter l'apparition de nouveaux extremas
c
         rmx=max(yp0,yp1)
         rmn=min(yp0,yp1)
         ppz=max(min(ppz,rmx),rmn)
c
         pres(n)=ppz
c
 15   continue
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine intert(tj1,t0)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                              c
c  sous-programme d'interpolation de la temperature            c
c  interpolateur de Lagrange: bicubique sur l'horizontale      c
c                             lineaire sur la verticale        c
c  (ecriture valable pour une grille horizontale reguliere )   c
c                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  parametres du modele de transport
c
      parameter(nlon = 180, nlat = 91, niv = 137)
c
      parameter(lpro=nlon*nlat*niv)
      parameter(lola=nlon*nlat)
c
      common/grille/dlon(0:nlon+2),ddlon,dlat(0:nlat+1),ddlat
     +             ,slon(nlon),clon(nlon),slat(nlat),clat(nlat)
     +             ,slonr(nlon,nlat),clonr(nlon,nlat)
     +             ,slatr(nlon,nlat),clatr(nlon,nlat)
     +             ,xcoef(-1:2),ycoef(-1:2),ip(nlon)
c
      common/poids/ alfxm1(lpro),alfxp0(lpro),alfxp1(lpro),
     +              alfxp2(lpro),alfym1(lpro),alfyp0(lpro),
     +              alfyp1(lpro),alfyp2(lpro),alfzm1(lpro),
     +              alfzp0(lpro),alfzp1(lpro),alfzp2(lpro)
c
      common/indices/ ig0(lpro),kg0(lpro),lg0(lpro)
c
      dimension y(lpro,-1:2),z(lpro,-1:2)
      dimension x(0:nlon+2,0:nlat+1,0:niv+1)
c
      real tj1(nlon,nlat,niv)
      real t0(nlon,nlat,niv)
      real tempo(lpro)
c
c  extension de la grille pour les interpolations cubiques
c
      do l = 1,niv
         do k = 1,nlat
            do i = 1,nlon
               x(i,k,l)=tj1(i,k,l)
            end do
         end do
      end do 
c
c  direction meridienne: symetrie par rapport au pole nord (k=0)
c                      + symetrie par rapport au pole sud  (k=nlat+1)
c
      do l = 1,niv
         do i = 1,nlon
            x(i,     0,l)=x(ip(i),     2,l)
            x(i,nlat+1,l)=x(ip(i),nlat-1,l)
         end do
      end do
c
c  direction zonale (periodicite)
c
      do l = 1,niv
         do k = 0,nlat+1
            x(0    ,k,l)=x(nlon,k,l)
            x(nlon+1,k,l)=x(1  ,k,l)
            x(nlon+2,k,l)=x(2  ,k,l)
         end do
      end do
c
c  direction verticale
c
      do i = 0,nlon+2
         do k = 0,nlat+1
            x(i,k,0    )=x(i,k,1)
            x(i,k,niv+1)=x(i,k,niv)
         end do
      end do
c
c  interpolations proprement dites
c
c  boucle sur les niveaux
c
      do 15 ls=-1,2
c
c  boucle sur les latitudes
c
         do 150 ks=-1,2
c
c  interpolation selon x
c
            do 1500 n=1,lpro
c
               i=ig0(n)
               k=kg0(n)+ks
               l=lg0(n)+ls
c
               xm1=x(i-1,k,l)
               xp0=x(i  ,k,l)
               xp1=x(i+1,k,l)
               xp2=x(i+2,k,l)
c
               yy=alfxm1(n)*xm1+alfxp0(n)*xp0
     +           +alfxp1(n)*xp1+alfxp2(n)*xp2
c
c  limitation pour eviter l'apparition de nouveaux extremas
c
               rmx=max(xp0,xp1)
               rmn=min(xp0,xp1)
               yy=max(min(yy,rmx),rmn)
c
               y(n,ks)=yy
c
 1500       continue
c
 150     continue
c
c  interpolation selon y
c
         do 151 n=1,lpro
c
            ym1=y(n,-1)
            yp0=y(n, 0)
            yp1=y(n, 1)
            yp2=y(n, 2)
c
            zz=alfym1(n)*ym1+alfyp0(n)*yp0
     +        +alfyp1(n)*yp1+alfyp2(n)*yp2
c
c  limitation pour eviter l'apparition de nouveaux extremas
c
            rmx=max(yp0,yp1)
            rmn=min(yp0,yp1)
            zz=max(min(zz,rmx),rmn)
c
            z(n,ls)=zz
c
 151     continue
c
 15   continue
c
c  interpolation selon z
c
      do 16 n=1,lpro
c
         zm1=z(n,-1)
         zp0=z(n, 0)
         zp1=z(n, 1)
         zp2=z(n, 2)
c
         zz=alfzm1(n)*zm1+alfzp0(n)*zp0
     +     +alfzp1(n)*zp1+alfzp2(n)*zp2
c
c  limitation pour eviter l'apparition de nouveaux extremas
c
         rmx=max(zp0,zp1)
         rmn=min(zp0,zp1)
         tempo(n)=max(min(zz,rmx),rmn)
c
 16   continue
c
      do l = 1,niv
         do k = 1,nlat
            do i = 1,nlon
               ikl = i+(k-1)*nlon+(l-1)*lola
               t0(i,k,l) = tempo(ikl)
            end do
         end do
      end do
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine corlat(x,y,ndim)
cccccccccccccccccccccccccccccccccccccccccccccccc
c                                              c
c  correction des latitudes phi>pi/2 et phi<0  c
c                                              c
cccccccccccccccccccccccccccccccccccccccccccccccc
c
c  common des constantes universelles
c
      common/cstes/ gg, ra, cp, omega, rter, xp00, rter2, deuomg
     +            , xpi, conv, xpi2, xpih, dlatlim
c
      dimension x(ndim),y(ndim)
c
c  la boucle 1 ramene x et y a des angles compris entre 0 et 2 pi.
c
      do 1 n=1,ndim
         x(n)=x(n)-xpi2*int(x(n)/xpi2)
         if (x(n) .lt. 0.) then
            x(n) = x(n) + xpi2
         end if
         y(n)=y(n)-xpi2*int(y(n)/xpi2)
 1    continue
c
c  in temp1 longitudinal angles less than pi get translated by pi,
c                        angles greater than pi get translated by -pi
c  ie. temp1 is rotated about the pole by pi.
c
c  temp2 reflects the latitudinal angle about the pole
c           (eg. 91 to 89 degrees).
c
      do 2 n=1,ndim
         temp0=xpih-abs(y(n))
         if (x(n) - xpi .lt. 0.) then
            temp1 = x(n) + xpi
         else
            temp1 = x(n) - xpi
         end if
         temp2=sign(xpi,y(n))-y(n)
         if (temp0 .lt. 0.) then
            x(n) = temp1
            y(n) = temp2
         end if
c        y(n)=max(y(n),0.)
 2    continue
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine corniv(z,p,ndim)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                       c
c  correction des pressions (ramene z entre niveaux pmb(1) et pmb(niv)  c
c                                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  parametres du modele de transport
c
      parameter(nlon = 180, nlat = 91, niv = 137)
c
c  definition de la grille du modele de transport
c
      common/grille/dlon(0:nlon+2),ddlon,dlat(0:nlat+1),ddlat
     +             ,slon(nlon),clon(nlon),slat(nlat),clat(nlat)
     +             ,slonr(nlon,nlat),clonr(nlon,nlat)
     +             ,slatr(nlon,nlat),clatr(nlon,nlat)
     +             ,xcoef(-1:2),ycoef(-1:2),ip(nlon)
c
      common/verniv/aaa(niv),bbb(niv)
c
c  common des constantes universelles
c
      common/cstes/ gg, ra, cp, omega, rter, xp00, rter2, deuomg
     +            , xpi, conv, xpi2, xpih, dlatlim
c
      dimension z(ndim),p(ndim)
c
      do 1 n = 1,ndim
         pmin=( aaa(  1)+bbb(  1)*p(n) )*100.
         pmax=( aaa(niv)+bbb(niv)*p(n) )*100.
         z(n)=min(pmax,max(pmin,z(n)))
 1    continue
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine posx(xx,indx,dx,ndim)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                      c
c  sous-programme de positionnement des longitudes xx  c
c  par rapport a celles du tableau dlon (regulier)     c
c                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  parametres du modele de transport
c
      parameter(nlon = 180, nlat = 91)
c
c  definition de la grille du modele de transport
c
      common/grille/dlon(0:nlon+2),ddlon,dlat(0:nlat+1),ddlat
     +             ,slon(nlon),clon(nlon),slat(nlat),clat(nlat)
     +             ,slonr(nlon,nlat),clonr(nlon,nlat)
     +             ,slatr(nlon,nlat),clatr(nlon,nlat)
     +             ,xcoef(-1:2),ycoef(-1:2),ip(nlon)
c
c  common des constantes universelles
c
      common/cstes/ gg, ra, cp, omega, rter, xp00, rter2, deuomg
     +            , xpi, conv, xpi2, xpih, dlatlim
c
      dimension xx(ndim)
      dimension indx(ndim),dx(ndim)
c
c  determination de la longitude dlon immediatement inferieure a xx
c
      do 1 n=1,ndim
         ix=int(xx(n)/ddlon)
         indx(n)=ix+1
         dx(n)=xx(n)-real(ix)*ddlon
 1    continue
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine posy(yy,indy,dy,ndim)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                     c
c  sous-programme de positionnement des latitudes yy  c
c  par rapport a celles du tableau dlat (regulier)    c
c                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  parametres du modele de transport
c
      parameter(nlon = 180, nlat = 91)
c
c  definition de la grille du modele de transport
c
      common/grille/dlon(0:nlon+2),ddlon,dlat(0:nlat+1),ddlat
     +             ,slon(nlon),clon(nlon),slat(nlat),clat(nlat)
     +             ,slonr(nlon,nlat),clonr(nlon,nlat)
     +             ,slatr(nlon,nlat),clatr(nlon,nlat)
     +             ,xcoef(-1:2),ycoef(-1:2),ip(nlon)
c
c  common des constantes universelles
c
      common/cstes/ gg, ra, cp, omega, rter, xp00, rter2, deuomg
     +            , xpi, conv, xpi2, xpih, dlatlim
c
      dimension yy(ndim)
      dimension indy(ndim),dy(ndim)
c
c  determination de la latitude dlat immediatement superieure a yy
c
      do 1 n=1,ndim
         iy=int((xpih-yy(n))/ddlat)
         indy(n)=iy+1
         dy(n)=xpih-real(iy)*ddlat-yy(n)
 1    continue
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine posz(zz,pp,indz,ndim)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                     c
c  sous-programme de positionnement des pressions zz  c
c                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  parametres du modele de transport
c
      parameter(nlon = 180, nlat = 91, niv = 137)
c
c  definition de la grille du modele de transport
c
      common/grille/dlon(0:nlon+2),ddlon,dlat(0:nlat+1),ddlat
     +             ,slon(nlon),clon(nlon),slat(nlat),clat(nlat)
     +             ,slonr(nlon,nlat),clonr(nlon,nlat)
     +             ,slatr(nlon,nlat),clatr(nlon,nlat)
     +             ,xcoef(-1:2),ycoef(-1:2),ip(nlon)
c
      common/verniv/aaa(niv),bbb(niv)
c
c  common des constantes universelles
c
      common/cstes/ gg, ra, cp, omega, rter, xp00, rter2, deuomg
     +            , xpi, conv, xpi2, xpih, dlatlim
c
      dimension zz(ndim),pp(ndim)
      dimension indz(ndim)
c
c  determination de la pression pres immediatement inferieure a zz
c
      do 1 n=1,ndim
         indz(n)=1
 1    continue
c
      do 2 l = 2,niv-1
         do 20 n = 1,ndim
            pres=( aaa(l)+bbb(l)*pp(n) )*100.
            if (zz(n) - pres .ge. 0.) then
               indz(n) = l
            end if
 20      continue
 2    continue
c
      return
      end
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine poidx(dx,alfxm1,alfxp0,alfxp1,alfxp2,ndim)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  sous-programme de calcul des poids pour l'interpolateur de lagrange
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  parametres du modele de transport
c
      parameter(nlon = 180, nlat = 91)
c
c  definition de la grille du modele de transport
c
      common/grille/dlon(0:nlon+2),ddlon,dlat(0:nlat+1),ddlat
     +             ,slon(nlon),clon(nlon),slat(nlat),clat(nlat)
     +             ,slonr(nlon,nlat),clonr(nlon,nlat)
     +             ,slatr(nlon,nlat),clatr(nlon,nlat)
     +             ,xcoef(-1:2),ycoef(-1:2),ip(nlon)
c
c  common des constantes universelles
c
      common/cstes/ gg, ra, cp, omega, rter, xp00, rter2, deuomg
     +            , xpi, conv, xpi2, xpih, dlatlim
c
      dimension dx(ndim)
      dimension alfxm1(ndim),alfxp0(ndim),alfxp1(ndim),alfxp2(ndim)
c
      do 1 n=1,ndim
         dxm1=dx(n)+ddlon
         dxp0=dx(n)
         dxp1=dx(n)-ddlon
         dxp2=dx(n)-ddlon*2.
         alfxm1(n)=xcoef(-1)*dxp0*dxp1*dxp2
         alfxp0(n)=xcoef( 0)*dxm1*dxp1*dxp2
         alfxp1(n)=xcoef( 1)*dxm1*dxp0*dxp2
         alfxp2(n)=xcoef( 2)*dxm1*dxp0*dxp1
 1    continue
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine poidy(dy,alfym1,alfyp0,alfyp1,alfyp2,ndim)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  sous-programme de calcul des poids pour l'interpolateur de lagrange
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  parametres du modele de transport
c
      parameter(nlon = 180, nlat = 91)
c
c  definition de la grille du modele de transport
c
      common/grille/dlon(0:nlon+2),ddlon,dlat(0:nlat+1),ddlat
     +             ,slon(nlon),clon(nlon),slat(nlat),clat(nlat)
     +             ,slonr(nlon,nlat),clonr(nlon,nlat)
     +             ,slatr(nlon,nlat),clatr(nlon,nlat)
     +             ,xcoef(-1:2),ycoef(-1:2),ip(nlon)
c
c  common des constantes universelles
c
      common/cstes/ gg, ra, cp, omega, rter, xp00, rter2, deuomg
     +            , xpi, conv, xpi2, xpih, dlatlim
c
      dimension dy(ndim)
      dimension alfym1(ndim),alfyp0(ndim),alfyp1(ndim),alfyp2(ndim)
c
      do 1 n=1,ndim
         dym1=dy(n)+ddlat
         dyp0=dy(n)
         dyp1=dy(n)-ddlat
         dyp2=dy(n)-ddlat*2.
         alfym1(n)=ycoef(-1)*dyp0*dyp1*dyp2
         alfyp0(n)=ycoef( 0)*dym1*dyp1*dyp2
         alfyp1(n)=ycoef( 1)*dym1*dyp0*dyp2
         alfyp2(n)=ycoef( 2)*dym1*dyp0*dyp1
 1    continue
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine poidz(zz,pp,indz,alfzm1,alfzp0,alfzp1,alfzp2,ndim)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  sous-programme de calcul des poids pour l'interpolateur de lagrange
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  parametres du modele de transport
c
      parameter(nlon = 180, nlat = 91, niv = 137)
c
c  definition de la grille du modele de transport
c
      common/grille/dlon(0:nlon+2),ddlon,dlat(0:nlat+1),ddlat
     +             ,slon(nlon),clon(nlon),slat(nlat),clat(nlat)
     +             ,slonr(nlon,nlat),clonr(nlon,nlat)
     +             ,slatr(nlon,nlat),clatr(nlon,nlat)
     +             ,xcoef(-1:2),ycoef(-1:2),ip(nlon)
c
      common/verniv/aaa(niv),bbb(niv)
c
c  common des constantes universelles
c
      common/cstes/ gg, ra, cp, omega, rter, xp00, rter2, deuomg
     +            , xpi, conv, xpi2, xpih, dlatlim
c
      dimension zz(ndim),pp(ndim),indz(ndim)
      dimension alfzm1(ndim),alfzp0(ndim),alfzp1(ndim),alfzp2(ndim)
c
      do 1 n = 1,ndim
         ll = indz(n)
         if (ll .gt. 1) then
            prem1 = (aaa(ll-1) + bbb(ll-1)*pp(n))*100.
         else
            prem1 = 5.e-3*100.
         end if
         prep0 = (aaa(ll  ) + bbb(ll  )*pp(n))*100.
         prep1 = (aaa(ll+1) + bbb(ll+1)*pp(n))*100.
         if (ll .lt. (niv-1)) then
            prep2 = (aaa(ll+2) + bbb(ll+2)*pp(n))*100.
         else
            prep2 = 1200.*100.
         end if
c
         dzm1=zz(n)-prem1
         dzp0=zz(n)-prep0
         dzp1=zz(n)-prep1
         dzp2=zz(n)-prep2
c
         dzp2m1=prep2-prem1
         dzp2p0=prep2-prep0
         dzp2p1=prep2-prep1
         dzp1m1=prep1-prem1
         dzp1p0=prep1-prep0
         dzp0m1=prep0-prem1
c
         alfzm1(n)=-(dzp0*dzp1*dzp2)/(dzp0m1*dzp1m1*dzp2m1)
         alfzp0(n)= (dzm1*dzp1*dzp2)/(dzp0m1*dzp1p0*dzp2p0)
         alfzp1(n)=-(dzm1*dzp0*dzp2)/(dzp1m1*dzp1p0*dzp2p1)
         alfzp2(n)= (dzm1*dzp0*dzp1)/(dzp2m1*dzp2p0*dzp2p1)
c
 1    continue
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine inter(qj1,jc)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                              c
c  sous-programme d'interpolation pour les champs de traceurs  c
c  interpolateur de Lagrange: bicubique sur l'horizontale      c
c                             lineaire sur la verticale        c
c  (ecriture valable pour une grille horizontale reguliere )   c
c                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  parametres du modele de transport
c
      parameter(nlon = 180, nlat = 91, niv = 137)
      parameter(nbcon = 43)
c
      parameter(lpro=nlon*nlat*niv)
      parameter(lola=nlon*nlat)
c
      common/grille/dlon(0:nlon+2),ddlon,dlat(0:nlat+1),ddlat
     +             ,slon(nlon),clon(nlon),slat(nlat),clat(nlat)
     +             ,slonr(nlon,nlat),clonr(nlon,nlat)
     +             ,slatr(nlon,nlat),clatr(nlon,nlat)
     +             ,xcoef(-1:2),ycoef(-1:2),ip(nlon)
c
      common/poids/ alfxm1(lpro),alfxp0(lpro),alfxp1(lpro),
     +              alfxp2(lpro),alfym1(lpro),alfyp0(lpro),
     +              alfyp1(lpro),alfyp2(lpro),alfzm1(lpro),
     +              alfzp0(lpro),alfzp1(lpro),alfzp2(lpro)
c
      common/indices/ ig0(lpro),kg0(lpro),lg0(lpro)
c
      dimension y(lpro,-1:2),z(lpro,-1:2)
      dimension x(0:nlon+2,0:nlat+1,0:niv+1)
c
      real qj1(nlon,nlat,niv,nbcon)
      real tempo(lpro)
c
c  extension de la grille pour les interpolations cubiques
c
      do l = 1,niv
         do k = 1,nlat
            do i = 1,nlon
               x(i,k,l) = qj1(i,k,l,jc)
            end do
         end do
      end do
c
c  direction meridienne: symetrie par rapport au pole nord (k=0)
c                      + symetrie par rapport au pole sud  (k=nlat+1)
c
      do l = 1,niv
         do i = 1,nlon
            x(i,     0,l)=x(ip(i),     2,l)
            x(i,nlat+1,l)=x(ip(i),nlat-1,l)
         end do
      end do
c
c  direction zonale (periodicite)
c
      do l = 1,niv
         do k = 0,nlat+1
            x(0    ,k,l)=x(nlon,k,l)
            x(nlon+1,k,l)=x(1  ,k,l)
            x(nlon+2,k,l)=x(2  ,k,l)
         end do
      end do
c
c  direction verticale
c
      do i = 0,nlon+2
         do k = 0,nlat+1
            x(i,k,0    )=x(i,k,1)
            x(i,k,niv+1)=x(i,k,niv)
         end do
      end do
c
c  interpolations proprement dites
c
c  boucle sur les niveaux
c
      do 15 ls=-1,2
c
c  boucle sur les latitudes
c
         do 150 ks=-1,2
c
c  interpolation selon x
c
            do 1500 n=1,lpro
c
               i=ig0(n)
               k=kg0(n)+ks
               l=lg0(n)+ls
c
               xm1=x(i-1,k,l)
               xp0=x(i  ,k,l)
               xp1=x(i+1,k,l)
               xp2=x(i+2,k,l)
c
               yy=alfxm1(n)*xm1+alfxp0(n)*xp0
     +           +alfxp1(n)*xp1+alfxp2(n)*xp2
c
c  limitation pour eviter l'apparition de nouveaux extremas
c
               rmx=max(xp0,xp1)
               rmn=min(xp0,xp1)
               yy=max(min(yy,rmx),rmn)
c
               y(n,ks)=yy
c
 1500       continue
c
 150     continue
c
c  interpolation selon y
c
         do 151 n=1,lpro
c
            ym1=y(n,-1)
            yp0=y(n, 0)
            yp1=y(n, 1)
            yp2=y(n, 2)
c
            zz=alfym1(n)*ym1+alfyp0(n)*yp0
     +        +alfyp1(n)*yp1+alfyp2(n)*yp2
c
c  limitation pour eviter l'apparition de nouveaux extremas
c
            rmx=max(yp0,yp1)
            rmn=min(yp0,yp1)
            zz=max(min(zz,rmx),rmn)
c
            z(n,ls)=zz
c
 151     continue
c
 15   continue
c
c  interpolation selon z
c
      do 16 n=1,lpro
c
         zm1=z(n,-1)
         zp0=z(n, 0)
         zp1=z(n, 1)
         zp2=z(n, 2)
c
         zz=alfzm1(n)*zm1+alfzp0(n)*zp0
     +     +alfzp1(n)*zp1+alfzp2(n)*zp2
c
c  limitation pour eviter l'apparition de nouveaux extremas
c
         rmx=max(zp0,zp1)
         rmn=min(zp0,zp1)
         tempo(n)=max(min(zz,rmx),rmn)
c
 16   continue
c
      do l = 1,niv
         do k = 1,nlat
            do i = 1,nlon
               ikl = i+(k-1)*nlon+(l-1)*lola
               qj1(i,k,l,jc) = tempo(ikl)
            end do
         end do
      end do
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine profils(namexp,
     $                   daynum,
     $                   ian,
     $                   imois,
     $                   ijour,
     $                   iheure,
     $                   imin,
     $                   alt,
     $                   trajpre,
     $                   trajtem,
     $                   trajlon,
     $                   trajlat,
     $                   qj1,
     $                   hc,
     $                   irapp)
c
      parameter(nlon = 180, nlat = 91, niv = 137, nivbas = 94)
      parameter(nrappmax = 12)
      parameter(nbcon = 43, ncm = 15)
      parameter(nsta = 13)
c
      real daynum
      real trajpre(nlon,nlat,niv,nrappmax)
      real trajtem(nlon,nlat,niv,nrappmax)
      real trajlon(nlon,nlat,niv,nrappmax)
      real trajlat(nlon,nlat,niv,nrappmax)
      real alt(nlon,nlat,niv)
      real qj1(nlon,nlat,niv,nbcon), long(nbcon,niv)
      real hc(nlon,nlat,niv,ncm), short(ncm,niv)
      real sza(nlon,nivbas), hnm(niv), theta(niv)
      real pmb(niv)
      real xlas(nsta), xlos(nsta)
c
      character*6  namexp
      character*28 lab2
      character*12 lab3d(nbcon)
      character*12 labsh(ncm)
      character*10 nommois(12)
      character*10 mois
      character*10 nomsta(nsta)
c
      data nomsta/ 'KIRUNA    ','OHP       ',
     $             'NYAALESUND','SODANKYLA ',
     $             'YAKUTSK   ','DDU       ',
     $             'HARESTUA  ','MARAMBIO  ',
     $             'SOUTH POLE','AIRESADOUR',
     $             'EUREKA    ','NIAMEY    ',
     $             'TERESINA  '/
c
      data xlas  / 67.80,  44.00,  78.93,  67.40,  62.03, -66.70,
     $             60.20, -64.27, -90.00,  43.70, 79.59, 13.48,
     $             -6.00/
      data xlos  / 20.40,  05.75,  11.95,  26.60, 129.62, 140.01,
     $             10.75, 303.28,  00.00,  00.25, 274.52, 2.16,
     $             318.00/
c
      data nommois/'  JANUARY ',' FEBRUARY ','   MARCH  ','   APRIL  ',
     +             '    MAY   ','   JUNE   ','   JULY   ','  AUGUST  ',
     +             ' SEPTEMBER','  OCTOBER ',' NOVEMBER ',' DECEMBER '/
c
      data lab3d/'N2O         ','CH4         ','H2O         ',
     +           'NOy         ','HNO3        ','N2O5        ',
     +           'Cly         ','Ox          ','CO          ',
     +           'OClO        ','Passive Ox  ','H2SO4       ',
     +           'HCl         ','ClONO2      ','HOCl        ',
     +           'Cl2         ','H2O2        ','ClNO2       ',
     +           'HBr         ','BrONO2      ','NOx         ',
     +           'HNO4        ','ClOx        ','BrOx        ',
     +           'Cl2O2       ','HOBr        ','BrCl        ',
     +           'CH2O        ','CH3O2       ','CH3O2H      ',
     +           'CFC-11      ','CFC-12      ','CFC-113     ',
     +           'CCl4        ','CH3CCl3*    ','CH3Cl       ',
     +           'HCFC-22*    ','CH3Br       ','H-1211*     ',
     +           'H-1301      ','Bry         ','CH2Br2*     ',
     +           'HNO3 GAS    '/
c
      data labsh/'O(1D)       ','OH          ','Cl          ',
     +           'O(3P)       ','O3          ','HO2         ',
     +           'NO2         ','NO          ','Br          ',
     +           'N           ','ClO         ','BrO         ',
     +           'NO3         ','H           ','CH3         '/
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     mise en forme de la date
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      mois = nommois(imois)
      write(lab2,13)ijour,mois,ian,iheure,imin,'UT'
c
 13   format(i2,a12,i6,2x,2i2,a2)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     boucle sur les stations
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do ista = 1,nsta
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     determination des indices de latitude et longitude
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      dlon = 360./real(nlon)
      dlat = 180./real(nlat - 1)
c
      ilon = nint(xlos(ista)/dlon) + 1
      if (ilon .gt. nlon) then
         ilon = 1
      end if
      ilat = nint((90. - xlas(ista))/dlat) + 1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calcul de l'angle zenithal pour chaque station
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call zenith(trajlon, trajlat, irapp, ilat, daynum, sza)
      szasta = sza(ilon,1)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ecriture d'informations utiles
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      write(59+ista,332) 'RUN ',namexp
      write(59+ista,*)
      write(59+ista,333)' REPROBUS JULIAN DAY = ',daynum
      write(59+ista,335) lab2
      write(59+ista,*)
      write(59+ista,115) nomsta(ista)
      write(59+ista,*)
      write(59+ista,112)'LONGITUDE = ',xlos(ista),
     $                  ' LATITUDE = ',xlas(ista)
      write(59+ista,114)'ILON      = ',ilon,
     $                  ' ILAT     = ',ilat
      write(59+ista,*)
      write(59+ista,113)'SOLAR ZENITH ANGLE = ',szasta
      write(59+ista,*)
c
 112  format(1x,a12,f8.3,a12,f8.3)
 113  format(1x,a21,f7.2)
 114  format(1x,a12,i8,a12,i8)
 115  format(1x,a10)
 332  format(1x,a4,a6)
 333  format(a23,f8.4)
 335  format(a28)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calcul de la concentration et de theta
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do iniv = 1,niv
         pmb(iniv)   =  trajpre(ilon,ilat,iniv,irapp)
         hnm(iniv)   =  trajpre(ilon,ilat,iniv,irapp)
     $                  /(trajtem(ilon,ilat,iniv,irapp)*1.38e-19)
         theta(iniv) =  trajtem(ilon,ilat,iniv,irapp)
     $                  *(1000.
     $                  /trajpre(ilon,ilat,iniv,irapp))**(2./7.)
      end do
      do iniv = 1,niv
         do is = 1,nbcon
            long(is,iniv)  = qj1(ilon,ilat,iniv,is)
         end do
         do is = 1,ncm
            short(is,iniv) = hc(ilon,ilat,iniv,is)
         end do
      end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ecriture proprement dite, en colonnes
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      write(59+ista,1000)'ALTITUDE    ',
     $                   'PRESSURE    ',
     $                   'TEMPERATURE ',
     $                   'THETA       ',
     $                   'DENSITY     ',
     $                   (lab3d(is),is = 1,1)
      write(59+ista,*)
      do iniv = 1,niv
         write(59+ista,1001) alt(ilon,ilat,iniv),
     $                       pmb(iniv),
     $                       trajtem(ilon,ilat,iniv,irapp),
     $                       theta(iniv),
     $                       hnm(iniv),
     $                       (long(is,iniv),is=1,1)
      end do
      write(59+ista,*)
      write(59+ista,1000)'PRESSURE    ',(lab3d(is),is = 2,6)
      write(59+ista,*)
      do iniv = 1,niv
         write(59+ista,1001) pmb(iniv),
     $                       (long(is,iniv),is=2,6)
      end do
      write(59+ista,*)
      write(59+ista,1000)'PRESSURE    ',(lab3d(is),is = 7,11)
      write(59+ista,*)
      do iniv = 1,niv
         write(59+ista,1001) pmb(iniv),
     $                       (long(is,iniv),is=7,11) 
      end do
      write(59+ista,*)
      write(59+ista,1000)'PRESSURE    ',(lab3d(is),is = 12,16)
      write(59+ista,*)
      do iniv = 1,niv
         write(59+ista,1001) pmb(iniv),
     $                       (long(is,iniv),is=12,16)
      end do
      write(59+ista,*)
      write(59+ista,1000)'PRESSURE    ',(lab3d(is),is = 17,21)
      write(59+ista,*)
      do iniv = 1,niv
         write(59+ista,1001) pmb(iniv),
     $                       (long(is,iniv),is=17,21)
      end do
      write(59+ista,*)
      write(59+ista,1000)'PRESSURE    ',(lab3d(is),is = 22,26)
      write(59+ista,*)
      do iniv = 1,niv
         write(59+ista,1001) pmb(iniv),
     $                       (long(is,iniv),is=22,26)
      end do
      write(59+ista,*)
      write(59+ista,1000)'PRESSURE    ',(lab3d(is),is = 27,31)
      write(59+ista,*)
      do iniv = 1,niv
         write(59+ista,1001) pmb(iniv),
     $                       (long(is,iniv),is=27,31)
      end do
      write(59+ista,*)
      write(59+ista,1000)'PRESSURE    ',(lab3d(is),is = 32,36)
      write(59+ista,*)
      do iniv = 1,niv
         write(59+ista,1001) pmb(iniv),
     $                       (long(is,iniv),is=32,36)
      end do
      write(59+ista,*)
      write(59+ista,1000)'PRESSURE    ',(lab3d(is),is = 37,41)
      write(59+ista,*)
      do iniv = 1,niv
         write(59+ista,1001) pmb(iniv),
     $                       (long(is,iniv),is=37,41)
      end do
      write(59+ista,*)
      write(59+ista,1000)'PRESSURE    ',(lab3d(is),is = 42,43)
      write(59+ista,*)
      do iniv = 1,niv
         write(59+ista,1001) pmb(iniv),
     $                       (long(is,iniv),is=42,43)
      end do
      write(59+ista,*)
      write(59+ista,1000)'PRESSURE    ',(labsh(is),is = 1,5)
      write(59+ista,*)
      do iniv = 1,niv
         write(59+ista,1001) pmb(iniv),
     $                       (short(is,iniv),is=1,5)
      end do
      write(59+ista,*)
      write(59+ista,1000)'PRESSURE    ',(labsh(is),is = 6,10)
      write(59+ista,*)
      do iniv = 1,niv
         write(59+ista,1001) pmb(iniv),
     $                       (short(is,iniv),is=6,10)
      end do
      write(59+ista,*)
      write(59+ista,1000)'PRESSURE    ',(labsh(is),is = 11,15)
      write(59+ista,*)
      do iniv = 1,niv
         write(59+ista,1001) pmb(iniv),
     $                       (short(is,iniv),is=11,15)
      end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     fin de boucle sur les stations
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      end do
c
 1000 format(2x,6a12)
 1001 format(6e12.4)
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine martyn(namexp, daynum, tj1, qj1, hc)
c
      parameter(nlon = 180, nlat = 91, niv = 137)
      parameter(lolani = nlon*nlat*niv)
      parameter(nbcon = 43, ncm = 15)
      parameter(nisotheta = 6)
      parameter(nlatbin   = 7)
c
      common /t21l30/ pmb(nlon,nlat,niv), xlat(nlat)
c
      real tj1(nlon,nlat,niv)
      real qj1(nlon,nlat,niv,nbcon)
      real hc(nlon,nlat,niv,ncm)
      real theta(nlon,nlat,niv)
      real isotheta(nisotheta)
      real xlatbin(nlatbin + 1)
      real cinf(nlon,nlat), csup(nlon,nlat)
      real tj1_iso(nlat,nisotheta)
      real qj1_iso(nlat,nisotheta,nbcon)
      real hc_iso(nlat,nisotheta,ncm)
      real tj1_bin(nlatbin,nisotheta)
      real qj1_bin(nlatbin,nisotheta,nbcon)
      real hc_bin(nlatbin,nisotheta,ncm)
      real nl
c
      integer index(nlon,nlat)
      character*6 namexp
c
      data isotheta/694., 550., 475., 430., 410., 380./
      data xlatbin/90., 80., 70., 60., 50., 40., 30., 20./
c
      rascp = 2./7.
      p0 = 1000.
c
c     calcul de theta a chaque point de grille
c
      do ilon = 1,lolani
         theta(ilon,1,1) = tj1(ilon,1,1)*((p0/pmb(ilon,1,1))**rascp)
      end do
c
      do i = 1,nisotheta
c
c        detection des niveaux encadrant isotheta, et calcul des poids
c
         do ilat = 1,nlat
            do ilon = 1,nlon
               do iniv = 2,niv
                  if(theta(ilon,ilat,iniv) .lt. isotheta(i)) then
                     index(ilon,ilat) = iniv
                     cinf(ilon,ilat) =
     $               (isotheta(i) - theta(ilon,ilat,iniv))
     $              /(theta(ilon,ilat,iniv-1)-theta(ilon,ilat,iniv))
                     csup(ilon,ilat) = 1. - cinf(ilon,ilat)
                     goto 1000
                  end if
               end do
 1000          continue
            end do
         end do
c
c        interpolation sur le niveau isotheta et moyenne zonale
c
c        temperature
c
         do ilat = 1,nlat
            tj1_iso(ilat,i) = 0.
            do ilon = 1,nlon
               iniv = index(ilon,ilat)
               tj1_iso(ilat,i) = tj1_iso(ilat,i)
     $                      +(cinf(ilon,ilat)*tj1(ilon,ilat,iniv-1)
     $                      + csup(ilon,ilat)*tj1(ilon,ilat,iniv) )
     $                      /real(nlon)
            end do
         end do
c
c        especes tranportees
c
         do is = 1,nbcon
            do ilat = 1,nlat
               qj1_iso(ilat,i,is) = 0.
               do ilon = 1,nlon
                  iniv = index(ilon,ilat)
                  qj1_iso(ilat,i,is) = qj1_iso(ilat,i,is)
     $                      +(cinf(ilon,ilat)*qj1(ilon,ilat,iniv-1,is)
     $                      + csup(ilon,ilat)*qj1(ilon,ilat,iniv  ,is))
     $                      /real(nlon)
               end do
            end do
         end do
c
c        especes a l'equilibre
c
         do is = 1,ncm
            do ilat = 1,nlat
               hc_iso(ilat,i,is) = 0.
               do ilon = 1,nlon
                  iniv = index(ilon,ilat)
                  hc_iso(ilat,i,is) = hc_iso(ilat,i,is)
     $                       +(cinf(ilon,ilat)*hc(ilon,ilat,iniv-1,is)
     $                       + csup(ilon,ilat)*hc(ilon,ilat,iniv  ,is))
     $                       /real(nlon)
               end do
            end do
         end do
c
c        moyenne sur les latitudes bins
c
         do ilatbin = 1,nlatbin
            nl = (xlatbin(ilatbin) - xlatbin(ilatbin + 1))
     $           /2. + 1.
c
c           temperature
c
            tj1_bin(ilatbin,i) = 0.
            do ilat = 1,nlat
               if (xlat(ilat) .le. xlatbin(ilatbin) .and.
     $             xlat(ilat) .ge. xlatbin(ilatbin + 1) ) then
                  tj1_bin(ilatbin,i) = tj1_bin(ilatbin,i)
     $                               + tj1_iso(ilat,i)/nl
               end if 
            end do
c
c           especes transportees
c
            do is = 1,nbcon
               qj1_bin(ilatbin,i,is) = 0.
               do ilat = 1,nlat
                  if (xlat(ilat) .le. xlatbin(ilatbin) .and.
     $                xlat(ilat) .ge. xlatbin(ilatbin + 1) ) then
                     qj1_bin(ilatbin,i,is) = 
     $                                  qj1_bin(ilatbin,i,is)
     $                                + qj1_iso(ilat,i,is)/nl
                  end if 
               end do
            end do
c
c           especes a l'equilibre
c
            do is = 1,ncm
               hc_bin(ilatbin,i,is) = 0.
               do ilat = 1,nlat
                  if (xlat(ilat) .le. xlatbin(ilatbin) .and.
     $                xlat(ilat) .ge. xlatbin(ilatbin + 1) ) then
                     hc_bin(ilatbin,i,is) = 
     $                                  hc_bin(ilatbin,i,is)
     $                                + hc_iso(ilat,i,is)/nl
                  end if 
               end do
            end do
         end do
c
      end do
c
c     ecriture au format commun
c
      ntp = 97
      day = daynum/86400. - 1.
      if (day .ge. 367.) then 
         day = day - 366.
      end if
c
      write(ntp,1012) (isotheta(i), i = 1,nisotheta)
      write(ntp,1012) (xlatbin(i), i = 1,nlatbin + 1)
      write(ntp,1010) day, namexp
      write(ntp,1012) ((tj1_bin(ilat,iniv), 
     $                             ilat = 1,nlatbin),
     $                             iniv = 1,nisotheta)
      write(ntp,1012) ((qj1_bin(ilat,iniv,8), 
     $                             ilat = 1,nlatbin),
     $                             iniv = 1,nisotheta)
      write(ntp,1012) ((qj1_bin(ilat,iniv,11), 
     $                             ilat = 1,nlatbin),
     $                             iniv = 1,nisotheta)
      write(ntp,1012) ((qj1_bin(ilat,iniv,1), 
     $                             ilat = 1,nlatbin),
     $                             iniv = 1,nisotheta)
      write(ntp,1012) ((qj1_bin(ilat,iniv,2), 
     $                             ilat = 1,nlatbin),
     $                             iniv = 1,nisotheta)
      write(ntp,1012) ((qj1_bin(ilat,iniv,3), 
     $                             ilat = 1,nlatbin),
     $                             iniv = 1,nisotheta)
      write(ntp,1012) ((qj1_bin(ilat,iniv,4), 
     $                             ilat = 1,nlatbin),
     $                             iniv = 1,nisotheta)
      write(ntp,1012) ((qj1_bin(ilat,iniv,5), 
     $                             ilat = 1,nlatbin),
     $                             iniv = 1,nisotheta)
      write(ntp,1012) ((qj1_bin(ilat,iniv,6), 
     $                             ilat = 1,nlatbin),
     $                             iniv = 1,nisotheta)
      write(ntp,1012) ((qj1_bin(ilat,iniv,21), 
     $                             ilat = 1,nlatbin),
     $                             iniv = 1,nisotheta)
      write(ntp,1012) ((qj1_bin(ilat,iniv,7), 
     $                             ilat = 1,nlatbin),
     $                             iniv = 1,nisotheta)
      write(ntp,1012) ((qj1_bin(ilat,iniv,13), 
     $                             ilat = 1,nlatbin),
     $                             iniv = 1,nisotheta)
      write(ntp,1012) ((qj1_bin(ilat,iniv,14), 
     $                             ilat = 1,nlatbin),
     $                             iniv = 1,nisotheta)
      write(ntp,1012) ((hc_bin(ilat,iniv,11), 
     $                             ilat = 1,nlatbin),
     $                             iniv = 1,nisotheta)
      write(ntp,1012) ((qj1_bin(ilat,iniv,25), 
     $                             ilat = 1,nlatbin),
     $                             iniv = 1,nisotheta)
      write(ntp,1012) ((qj1_bin(ilat,iniv,9), 
     $                             ilat = 1,nlatbin),
     $                             iniv = 1,nisotheta)
      write(ntp,1012) ((hc_bin(ilat,iniv,2), 
     $                             ilat = 1,nlatbin),
     $                             iniv = 1,nisotheta)
      write(ntp,1012) ((hc_bin(ilat,iniv,6), 
     $                             ilat = 1,nlatbin),
     $                             iniv = 1,nisotheta)
      write(ntp,1012) ((qj1_bin(ilat,iniv,17), 
     $                             ilat = 1,nlatbin),
     $                             iniv = 1,nisotheta)
      write(ntp,1012) ((hc_bin(ilat,iniv,12), 
     $                             ilat = 1,nlatbin),
     $                             iniv = 1,nisotheta)
      write(ntp,1012) ((qj1_bin(ilat,iniv,12), 
     $                             ilat = 1,nlatbin),
     $                             iniv = 1,nisotheta)
c
 1010 format(f10.4,a10)
 1012 format(1p,8e10.3)
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine internew(qj1,jc)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                              c
c  sous-programme d'interpolation pour les champs de traceurs  c
c  interpolateur de lagrange: bicubique sur l'horizontale      c
c  (ecriture valable pour une grille horizontale reguliere )   c
c                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  parametres du modele de transport
c
      parameter(nlon = 180, nlat = 91, niv = 137)
      parameter(nbcon = 43)
c
      parameter(lpro=nlon*nlat*niv)
      parameter(lola=nlon*nlat)
c
      common/grille/dlon(0:nlon+2),ddlon,dlat(0:nlat+1),ddlat
     +             ,slon(nlon),clon(nlon),slat(nlat),clat(nlat)
     +             ,slonr(nlon,nlat),clonr(nlon,nlat)
     +             ,slatr(nlon,nlat),clatr(nlon,nlat)
     +             ,xcoef(-1:2),ycoef(-1:2),ip(nlon)
c
      common/poids/ alfxm1(lpro),alfxp0(lpro),alfxp1(lpro),
     +              alfxp2(lpro),alfym1(lpro),alfyp0(lpro),
     +              alfyp1(lpro),alfyp2(lpro),alfzm1(lpro),
     +              alfzp0(lpro),alfzp1(lpro),alfzp2(lpro)
c
      common/indices/ ig0(lpro),kg0(lpro),lg0(lpro)
c
      real x(0:nlon+2,0:nlat+1,0:niv+1)
      real y(lpro,-1:2,-1:2),z(lpro,-1:2)
      real rmx(lpro), rmn(lpro)
c
      real qj1(nlon,nlat,niv,nbcon)
      real tempo(lpro)
c
c  extension de la grille pour les interpolations cubiques
c
      do l = 1,niv
         do k = 1,nlat
            do i = 1,nlon
               x(i,k,l)=qj1(i,k,l,jc)
            end do
         end do
      end do
c
c  direction meridienne: symetrie par rapport au pole nord (k=0)
c                      + symetrie par rapport au pole sud  (k=nlat+1)
c
      do l = 1,niv
         do i = 1,nlon
            x(i,     0,l)=x(ip(i),     2,l)
            x(i,nlat+1,l)=x(ip(i),nlat-1,l)
         end do
      end do
c
c  direction zonale (periodicite)
c
      do l = 1,niv
         do k = 0,nlat+1
            x(0    ,k,l)=x(nlon,k,l)
            x(nlon+1,k,l)=x(1  ,k,l)
            x(nlon+2,k,l)=x(2  ,k,l)
         end do
      end do
c
c  direction verticale
c
      do i = 0,nlon+2
         do k = 0,nlat+1
            x(i,k,0    )=x(i,k,1)
            x(i,k,niv+1)=x(i,k,niv)
         end do 
      end do
c
c  interpolations proprement dites
c
c  interpolation selon x                                                               
c                                                                                      
      do 2 n=1,lpro                                           
c                                                                                      
         i=ig0(n)                                             
         k=kg0(n)                                             
         l=lg0(n)                                             
c                                                                                      
         im1=i-1                                              
         ip0=i                                                
         ip1=i+1                                              
         ip2=i+2                                              
c                                                                                      
         km1=k-1                                              
         kp0=k                                                
         kp1=k+1                                              
         kp2=k+2                                              
c                                                                                      
         lm1=max(l-1,0)                                       
         lp0=l                                                
         lp1=l+1                                              
         lp2=min(l+2,niv+1)                                   
c                                                                                      
         xm1m1m1=x(im1,km1,lm1)                                
         xp0m1m1=x(ip0,km1,lm1)                                     
         xp1m1m1=x(ip1,km1,lm1)                                      
         xp2m1m1=x(ip2,km1,lm1)                                      
c                                                                                      
         xm1p0m1=x(im1,kp0,lm1)                                      
         xp0p0m1=x(ip0,kp0,lm1)                                     
         xp1p0m1=x(ip1,kp0,lm1)                                      
         xp2p0m1=x(ip2,kp0,lm1)                                      
c                                                                                      
         xm1p1m1=x(im1,kp1,lm1)                                      
         xp0p1m1=x(ip0,kp1,lm1)                                    
         xp1p1m1=x(ip1,kp1,lm1)                                     
         xp2p1m1=x(ip2,kp1,lm1)                                     
c                                                                                     
         xm1p2m1=x(im1,kp2,lm1)                                     
         xp0p2m1=x(ip0,kp2,lm1)                                    
         xp1p2m1=x(ip1,kp2,lm1)                                     
         xp2p2m1=x(ip2,kp2,lm1)                                     
c                                                                                     
         xm1m1p0=x(im1,km1,lp0)                                 
         xp0m1p0=x(ip0,km1,lp0)                                    
         xp1m1p0=x(ip1,km1,lp0)                                     
         xp2m1p0=x(ip2,km1,lp0)                                     
c                                                                                     
         xm1p0p0=x(im1,kp0,lp0)                                     
         xp0p0p0=x(ip0,kp0,lp0)                                    
         xp1p0p0=x(ip1,kp0,lp0)                                     
         xp2p0p0=x(ip2,kp0,lp0)                                     
c                                                                                     
         xm1p1p0=x(im1,kp1,lp0)                                     
         xp0p1p0=x(ip0,kp1,lp0)                                    
         xp1p1p0=x(ip1,kp1,lp0)                                     
         xp2p1p0=x(ip2,kp1,lp0)                                     
c                                                                                     
         xm1p2p0=x(im1,kp2,lp0)                                     
         xp0p2p0=x(ip0,kp2,lp0)                                    
         xp1p2p0=x(ip1,kp2,lp0)                                     
         xp2p2p0=x(ip2,kp2,lp0)                                     
c                                                                                     
         xm1m1p1=x(im1,km1,lp1)                                  
         xp0m1p1=x(ip0,km1,lp1)                                    
         xp1m1p1=x(ip1,km1,lp1)                                     
         xp2m1p1=x(ip2,km1,lp1)                                     
c                                                                                     
         xm1p0p1=x(im1,kp0,lp1)                                     
         xp0p0p1=x(ip0,kp0,lp1)                                    
         xp1p0p1=x(ip1,kp0,lp1)                                     
         xp2p0p1=x(ip2,kp0,lp1)                                     
c                                                                                     
         xm1p1p1=x(im1,kp1,lp1)                                     
         xp0p1p1=x(ip0,kp1,lp1)                                    
         xp1p1p1=x(ip1,kp1,lp1)                                     
         xp2p1p1=x(ip2,kp1,lp1)                                     
c                                                                                     
         xm1p2p1=x(im1,kp2,lp1)                                     
         xp0p2p1=x(ip0,kp2,lp1)                                    
         xp1p2p1=x(ip1,kp2,lp1)                                     
         xp2p2p1=x(ip2,kp2,lp1)                                     
c                                                                                     
         xm1m1p2=x(im1,km1,lp2)                                   
         xp0m1p2=x(ip0,km1,lp2)                                    
         xp1m1p2=x(ip1,km1,lp2)                                     
         xp2m1p2=x(ip2,km1,lp2)                                     
c                                                                                     
         xm1p0p2=x(im1,kp0,lp2)                                     
         xp0p0p2=x(ip0,kp0,lp2)                                    
         xp1p0p2=x(ip1,kp0,lp2)                                     
         xp2p0p2=x(ip2,kp0,lp2)                                     
c                                                                                     
         xm1p1p2=x(im1,kp1,lp2)                                     
         xp0p1p2=x(ip0,kp1,lp2)                                    
         xp1p1p2=x(ip1,kp1,lp2)                                     
         xp2p1p2=x(ip2,kp1,lp2)                                     
c                                                                                     
         xm1p2p2=x(im1,kp2,lp2)                                     
         xp0p2p2=x(ip0,kp2,lp2)                                    
         xp1p2p2=x(ip1,kp2,lp2)                                     
         xp2p2p2=x(ip2,kp2,lp2)                                     
c                                                                                     
         y(n,-1,-1)=alfxm1(n)*xm1m1m1+alfxp0(n)*xp0m1m1             
     +             +alfxp1(n)*xp1m1m1+alfxp2(n)*xp2m1m1              
         y(n, 0,-1)=alfxm1(n)*xm1p0m1+alfxp0(n)*xp0p0m1               
     +             +alfxp1(n)*xp1p0m1+alfxp2(n)*xp2p0m1                
         y(n, 1,-1)=alfxm1(n)*xm1p1m1+alfxp0(n)*xp0p1m1          
     +             +alfxp1(n)*xp1p1m1+alfxp2(n)*xp2p1m1           
         y(n, 2,-1)=alfxm1(n)*xm1p2m1+alfxp0(n)*xp0p2m1            
     +             +alfxp1(n)*xp1p2m1+alfxp2(n)*xp2p2m1             
c                                                                                     
         y(n,-1, 0)=alfxm1(n)*xm1m1p0+alfxp0(n)*xp0m1p0          
     +             +alfxp1(n)*xp1m1p0+alfxp2(n)*xp2m1p0           
         y(n, 0, 0)=alfxm1(n)*xm1p0p0+alfxp0(n)*xp0p0p0            
     +             +alfxp1(n)*xp1p0p0+alfxp2(n)*xp2p0p0           
         y(n, 1, 0)=alfxm1(n)*xm1p1p0+alfxp0(n)*xp0p1p0       
     +             +alfxp1(n)*xp1p1p0+alfxp2(n)*xp2p1p0           
         y(n, 2, 0)=alfxm1(n)*xm1p2p0+alfxp0(n)*xp0p2p0           
     +             +alfxp1(n)*xp1p2p0+alfxp2(n)*xp2p2p0           
c                                                                                     
         y(n,-1, 1)=alfxm1(n)*xm1m1p1+alfxp0(n)*xp0m1p1        
     +             +alfxp1(n)*xp1m1p1+alfxp2(n)*xp2m1p1           
         y(n, 0, 1)=alfxm1(n)*xm1p0p1+alfxp0(n)*xp0p0p1         
     +             +alfxp1(n)*xp1p0p1+alfxp2(n)*xp2p0p1           
         y(n, 1, 1)=alfxm1(n)*xm1p1p1+alfxp0(n)*xp0p1p1          
     +             +alfxp1(n)*xp1p1p1+alfxp2(n)*xp2p1p1           
         y(n, 2, 1)=alfxm1(n)*xm1p2p1+alfxp0(n)*xp0p2p1           
     +             +alfxp1(n)*xp1p2p1+alfxp2(n)*xp2p2p1           
c                                                                                     
         y(n,-1, 2)=alfxm1(n)*xm1m1p2+alfxp0(n)*xp0m1p2            
     +             +alfxp1(n)*xp1m1p2+alfxp2(n)*xp2m1p2             
         y(n, 0, 2)=alfxm1(n)*xm1p0p2+alfxp0(n)*xp0p0p2              
     +             +alfxp1(n)*xp1p0p2+alfxp2(n)*xp2p0p2           
         y(n, 1, 2)=alfxm1(n)*xm1p1p2+alfxp0(n)*xp0p1p2         
     +             +alfxp1(n)*xp1p1p2+alfxp2(n)*xp2p1p2           
         y(n, 2, 2)=alfxm1(n)*xm1p2p2+alfxp0(n)*xp0p2p2           
     +             +alfxp1(n)*xp1p2p2+alfxp2(n)*xp2p2p2           
c                                                                                     
c     limitation pour eviter l'apparition de nouveaux extrema
c                                                                                     
         rp0p0mx=max(xp0p0p0,xp0p1p0)                           
         rp0p0mn=min(xp0p0p0,xp0p1p0)                            
c                                                                                     
         rp1p0mx=max(xp1p0p0,xp1p1p0)                             
         rp1p0mn=min(xp1p0p0,xp1p1p0)                         
c                                                                                     
         rp0p1mx=max(xp0p0p1,xp0p1p1)                         
         rp0p1mn=min(xp0p0p1,xp0p1p1)                         
c                                                                                     
         rp1p1mx=max(xp1p0p1,xp1p1p1)                         
         rp1p1mn=min(xp1p0p1,xp1p1p1)                         
c                                                                                     
         rp0mx=max(rp0p0mx,rp1p0mx)                           
         rp0mn=min(rp0p0mn,rp1p0mn)                           
c                                                                                     
         rp1mx=max(rp0p1mx,rp1p1mx)                           
         rp1mn=min(rp0p1mn,rp1p1mn)                           
c                                                                 
         rmx(n)=max(rp0mx,rp1mx)                              
         rmn(n)=min(rp0mn,rp1mn)                              
c                                                                                     
c                                                                                     
 2    continue                                                    
c                                                                                     
c  interpolation selon y                                                              
c                                                                                     
      do 3 n=1,lpro                                               
         z(n,-1)=alfym1(n)*y(n,-1,-1)+alfyp0(n)*y(n, 0,-1)        
     +          +alfyp1(n)*y(n, 1,-1)+alfyp2(n)*y(n, 2,-1)        
         z(n, 0)=alfym1(n)*y(n,-1, 0)+alfyp0(n)*y(n, 0, 0)        
     +          +alfyp1(n)*y(n, 1, 0)+alfyp2(n)*y(n, 2, 0)        
         z(n, 1)=alfym1(n)*y(n,-1, 1)+alfyp0(n)*y(n, 0, 1)        
     +          +alfyp1(n)*y(n, 1, 1)+alfyp2(n)*y(n, 2, 1)        
         z(n, 2)=alfym1(n)*y(n,-1, 2)+alfyp0(n)*y(n, 0, 2)        
     +          +alfyp1(n)*y(n, 1, 2)+alfyp2(n)*y(n, 2, 2)        
 3    continue                                                    
c                                                                 
c   interpolation selon z                                                             
c                                                                                     
      do 4 n=1,lpro                                               
         tempo(n)=alfzm1(n)*z(n,-1)+alfzp0(n)*z(n, 0)          
     +           +alfzp1(n)*z(n, 1)+alfzp2(n)*z(n, 2)          
 4    continue                                                    
c                                                                 
c     limitation pour eviter l'apparition de nouveaux extrema
c                                                                 
      do 5 n=1,lpro                                             
         tempo(n)=max(min(tempo(n),rmx(n)),rmn(n))    
 5    continue                                                  
c                                                                 
      do l = 1,niv
         do k = 1,nlat
            do i = 1,nlon
               ikl = i+(k-1)*nlon+(l-1)*lola
               qj1(i,k,l,jc) = tempo(ikl)
            end do
         end do
      end do
c
      return                                                      
      end                                                        
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ecmwfh2o (daynum, qj1,
     $                     nivh2oecmwf, ntime, nintecmwf) 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     traitement de la vapeur d'eau tropospherique:
c     analyses ecmwf a partir du niveau nivh2oecmwf
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      parameter(nlon = 180, nlat = 91, niv = 137)
      parameter(lolani = nlon*nlat*niv)
      parameter(nbcon = 43)
c
      real qj1(nlon,nlat,niv,nbcon)
c
      common/forcj1/uj1(nlon,nlat,niv),vj1(nlon,nlat,niv)
     +             ,wj1(nlon,nlat,niv),tj1(nlon,nlat,niv)
     +             ,hj1(nlon,nlat,niv),pj1(nlon,nlat)
c
      common/forcm1/um1(nlon,nlat,niv),vm1(nlon,nlat,niv)
     +             ,wm1(nlon,nlat,niv),tm1(nlon,nlat,niv)
     +             ,hm1(nlon,nlat,niv),pm1(nlon,nlat),daym1
c
      common/forcp1/up1(nlon,nlat,niv),vp1(nlon,nlat,niv)
     +             ,wp1(nlon,nlat,niv),tp1(nlon,nlat,niv)
     +             ,hp1(nlon,nlat,niv),pp1(nlon,nlat),dayp1
c
      call instant1(daynum,hj1,hm1,daym1,hp1,dayp1,lolani,
     $              ntime,nintecmwf)
c
      do iniv = nivh2oecmwf,niv
         do ilat = 1,nlat
            do ilon = 1,nlon
               qj1(ilon,ilat,iniv,3) = hj1(ilon,ilat,iniv)*28.97/18.
            end do
         end do
      end do
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine avion(namexp, day, ian, imois, ijour,
     $                 iheure, imin, trajlon, trajlat, 
     $                 trajpre, trajtem, irapp, qj1, hc)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      parameter(nlon = 180, nlat = 91, niv = 137)
      parameter(nbcon = 43, ncm = 15)
      parameter(nrappmax = 12)
      parameter(rter = 6371., pi = 3.141592)
c
      real dayavion, xlatavion, xlonavion
      real trajlon(nlon,nlat,niv,nrappmax)
      real trajlat(nlon,nlat,niv,nrappmax)
      real trajpre(nlon,nlat,niv,nrappmax)
      real trajtem(nlon,nlat,niv,nrappmax)
      real qj1(nlon,nlat,niv,nbcon), long(nbcon,niv)
      real hc(nlon,nlat,niv,ncm), short(nbcon,niv)
      real distcurv(niv)
      real pmb(niv), hnm(niv), theta(niv)
c
      integer indlon(niv), indlat(niv)
c
      character*6  namexp
      character*10 nommois(12)
      character*10 mois
      character*28 lab2
      character*12 lab3d(nbcon)
      character*12 labsh(ncm)
c
      data nommois/'  JANUARY ',' FEBRUARY ','   MARCH  ','   APRIL  ',
     +             '    MAY   ','   JUNE   ','   JULY   ','  AUGUST  ',
     +             ' SEPTEMBER','  OCTOBER ',' NOVEMBER ',' DECEMBER '/
c
      data lab3d/'N2O         ','CH4         ','H2O         ',
     +           'NOy         ','HNO3        ','N2O5        ',
     +           'Cly         ','Ox          ','CO          ',
     +           'OClO        ','Passive Ox  ','H2SO4       ',
     +           'HCl         ','ClONO2      ','HOCl        ',
     +           'Cl2         ','H2O2        ','ClNO2       ',
     +           'HBr         ','BrONO2      ','NOx         ',
     +           'HNO4        ','ClOx        ','BrOx        ',
     +           'Cl2O2       ','HOBr        ','BrCl        ',
     +           'CH2O        ','CH3O2       ','CH3O2H      ',
     +           'CFC-11      ','CFC-12      ','CFC-113     ',
     +           'CCl4        ','CH3CCl3*    ','CH3Cl       ',
     +           'HCFC-22*    ','CH3Br       ','H-1211*     ',
     +           'H-1301      ','Bry         ','CH2Br2*     ',
     +           'HNO3 GAS    '/
c
      data labsh/'O(1D)       ','OH          ','Cl          ',
     +           'O(3P)       ','O3          ','HO2         ',
     +           'NO2         ','NO          ','Br          ',
     +           'N           ','ClO         ','BrO         ',
     +           'NO3         ','H           ','CH3         '/
c
      rewind(79)
c
      do iread = 1,1000
         read(79,*,end = 20) dayavion, secavion, 
     $                       xlatavion, xlonavion,
     $                       xpreavion
         dayavion = dayavion - 1.
         delta  = abs(dayavion - day)
         if (delta .le. 1.e-3) goto 10
      end do
 10   continue
c
      write(6,*) 'coincidence dc8 pour xlat = ',xlatavion,
     $           ' et xlon = ',xlonavion
c
      xavion = rter*cos(xlatavion*pi/180.)*cos(xlonavion*pi/180.)
      yavion = rter*cos(xlatavion*pi/180.)*sin(xlonavion*pi/180.)
      zavion = rter*sin(xlatavion*pi/180.)
c
      do iniv = 1,niv
         diststra = 20000.
         do ilat = 1,nlat
            do ilon = 1,nlon
               xtheta = trajlon(ilon,ilat,iniv,irapp)*pi/180.
               xphi   = trajlat(ilon,ilat,iniv,irapp)*pi/180.
               x = rter*cos(xphi)*cos(xtheta)
               y = rter*cos(xphi)*sin(xtheta)
               z = rter*sin(xphi)
               d = sqrt((xavion - x)**2 + (yavion - y)**2
     $                 +(zavion - z)**2)
               if (d .lt. diststra) then
                  diststra = d
                  indlon(iniv)   = ilon
                  indlat(iniv)   = ilat
                  distcurv(iniv) = 2.*rter*asin(diststra/(2.*rter)) 
               end if 
            end do
         end do
c        write(6,*) 'straight minimal distance = ',diststra
c        write(6,*) 'curve    minimal distance = ',distcurv
c        write(6,*) 'pour lon = ',trajlon(indlon,indlat,1,irapp)
c        write(6,*) 'pour lat = ',trajlat(indlon,indlat,1,irapp)
c        write(6,*) 'indlon = ',indlon,' indlat = ',indlat
      end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     mise en forme de la date
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      mois = nommois(imois)
      write(lab2,13)ijour,mois,ian,iheure,imin,'UT'
 13   format(i2,a12,i6,2x,2i2,a2)
      print*,lab2
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ecriture d'informations utiles
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      write(80,332) 'RUN ',namexp
      write(80,*)
      write(80,333)' REPROBUS JULIAN DAY = ',day
      write(80,*)
      write(80,335) lab2
      write(80,*)
      write(80,334)' UT SECONDS = ',int(secavion)
      write(80,*)
      write(80,112)'LONGITUDE DC8 = ',xlonavion,
     $                  ' LATITUDE DC8 = ',xlatavion
      write(80,*)
c
 112  format(1x,a16,f8.3,a16,f8.3)
 332  format(1x,a4,a6)
 333  format(a23,f8.4)
 334  format(a14,i6)
 335  format(1x,a28)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calcul de la concentration et de theta
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do iniv = 1,niv
         pmb(iniv)   =  
     $      trajpre(indlon(iniv),indlat(iniv),iniv,irapp)
         hnm(iniv)   =  
     $      trajpre(indlon(iniv),indlat(iniv),iniv,irapp)
     $    /(trajtem(indlon(iniv),indlat(iniv),iniv,irapp)*1.38e-19)
         theta(iniv) = 
     $      trajtem(indlon(iniv),indlat(iniv),iniv,irapp)*(1000.
     $     /trajpre(indlon(iniv),indlat(iniv),iniv,irapp))**(2./7.)
      end do
c
      do iniv = 1,niv
         do is = 1,nbcon
            long(is,iniv)  = 
     $         qj1(indlon(iniv),indlat(iniv),iniv,is)
         end do
         do is = 1,ncm
            short(is,iniv) = 
     $         hc(indlon(iniv),indlat(iniv),iniv,is)
         end do
      end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ecriture proprement dite, en colonnes
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      write(80,1000)'DIST.(KM)   ',
     $              'PRESS.(HPA) ',
     $              'TEMPERATURE ',
     $              'THETA       ',
     $              'DENSITY     ',
     $              (lab3d(is),is = 1,1)
      write(80,*)
      do iniv = 1,niv
         write(80,1002) distcurv(iniv), pmb(iniv),
     $            trajtem(indlon(iniv),indlat(iniv),iniv,irapp),
     $                  theta(iniv),
     $                  hnm(iniv),
     $                  (long(is,iniv),is=1,1)
      end do
      write(80,*)
      write(80,1000)'PRESS.(HPA) ',(lab3d(is),is = 2,6)
      write(80,*)
      do iniv = 1,niv
         write(80,1001) pmb(iniv),
     $                       (long(is,iniv),is= 2,6)
      end do
      write(80,*)
      write(80,1000)'PRESS.(HPA) ',(lab3d(is),is = 7,11)
      write(80,*)
      do iniv = 1,niv
         write(80,1001) pmb(iniv),
     $                       (long(is,iniv),is=7,11)
      end do
      write(80,*)
      write(80,1000)'PRESS.(HPA) ',(lab3d(is),is = 12,16)
      write(80,*)
      do iniv = 1,niv
         write(80,1001) pmb(iniv),
     $                       (long(is,iniv),is= 12,16)
      end do
      write(80,*)
      write(80,1000)'PRESS.(HPA) ',(lab3d(is),is = 17,21)
      write(80,*)
      do iniv = 1,niv
         write(80,1001) pmb(iniv),
     $                       (long(is,iniv),is= 17,21)
      end do
      write(80,*)
      write(80,1000)'PRESS.(HPA) ',(lab3d(is),is = 22,26)
      write(80,*)
      do iniv = 1,niv
         write(80,1001) pmb(iniv),
     $                       (long(is,iniv),is= 22,26)
      end do
      write(80,*)
      write(80,1000)'PRESS.(HPA) ',(lab3d(is),is = 27,31)
      write(80,*)
      do iniv = 1,niv
         write(80,1001) pmb(iniv),
     $                       (long(is,iniv),is= 27,31)
      end do
      write(80,*)
      write(80,1000)'PRESS.(HPA) ',(lab3d(is),is = 32,36)
      write(80,*)
      do iniv = 1,niv
         write(80,1001) pmb(iniv),
     $                       (long(is,iniv),is= 32,36)
      end do
      write(80,*)
      write(80,1000)'PRESS.(HPA) ',(lab3d(is),is = 37,41)
      write(80,*)
      do iniv = 1,niv
         write(80,1001) pmb(iniv),
     $                       (long(is,iniv),is= 37,41)
      end do
      write(80,*)
      write(80,1000)'PRESS.(HPA) ',(lab3d(is),is = 42,43)
      write(80,*)
      do iniv = 1,niv
         write(80,1001) pmb(iniv),
     $                       (long(is,iniv),is= 42,43)
      end do
      write(80,*)
      write(80,1000)'PRESS.(HPA) ',(labsh(is),is = 1,5)
      write(80,*)
      do iniv = 1,niv
         write(80,1001) pmb(iniv),
     $                       (short(is,iniv),is=  1, 5)
      end do
      write(80,*)
      write(80,1000)'PRESS.(HPA) ',(labsh(is),is = 6,10)
      write(80,*)
      do iniv = 1,niv
         write(80,1001) pmb(iniv),
     $                       (short(is,iniv),is= 6,10)
      end do
      write(80,*)
      write(80,1000)'PRESS.(HPA) ',(labsh(is),is = 11,15)
      write(80,*)
      do iniv = 1,niv
         write(80,1001) pmb(iniv),
     $                       (short(is,iniv),is= 11,15)
      end do
      write(80,*)
c
 1000 format(2x,6a12)
 1001 format(6e12.4)
 1002 format(f12.4,5e12.4)
c
 20   return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine odin(namexp, day, ian, imois, ijour,
     $                iheure, imin, alt, sza3d, trajlon, trajlat, 
     $                trajpre, trajtem, irapp, qj1, hc)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      parameter(nlon = 180, nlat = 91, niv = 137, nivbas = 94)
      parameter(nbcon = 43, ncm = 15)
      parameter(nrappmax = 24)
      parameter(rter = 6371., pi = 3.141592)
c
      real dayavion, xlatavion, xlonavion
      real alt(nlon,nlat,niv), sza3d(nlon,nlat,nivbas)
      real trajlon(nlon,nlat,niv,nrappmax)
      real trajlat(nlon,nlat,niv,nrappmax)
      real trajpre(nlon,nlat,niv,nrappmax)
      real trajtem(nlon,nlat,niv,nrappmax)
      real qj1(nlon,nlat,niv,nbcon), long(nbcon,niv)
      real hc(nlon,nlat,niv,ncm), short(nbcon,niv)
      real distcurv(niv), xlat(niv), xlon(niv)
      real pmb(niv), hnm(niv), theta(niv), altsmr(niv), sza(niv)
c
      integer indlon(niv), indlat(niv)
c
      character*6  namexp
      character*10 nommois(12)
      character*10 mois
      character*28 lab2
      character*12 lab3d(nbcon)
      character*12 labsh(ncm)
c
      data nommois/'  JANUARY ',' FEBRUARY ','   MARCH  ','   APRIL  ',
     +             '    MAY   ','   JUNE   ','   JULY   ','  AUGUST  ',
     +             ' SEPTEMBER','  OCTOBER ',' NOVEMBER ',' DECEMBER '/
c
      data lab3d/'N2O         ','CH4         ','H2O         ',
     +           'NOy         ','HNO3        ','N2O5        ',
     +           'Cly         ','Ox          ','CO          ',
     +           'OClO        ','Passive Ox  ','H2SO4       ',
     +           'HCl         ','ClONO2      ','HOCl        ',
     +           'Cl2         ','H2O2        ','ClNO2       ',
     +           'HBr         ','BrONO2      ','NOx         ',
     +           'HNO4        ','ClOx        ','BrOx        ',
     +           'Cl2O2       ','HOBr        ','BrCl        ',
     +           'CH2O        ','CH3O2       ','CH3O2H      ',
     +           'CFC-11      ','CFC-12      ','CFC-113     ',
     +           'CCl4        ','CH3CCl3*    ','CH3Cl       ',
     +           'HCFC-22*    ','CH3Br       ','H-1211*     ',
     +           'H-1301      ','Bry         ','CH2Br2*     ',
     +           'HNO3 GAS    '/
c
      data labsh/'O(1D)       ','OH          ','Cl          ',
     +           'O(3P)       ','O3          ','HO2         ',
     +           'NO2         ','NO          ','Br          ',
     +           'N           ','ClO         ','BrO         ',
     +           'NO3         ','H           ','CH3         '/
c
      rewind(79)
c
 30   read(79,*,end = 20) dayavion, xlatavion, xlonavion
      delta  = abs(dayavion - day)
      if (delta .le. 5.2e-3) then
         goto 10
      else
         goto 30
      end if
 10   continue
c
      write(6,*) 'coincidence odin pour xlat = ',xlatavion,
     $           ' et xlon = ',xlonavion
c
      xavion = rter*cos(xlatavion*pi/180.)*cos(xlonavion*pi/180.)
      yavion = rter*cos(xlatavion*pi/180.)*sin(xlonavion*pi/180.)
      zavion = rter*sin(xlatavion*pi/180.)
c
      do iniv = 1,niv
         diststra = 20000.
         do ilat = 1,nlat
            do ilon = 1,nlon
               xtheta = trajlon(ilon,ilat,iniv,irapp)*pi/180.
               xphi   = trajlat(ilon,ilat,iniv,irapp)*pi/180.
               x = rter*cos(xphi)*cos(xtheta)
               y = rter*cos(xphi)*sin(xtheta)
               z = rter*sin(xphi)
               d = sqrt((xavion - x)**2 + (yavion - y)**2
     $                 +(zavion - z)**2)
               if (d .lt. diststra) then
                  diststra = d
                  indlon(iniv)   = ilon
                  indlat(iniv)   = ilat
                  distcurv(iniv) = 2.*rter*asin(diststra/(2.*rter)) 
                  xlat(iniv)     = trajlat(ilon,ilat,iniv,irapp)
                  xlon(iniv)     = trajlon(ilon,ilat,iniv,irapp)
               end if 
            end do
         end do
      end do
c     do iniv = 1,niv
c        write(6,*) 'iniv = ',iniv
c        write(6,*) 'curve    minimal distance = ',distcurv(iniv)
c        write(6,*) 'pour lon = ',
c    $               trajlon(indlon(iniv),indlat(iniv),iniv,irapp)
c        write(6,*) 'pour lat = ',
c    $               trajlat(indlon(iniv),indlat(iniv),iniv,irapp)
c        write(6,*) 'indlon = ',indlon(iniv),
c    $              ' indlat = ',indlat(iniv)
c     end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     mise en forme de la date
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      mois = nommois(imois)
      write(lab2,13)ijour,mois,ian,iheure,imin,'UT'
 13   format(i2,a12,i6,2x,2i2,a2)
      print*,lab2
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ecriture d'informations utiles
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      write(80,337) '******************************'
      write(80,332) 'RUN ',namexp
      write(80,337) '******************************'
      write(80,*)
      write(80,333)' REPROBUS JULIAN DAY = ',day
      write(80,333)' ODIN/SMR JULIAN DAY = ',dayavion
      write(80,336)' DIFFERENCE          = ',abs(day - dayavion)
     $                                       *24.*60.,' MN'
      write(80,*)
      write(80,335) lab2
      write(80,*)
      write(80,112)'LONGITUDE SMR = ',xlonavion,
     $                  ' LATITUDE SMR = ',xlatavion
      write(80,*)
c
 112  format(1x,a16,f8.3,a16,f8.3)
 332  format(1x,a4,a6)
 333  format(a23,f8.4)
 334  format(a14,i6)
 335  format(1x,a28)
 336  format(a23,f8.4,a3)
 337  format(a30)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calcul de la concentration et de theta
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do iniv = 1,niv
         pmb(iniv)   =  
     $      trajpre(indlon(iniv),indlat(iniv),iniv,irapp)
         hnm(iniv)   =  
     $      trajpre(indlon(iniv),indlat(iniv),iniv,irapp)
     $    /(trajtem(indlon(iniv),indlat(iniv),iniv,irapp)*1.38e-19)
         theta(iniv) = 
     $      trajtem(indlon(iniv),indlat(iniv),iniv,irapp)*(1000.
     $     /trajpre(indlon(iniv),indlat(iniv),iniv,irapp))**(2./7.)
         altsmr(iniv) = alt(indlon(iniv),indlat(iniv),iniv)
      end do
      do iniv = 1,nivbas
         sza(iniv) = sza3d(indlon(iniv),indlat(iniv),iniv)
      end do
      do iniv = nivbas+1,niv
         sza(iniv) = sza(nivbas)
      end do
c
      do iniv = 1,niv
         do is = 1,nbcon
            long(is,iniv)  = 
     $         qj1(indlon(iniv),indlat(iniv),iniv,is)
         end do
         do is = 1,ncm
            short(is,iniv) = 
     $         hc(indlon(iniv),indlat(iniv),iniv,is)
         end do
      end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ecriture proprement dite, en colonnes
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      write(80,1003)'LATITUDE    ',
     $              'LONGITUDE   ',
     $              'DIST.(KM)   ',
     $              'PRESS.(HPA) ',
     $              'TEMPERATURE ',
     $              'THETA       ',
     $              'DENSITY     ',
     $              'SZA (DEG)   ' 
      write(80,*)
      do iniv = 1,niv
         write(80,1005) xlat(iniv), xlon(iniv),
     $            distcurv(iniv), pmb(iniv),
     $            trajtem(indlon(iniv),indlat(iniv),iniv,irapp),
     $            theta(iniv), hnm(iniv), sza(iniv)
      end do
      write(80,*)
      write(80,1000)'ALTI. (KM)  ',(lab3d(is),is = 1,7)
      write(80,*)
      do iniv = 1,niv
         write(80,1001) altsmr(iniv),
     $                       (long(is,iniv),is= 1,7)
      end do
      write(80,*)
      write(80,1000)'ALTI. (KM)  ',(lab3d(is),is = 8,14)
      write(80,*)
      do iniv = 1,niv
         write(80,1001) altsmr(iniv),
     $                       (long(is,iniv),is= 8,14)
      end do
      write(80,*)
      write(80,1000)'ALTI. (KM)  ',(lab3d(is),is = 15,21)
      write(80,*)
      do iniv = 1,niv
         write(80,1001) altsmr(iniv),
     $                       (long(is,iniv),is= 15,21)
      end do
      write(80,*)
      write(80,1000)'ALTI. (KM)  ',(lab3d(is),is = 22,28)
      write(80,*)
      do iniv = 1,niv
         write(80,1001) altsmr(iniv),
     $                       (long(is,iniv),is= 22,28)
      end do
      write(80,*)
      write(80,1000)'ALTI. (KM)  ',(lab3d(is),is = 29,35)
      write(80,*)
      do iniv = 1,niv
         write(80,1001) altsmr(iniv),
     $                       (long(is,iniv),is= 29,35)
      end do
      write(80,*)
      write(80,1000)'ALTI. (KM)  ',(lab3d(is),is = 36,42)
      write(80,*)
      do iniv = 1,niv
         write(80,1001) altsmr(iniv),
     $                       (long(is,iniv),is= 36,42)
      end do
      write(80,*)
      write(80,1000)'ALTI. (KM)  ',(lab3d(is),is = 43,43)
      write(80,*)
      do iniv = 1,niv
         write(80,1001) altsmr(iniv),
     $                       (long(is,iniv),is= 43,43)
      end do
      write(80,*)
      write(80,1000)'ALTI. (KM)  ',(labsh(is),is = 1,7)
      write(80,*)
      do iniv = 1,niv
         write(80,1001) altsmr(iniv),
     $                       (short(is,iniv),is= 1,7)
      end do
      write(80,*)
      write(80,1000)'ALTI. (KM)  ',(labsh(is),is = 8,14)
      write(80,*)
      do iniv = 1,niv
         write(80,1001) altsmr(iniv),
     $                       (short(is,iniv),is= 8,14)
      end do
      write(80,*)
      write(80,1000)'ALTI. (KM)  ',(labsh(is),is = 15,15)
      write(80,*)
      do iniv = 1,niv
         write(80,1001) altsmr(iniv),
     $                       (short(is,iniv),is= 15,15)
      end do
      write(80,*)
c
      goto 30
c
 1000 format(2x,8a12)
 1001 format(8e12.4)
 1003 format(2x,8a12)
 1005 format(3f12.4,4e12.4,f12.4)
c
 20   return
      end
