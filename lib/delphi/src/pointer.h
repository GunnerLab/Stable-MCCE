c**************************************************************
c     integer*8 i_iepsmp,i_idebmap,i_ioff,i_phimap,i_phimap1,i_phimap2,
c    &i_phimap3,i_db,i_idpos,i_sf1,i_sf2,i_qmap1,i_qmap2,i_debmap1,
c    &i_debmap2,i_bndx1,i_bndx2,i_bndx3,i_bndx4,i_ibndx,i_ibndy,
c    &i_ibndz,i_neps,i_keps,i_cgrid,i_spdiv,i_sen,i_spot,i_sqs,
c    &i_iepsv,i_rfield,i_r0,i_r02,i_rs2,i_expos,i_pls,i_ast,i_ast2,
c    &i_ibnd,i_ibgrd,i_bndeps,i_cbn1,i_cbn2,i_cbal,i_icbn,i_iab1,
c    &i_iab2,i_icume,i_iexpos,i_atndx,i_scspos,i_atsurf,i_atpos,i_xn2,
c    &i_rad3,i_chrgv4,i_atinf,i_iatmmed,i_medeps,i_dataobject,
c    &i_datadistr,i_nqgrdtonqass,i_limobject,i_iatmobj,i_coi,
c    &i_internal,i_atmeps,i_tmpmap,i_limgunit,i_sout,i_gchrgtmp,
c    &i_cgbp,i_atmcrg,i_chgpos,i_atmforce,i_polariz,i_sitephi,i_qval,
c    &i_gchrg,i_gchrgp,i_iqpos,i_chrgv2,i_cqs,i_vert,i_vindx,i_vtemp,
c    &i_scsnor,i_vnorm,i_nsel,i_vtlen,i_vtlst,i_tmlst,i_vtpnt,i_schrg,
c    &i_crgatn,i_gchrgd,i_gchrg2,i_gval,i_phimap4

c**************************************************************

	pointer (i_iepsmp,iepsmp),(i_idebmap,idebmap),(i_ioff,ioff)
	pointer (i_phimap,phimap),(i_phimap1,phimap1)
	pointer (i_phimap2,phimap2),(i_phimap3,phimap3)
	pointer (i_db,db),(i_idpos,idpos),(i_sf1,sf1),(i_sf2,sf2)
	pointer (i_qmap1,qmap1),(i_qmap2,qmap2)
	pointer (i_debmap1,debmap1),(i_debmap2,debmap2)
	pointer (i_bndx1,bndx1),(i_bndx2,bndx2),(i_bndx3,bndx3)
	pointer (i_bndx4,bndx4)
	pointer (i_ibndx,ibndx),(i_ibndy,ibndy),(i_ibndz,ibndz)
	pointer (i_neps,neps),(i_keps,keps)
	pointer (i_cgrid,cgrid),(i_spdiv,spdiv),(i_sen,sen)
	pointer (i_spot,spot),(i_sqs,sqs),(i_iepsv,iepsv)
	pointer (i_rfield,rfield)

	pointer (i_r0,r0),(i_r02,r02),(i_rs2,rs2)
	pointer (i_expos,expos),(i_pls,pls),(i_ast,ast),(i_ast2,ast2)
	pointer (i_ibnd,ibnd),(i_ibgrd,ibgrd),(i_bndeps,bndeps)
	pointer (i_cbn1,cbn1),(i_cbn2,cbn2),(i_cbal,cbal),(i_icbn,icbn)
	pointer (i_iab1,iab1),(i_iab2,iab2),(i_icume,icume)
	pointer (i_iexpos,iexpos)
	pointer (i_atndx,atndx),(i_scspos,scspos),(i_atsurf,atsurf)

	pointer (i_atpos,atpos),(i_xn2,xn2),(i_rad3,rad3)
	pointer (i_chrgv4,chrgv4),(i_atinf,atinf)
c added by walter
        pointer (i_iatmmed,iatmmed) 
        pointer (i_medeps,medeps)
        pointer (i_dataobject,dataobject)
        pointer (i_datadistr,datadistr)
        pointer (i_nqgrdtonqass,nqgrdtonqass)
        pointer (i_limobject,limobject)
        pointer (i_iatmobj,iatmobj)
        pointer (i_coi,coi)
        pointer (i_internal,internal)
        pointer (i_atmeps,atmeps)
        pointer (i_tmpmap,tmpmap)
        pointer (i_limgunit,limgunit)
        pointer (i_sout,sout)
        pointer (i_gchrgtmp,gchrgtmp)
        pointer (i_cgbp,cgbp)
        pointer (i_atmcrg,atmcrg)
        pointer (i_chgpos,chgpos)
        pointer (i_atmforce,atmforce)
        pointer (i_polariz,polariz)
        pointer (i_sitephi,sitephi)
		pointer (i_epsmap,epsmap)
		pointer (i_debmap,debmap),(i_qmap,qmap)
c end addition
	pointer (i_qval,qval),(i_gchrg,gchrg),(i_gchrgp,gchrgp)
	pointer (i_iqpos,iqpos)
	pointer (i_chrgv2,chrgv2),(i_cqs,cqs)

	pointer (i_vert,vert),(i_vindx,vindx),(i_vtemp,vtemp)
	pointer (i_scsnor,scsnor),(i_vnorm,vnorm)
	pointer (i_nsel,nsel),(i_vtlen,vtlen),(i_vtlst,vtlst)
	pointer (i_tmlst,tmlst),(i_vtpnt,vtpnt)
	pointer (i_schrg,schrg),(i_crgatn,crgatn)
	pointer (i_gchrgd,gchrgd),(i_gchrg2,gchrg2),(i_gval,gval)
       pointer (i_phimap4,phimap4)
C COMPILER-based  DIFFERENCES
c       use this for 32bits architectures
c      integer*4 pntr(102)
c       use this for 64bits architectures
 	integer*8 pntr(102)
c       use also the following line when using pgi compiler for Linux 64bits machine
         integer*8 memalloc


	equivalence (i_iepsmp,pntr(1)),(i_idebmap,pntr(2))
	equivalence (i_ioff,pntr(3))
	equivalence (i_phimap,pntr(4)),(i_phimap1,pntr(5))
	equivalence (i_phimap2,pntr(6)),(i_phimap3,pntr(7))
	equivalence (i_db,pntr(8)),(i_idpos,pntr(9))
	equivalence (i_sf1,pntr(10)),(i_sf2,pntr(11))
	equivalence (i_qmap1,pntr(12)),(i_qmap2,pntr(13))
	equivalence (i_debmap1,pntr(14)),(i_debmap2,pntr(15))
	equivalence (i_bndx1,pntr(16)),(i_bndx2,pntr(17))
	equivalence (i_bndx3,pntr(18)),(i_bndx4,pntr(19))
	equivalence (i_ibndx,pntr(20)),(i_ibndy,pntr(21))
	equivalence (i_ibndz,pntr(22))
	equivalence (i_neps,pntr(23)),(i_keps,pntr(24))
	equivalence (i_cgrid,pntr(25))

	equivalence (i_r0,pntr(26)),(i_r02,pntr(27)),(i_rs2,pntr(28))
	equivalence (i_expos,pntr(29)),(i_pls,pntr(30))
	equivalence (i_ast,pntr(31)),(i_ast2,pntr(32))
	equivalence (i_ibnd,pntr(33)),(i_ibgrd,pntr(34))
	equivalence (i_bndeps,pntr(35))
	equivalence (i_cbn1,pntr(36)),(i_cbn2,pntr(37))
	equivalence (i_cbal,pntr(38)),(i_icbn,pntr(39))
	equivalence (i_iab1,pntr(40)),(i_iab2,pntr(41))
	equivalence (i_icume,pntr(42)),(i_iexpos,pntr(43))
	equivalence (i_atndx,pntr(44)),(i_scspos,pntr(45))
	equivalence (i_atsurf,pntr(46))
c added by walter
        equivalence (i_atmforce,pntr(47))
        equivalence (i_polariz,pntr(48))
        equivalence (i_sitephi,pntr(49))
		equivalence (i_epsmap,pntr(50))
		equivalence (i_debmap,pntr(101)),(i_qmap,pntr(102))
c end addition
	equivalence (i_spdiv,pntr(51)),(i_sen,pntr(52))
	equivalence (i_spot,pntr(53)),(i_iepsv,pntr(54))
	equivalence (i_rfield,pntr(55)),(i_sqs,pntr(56))

	equivalence (i_atpos,pntr(57)),(i_xn2,pntr(58))
	equivalence (i_rad3,pntr(59)),(i_chrgv4,pntr(60))
	equivalence (i_atinf,pntr(61))
	equivalence (i_qval,pntr(62)),(i_gchrg,pntr(63))
	equivalence (i_gchrgp,pntr(64)),(i_iqpos,pntr(65))
	equivalence (i_chrgv2,pntr(66)),(i_cqs,pntr(67))

	equivalence (i_vert,pntr(68)),(i_vindx,pntr(69))
	equivalence (i_scsnor,pntr(70)),(i_vnorm,pntr(71))
	equivalence (i_vtemp,pntr(72))
	equivalence (i_nsel,pntr(73)),(i_vtlen,pntr(74))
	equivalence (i_vtlst,pntr(75)),(i_tmlst,pntr(76))
	equivalence (i_vtpnt,pntr(77))
	equivalence (i_schrg,pntr(78)),(i_crgatn,pntr(79))
	equivalence (i_gchrgd,pntr(80)),(i_gchrg2,pntr(81))
	equivalence (i_gval,pntr(82))
c added by walter
        equivalence (i_iatmmed,pntr(83))
        equivalence (i_medeps,pntr(84))
        equivalence (i_dataobject,pntr(85))
        equivalence (i_nqgrdtonqass,pntr(86))
        equivalence (i_limobject,pntr(87))
        equivalence (i_datadistr,pntr(88))
        equivalence (i_iatmobj,pntr(89))
        equivalence (i_coi,pntr(90))
        equivalence (i_internal,pntr(91))
        equivalence (i_atmeps,pntr(92))
        equivalence (i_tmpmap,pntr(94))
        equivalence (i_sout,pntr(95))
        equivalence (i_gchrgtmp,pntr(96))
        equivalence (i_cgbp,pntr(97))
        equivalence (i_atmcrg,pntr(98))
        equivalence (i_chgpos,pntr(99))
        equivalence (i_limgunit,pntr(100))
c end addition
        equivalence (i_phimap4,pntr(93))

	common /pointr/ pntr