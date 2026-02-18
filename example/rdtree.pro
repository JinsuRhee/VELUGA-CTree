FUNCTION rdtree_mkarray, n, tag

	CASE tag OF
		1L: RETURN, LONARR(n)
		2L: RETURN, LON64ARR(n)
		3L: RETURN, FLTARR(n)
		4L: RETURN, DBLARR(n)
	ENDCASE
END

FUNCTION rdtree_rdkey, fname
	OPENR, 1, fname

	bidtag	= 0L
	READU, 1, bidtag

	nkey 	= 0LL
	READU, 1, nkey

	key 	= rdtree_mkarray(nkey, bidtag)

	keyval 	= 0LL
	READU, 1, keyval

	READU, 1, key


	CLOSE, 1 


	key(0)	= keyval


	RETURN, key
END

FUNCTION rdtree_rdtree, fname
	OPENR, 1, fname

	gidtag 	= 0L
	snaptag = 0L
	bidtag 	= 0L
	mertag	= 0L

	READU, 1, gidtag
	READU, 1, snaptag
	READU, 1, bidtag
	READU, 1, mertag

	;ntree 	= 0LL
	;READU, 1, ntree

	lind	= 0LL
	READU, 1, lind

	tree 	= PTRARR(lind+1LL)

	FOR i=0L, lind DO BEGIN
		nbranch 	= 0L
		nmerge 		= 0L
		;fid 		= 0L
		fid 		= rdtree(1, bidtag)
		tstat		= 0L
		READU, 1, nbranch

		IF nbranch LE 0L THEN BEGIN
			tree(i)	= PTR_NEW({stat:-1L})
			CONTINUE
		ENDIF

		READU, 1, nmerge
		READU, 1, fid
		READU, 1, tstat

		idtmp 		= rdtree_mkarray(nbranch, gidtag)
		snaptmp		= rdtree_mkarray(nbranch, snaptag)
		;pidtmp		= rdtree_mkarray(nbranch, gidtag)
		;psnaptmp	= rdtree_mkarray(nbranch, snaptag)
		mertmp 		= rdtree_mkarray(nbranch, mertag)
		
		READU, 1, idtmp
		READU, 1, snaptmp
		;READU, 1, pidtmp
		;READU, 1, psnaptmp
		READU, 1, mertmp

		IF nmerge GE 1L THEN BEGIN
			midtmp 			= rdtree_mkarray(nmerge, gidtag)
			msnaptmp 		= rdtree_mkarray(nmerge, snaptag)
			mmertmp 		= rdtree_mkarray(nmerge, mertag)
			mbidtmp 		= rdtree_mkarray(nmerge, bidtag)

			READU, 1, midtmp
			READU, 1, msnaptmp
			READU, 1, mmertmp
			READU, 1, mbidtmp
		ENDIF ELSE BEGIN
			midtmp 			= rdtree_mkarray(1, gidtag)
			msnaptmp 		= rdtree_mkarray(1, snaptag)
			mmertmp 		= rdtree_mkarray(1, mertag)
			mbidtmp 		= rdtree_mkarray(1, bidtag)

		ENDELSE

		tmp 	= {br_len:nbranch, father_ID:fid, n_mergebr:nmerge, $
			id:idtmp, snap:snaptmp, p_id:pidtmp, p_snap:psnaptmp, merit:mertmp, $
			m_id:midtmp, m_snap:msnaptmp, m_merit:mmertmp, m_bid:mbidtmp, $
			stat:1L}
		tree(i)	= PTR_NEW(tmp)
	ENDFOR

	CLOSE, 1

	tree 	= tree(0L:lind+1L)
	RETURN, tree
END

PRO rdtree

	treename		= './ctree_key.dat'
	keyname			= './ctree_tree.dat'
	varname 		= './ctree.sav'

	tree_key		= rdtree_rdkey(treename)
	tree_data		= rdtree_rdtree(keyname)


	SAVE, filename=varname, tree_key, tree_data
	

	
	STOP
END
