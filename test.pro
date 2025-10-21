PRO test

	 treeset = {$
                        tag_num : 'NumDesc', $
                        tag_off : 'DescOffsets', $
                        tag_result      : 'Descendants', $
                        tag_npart       : 'DescNpart', $
                        tag_merit       : 'Merits', $
                        tag_nlink       : 'Nsteps_search_new_links' $
                        }

		fname 	= '/data100/jinsu/NG/parentDM/VELOCIraptor/Halo/tree/tfout/tree.snapshot_0098VELOCIraptor.tree'

        fid     = H5F_OPEN(fname)

        did     = H5D_OPEN(fid, treeset.tag_num)
        num = H5D_READ(did) & H5D_CLOSE, did

        did     = H5D_OPEN(fid, treeset.tag_off)
        off = H5D_READ(did) & H5D_CLOSE, did

        did     = H5D_OPEN(fid, treeset.tag_result)
        res = H5D_READ(did) & H5D_CLOSE, did

        did     = H5D_OPEN(fid, treeset.tag_merit)
        mer = H5D_READ(did) & H5D_CLOSE, did

        did     = H5D_OPEN(fid, treeset.tag_npart)
        npart = H5D_READ(did) & H5D_CLOSE, did

        did     = H5A_OPEN_NAME(fid, treeset.tag_nlink)
        nlink = H5A_READ(did) & H5A_CLOSE, did

        did     = H5D_OPEN(fid, 'ID')
        id = H5D_READ(did) & H5D_CLOSE, did

        H5F_CLOSE, fid

        dum     = {num:num, off:off, res:res, merit:mer, npart:npart, nlink:nlink, id:id}
	STOP
END
