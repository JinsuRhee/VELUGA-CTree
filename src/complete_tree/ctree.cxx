
#include "global/io.h"
#include "complete_tree/ctree.h"
#include "utilities/utilities.h"

namespace Ctree{

	//-----
	// ETC
	//-----
	void expandbr(const vctree_set::Settings& vh, ControlArray& data, CT_I32 ind, Tree::TreeArray& tree, Tree::TreeKeyArray& key, CT_I32 id_to_link, CT_I32 snap_to_link, CT_double merit_to_link){
		bool istree = Tree::istree(key, data[ind].snap, data[ind].id);


		if(!istree){
			LOG()<<" why this: should not be happened";
			u_stop();

			makenewbr(vh, data, ind, data[ind].snap, data[ind].id, tree, key);
		}

		//Tree::TreeSt tree0 	= Tree::gettree(tree, key, data[ind].snap, data[ind].id);

		Tree::treeinput(tree, key, data[ind].snap, data[ind].id, snap_to_link, id_to_link, merit_to_link);
	}


	void makenewbr(const vctree_set::Settings& vh, ControlArray& data, CT_I32 ind, CT_snap snap0, CT_ID id0, Tree::TreeArray& tree, Tree::TreeKeyArray& key){

		Tree::TreeSt newtree;

		Tree::treevecinput<Tree::Tree_GID>(newtree.id, 0, id0);
		Tree::treevecinput<Tree::Tree_Snap>(newtree.snap, 0, snap0); 
		Tree::treevecinput<Tree::Tree_merit>(newtree.p_merit, 0, -1.);
		newtree.endind = 0;

		tree[0].lind ++;

		if((Tree::Tree_BID) tree.size() <= tree[0].lind){
			tree.resize( (Tree::Tree_BID) tree.size() + vh.ctree_nstep);
		}
		
		tree[tree[0].lind] 	= std::move(newtree);


		Tree::Tree_I64 keyval;
		keyval 	= snap0 + key[0].key * id0;

		if((Tree::Tree_BID) key.size() <= keyval){
			key.resize( (Tree::Tree_BID) key.size() + vh.ctree_nstep);
		}

		key[keyval].ind 	= tree[0].lind;

		data[ind].id 	= id0;
		data[ind].snap 	= snap0;

	}




	void get_merit(std::vector<CT_PID>& pid_g, std::vector<CT_I32>& gid_g, 
		std::vector<CT_PID>& pid_s, std::vector<CT_I32>& gid_s, 
		std::vector<CT_I32>& hash, std::vector<CT_I32>& hash_next, 
		std::vector<CT_I32>& npart_g, std::vector<CT_I32>& npart_s, 
		std::vector<CT_I32>& match_id, std::vector<CT_double>& match_merit){

		// local variables
		CT_I32 n_pg = pid_g.size();
		//CT_I32 n_ps = pid_s.size();
		CT_I32 n_g 	= npart_g.size();
		CT_I32 n_s 	= npart_s.size();
		CT_I32 n_dn = pid_s.size();
		//CT_I32 n_hashth	= hash.size();


		// Make Index Table
		std::vector<CT_I32> tbl_g(n_g, CT_I32{0});
		std::vector<CT_I32> tbl_s(n_s, CT_I32{0});
		CT_I32 nn_g = 0;
		CT_I32 nn_s = 0;

		for(CT_I32 i=0; i<n_g; i++){
			if(npart_g[i] > 0){
				tbl_g[i] = nn_g;
				nn_g ++;
			}
		}

		for(CT_I32 i=0; i<n_s; i++){
			if(npart_s[i] > 0){
				tbl_s[i] = nn_s;
				nn_s ++;
			}
		}

		//nn_g --;
		//nn_s --;

		// merit allocation
		
		MeritMatrix<CT_double> merit(nn_g , nn_s, CT_double{0});

		// Hash
		CT_I32 i0, ind;

		for(CT_I32 i=0; i<n_pg; i++){
			ind 	= std::abs(pid_g[i]) % n_dn;

			if(ind < 0) ind = 0;
			

			
			i0 	= hash[ind];

			if(i0 < 0) continue;

			while(true){
				if(pid_g[i]==pid_s[i0]){
					merit(tbl_g[gid_g[i]], tbl_s[gid_s[i0]]) += 1.;
					break;
				} else{
					if(hash_next[i0] >= 0){
						i0 = hash_next[i0];
					} else{
						break;
					}
				}
			}		
		}


		for(CT_I32 i=0; i<nn_g; i++){
			for(CT_I32 j=0; j<nn_s; j++){
				merit(i,j) 	= merit(i,j) * merit(i,j);

				if(npart_g[i]>0){
					merit(tbl_g[i],j) /= npart_g[i];
				}

				if(npart_s[j]>0){
					merit(i,tbl_s[j]) /= npart_s[j];
				}
			}
		}


		ind = 0;
		std::vector<CT_I32> idtoind(n_g, CT_I32{0});

		for(CT_I32 i=0; i<n_g; i++){
			if(npart_g[i] > 0){
				idtoind[i] = ind;
				ind += 1;
			}
		}



		CT_double dum;
		CT_I32 dumid;
		for(CT_I32 i=0; i<n_g; i++){
			if(npart_g[i] <= 0) continue;

			dum = -1.;
			dumid = -1;


			for(CT_I32 j=0; j<n_s; j++){
				if(npart_s[j] <= 0) continue;
				if(merit(tbl_g[i], tbl_s[j]) > dum && merit(tbl_g[i], tbl_s[j]) > 0){
					dum = merit(tbl_g[i], tbl_s[j]);
					dumid = j;
				}
			}

			if(dumid > 0){
				match_id[idtoind[i]] = dumid;
				match_merit[idtoind[i]] = dum;
			}

		}

	}

	// for simple merit
	CT_double get_merit2(std::vector<CT_PID>& pid0, std::vector<CT_PID>& pid, std::vector<CT_double>& w0, std::vector<CT_double>& w1, CT_I32 merittype){

		CT_I32 n_raw 	= pid0.size();
		CT_I32 n_match 	= pid.size();
		CT_double share = 0.;
		CT_double weight1 = 0.;
		CT_double weight2 = 0.;

		std::vector<CT_I32> hash(n_raw, -1);
		std::vector<CT_I32> hash_next(n_raw, -1);

		// Make Hash
		CT_I32 ind, i0;
		for(CT_I32 i=0; i<n_raw; i++){
			ind 	= std::abs(pid0[i]) % n_raw;

			if(ind < 0) ind = 0;

			if(hash[ind] < 0){
				hash[ind] = i;
			} else{
				i0 	= hash[ind];
				while(true){
					if(hash_next[i0] < 0){
						hash_next[i0] = i;
						break;
					} else{
						i0 	= hash_next[i0];
					}
				}
			}
		}

#ifdef CTREE_USE_OMP
		#pragma omp parallel for default(none) \
			private(ind, i0) \
    		shared(n_match, pid, pid0, n_raw, hash, hash_next, w0, w1) \
    		reduction(+:share, weight1, weight2)
#endif
		for(CT_I32 i=0; i<n_match; i++){
			ind 	= std::abs(pid[i]) % n_raw;

			if(ind < 0) ind = 0;
			i0 	= hash[ind];

			if(i0 < 0) continue;

			while(true){
				if(pid[i] == pid0[i0]){
					share += 1.;
					weight1 	+= w0[i0];
					weight2		+= w1[i];
					break;
				}else{
					if(hash_next[i0] >= 0){
						i0 	= hash_next[i0];
					}else{
						break;
					}
				}
			}
		}


		if(merittype == 1) share 	= (share*share) / ((CT_double) n_raw) / ( (CT_double) n_match );
		if(merittype == 2) share 	= share / ( (CT_double) n_raw );
		if(merittype == 3) share 	= (share*share) / ((CT_double) n_raw) / ( (CT_double) n_match ) * weight1 * weight2;

		return share;
	}

	std::vector<CT_I32> get_control(ControlArray& data, CT_I32 type){
		std::vector<CT_I32> cut(data[0].last_ind+1);
		CT_I32 ncut = 0;
		for(CT_I32 i=0; i<data[0].last_ind+1; i++){

			bool is = false;

			if(type == 0){
				is 	= data[i].stat == 0;
			} else if(type == 1){
				is 	= (data[i].stat == 0) && (data[i].n_ptcl < 0);
			}

			if(is){
				cut[ncut] 	= i;
				ncut ++;
			}
		}

		cut.resize(ncut);
		return cut;
	}

	// Reallocate PID Array
	void PIDReallocate(const vctree_set::Settings& vh, PIDArray& pid, CT_I32 ind){
		if((CT_I32) pid.size() > ind) return;
		CT_I32 stepsize = vh.ctree_npid;
		CT_I32 nn = 1;
		while(true){
			if(stepsize > ind) break;
			stepsize += vh.ctree_npid*nn;
			nn ++;
		}
		pid.resize(stepsize);
	}

	// merit related
	std::vector<CT_double> get_weight(const vctree_set::Settings& vh, std::vector<CT_PID> pid){
		CT_I32 npart = pid.size();
		std::vector<CT_double> weight(npart);
		
		CT_double dind;
		CT_double factor = (0.5772156649 + std::log( (CT_double) npart ));

		if(vh.ctree_weighttype == 1){

#ifdef CTREE_USE_OMP
			#pragma omp parallel for default(none) \
    			shared(factor, weight, npart) \
    			private(dind)
#endif
			for(CT_I32 i=0; i<(CT_I32) npart; i++){
				dind 	= npart - i;
				weight[i] 	= dind / ((CT_double) npart) / factor;
			}
		}else{
			LOG()<<" Not implemented yet";
			u_stop();
		}

		return weight;
	}

	// Find index in snapshot array
	CT_I32 wheresnap(IO::snapinfo& sinfo, CT_I32 snap_curr){
		CT_I32 ind=0;
		for(CT_I32 j=0; j<(CT_I32) sinfo.size() + 1; j++){
			if(sinfo[j].snum == snap_curr){
				ind = j;
				break;
			}
		}
		return ind;
	}

	// Extract Core Particles
	PIDArray get_coreptcl(const vctree_set::Settings& vh, PIDArray& pid){

		CT_I32 npid = pid.size();

		std::sort(pid.begin(), pid.end(), [](const PIDSt& a, const PIDSt& b) {
            if (a.pid != b.pid) return a.pid < b.pid;
            return a.weight > b.weight; 
        });

		std::vector<CT_I32> uind(npid);

		CT_I32 uind_n = 1;
		CT_PID uind_ptdum = pid[0].pid;
		uind[0] 	= 0;

		for(CT_I32 i=1; i<npid; i++){

			if(pid[i].pid > uind_ptdum){
				uind[uind_n]	= i;
				uind_ptdum 		= pid[i].pid;
				uind_n ++;
			}
		}
		uind.resize(uind_n);

		// step
		std::vector<CT_I32> numid(uind.size());
		//numid[0]	= uind[0] + 1;
		for(CT_I32 i=0; i<(CT_I32) uind.size()-1; i++){
			numid[i]	= uind[i+1] - uind[i];
		}
		numid[uind.size()-1] 	= pid.size() - uind[uind.size()-1];

		// resize particle array for exsiting ones for multiple snapshots
		std::vector<CT_I32> ucut(npid);
		CT_I32 uncut = 0;
		CT_I32 n_step_bw0 = vh.ctree_n_step_n;
		while(true){
			uncut = 0;
			for(CT_I32 i=0; i<uind_n; i++){
				if(numid[i] >= n_step_bw0){
					ucut[uncut]	= uind[i];
					uncut ++;
				}
			}
			if( ((CT_double) uncut) / ((CT_double) uind_n) > vh.ctree_minfrac) break;
			n_step_bw0 --;
		}

		ucut.resize(uncut);

		// input
		PIDArray pid2(uncut);

		for(CT_I32 i=0; i<uncut; i++){
			pid2[i] 		= pid[ucut[i]];
		}


		pid2.resize(uncut);

		pid2[0].n_con 	= vh.ctree_n_step_dn * (n_step_bw0+1);
		return pid2;
	}

	// Free
	void ctfree(const vctree_set::Settings& vh, ControlArray& data, CT_I32 ind, CT_I32 s_end, CT_I32 id_end, CT_I32 snap0){

		//dkey[ data[ind].snap + dkey[0]*data[ind].id ] = -1;
		//dkey[ data[ind].snap0 + dkey[0]*data[ind].id0 ] = -1;

		CT_I32 ncut;
		data[ind].detstat 	= -1;

		data[ind].Free_plist();

		std::vector<CT_I32> cut;

		data[ind].n_ptcl 	= -1;
		if(s_end<0){
			//dkey[ data[ind].snap + dkey[0]*data[ind].id ] = -1;
			data[ind].id = -1;
			data[ind].snap = -1;

			for(CT_I32 i=0; i<vh.ctree_n_search; i++){
				data[ind].list[i].merit = -1.;
				data[ind].list[i].id 	= -1;
				data[ind].list[i].snap 	= -1;
			}
			data[ind].list_n 		= 0;
		}else if(s_end > snap0){
			data[ind].id 	= id_end;
			data[ind].snap 	= s_end;
			//dkey[ data[ind].snap + dkey[0]*data[ind].id ] = ind;

			ncut = 0;
			for(CT_I32 i=0; i<vh.ctree_n_search; i++){
				if(data[ind].list[i].snap < s_end && data[ind].list[i].snap > 0){
					ncut ++;
					cut.push_back(i);
				}
			}

			if(ncut == 0){
				for(CT_I32 i=0; i<vh.ctree_n_search; i++){
					data[ind].list[i].merit = -1.;
					data[ind].list[i].id 	= -1;
					data[ind].list[i].snap 	= -1;
				}
				data[ind].list_n 		= 0;
			}else{
				for(CT_I32 i=0; i<ncut; i++){
					data[ind].list[i].merit = data[ind].list[cut[i]].merit;
					data[ind].list[i].id 	= data[ind].list[cut[i]].id;
					data[ind].list[i].snap 	= data[ind].list[cut[i]].snap;
				}
				data[ind].list_n	= ncut;

				if(ncut <= vh.ctree_n_search-1){
					for(CT_I32 i=ncut; i<vh.ctree_n_search; i++){
						data[ind].list[i].merit = -1.;
						data[ind].list[i].id 	= -1;
						data[ind].list[i].snap 	= -1;
					}
				}
			}
			
		}else if(s_end <= snap0 && s_end > 0){
			data[ind].id 	= id_end;
			data[ind].snap 	= s_end;
			//dkey[ data[ind].snap + dkey[0]*data[ind].id ] = ind;

			for(CT_I32 i=0; i<vh.ctree_n_search; i++){
				data[ind].list[i].merit = -1.;
				data[ind].list[i].id 	= -1;
				data[ind].list[i].snap 	= -1;
			}
			data[ind].list_n 		= 0;
		}
	}

	// Reallocate Control Array
	void reallocate(const vctree_set::Settings& vh, ControlArray& data, CT_I32 nn){
		CT_I32 old_n = data.size();
		data.reserve(old_n + nn);
		for(CT_I32 i=0; i<nn; i++) data.emplace_back(vh);
	}

	 
	//-----
	// Input Galaxy to control array
	//-----
	void inputgal(const vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& key, ControlArray& data, ControlKey& dkey, IO_dtype::GalArray& gal){

		CT_I32 i0, i1, dsize, nstep;

		i0	= data[0].last_ind + 1;
		i1 	= i0 + ( (CT_I32) gal.size() ) - 1;
		dsize 	= data.size();

		if(i1 > dsize-1){
			nstep 	= vh.ctree_nstep;
			while(true){
				if(nstep > i1-dsize) break;
				nstep += vh.ctree_nstep;
			}
			reallocate(vh, data, nstep);
		}

		// OMP HERE
		for(CT_I32 i=i0; i<i1+1; i++){
			data[i].id0 	= gal[i-i0].id;
			data[i].snap0 	= gal[i-i0].snap;

			if(Tree::istree(key, data[i].snap0, data[i].id0)){
				Tree::TreeSt tree0 = Tree::gettree(tree, key, data[i].snap0, data[i].id0);
				data[i].id 		= tree0.id[0];
				data[i].snap 	= tree0.snap[0];
			}else{
				// no branch exists
				makenewbr(vh, data, i, data[i].snap0, data[i].id0, tree, key);
				//data[i].id 		= data[i].id0;
				//data[i].snap 	= data[i].snap0;
			}

			dkey[ data[i].snap0 + dkey[0] * data[i].id0 ]	= i;


			//dkey[ data[i].snap + dkey[0] * data[i].id ]	= i;

		}
		data[0].last_ind 	= i1;

		
		
	}

	//-----
	// TREE CLASSIFICATION
	//-----
	void classify(const vctree_set::Settings& vh, ControlArray& data, IO::snapinfo& sinfo, CT_I32 snap_curr, ctree_num& number){
		//Log Initialize
		number.T = 0;
		number.B = 0;
		number.C = 0;

		//OMP or MPI HERE
		CT_I32 ind = 0;
		CT_I32 ind0 = 0;
		CT_I32 ind1 = 0;

		ind0	= wheresnap(sinfo, snap_curr);
		for(CT_I32 i=0; i<data[0].last_ind + 1; i++){
			if(data[i].id0 <= 0){
				data[i].stat = 2; // dummy
				continue;
			}

			if(data[i].stat == -1){
				number.B ++;
				continue;
			}
			// tree exists
			if(data[i].snap <= snap_curr){
				data[i].stat = 1;
				number.T ++;
				continue;
			}
			// tree length check
			ind 	= wheresnap(sinfo, data[i].snap);

			if(data[i].list[0].snap > 0){
				ind1 	= wheresnap(sinfo, data[i].list[0].snap);
			}else{
				ind1 	= ind0 + vh.ctree_n_search*2;
			}

			// determine
			if((ind - ind0) <= vh.ctree_n_search || (ind1-ind0) < vh.ctree_n_search){
				data[i].stat = 0;
				number.C ++;
			} else{
				data[i].stat = -1;
				ctfree(vh, data, i, -1, -1, snap_curr);
				number.B ++;
			}
		} 
	}

	//-----
	// Collect PID
	//-----
	PIDArray collectpidalongbranch(const vctree_set::Settings& vh, std::vector<CT_snap>& slist, std::vector<CT_ID>& idlist){

		PIDArray pid;
		std::vector<CT_double> weight;
		IO_dtype::GalArray gal0;


		gal0 	= IO::r_gal(vh, slist[0], idlist[0], true);

		CT_I32 nn = gal0[0].pid.size();

		// initial allocation
		pid.resize(nn*vh.ctree_n_step_n*2);

		// input
		CT_I32 ind1 	= 0;
		CT_I32 ind0 	= 0;
		CT_I32 loop_n 	= 0;
		CT_I32 loop_ind = 0;
		//CT_I32 stepsize = 0;

		while(true){
			if(loop_ind >= (CT_I32) slist.size()) break;
			weight.clear();

			gal0 	= IO::r_gal(vh, slist[loop_ind], idlist[loop_ind], true);

			ind1 	= ind0 + gal0[0].pid.size()-1;
			weight 	= get_weight(vh, gal0[0].pid);

			PIDReallocate(vh, pid, ind1);


			//if(ind1 >= (CT_I32) pid.size()){
			//	stepsize 	= pid.size() + ind1 + vh.ctree_npid;
			//	while(true){
			//		if(ind1 < stepsize) break;
			//		stepsize += vh.ctree_npid;
			//	}
			//	pid.resize( stepsize );
			//}

			// OMP?
#ifdef CTREE_USE_OMP
			#pragma omp parallel for default(none) \
    			shared(ind0, ind1, pid, gal0, weight)
#endif
			for(CT_I32 i=ind0; i<ind1+1; i++){
				pid[i].pid 	= gal0[0].pid[i-ind0];
				pid[i].weight 	= weight[i-ind0];
			}

			ind0 	= ind1 + 1;
			loop_n ++;
			loop_ind 	+= vh.ctree_n_step_dn;
			if(loop_n >= vh.ctree_n_step_n || loop_ind >= (CT_I32) slist.size()) break;
		}
		// remove particles appearing multiple times
		pid.resize(ind0);

		PIDArray pid2;
		pid2 	= get_coreptcl(vh, pid);


//		std::sort(pid.begin(), pid.end(), [](const PIDSt& a, const PIDSt& b) {
//            if (a.pid != b.pid) return a.pid < b.pid;
//            return a.weight > b.weight; 
//        });
//
//		std::vector<CT_I32> uind(ind0);
//
//		CT_I32 uind_n = 1;
//		CT_PID uind_ptdum = pid[0].pid;
//		uind[0] 	= 0;
//
//		for(CT_I32 i=1; i<ind0; i++){
//			if(pid[i].pid > uind_ptdum){
//				uind[uind_n]	= i;
//				uind_ptdum 		= pid[i].pid;
//				uind_n ++;
//			}
//		}
//		uind.resize(uind_n);
//
//		// step
//		std::vector<CT_I32> numid(uind.size());
//		//numid[0]	= uind[0] + 1;
//		for(CT_I32 i=0; i<uind.size()-1; i++){
//			numid[i]	= uind[i+1] - uind[i];
//		}
//		numid[uind.size()-1] 	= pid.size() - uind[uind.size()-1];
//
//		// resize particle array for exsiting ones for multiple snapshots
//		std::vector<CT_I32> ucut(ind0);
//		CT_I32 uncut = 0;
//		CT_I32 n_step_bw0 = vh.ctree_n_step_n;
//		while(true){
//			uncut = 0;
//			for(CT_I32 i=0; i<uind_n; i++){
//				if(numid[i] >= n_step_bw0){
//					ucut[uncut]	= uind[i];
//					uncut ++;
//				}
//			}
//			if( ((CT_double) uncut) / ((CT_double) uind_n) > 0.25) break;
//			n_step_bw0 --;
//		}
//
//		ucut.resize(uncut);
//
//		// input
//		PIDArray pid2(uncut);
//
//		for(CT_I32 i=0; i<uncut; i++){
//			pid2[i] 		= pid[ucut[i]];
//		}
//
//
//		pid2.resize(uncut);

		return pid2;
		
	}

	PIDArray collectpid(const vctree_set::Settings& vh, ControlArray& data, Tree::TreeArray& tree, Tree::TreeKeyArray& key){

		// Gather Target
		//int myrank = mpi_rank();

		std::vector<CT_I32> cut;
		CT_I32 ncut = 0;
		

		cut 	= get_control(data, 1);
		ncut 	= cut.size();
		//for(CT_I32 i=0; i<data[0].last_ind+1; i++){
		//	if(data[i].stat == 0 && data[i].n_ptcl < 0){
		//		ncut ++;
		//		cut.push_back(i);
		//	}
		//}

		// Input particle list to the control array
		std::vector<CT_PID> pid0;
		PIDArray cpid, cpid2;
		Tree::TreeSt tree0;
		std::vector<CT_snap> t_slist;
		std::vector<CT_ID> t_idlist;
		std::vector<CT_double> weight;
		IO_dtype::GalArray gal0;

		CT_I32 ntcut;

		//t_slist.reserve(vh.ctree_n_search);
		//t_idlist.reserve(vh.ctree_n_search);
		t_slist.resize(vh.ctree_n_search);
		t_idlist.resize(vh.ctree_n_search);

		if(ncut > 0){

#ifdef CTREE_USE_MPI
			int rank = 0, size = 1;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	    	MPI_Comm_size(MPI_COMM_WORLD, &size);

			for(CT_I32 ind : cut){
				int owner = ind % size;
				if(owner != rank) continue;

				pid0.clear();
				cpid.clear();

				gal0 	= IO::r_gal(vh, data[ind].snap, data[ind].id, true);


				if(!istree(key, data[ind].snap, data[ind].id)){ 
					// no treee for this galaxy (read from a single snapshot)
					weight 	= get_weight(vh, gal0[0].pid);
					cpid2.resize(gal0[0].pid.size());
					for(CT_I32 k=0; k<(CT_I32) cpid2.size(); k++){
						cpid2[k].pid 	= gal0[0].pid[k];
						cpid2[k].weight 	= weight[k];

					}

					cpid 	= get_coreptcl(vh, cpid2);

					//pid0 = std::move(gal0[0].pid);

				} else {
					// this galaxy has a branch. Collect particles along its branch
					tree0 = Tree::gettree(tree, key, data[ind].snap, data[ind].id);

					ntcut = 0;
					t_slist.clear();
					t_idlist.clear();

					for(CT_I32 i=0; i<vh.ctree_n_search; i++){
						if(tree0.snap[i] >= data[ind].snap){
							ntcut ++;
							t_slist.push_back(tree0.snap[i]);
							t_idlist.push_back(tree0.id[i]);
						}
					}


					if(ntcut == 0){
						LOG()<<"Weird branch: Snap = "<<data[ind].snap0<<" / ID = "<<data[ind].id0;
						u_stop();
					}

					cpid 	= collectpidalongbranch(vh, t_slist, t_idlist);
				}


				// store to control array
				data[ind].p_list 	= std::move(cpid);
				data[ind].n_ptcl 	= data[ind].p_list.size();
				data[ind].pos[0]	= gal0[0].xc;
				data[ind].pos[1]	= gal0[0].yc;
				data[ind].pos[2]	= gal0[0].zc;

				data[ind].vel[0]	= gal0[0].vxc;
				data[ind].vel[1]	= gal0[0].vyc;
				data[ind].vel[2]	= gal0[0].vzc;

			}

			MPI_Barrier(MPI_COMM_WORLD);

			//Synchronize
			

    		for (CT_I32 ind : cut){
        		int owner = ind % size;

        		std::vector<std::uint8_t> blob_t;
        		if (rank == owner) blob_t = serialize(data[ind]);
        		bcast_blob_from_owner(owner, blob_t);
        		if (rank != owner) deserialize(blob_t, data[ind]);

	       	}
    		MPI_Barrier(MPI_COMM_WORLD);
#else 
			for(CT_I32 ind : cut){
				pid0.clear();
				cpid.clear();

				gal0 	= IO::r_gal(vh, data[ind].snap, data[ind].id, true);


				if(!istree(key, data[ind].snap, data[ind].id)){ 
					// no treee for this galaxy (read from a single snapshot)
					weight 	= get_weight(vh, gal0[0].pid);
					cpid2.resize(gal0[0].pid.size());
					for(CT_I32 k=0; k<(CT_I32) cpid2.size(); k++){
						cpid2[k].pid 	= gal0[0].pid[k];
						cpid2[k].weight 	= weight[k];

					}

					cpid 	= get_coreptcl(vh, cpid2);

					//pid0 = std::move(gal0[0].pid);

				} else {
					// this galaxy has a branch. Collect particles along its branch
					tree0 = Tree::gettree(tree, key, data[ind].snap, data[ind].id);

					ntcut = 0;
					t_slist.clear();
					t_idlist.clear();

					for(CT_I32 i=0; i<vh.ctree_n_search; i++){
						if(tree0.snap[i] >= data[ind].snap){
							ntcut ++;
							t_slist.push_back(tree0.snap[i]);
							t_idlist.push_back(tree0.id[i]);
						}
					}


					if(ntcut == 0){
						LOG()<<"Weird branch: Snap = "<<data[ind].snap0<<" / ID = "<<data[ind].id0;
						u_stop();
					}

					cpid 	= collectpidalongbranch(vh, t_slist, t_idlist);
				}


				// store to control array
				data[ind].p_list 	= std::move(cpid);
				data[ind].n_ptcl 	= data[ind].p_list.size();
				data[ind].pos[0]	= gal0[0].xc;
				data[ind].pos[1]	= gal0[0].yc;
				data[ind].pos[2]	= gal0[0].zc;

				data[ind].vel[0]	= gal0[0].vxc;
				data[ind].vel[1]	= gal0[0].vyc;
				data[ind].vel[2]	= gal0[0].vzc;
			}
#endif
		}



		// merge all particles
		std::vector<CT_I32> merge_cut;
		CT_I32 nmerge = 0;

		for(CT_I32 i=0; i<data[0].last_ind+1; i++){
			if(data[i].n_ptcl >= 0){
				merge_cut.push_back(i);
				nmerge ++;
			}
		}

		// no particles
		if(nmerge == 0){
			PIDArray pidout(1);
			pidout[0].pid = 0;
			pidout[0].gid = -1;
			pidout[0].weight = 0.;
			return pidout;
		}

		// Target Control Array
		std::vector<CT_I32> out_cut;
		CT_I32 nptcl = 0;



		for(CT_I32 i=0; i<data[0].last_ind+1; i++){
			if(data[i].stat == 0){
				out_cut.push_back(i);
				nptcl += data[i].n_ptcl;

			}
		}

		if(nptcl == 0){
			PIDArray pidout(1);
			pidout[0].pid = 0;
			pidout[0].gid = -1;
			pidout[0].weight = 0.;
			return pidout;
		}

		// Merge & Out
		PIDArray pidout(nptcl);

		CT_I32 i0, i1, maxgid;
		i0 = 0;
		i1 = 0;
		maxgid = -1;
		for(CT_I32 ind : out_cut){
			i1 	= i0 + data[ind].n_ptcl - 1;

			PIDReallocate(vh, pidout, i1);

			for(CT_I32 k=i0; k<i1+1; k++){
				pidout[k]	= data[ind].p_list[k-i0];

				pidout[k].gid 	= ind;
				if(ind >= maxgid) maxgid = ind;
			}
			i0 	= i1 + 1;
		}

		pidout.resize(i1+1);
		pidout[0].maxgid = maxgid;

		return pidout;
	}

	//-----
	// Read Snapshot particles
	//-----
	SnapPT readsnap(const vctree_set::Settings& vh, ControlArray& data, CT_I32 snap_curr){

		// Target control array
		SnapPT spt;
		std::vector<CT_I32> cut;
		CT_I32 ncut = 0;

		cut 	= get_control(data, 0);
		ncut 	= cut.size();
		
		// all tree ends; return null
		if(ncut==0){
			spt.pid.resize(1);
			spt.pid[0].pid 	= 0;
			spt.pid[0].gid 	= -1;
			spt.gid.push_back(-1);
			spt.hash.push_back(1);
			spt.hash_next.push_back(1);
			return spt;
		}
		
		// collect all particles of galaxies at this snapshot
		IO_dtype::GalArray gal = IO::r_gal(vh, snap_curr, -1, true);

		CT_I32 nptcl = 0;
		for(IO_dtype::GalSt g:gal){
			nptcl += g.npart;
		}

		//PIDArray pid(nptcl);
		spt.pid.resize(nptcl);
		spt.gid.resize(gal.size());
		spt.n_ptcl.resize(gal.size());


		CT_I32 i0, i1;
		i0 	= 0;
		i1 	= 0;
		//CT_PID maxpid=-1;
		CT_I32 maxgid=-1;

		for(CT_I32 i=0; i<(CT_I32) gal.size(); i++){
		//for(IO_dtype::GalSt g:gal){

			i1 	= i0 + gal[i].npart  - 1;
			for(CT_I32 j=0; j<(CT_I32) gal[i].npart; j++){
				spt.pid[j+i0].pid 	= gal[i].pid[j];
				spt.pid[j+i0].gid 	= gal[i].id;

				if(gal[i].id >= maxgid) maxgid = gal[i].id;
				//if(g.pid[i]>=maxpid) maxpid = g.pid[i];
			}
			spt.n_ptcl[i]	= gal[i].npart;
			spt.gid[i] 		= gal[i].id;
			i0 	= i1 + 1;
		}

		spt.maxgid 	= maxgid;
		//spt.maxpid 	= maxpid;

		//spt.pid 	= std::move(pid);

		// Make Hash
		CT_I32 dN = spt.pid.size();
		spt.hash.resize(dN);
		spt.hash_next.resize(dN);
		std::vector<CT_I32> hash_last(dN);

		for(CT_I32 i=0; i<dN; i++){
			spt.hash[i] = -1;
			spt.hash_next[i] = -1;
			hash_last[i] = -1;
		}

	
		CT_I64 pidtmp, ind;
		//noptcl 	= -92233720368547758;


		for(CT_I32 i=0; i<dN; i++){

			pidtmp 	= spt.pid[i].pid;

			//if(pidtmp < noptcl) continue;

			ind 	= std::abs(pidtmp) % dN;
			if(ind < 0) ind = 0;

			if(spt.hash[ind] < 0){
				spt.hash[ind] = i;
			} else{
				i0 	= spt.hash[ind];

				if(hash_last[i0] < 0){
					spt.hash_next[i0] = i;
					hash_last[i0] = i;
				} else{
					spt.hash_next[hash_last[i0]] = i;
					hash_last[i0] = i;
				}
			}
			
		}

//for(CT_I32 i=0; i<100; i++){
//	LOG()<<" --- "<<i<<" / "<<spt.hash[i]<<" / "<<spt.hash_next[i];
//}
		return spt;
	}

	//-----
	// Compute Merit
	//-----
	void commerit(const vctree_set::Settings& vh, ControlArray& data, PIDArray& pid, SnapPT& pid0, CT_I32 snap_curr){

		CT_I32 dummy = vh.snapi;
		dummy = dummy + 1;

		std::vector<CT_PID> pid_g, pid_s;
		std::vector<CT_ID> gid_g, gid_s;
		std::vector<CT_I32> hash, hash_next;
		CT_I32 n0, ncut;
		// No particles at this snapshot
		if(pid0.maxgid < 0 || pid[0].maxgid <0){
			std::vector<CT_I32> cut = get_control(data, 0);
			ncut = cut.size();


			if(ncut == 0) return;

			for(CT_I32 i=0; i<ncut; i++){
				n0 	= data[cut[i]].list_n;

				if(n0 >= vh.ctree_n_search){
					LOG()<<"indexing error occurs here";
					u_stop();
				}

				data[cut[i]].list[n0].merit 	= 0.;
				data[cut[i]].list[n0].id 		= -1;
				data[cut[i]].list[n0].snap 		= snap_curr;
				data[cut[i]].list_n ++;
			}
			return;
		}

		// Allocation for merit array
		std::vector<CT_I32> npart_g(data[0].last_ind+1, CT_I32{0}); //npart_g(pid[0].maxgid+1, CT_I32{0});
		std::vector<CT_I32> npart_s(pid0.maxgid+1, CT_I32{0});

		std::vector<CT_I32> cut = get_control(data, 0);
		ncut = cut.size();


		//extract
		pid_g.resize(pid.size());
		gid_g.resize(pid.size());

		pid_s.resize(pid0.pid.size());
		gid_s.resize(pid0.pid.size());
		hash.resize(pid0.pid.size());
		hash_next.resize(pid0.pid.size());

		for(CT_I32 i=0; i< (CT_I32) pid.size(); i++){
			pid_g[i]	= pid[i].pid;
			gid_g[i]	= pid[i].gid;
		}
		for(CT_I32 i=0; i<ncut; i++){
			npart_g[cut[i]] = data[cut[i]].n_ptcl;
		}


		for(CT_I32 i=0; i< (CT_I32) pid0.pid.size(); i++){
			pid_s[i]	= pid0.pid[i].pid;
			gid_s[i]	= pid0.pid[i].gid;
			hash[i]		= pid0.hash[i];
			hash_next[i]= pid0.hash_next[i];
			
		}
		for(CT_I32 i=0; i<(CT_I32) pid0.gid.size(); i++){
			npart_s[pid0.gid[i]] = pid0.n_ptcl[i];
		}

		// Merit Computation
		std::vector<CT_I32> match_id(ncut, CT_I32{-1});
		std::vector<CT_double> match_merit(ncut, CT_double{-1.});



		get_merit(pid_g, gid_g, pid_s, gid_s, hash, hash_next, npart_g, npart_s, match_id, match_merit);

		// input
		for(CT_I32 i=0; i<ncut; i++){

			n0 	= data[cut[i]].list_n;
			if(n0 >= vh.ctree_n_search){
					LOG()<<"indexing error occurs here";
					u_stop();
			}

			data[cut[i]].list[n0].merit 	= match_merit[i];
			data[cut[i]].list[n0].id 		= match_id[i];
			data[cut[i]].list[n0].snap 		= snap_curr;
			data[cut[i]].list_n ++;
		}

	}


	//-----
	// Link Branch
	//-----
	CT_double brcompare(const vctree_set::Settings& vh, CT_I32 s0, CT_I32 id0, std::vector<CT_I32>& slist, std::vector<CT_I32>& idlist){

		IO_dtype::GalArray gal0 	= IO::r_gal(vh, s0, id0, true);

		std::vector<CT_PID> pid0 = std::move(gal0[0].pid);
		std::vector<CT_double> pweight0 = get_weight(vh, pid0);

		CT_double n_occ;
		std::vector<CT_PID> pid1;
		std::vector<CT_double> pweight1;
		if( (CT_I32) slist.size() == 1 ){
			IO_dtype::GalArray gal1 	= IO::r_gal(vh, slist[0], idlist[0], true);
			pid1 = std::move(gal1[0].pid);
			pweight1 = get_weight(vh, pid1);
			n_occ = 1.;
		}else{
			PIDArray cpid = collectpidalongbranch(vh, slist, idlist);
			pid1.resize(cpid.size());
#ifdef CTREE_USE_OMP
			#pragma omp parallel for default(none) \
    			shared(pid1, cpid)
#endif
			for(CT_I32 i=0; i< (CT_I32) pid1.size(); i++){
				pid1[i] 	= cpid[i].pid;
			}
			
			pweight1 = get_weight(vh, pid1);

			n_occ 	= (CT_double) cpid[0].n_con;
		}

		CT_I32 merittype = 1;	// pre-selected to be 1
		CT_double meritdum = get_merit2(pid0, pid1, pweight0, pweight1, merittype);

		//-----
		// Factor calculation
		//	- should be 1 if all particles appear during the selected part of the branch
		// 	- higher occurance should have a stronger weight

		CT_double factor = ((CT_double) 1.) / std::pow( ((CT_double) vh.ctree_n_step_n * vh.ctree_n_step_dn),2.) * std::pow( ((CT_double )n_occ * vh.ctree_n_step_dn),2.);

		return factor*meritdum;
		
	}

	void linkbr(const vctree_set::Settings& vh, ControlArray& data, ControlKey& dkey, CT_I32 ind, IO::snapinfo& sinfo, Tree::TreeArray& tree, Tree::TreeKeyArray& key, CT_I32 id_to_link, CT_I32 snap_to_link, CT_double merit_to_link, CT_I32 snap_curr){

		bool istree = Tree::istree(key, data[ind].snap, data[ind].id);
		if(!istree){
			LOG()<<" why this: should not be happened"<<" / "<<ind;
			u_stop();

			makenewbr(vh, data, ind, data[ind].snap, data[ind].id, tree, key);
		}

		CT_I32 dummy = wheresnap(sinfo, snap_curr);
		dummy += 1;

		istree	= Tree::istree(key, snap_to_link, id_to_link);
		if(!istree){
			LOG()<< "Branch is stolen somehow";
			u_stop();
		}

		Tree::TreeSt tmp_tree 		= Tree::gettree(tree, key, data[ind].snap, data[ind].id);
		Tree::TreeSt tmp_tree_toc	= Tree::gettree(tree, key, snap_to_link, id_to_link);

		Tree::Tree_BID org_bid 	= Tree::getkey(key, data[ind].snap, data[ind].id);
		Tree::Tree_BID com_bid 	= Tree::getkey(key, snap_to_link, id_to_link);


		CT_I32 dind = dkey[tmp_tree_toc.snap[tmp_tree_toc.endind] + dkey[0] * tmp_tree_toc.id[tmp_tree_toc.endind]];
		if(dind>0 && (data[dind].id0 != tmp_tree_toc.id[tmp_tree_toc.endind] || data[dind].snap0 != tmp_tree_toc.snap[tmp_tree_toc.endind])){
			LOG()<<"?? : "<<dind<<" / "<<com_bid;
			LOG()<<data[dind].id0<<" / "<<data[dind].snap0;
			LOG()<<tmp_tree_toc.id[tmp_tree_toc.endind]<<" / "<<tmp_tree_toc.snap[tmp_tree_toc.endind];
			u_stop();
		}

		// Clean Merge
		if(tmp_tree_toc.snap[tmp_tree_toc.endind] < tmp_tree.snap[0]){ // end before the start: complete merge


			//if( !(tmp_tree_toc.id[tmp_tree_toc.endind] == id_to_link && tmp_tree_toc.snap[tmp_tree_toc.endind] == snap_to_link )){
			//	LOG()<<"Just check whether this is possible";
			//	u_stop();
			//}

			
			tmp_tree_toc.p_merit[0] 	= merit_to_link;

			// free the com- branch
			Tree::treefree(tree, key, snap_to_link, id_to_link);
			
			// input new info
			for(CT_I32 i=tmp_tree_toc.endind; i>=0; i--){
				Tree::treeinput(tree, key, tmp_tree.snap[0], tmp_tree.id[0], tmp_tree_toc.snap[i], tmp_tree_toc.id[i], tmp_tree_toc.p_merit[i]);
			}
			// Control free
			ctfree(vh, data, ind, tmp_tree_toc.snap[0], tmp_tree_toc.id[0], snap_curr);

			if(dind>=0){
				ctfree(vh, data, dind, -1, -1, snap_curr);
				data[dind].stat 	= -1;
			}
			
			return;
		}

		// Dirty Merge (two branches both survive simultaneously)
		std::vector<CT_I32> brorg_id(tmp_tree.endind+1), brorg_snap(tmp_tree.endind+1);
		std::vector<CT_I32> brcom_id(tmp_tree_toc.endind+1), brcom_snap(tmp_tree_toc.endind+1);
		CT_I32 org_n = 0;
		CT_I32 com_n = 0;

		//// Gather snap & ID
		for(CT_I32 i=0; i<tmp_tree.endind+1; i++){
			if(tmp_tree.snap[i] > snap_curr + vh.ctree_n_step_dn){
				brorg_id[org_n] 	= tmp_tree.id[i];
				brorg_snap[org_n]  	= tmp_tree.snap[i];
				org_n ++;
			}
		}
		if(org_n == 0){
			for(CT_I32 i=0; i<tmp_tree.endind+1; i++){
				if(tmp_tree.snap[i] > snap_curr){
					brorg_id[org_n] 	= tmp_tree.id[i];
					brorg_snap[org_n]  	= tmp_tree.snap[i];
					org_n ++;
				}
			}
		}

		brorg_id.resize(org_n);
		brorg_snap.resize(org_n);

		for(CT_I32 i=0; i<tmp_tree_toc.endind+1; i++){
			if(tmp_tree_toc.snap[i] > snap_curr + vh.ctree_n_step_dn){
				brcom_id[com_n]		= tmp_tree_toc.id[i];
				brcom_snap[com_n]	= tmp_tree_toc.snap[i];
				com_n ++;
			}
		}
		if(com_n == 0){
			for(CT_I32 i=0; i<tmp_tree_toc.endind+1; i++){
				if(tmp_tree_toc.snap[i] > snap_curr){
					brcom_id[com_n]		= tmp_tree_toc.id[i];
					brcom_snap[com_n]	= tmp_tree_toc.snap[i];
					com_n ++;
				}
			}
		}
		brcom_id.resize(com_n);
		brcom_snap.resize(com_n);

		CT_double merit_com = brcompare(vh, snap_to_link, id_to_link, brcom_snap, brcom_id);
		CT_double merit_org = brcompare(vh, snap_to_link, id_to_link, brorg_snap, brorg_id);

		// existing branch is better
		//Tree::Tree_I64 keyval, keyval2;
		if(merit_com > merit_org){
			data[ind].stat = -1;

			Tree::TreeSt& tree0 = tree[org_bid];
			tree0.stat = -2;
			tree0.frag_bid 	= com_bid;

			//Tree::modifytree(tree, key, data[ind].snap0, data[ind].id0, snap_to_link);
			ctfree(vh, data, ind, -1, -1, snap_curr);
			return;
		}

		// New branch is better
		//// adjust free
		Tree::modifytree(tree, key, snap_to_link, id_to_link, snap_to_link);
		//ctfree(vh, data, dind, )


		//// move tree from com to org
		for(CT_I32 i=tmp_tree_toc.endind; i>=0; i--){
			if(tmp_tree_toc.snap[i] > snap_to_link) continue;
			Tree::treeinput(tree, key, tmp_tree.snap[0], tmp_tree.id[0], tmp_tree_toc.snap[i], tmp_tree_toc.id[i], tmp_tree_toc.p_merit[i]);
		}


		//// control
		ctfree(vh, data, ind, tmp_tree_toc.snap[0], tmp_tree_toc.id[0], snap_curr);

		if(dind<0) return; // no control array

		//// Check whether data[dind] to close or not
		if(!Tree::istree(key, data[dind].snap0, data[dind].id0)){ // one snapshot branch case
			data[dind].stat = -1;
			ctfree(vh, data, dind, -1, -1, snap_curr);
			return;
		} else{
			Tree::TreeSt& modtree = tree[com_bid];

			if(modtree.snap[0] <= snap_curr){
				data[dind].stat = 1;
				ctfree(vh, data, dind, modtree.snap[0], modtree.id[0], snap_curr);
			} else{
				CT_I32 index0, index1;

				index0	= wheresnap(sinfo, snap_curr);
				index1 	= wheresnap(sinfo, modtree.snap[0]);
				if(std::abs(index0-index1) > vh.ctree_n_search){
					//modtree.snap[0] > snap_curr + vh.ctree_n_search){
					data[dind].stat = -1;
					ctfree(vh, data, dind, -1, -1, snap_curr);
					Tree::TreeSt& tree0 = tree[com_bid];
					tree0.stat = -2;
					tree0.frag_bid 	= org_bid;
				}else{
					data[dind].stat = 0;
					ctfree(vh, data, dind, modtree.snap[0], modtree.id[0], snap_curr);
				}
			}
			return;
		}

	}

	void link(const vctree_set::Settings& vh, ControlArray& data, ControlKey& dkey, Tree::TreeArray& tree, Tree::TreeKeyArray& key, IO::snapinfo& sinfo, CT_I32 snap_curr){



		CT_I32 ind, ind0, ind1, snap_int_cut;

		ind1 	= wheresnap(sinfo, snap_curr);
		ind0 	= wheresnap(sinfo, (CT_I32) vh.snapi);
		snap_int_cut 	= ind1-ind0-1;
		if(snap_int_cut >= vh.ctree_n_search) snap_int_cut = vh.ctree_n_search;
		//ind1 	= snap_curr;
		//ind0 	= vh.snapi;
		//snap_int_cut	= ind1-ind0-1;
		//if(snap_int_cut >= vh.ctree_n_search) snap_int_cut = vh.ctree_n_search;
		//----- Extract Target Control whose list is fully filled
		std::vector<CT_I32> cut(data[0].last_ind+1);
		CT_I32 ncut = 0;
		for(CT_I32 i=0; i<data[0].last_ind+1; i++){
			if(data[i].stat == 0 && data[i].list_n >= snap_int_cut && data[i].list_n > 0){
				cut[ncut] 	= i;
				ncut ++;
			}
		}

		if(ncut == 0) return;
		cut.resize(ncut);

		//----- Extract next points of the extracted controls
		NextArray next_point(ncut);

		CT_double dum_merit;
		

#ifdef CTREE_USE_OMP
		#pragma omp parallel for default(none) \
			private(ind, dum_merit) shared(ncut, cut, data, vh, next_point)
#endif
		for(CT_I32 i=0; i<ncut; i++){
			ind 	= cut[i];

			dum_merit = -1.;
			for(CT_I32 j=0; j<vh.ctree_n_search; j++){

				if(data[ind].list[j].merit >= dum_merit){
					dum_merit 	= data[ind].list[j].merit;

					next_point[i].id 	= data[ind].list[j].id;
					next_point[i].snap 	= data[ind].list[j].snap;
					next_point[i].merit = data[ind].list[j].merit;

				}
			}
		}

//#endif

		// Gathering all checkpoints
		std::vector<CT_I32> data_ind(data[0].last_ind+1);
		CT_I32 nind1 = 0;
		CT_I32 nall = 0;

		// include not finished branch
		for(CT_I32 i=0; i<data[0].last_ind+1; i++){
			if(data[i].list_n >= 1){
				data_ind[nind1] 	= i;
				nind1 ++;
				nall 	+= data[i].list_n;
			}
		}
		data_ind.resize(nind1);

		CheckArray checkarr(nall);


		CT_I32 i0 = 0;
		CT_I32 i1;

		for(CT_I32 i=0; i<nind1; i++){
			i1 	= i0 + data[data_ind[i]].list_n - 1;

			for(CT_I32 j=i0; j<i1+1; j++){
				checkarr[j].merit 	= data[data_ind[i]].list[j-i0].merit;
				checkarr[j].id 		= data[data_ind[i]].list[j-i0].id;
				checkarr[j].snap 	= data[data_ind[i]].list[j-i0].snap;
				checkarr[j].id0 	= data[data_ind[i]].id0;
				checkarr[j].snap0	= data[data_ind[i]].snap0;
				checkarr[j].ind 	= data_ind[i];
			}

			i0 	= i1 + 1;
		}


		// Decide Connectivity
		// Islink
		//	-1 	: low meirt or no link
		//  -2  : other branch has higher merit (fragmented)

		std::vector<CT_I32> islink(ncut, CT_I32{1});
		std::vector<CT_I32> ischeck(nall, CT_I32{0});
		CheckArray checkcon(nall);

		CT_I32 nischeck;
		bool istree;
		Tree::TreeSt tree0;

		std::vector<CT_snap> tt_snap;
		std::vector<CT_ID> tt_id;
		CT_I32 ntt;

		CT_double this_merit, other_merit;


#ifdef CTREE_USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		int rank = 0, size = 1;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
		for(CT_I32 i=0; i<ncut; i++){
			ind 	= cut[i];

#ifdef CTREE_USE_MPI
			int owner = i % size;
			if(owner != rank) continue;
#endif

			if(next_point[i].merit < vh.meritlimit){
				islink[i] = -1;
				continue;
			}
			nischeck = 0;
			for(CT_I32 j=0; j<nall; j++){
				if(checkarr[j].snap == next_point[i].snap && checkarr[j].id == next_point[i].id && checkarr[j].merit > vh.meritlimit){
					ischeck[nischeck] = j;
					nischeck ++;
				}
			}

			if(nischeck == 0){ // no link
				islink[i] = -1;
			} else if(nischeck >= 2){ // connection is overapped

				// merit comparison
				for(CT_I32 j=0; j<nischeck; j++){
					checkcon[j] 	= checkarr[ischeck[j]];

					istree 	= Tree::istree(key, checkcon[j].snap0, checkcon[j].id0);
					if(!istree){
						LOG()<< " no tree !?";
						LOG()<<" / "<<checkcon[j].snap0<<" / "<<checkcon[j].id0;
						u_stop();
					}

					tree0 	= Tree::gettree(tree, key, checkcon[j].snap0, checkcon[j].id0);

					tt_snap.resize(0);
					tt_id.resize(0);
					ntt = 0;
					for(CT_I32 k=0; k<tree0.endind+1; k++){
						if(tree0.snap[k] > checkcon[j].snap + vh.ctree_n_step_dn){
							tt_snap.push_back(tree0.snap[k]);
							tt_id.push_back(tree0.id[k]);
							ntt ++;
						}
					}

					if(ntt == 0){
						for(CT_I32 k=0; k<tree0.endind+1; k++){
							if(tree0.snap[k] > checkcon[j].snap){
								tt_snap.push_back(tree0.snap[k]);
								tt_id.push_back(tree0.id[k]);
								ntt ++;
							}
						}
					}
					

					checkcon[j].merit 	= brcompare(vh, checkcon[j].snap, checkcon[j].id, tt_snap, tt_id);

					
				}

				// merit comparison2
				this_merit = -1.;
				other_merit = -1.;

				for(CT_I32 j=0; j<nischeck; j++){

					if(checkcon[j].id0 == data[ind].id0 && checkcon[j].snap0 == data[ind].snap0){
						this_merit 	= checkcon[j].merit;
					}else{
						if(checkcon[j].merit > other_merit){
							other_merit 	= checkcon[j].merit;
						}
					}
				}


				if(this_merit < 0){
					LOG()<<"Weird merit matching!?";
					u_stop();
				}

				// determine
				if(this_merit < other_merit){
					islink[i] = -2;

					//next_point[i].id 	= -1;
					//next_point[i].snap 	= -1;
					continue;
				}
			} else{
				islink[i] = 1;
			}
		}

#ifdef CTREE_USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);

		//----- Synchronize islink
		for(CT_I32 i=0; i<ncut; i++){
			int owner = i % size;
			
			std::vector<std::uint8_t> blob_t;
			if (rank == owner) blob_t = serialize(islink[i]);
			bcast_blob_from_owner(owner, blob_t);
	        if (rank != owner) deserialize(blob_t, islink[i]);

		}
#endif

		// Close controls with islink < 0
		Tree::Tree_I64 keyval;
#ifdef CTREE_USE_OMP
		#pragma omp parallel for default(none) \
			private(tree0, keyval) \
			shared(ncut, islink, data, cut, next_point, snap_curr, vh, key, tree)
#endif
		for(CT_I32 i=0; i<ncut; i++){
			if(islink[i] < 0){

				data[cut[i]].stat = -1;

				keyval 	= data[cut[i]].snap0 + key[0].key * data[cut[i]].id0;
				Tree::TreeSt& tree0 = tree[keyval];



				if(islink[i] == -1){
					tree0.stat = -1;
				}else if(islink[i] == -2){
					tree0.stat 	= -2;

					keyval = next_point[i].snap + key[0].key * next_point[i].id;
					tree0.frag_bid 	= key[keyval].ind;

				}


				next_point[i].id 	= -1;
				next_point[i].snap 	= -1;

				ctfree(vh, data, cut[i], -1, -1, snap_curr);

			}
		}

//456456
		// Link to a next checkpoint
		CT_I32 snap_to_link;
		CT_I32 id_to_link;
		CT_double merit_to_link;

		// parallelizatin here..
		for(CT_I32 i=0; i<ncut; i++){
			if(islink[i] < 0) continue;

			ind 	= cut[i];
//if(myrank == 0){
//	if(!Tree::istree(key, data[461].snap0, data[461].id0)){
//		CT_I32 k1 = Tree::getkey(key, data[461].snap0, data[461].id0);
//		CT_I32 k2 = Tree::getkey(key, data[461].snap, data[461].id);
//
//		LOG()<<"data ind: "<<ind<<" / "<<oldind;
//		LOG()<<"key check: "<<k1<<" / "<<k2;
//		LOG()<<"data info: "<<data[461].snap0<<" / "<<data[461].id0<<" / "<<data[461].stat;
//		LOG()<<"data info: "<<data[461].snap<<" / "<<data[461].id<<" / "<<data[461].stat;
//	}
//
//	if(!Tree::istree(key, data[461].snap, data[461].id)){
//		LOG()<<" -- "<<" / "<<snap_curr;
//		CT_I32 k1 = Tree::getkey(key, data[461].snap0, data[461].id0);
//		CT_I32 k2 = Tree::getkey(key, data[461].snap, data[461].id);
//
//		LOG()<<"data ind: "<<ind<<" / "<<oldind;
//		LOG()<<"key check: "<<k1<<" / "<<k2;
//		LOG()<<"data info: "<<data[461].snap0<<" / "<<data[461].id0<<" / "<<data[461].stat;
//		LOG()<<"data info: "<<data[461].snap<<" / "<<data[461].id<<" / "<<data[461].stat;
//	}
//}
			if(data[ind].stat == -1) continue;		// already linked to another

			snap_to_link 	= next_point[i].snap;
			id_to_link 		= next_point[i].id;
			merit_to_link	= next_point[i].merit;
 
			istree 	= Tree::istree(key, snap_to_link, id_to_link);

			if(!istree){ // notree
				expandbr(vh, data, ind, tree, key, id_to_link, snap_to_link, merit_to_link);

				if(data[ind].stat == -1){
					LOG()<<"Why?";
					u_stop();
				}else{
					ctfree(vh, data, ind, snap_to_link, id_to_link, snap_curr);
				}

			}else{ // tree exist (compare two branches)
				linkbr(vh, data, dkey, ind, sinfo, tree, key, id_to_link, snap_to_link, merit_to_link, snap_curr);
			}



		}



	}

	//
	void addgal(const vctree_set::Settings& vh, ControlArray& data, ControlKey& dkey, Tree::TreeArray& tree, Tree::TreeKeyArray& key, CT_I32 snap_curr){

		IO_dtype::GalArray gal0 = IO::r_gal(vh, snap_curr, -1, false);
		IO_dtype::GalArray gal(gal0.size());

		CT_I32 nn = 0;

		for(CT_I32 i=0; i<(CT_I32) gal0.size(); i++){


			if(!Tree::istree( key, gal0[i].snap, gal0[i].id )){
				gal[nn] 	= gal0[i];
				nn ++;

				continue;
			} else{

				Tree::TreeSt tree0 = Tree::gettree(tree, key, gal0[i].snap, gal0[i].id);


				if( dkey[ tree0.snap[tree0.endind] + dkey[0]*tree0.id[tree0.endind] ] < 0){
					gal[nn]		= gal0[i];
					nn ++;
					continue;
				}else{
					continue;
				}	
			}
		}

		gal.resize(nn);

		inputgal(vh, tree, key, data, dkey, gal);

	}

	// finalize
	void finalize(const vctree_set::Settings& vh, ControlArray& data, ControlKey& dkey, Tree::TreeArray& tree, Tree::TreeKeyArray& key, IO::snapinfo& sinfo, CT_I32 snap_curr, ctree_num& number){

		std::vector<CT_I32> cut 	= get_control(data, 0);
		CT_I32 ncut = cut.size();

		if(ncut == 0) return;

		CT_I32 count_n = 0;
		while(true){
			link(vh, data, dkey, tree, key, sinfo, snap_curr);
			classify(vh, data, sinfo, snap_curr, number);

			bool isout = false;
			for(CT_I32 i=0; i<data[0].last_ind+1; i++){
				if(data[i].list_n > 0) isout = true;
			}

			count_n ++;
			if(isout) break;
			if(count_n>10) break;
		}
	}

	// Main Part

	void main(const vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& key){
// Check
// 	- different merit quantification but compared simultaneously (is this fine?)

		int myrank  = mpi_rank();
		double dt_classify, dt_collectpid, dt_readsnap, dt_commerit, dt_link, dt_addgal;

		//-----
		//tree.resize(tree[0].endind*2);

		//-----
		// Get Basic Information in each snapshot
		//-----
		if(myrank==0)LOG() <<"    Ctree) Reading snapshot info with "<<vh.simtype<<" format";
		IO::snapinfo sinfo = IO::get_snapinfo(vh);

		//-----
		// Branch Controller
		//-----
		IO_dtype::GalArray gal = IO::r_gal(vh, vh.snapf, -1, false);
		ControlArray data = allocate(vh, gal.size());
		ControlKey dkey(key.size(), -1);
		dkey[0] 	= key[0].key;
		// dkey is only called by snap0 and id0


		if(myrank==0)LOG() <<"    Ctree) Initial allocating for "<<gal.size()<<" galaxies";
		inputgal(vh, tree, key, data, dkey, gal);

		//-----
		// Main Loop
		//-----
		ctree_num number;
		if(myrank==0)LOG() <<"    Ctree) Main loop starts";
		for(CT_I32 i=sinfo.size()-1; i>=0; i--){
			if(sinfo[i].snum<0) continue;

			if(myrank==0)LOG() <<"      TREE CONNECTION AT SNAP "<<i4(sinfo[i].snum)<<" for "<<i6(data[0].last_ind+1)<<" galaxies";

			if(myrank==0){
				LOG()<<"        Memory report";
				LOG()<<"			"<<how_big<Tree::TreeSt>(tree)<<" GB for tree";
				LOG()<<"			"<<how_big<Tree::TreeKey>(key)<<" GB for key";
				LOG()<<"			"<<how_big<ControlSt>(data)<<" GB for data";
				LOG()<<"			"<<how_big<CT_I32>(dkey)<<" GB for dkey";
				LOG()<<" ";
				LOG()<<"		"<<tree[0].lind<<" / "<<tree.size();

			}
			//----- Classify Tree

			auto t0 = std::chrono::steady_clock::now();
			classify(vh, data, sinfo, sinfo[i].snum, number);
			if(myrank == 0){
				auto t1 = std::chrono::steady_clock::now();
				dt_classify = std::chrono::duration<double>(t1 - t0).count();
			}



			if(myrank==0)LOG() <<"        With Tree = "<< number.T;
			if(myrank==0)LOG() <<"        To be connected = "<< number.C;
			if(myrank==0)LOG() <<"        Broken = "<< number.B;
			if(number.T == data[0].last_ind+1){
				if(myrank == 0) LOG() <<"          SKIP this snapshot because all galaxies have tree";
				continue;
			}
#ifdef CTREE_USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif

			//----- Collect PID
			t0 = std::chrono::steady_clock::now();
			PIDArray pid 	= collectpid(vh, data, tree, key);
			if(myrank == 0){
				auto t1 = std::chrono::steady_clock::now();
				dt_collectpid = std::chrono::duration<double>(t1 - t0).count();
			}

			//----- Read particles at this snapshot
			t0 = std::chrono::steady_clock::now();
			SnapPT pid0 	= readsnap(vh, data, sinfo[i].snum);
			if(myrank == 0){
				auto t1 = std::chrono::steady_clock::now();
				dt_readsnap = std::chrono::duration<double>(t1 - t0).count();
			}

			//----- Compute Merit
			t0 = std::chrono::steady_clock::now();
			commerit(vh, data, pid, pid0, sinfo[i].snum);
			if(myrank == 0){
				auto t1 = std::chrono::steady_clock::now();
				dt_commerit = std::chrono::duration<double>(t1 - t0).count();
			}

			//----- Link
			t0 = std::chrono::steady_clock::now();
			link(vh, data, dkey, tree, key, sinfo, sinfo[i].snum);
			if(myrank == 0){
				auto t1 = std::chrono::steady_clock::now();
				dt_link = std::chrono::duration<double>(t1 - t0).count();
			}


			//----- Add Galaxies starting from this snapshot
			t0 = std::chrono::steady_clock::now();
			addgal(vh, data, dkey, tree, key, sinfo[i].snum);
			if(myrank == 0){
				auto t1 = std::chrono::steady_clock::now();
				dt_addgal = std::chrono::duration<double>(t1 - t0).count();
			}

			//----- Time Log
			if(myrank == 0){
				LOG()<<"        Time report";
				LOG()<<"          Classify Galaxies in   "<<dt_classify;
				LOG()<<"          Collect PID in         "<<dt_collectpid;
				LOG()<<"          Read Snapshot ptcls in "<<dt_readsnap;
				LOG()<<"          Compute merit in       "<<dt_commerit;
				LOG()<<"          Link branch in         "<<dt_link;
				LOG()<<"          Add new galaxies in    "<<dt_addgal;
			}

#ifdef CTREE_USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);	
#endif
		}

		classify(vh, data, sinfo, vh.snapi, number);
		finalize(vh, data, dkey, tree, key, sinfo, vh.snapi, number);
		LOG()<<" finialize";
		LOG()<<" merging ?";

		//std::vector<IO_VR::VRT_Snap> snaplist;
		//IO_VR::VRT_Snap maxsnap = -1;

		//for(IO_VR::VRT_Snap i=vh.snapi; i<vh.snapf+1; i++){
		//	if(is_snap(vh, i)){
		//		snaplist.push_back(i);
		//		if(i > maxsnap) maxsnap = i;
		//	}
		//}


		//auto snapinfo 	= IO::g_snapinfo(vh);
		//IO_VR::GalArray gal = IO_VR::r_gal(vh, 100, -1, false);
	}

	    //----
    // MPI Helper
    //-----
#ifdef CTREE_USE_MPI
    //---- raw POD append / read ----
	template<typename T>
	void append_pod(std::vector<std::uint8_t>& buf, const T& v) {
	    static_assert(std::is_trivially_copyable<T>::value, "POD only");
	    const auto* p = reinterpret_cast<const std::uint8_t*>(&v);
	    buf.insert(buf.end(), p, p + sizeof(T));
	}

	template<typename T>
	T read_pod(const std::uint8_t*& p, const std::uint8_t* end) {
	    static_assert(std::is_trivially_copyable<T>::value, "POD only");
	    if (size_t(end - p) < sizeof(T)) throw std::runtime_error("deserialize: buffer underflow");
	    T v;
	    std::memcpy(&v, p, sizeof(T));
	    p += sizeof(T);
	    return v;
	}

	template<typename T>
	void append_vec(std::vector<std::uint8_t>& buf, const std::vector<T>& v) {
	    std::uint64_t n = static_cast<std::uint64_t>(v.size());
	    append_pod(buf, n);
	    if (n) {
	        static_assert(std::is_trivially_copyable<T>::value, "vector element must be POD");
	        const auto* p = reinterpret_cast<const std::uint8_t*>(v.data());
	        buf.insert(buf.end(), p, p + n * sizeof(T));
	    }
	}

	template<typename T>
	std::vector<T> read_vec(const std::uint8_t*& p, const std::uint8_t* end) {
	    std::uint64_t n = read_pod<std::uint64_t>(p, end);
	    std::vector<T> v;
	    if (n) {
	        if (size_t(end - p) < n * sizeof(T)) throw std::runtime_error("deserialize vec: buffer underflow");
	        v.resize(static_cast<size_t>(n));
	        std::memcpy(v.data(), p, n * sizeof(T));
	        p += n * sizeof(T);
	    }
	    return v;
	}	

    // Serialize & Deserialize for Control Array
	std::vector<std::uint8_t> serialize(const ControlSt& x) {
	    std::vector<std::uint8_t> buf;

	    append_pod<CT_ID>(buf, x.id0);
	    append_pod<CT_snap>(buf, x.snap0);
	    append_pod<CT_ID>(buf, x.id);
	    append_pod<CT_snap>(buf, x.snap);
	    
	    //append_vec<CT_double>(buf, x.pos);
	    //append_vec<CT_double>(buf, x.vel);
	    
	    append_pod<CT_I32>(buf, x.stat);
	    append_pod<CT_I32>(buf, x.detstat);

	    append_vec<PIDSt>(buf, x.p_list);
		append_pod<CT_I32>(buf, x.n_ptcl);

		append_vec<ListSt>(buf, x.list);
		append_pod<CT_I32>(buf, x.list_n);

		append_pod<CT_I32>(buf, x.last_ind);
	    return buf;


	}

	void deserialize(const std::vector<std::uint8_t>& buf, ControlSt& out) {
	    const std::uint8_t* p   = buf.data();
	    const std::uint8_t* end = p + buf.size();

	    out.id0 		= read_pod<CT_ID>(p, end);
	    out.snap0 		= read_pod<CT_snap>(p, end);
	    out.id 			= read_pod<CT_ID>(p, end);
	    out.snap 		= read_pod<CT_snap>(p, end);

	    //out.pos 		= read_vec<CT_double>(p, end);
	    //out.vel 		= read_vec<CT_double>(p, end);

	    out.stat 		= read_pod<CT_I32>(p, end);
	    out.detstat		= read_pod<CT_I32>(p, end);

	    out.p_list 		= read_vec<PIDSt>(p, end);
	    out.n_ptcl 		= read_pod<CT_I32>(p, end);

	    out.list 		= read_vec<ListSt>(p, end);
	    out.list_n 		= read_pod<CT_I32>(p, end);

	    out.last_ind 	= read_pod<CT_I32>(p, end);

	    if (p != end) throw std::runtime_error("TFSt deserialize: trailing bytes");
	}

	// Serialize & Deserialize for Next Array
	std::vector<std::uint8_t> serialize(const NextSt& x) {
	    std::vector<std::uint8_t> buf;

	    append_pod<CT_ID>(buf, x.id);
	    append_pod<CT_snap>(buf, x.snap);
	    append_pod<CT_double>(buf, x.merit);
	    return buf;
	}

	void deserialize(const std::vector<std::uint8_t>& buf, NextSt& out) {
	    const std::uint8_t* p   = buf.data();
	    const std::uint8_t* end = p + buf.size();

	    out.id 			= read_pod<CT_ID>(p, end);
	    out.snap 		= read_pod<CT_snap>(p, end);
	    out.merit		= read_pod<CT_double>(p, end);

	    if (p != end) throw std::runtime_error("TFSt deserialize: trailing bytes");
	}

	// Serialize & Deserialize for Islink Array
	std::vector<std::uint8_t> serialize(const CT_I32& x) {
		std::vector<std::uint8_t> buf;
		append_pod<CT_I32>(buf, x);
		return buf;
	}
	void deserialize(const std::vector<std::uint8_t>& buf, CT_I32& out) {
	    const std::uint8_t* p   = buf.data();
	    const std::uint8_t* end = p + buf.size();

	    out 	= read_pod<CT_I32>(p, end);

	    if (p != end) throw std::runtime_error("TFSt deserialize: trailing bytes");
	}


	// Broad Cast Helper
	void bcast_blob_from_owner(int owner, std::vector<std::uint8_t>& blob) {
	    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	    int len = (rank == owner) ? static_cast<int>(blob.size()) : 0;
	    MPI_Bcast(&len, 1, MPI_INT, owner, MPI_COMM_WORLD);

	    if (rank != owner) blob.resize(len);
	    if (len > 0) {
	        MPI_Bcast(blob.data(), len, MPI_BYTE, owner, MPI_COMM_WORLD);
	    }
	}

#endif
}