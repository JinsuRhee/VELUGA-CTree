
#include "global/io.h"
#include "complete_tree/ctree.h"
#include "utilities/utilities.h"

namespace Ctree{

	//-----
	// ETC
	//-----
	void expandbr(const vctree_set::Settings& vh, ControlArray& data, CT_I32 ind, Tree::TreeArray& tree, Tree::TreeKeyArray& key, CT_I32 id_to_link, CT_I32 snap_to_link, CT_Merit merit_to_link){
		
		// this should be snap0 & id0 in mpi parallelization
		bool istree = Tree::istree(key, data[ind].snap0, data[ind].id0);


		if(!istree){
			int rank = 0, size = 1;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	    	MPI_Comm_size(MPI_COMM_WORLD, &size);

	    	for(int i=0; i<size; i++){
	    		if(rank == i){
	    			LOG()<<" why this: should not be happened : "<<rank<<" / "<<ind<<" / "<<data[ind].snap<<" / "<<data[ind].id;
	    		}
	    		MPI_Barrier(MPI_COMM_WORLD);
	    	}
			
			u_stop();

			makenewbr(vh, data, ind, data[ind].snap, data[ind].id, tree, key);
		}

		//Tree::TreeSt tree0 	= Tree::gettree(tree, key, data[ind].snap, data[ind].id);

		Tree::treeinput(tree, key, data[ind].snap0, data[ind].id0, snap_to_link, id_to_link, merit_to_link);
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


		Tree::Tree_BID keyval;
		//keyval 	= snap0 + key[0].key * id0;
		keyval 	= snap0 + key[0] * id0;

		if((Tree::Tree_BID) key.size() <= keyval){
			key.resize( (Tree::Tree_BID) key.size() + vh.ctree_nstep);
		}

		//key[keyval].ind 	= tree[0].lind;
		key[keyval] 	= tree[0].lind;

		data[ind].id 	= id0;
		data[ind].snap 	= snap0;



	}




	void get_merit(std::vector<CT_PID>& pid_g, std::vector<CT_I32>& gid_g, 
		std::vector<CT_PID>& pid_s, std::vector<CT_I32>& gid_s, 
		std::vector<CT_I32>& hash, std::vector<CT_I32>& hash_next, 
		std::vector<CT_I32>& npart_g, std::vector<CT_I32>& npart_s, 
		std::vector<CT_I32>& match_id, std::vector<CT_Merit>& match_merit){

		// local variables
		CT_I32 n_pg = pid_g.size();
		//CT_I32 n_ps = pid_s.size();
		CT_I32 n_g 	= npart_g.size();
		CT_I32 n_s 	= npart_s.size();
		CT_I32 n_dn = pid_s.size();
		//CT_I32 n_hashth	= hash.size();


		// Make Index Table
LOG()<<" n_g & n_s : "<<n_g<<" / "<<n_s;
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
LOG()<<" nn_g & nn_s : "<<nn_g<<" / "<<nn_s;
		MeritMatrix<CT_Merit> merit(nn_g , nn_s, CT_Merit{0});

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
LOG()<<" n_g: "<<n_g;
		std::vector<CT_I32> idtoind(n_g, CT_I32{0});

		for(CT_I32 i=0; i<n_g; i++){
			if(npart_g[i] > 0){
				idtoind[i] = ind;
				ind += 1;
			}
		}



		CT_Merit dum;
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
	CT_Merit get_merit2(std::vector<CT_PID>& pid0, std::vector<CT_PID>& pid, std::vector<CT_Merit>& w0, std::vector<CT_Merit>& w1, CT_I32 merittype){

		CT_I32 n_raw 	= pid0.size();
		CT_I32 n_match 	= pid.size();
		CT_Merit share = 0.;
		CT_Merit weight1 = 0.;
		CT_Merit weight2 = 0.;

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


		if(merittype == 1) share 	= (share*share) / ((CT_Merit) n_raw) / ( (CT_Merit) n_match );
		if(merittype == 2) share 	= share / ( (CT_Merit) n_raw );
		if(merittype == 3) share 	= (share*share) / ((CT_Merit) n_raw) / ( (CT_Merit) n_match ) * weight1 * weight2;

		return share;
	}

// for memory efficient ver
	MeritSt get_merit3(SnapPT& pid0, PIDArray& pid){

		CT_I32 np_snap = pid0.pid.size();
		CT_I32 np_data = pid.size();
		CT_I32 ind, i0;

		std::vector<CT_Merit> merit(pid0.maxgid+1, CT_Merit{0.});

		for(CT_I32 i=0; i<np_data; i++){
			ind 	= std::abs(pid[i].pid) % np_snap;

			if(ind < 0) ind = 0;

			i0 	= pid0.hash[ind];

			if(i0 < 0) continue;

			while(true){
				if(pid[i].pid == pid0.pid[i0].pid){
					merit[pid0.pid[i0].gid] += 1.;
					break;
				} else {
					if(pid0.hash_next[i0] >= 0){
						i0 	= pid0.hash_next[i0];
					} else{
						break;
					}
				}
			}
		}



#ifdef CTREE_USE_OMP
		#pragma omp parallel for default(none) \
			shared(np_snap, pid0, merit, np_data)
#endif
		for(CT_I32 i=0; i<pid0.maxgid+1; i++){
			if(pid0.n_ptcl2[i] <= 0) continue;
			merit[i] 	= merit[i] * merit[i] / pid0.n_ptcl2[i] / ((CT_Merit) np_data);
		}

		//for(CT_I32 i=0; i<np_snap; i++){
		//	if(pid0.n_ptcl[i]<=0) continue;
//
//		//	merit[pid0.pid[i].gid]	= merit[pid0.pid[i].gid] * merit[pid0.pid[i].gid] / pid0.n_ptcl[i];
		//}

		MeritSt meritcom;
		CT_Merit dum = -1;
		for(CT_I32 i=0; i<pid0.maxgid+1; i++){
			if(pid0.n_ptcl2[i]<=0) continue;

			if(merit[i] > dum){
				dum 	= merit[i];
				meritcom.id 	= i;
				meritcom.merit 	= dum;
			}
		}

		return meritcom;
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
			} else if(type == -1){
				is  = (data[i].stat >= 0);
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
	std::vector<CT_Merit> get_weight(const vctree_set::Settings& vh, std::vector<CT_PID> pid){
		CT_I32 npart = pid.size();
		std::vector<CT_Merit> weight(npart);
		
		CT_Merit dind;
		CT_Merit factor = (0.5772156649 + std::log( (CT_Merit) npart ));

		if(vh.ctree_weighttype == 1){

#ifdef CTREE_USE_OMP
			#pragma omp parallel for default(none) \
    			shared(factor, weight, npart) \
    			private(dind)
#endif
			for(CT_I32 i=0; i<(CT_I32) npart; i++){
				dind 	= npart - i;
				weight[i] 	= dind / ((CT_Merit) npart) / factor;
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
			if( ((CT_Merit) uncut) / ((CT_Merit) uind_n) > vh.ctree_minfrac) break;
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
		std::vector<CT_Merit> weight;
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

	PIDArray collectpidalongbranch2(const vctree_set::Settings& vh, CollectArray& CArr){

		PIDArray pid;
		std::vector<CT_Merit> weight;
		IO_dtype::GalArray gal0;

int myrank = mpi_rank();
if(myrank == 0 && CArr[0].snap>100) LOG()<<" here ";
		gal0 	= IO::r_gal(vh, CArr[0].snap, CArr[0].id, true);

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
			if(loop_ind >= (CT_I32) CArr[0].lind + 1) break;
			weight.clear();

int myrank = mpi_rank();
if(myrank == 0 && CArr[loop_ind].snap>100) LOG()<<" here ";
			gal0 	= IO::r_gal(vh, CArr[loop_ind].snap, CArr[loop_ind].id, true);

			ind1 	= ind0 + gal0[0].pid.size()-1;
			weight 	= get_weight(vh, gal0[0].pid);

			PIDReallocate(vh, pid, ind1);

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
			if(loop_n >= vh.ctree_n_step_n || loop_ind >= (CT_I32) CArr[0].lind + 1) break;
		}
		// remove particles appearing multiple times
		pid.resize(ind0);

		PIDArray pid2;
		pid2 	= get_coreptcl(vh, pid);


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
		std::vector<CT_Merit> weight;
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


				if(!Tree::istree(key, data[ind].snap, data[ind].id)){ 
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
        		if (rank == owner) blob_t = serialize_d1(data[ind]);
        		bcast_blob_from_owner(owner, blob_t);
        		if (rank != owner) deserialize_d1(blob_t, data[ind]);

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
			if(gal[i].id >= maxgid) maxgid = gal[i].id;
		}
		spt.maxgid 	= maxgid;


		std::vector<CT_I32> n_ptcl2(maxgid+1, CT_I32{0});


		for(CT_I32 i=0; i<(CT_I32) gal.size(); i++){
		//for(IO_dtype::GalSt g:gal){

			i1 	= i0 + gal[i].npart  - 1;
			for(CT_I32 j=0; j<(CT_I32) gal[i].npart; j++){
				spt.pid[j+i0].pid 	= gal[i].pid[j];
				spt.pid[j+i0].gid 	= gal[i].id;

				//if(gal[i].id >= maxgid) maxgid = gal[i].id;
				//if(g.pid[i]>=maxpid) maxpid = g.pid[i];
			}
			spt.n_ptcl[i]	= gal[i].npart;
			spt.gid[i] 		= gal[i].id;

			n_ptcl2[gal[i].id]	= gal[i].npart;

			i0 	= i1 + 1;
		}

		spt.n_ptcl2 	= std::move(n_ptcl2);
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
		std::vector<CT_Merit> match_merit(ncut, CT_Merit{-1.});



		get_merit(pid_g, gid_g, pid_s, gid_s, hash, hash_next, npart_g, npart_s, match_id, match_merit);


int myrank = mpi_rank();

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

if(cut[i]<2000 && myrank == 0){
	LOG()<<cut[i]<<" / "<<data[cut[i]].id0<<" / "<<data[cut[i]].snap0<<" / "<<match_merit[i]<<" / "<<match_id[i];
}
		}

	}

	void commerit2(const vctree_set::Settings& vh, ControlArray& data, Tree::TreeArray& tree, Tree::TreeKeyArray& key, SnapPT& pid0, CT_I32 snap_curr){
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
		std::vector<CT_snap> t_slist;
		std::vector<CT_ID> t_idlist;
		PIDArray cpid, cpid2;
		IO_dtype::GalArray gal0;
		Tree::TreeSt tree0;
		CT_I32 ntcut;
		CollectArray CArr(vh.ctree_n_search+1);

							

		//t_slist.reserve(vh.ctree_n_search);
		//t_idlist.reserve(vh.ctree_n_search);
		//t_slist.resize();
		//t_idlist.resize(vh.ctree_n_search);

		if(ncut > 0){

#ifdef CTREE_USE_MPI
			int rank = 0, size = 1;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	    	MPI_Comm_size(MPI_COMM_WORLD, &size);

			for(CT_I32 ind : cut){
				int owner = ind % size;
				if(owner != rank) continue;

				
				cpid.clear();


int myrank = mpi_rank();
if(myrank == 0 && data[ind].snap>100) LOG()<<" here : "<<ind;
				gal0 	= IO::r_gal(vh, data[ind].snap, data[ind].id, true);


				if(!Tree::istree(key, data[ind].snap, data[ind].id)){ 
					// no treee for this galaxy (read from a single snapshot)
					//weight 	= get_weight(vh, gal0[0].pid);
					cpid2.resize(gal0[0].pid.size());
					for(CT_I32 k=0; k<(CT_I32) cpid2.size(); k++){
						cpid2[k].pid 	= gal0[0].pid[k];
						//cpid2[k].weight 	= weight[k];

					}

					cpid 	= get_coreptcl(vh, cpid2);

					//pid0 = std::move(gal0[0].pid);

				} else {
					// this galaxy has a branch. Collect particles along its branch
					tree0 = Tree::gettree(tree, key, data[ind].snap, data[ind].id);


					ntcut = 0;
					CArr[0].lind = -1;
					//t_slist.clear();
					//t_idlist.clear();

					for(CT_I32 i=0; i<tree0.endind+1; i++){
						if(tree0.snap[i] >= data[ind].snap){
							CArr[ntcut].id 	= tree0.id[i];
							CArr[ntcut].snap= tree0.snap[i];
							CArr[0].lind 	= ntcut;
							ntcut ++;
							//t_slist.push_back(tree0.snap[i]);
							//t_idlist.push_back(tree0.id[i]);
							if(ntcut >= vh.ctree_n_search) break;
						}
					}


					if(ntcut == 0){
						LOG()<<"Weird branch: Snap = "<<data[ind].snap0<<" / ID = "<<data[ind].id0;
						u_stop();
					}
					//cpid 	= collectpidalongbranch(vh, t_slist, t_idlist);
					cpid 	= collectpidalongbranch2(vh, CArr);
				}


				// store to control array
				data[ind].p_list 	= std::move(cpid);
				data[ind].n_ptcl 	= data[ind].p_list.size();
				//data[ind].pos[0]	= gal0[0].xc;
				//data[ind].pos[1]	= gal0[0].yc;
				//data[ind].pos[2]	= gal0[0].zc;

				//data[ind].vel[0]	= gal0[0].vxc;
				//data[ind].vel[1]	= gal0[0].vyc;
				//data[ind].vel[2]	= gal0[0].vzc;

			}

			MPI_Barrier(MPI_COMM_WORLD);

			//Synchronize
			

    		for (CT_I32 ind : cut){
        		int owner = ind % size;

        		std::vector<std::uint8_t> blob_t;
        		if (rank == owner) blob_t = serialize_d1(data[ind]);
        		bcast_blob_from_owner(owner, blob_t);
        		if (rank != owner) deserialize_d1(blob_t, data[ind]);

	       	}
    		MPI_Barrier(MPI_COMM_WORLD);
#else 
			for(CT_I32 ind : cut){
				
				cpid.clear();

				gal0 	= IO::r_gal(vh, data[ind].snap, data[ind].id, true);


				if(!istree(key, data[ind].snap, data[ind].id)){ 
					// no treee for this galaxy (read from a single snapshot)
					//weight 	= get_weight(vh, gal0[0].pid);
					cpid2.resize(gal0[0].pid.size());
					for(CT_I32 k=0; k<(CT_I32) cpid2.size(); k++){
						cpid2[k].pid 	= gal0[0].pid[k];
						//cpid2[k].weight 	= weight[k];

					}

					cpid 	= get_coreptcl(vh, cpid2);

					//pid0 = std::move(gal0[0].pid);

				} else {
					// this galaxy has a branch. Collect particles along its branch
					tree0 = Tree::gettree(tree, key, data[ind].snap, data[ind].id);

					ntcut = 0;
					CArr[0].lind = -1;
					//t_slist.clear();
					//t_idlist.clear();


					for(CT_I32 i=0; i<vh.ctree_n_search; i++){
						if(tree0.snap[i] >= data[ind].snap){
							CArr[ntcut].id 	= tree0.id[i];
							CArr[ntcut].snap= tree0.snap[i];
							CArr[0].lind 	= ntcut;
							ntcut ++;

							//ntcut ++;
							//t_slist.push_back(tree0.snap[i]);
							//t_idlist.push_back(tree0.id[i]);
						}
					}


					if(ntcut == 0){
						LOG()<<"Weird branch: Snap = "<<data[ind].snap0<<" / ID = "<<data[ind].id0;
						u_stop();
					}

					//cpid 	= collectpidalongbranch(vh, t_slist, t_idlist);
					cpid 	= collectpidalongbranch2(vh, CArr);
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



		// Compute Merit

		std::vector<CT_I32> cut2 = get_control(data, 0);
		ncut = cut2.size();

		CT_I32 n0;

#ifdef CTREE_USE_MPI
		int rank = 0, size = 1;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	    MPI_Comm_size(MPI_COMM_WORLD, &size);
	    for(CT_I32 ind : cut2){
			int owner = ind % size;
			if(owner != rank) continue;

			MeritSt meritcom = get_merit3(pid0, data[ind].p_list);

			n0 	= data[ind].list_n;

			data[ind].list[n0].merit 	= meritcom.merit;
			data[ind].list[n0].id		= meritcom.id;
			data[ind].list[n0].snap 	= snap_curr;
			data[ind].list_n ++;
		}

		MPI_Barrier(MPI_COMM_WORLD);

		//Synchronize
    	for (CT_I32 ind : cut2){
        	int owner = ind % size;

        	std::vector<std::uint8_t> blob_t;
        	if (rank == owner) blob_t = serialize_d2(data[ind]);
        	bcast_blob_from_owner(owner, blob_t);
        	if (rank != owner) deserialize_d2(blob_t, data[ind]);

	    }
    	MPI_Barrier(MPI_COMM_WORLD);
#else
		for(CT_I32 ind : cut2){
			MeritSt meritcom = get_merit3(pid0, data[ind].p_list);

			n0 	= data[ind].list_n;

			data[ind].list[n0].merit 	= meritcom.merit;
			data[ind].list[n0].id		= meritcom.id;
			data[ind].list[n0].snap 	= snap_curr;
			data[ind].list_n ++;

		}
#endif

	}

	//-----
	// Link Branch
	//-----
	CT_Merit brcompare(const vctree_set::Settings& vh, CT_I32 s0, CT_I32 id0, std::vector<CT_I32>& slist, std::vector<CT_I32>& idlist){

		IO_dtype::GalArray gal0 	= IO::r_gal(vh, s0, id0, true);

		std::vector<CT_PID> pid0 = std::move(gal0[0].pid);
		std::vector<CT_Merit> pweight0 = get_weight(vh, pid0);

		CT_Merit n_occ;
		std::vector<CT_PID> pid1;
		std::vector<CT_Merit> pweight1;
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

			n_occ 	= (CT_Merit) cpid[0].n_con;
		}

		CT_I32 merittype = 1;	// pre-selected to be 1
		CT_Merit meritdum = get_merit2(pid0, pid1, pweight0, pweight1, merittype);

		//-----
		// Factor calculation
		//	- should be 1 if all particles appear during the selected part of the branch
		// 	- higher occurance should have a stronger weight

		CT_Merit factor = ((CT_Merit) 1.) / std::pow( ((CT_Merit) vh.ctree_n_step_n * vh.ctree_n_step_dn),2.) * std::pow( ((CT_Merit )n_occ * vh.ctree_n_step_dn),2.);

		return factor*meritdum;
		
	}

	CT_Merit brcompare2(const vctree_set::Settings& vh, CT_I32 s0, CT_I32 id0, CollectArray& CArr){

int myrank = mpi_rank();
if(myrank == 0 && s0>100) LOG()<<" here ";

		IO_dtype::GalArray gal0 	= IO::r_gal(vh, s0, id0, true);

		std::vector<CT_PID> pid0 = std::move(gal0[0].pid);
		std::vector<CT_Merit> pweight0 = get_weight(vh, pid0);

		CT_Merit n_occ;
		std::vector<CT_PID> pid1;
		std::vector<CT_Merit> pweight1;

		if( CArr[0].lind == 0 ){
if(myrank == 0 && CArr[0].snap>100) LOG()<<" here ";
			IO_dtype::GalArray gal1 	= IO::r_gal(vh, CArr[0].snap, CArr[0].id, true);
			pid1 = std::move(gal1[0].pid);
			pweight1 = get_weight(vh, pid1);
			n_occ = 1.;
		}else{
			PIDArray cpid = collectpidalongbranch2(vh, CArr);
			pid1.resize(cpid.size());
#ifdef CTREE_USE_OMP
			#pragma omp parallel for default(none) \
    			shared(pid1, cpid)
#endif
			for(CT_I32 i=0; i< (CT_I32) pid1.size(); i++){
				pid1[i] 	= cpid[i].pid;
			}
			
			pweight1 = get_weight(vh, pid1);

			n_occ 	= (CT_Merit) cpid[0].n_con;
		}

		CT_I32 merittype = 1;	// pre-selected to be 1
		CT_Merit meritdum = get_merit2(pid0, pid1, pweight0, pweight1, merittype);

		//-----
		// Factor calculation
		//	- should be 1 if all particles appear during the selected part of the branch
		// 	- higher occurance should have a stronger weight

		CT_Merit factor = ((CT_Merit) 1.) / std::pow( ((CT_Merit) vh.ctree_n_step_n * vh.ctree_n_step_dn),2.) * std::pow( ((CT_Merit )n_occ * vh.ctree_n_step_dn),2.);

		return factor*meritdum;
		
	}

	void linkbr(const vctree_set::Settings& vh, ControlArray& data, ControlKey& dkey, CT_I32 ind, IO::snapinfo& sinfo, Tree::TreeArray& tree, Tree::TreeKeyArray& key, CT_I32 id_to_link, CT_I32 snap_to_link, CT_Merit merit_to_link, CT_I32 snap_curr){

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

		CT_Merit merit_com = brcompare(vh, snap_to_link, id_to_link, brcom_snap, brcom_id);
		CT_Merit merit_org = brcompare(vh, snap_to_link, id_to_link, brorg_snap, brorg_id);

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

double howlong1, howlong2, howlong3, howlong4, howlong5;
int myrank = mpi_rank();
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

		CT_Merit dum_merit;
		
auto t0 = std::chrono::steady_clock::now();
	
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

if(myrank == 0){
	auto t1 = std::chrono::steady_clock::now();
	howlong1 	= std::chrono::duration<double>(t1 - t0).count();
}


t0 = std::chrono::steady_clock::now();
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

if(myrank == 0){
	auto t1 = std::chrono::steady_clock::now();
	howlong2 	= std::chrono::duration<double>(t1 - t0).count();
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

		CollectArray CArr(vh.snapf + 1);

		CT_I32 ntt;

		CT_Merit this_merit, other_merit;

t0 = std::chrono::steady_clock::now();

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

					//tt_snap.resize(0);
					//tt_id.resize(0);

					CArr[0].lind = -1;

					ntt = 0;
					for(CT_I32 k=0; k<tree0.endind+1; k++){
						if(tree0.snap[k] > checkcon[j].snap + vh.ctree_n_step_dn){
							//tt_snap.push_back(tree0.snap[k]);
							//tt_id.push_back(tree0.id[k]);
							//ntt ++;

							CArr[ntt].snap 	= tree0.snap[k];
							CArr[ntt].id 	= tree0.id[k];
							CArr[0].lind 	= ntt;
							ntt ++;
						}
					}

					if(ntt == 0){
						for(CT_I32 k=0; k<tree0.endind+1; k++){
							if(tree0.snap[k] > checkcon[j].snap){
								//tt_snap.push_back(tree0.snap[k]);
								//tt_id.push_back(tree0.id[k]);
								//ntt ++;

								CArr[ntt].snap 	= tree0.snap[k];
								CArr[ntt].id 	= tree0.id[k];
								CArr[0].lind 	= ntt;
								ntt ++;
							}
						}
					}
					

					//checkcon[j].merit 	= brcompare(vh, checkcon[j].snap, checkcon[j].id, tt_snap, tt_id);
					checkcon[j].merit 	= brcompare2(vh, checkcon[j].snap, checkcon[j].id, CArr);

					
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
if(myrank == 0){
	auto t1 = std::chrono::steady_clock::now();
	howlong3 	= std::chrono::duration<double>(t1 - t0).count();
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
		Tree::Tree_BID keyval;
t0 = std::chrono::steady_clock::now();
#ifdef CTREE_USE_OMP
		#pragma omp parallel for default(none) \
			private(tree0, keyval) \
			shared(ncut, islink, data, cut, next_point, snap_curr, vh, key, tree)
#endif
		for(CT_I32 i=0; i<ncut; i++){
			if(islink[i] < 0){

				data[cut[i]].stat = -1;

				keyval 	= Tree::getkey(key, data[cut[i]].snap0, data[cut[i]].id0);

				//data[cut[i]].snap0 + key[0].key * data[cut[i]].id0;

				Tree::TreeSt& tree0 = tree[keyval];



				if(islink[i] == -1){
					tree0.stat = -1;
				}else if(islink[i] == -2){
					tree0.stat 	= -2;

					//keyval = next_point[i].snap + key[0].key * next_point[i].id;
					keyval 	= Tree::getkey(key, next_point[i].snap, next_point[i].id);
					tree0.frag_bid 	= keyval;//key[keyval].ind;

				}


				next_point[i].id 	= -1;
				next_point[i].snap 	= -1;

				ctfree(vh, data, cut[i], -1, -1, snap_curr);

			}
		}

if(myrank == 0){
	auto t1 = std::chrono::steady_clock::now();
	howlong4 	= std::chrono::duration<double>(t1 - t0).count();
}


//456456
		// Link to a next checkpoint
		CT_I32 snap_to_link;
		CT_I32 id_to_link;
		CT_Merit merit_to_link;
		
		ind = 0;
t0 = std::chrono::steady_clock::now();
	

#ifdef CTREE_USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);

		LinkJob thisjob;
		;
		CT_I32 rank_index = rank;
		std::vector<CT_I32> keytracer(tree.size(), {-1});
		MPI_Datatype LINKJOB_T = make_linkjob_type();
		MPI_Datatype NEXT_T = make_next_type();

		JobArray job_ind, job_dind, job_tind, job_tind2;

		int jobdone = 0;
CT_I32 oldncut	= ncut;
CT_I32 oldnn 	= next_point.size();
		while(true){
			init_job(thisjob);
			

			if(rank_index < ncut){
				ind 	= cut[rank_index];

				if(!(data[ind].stat == -1 || islink[rank_index] < 0)){
					snap_to_link 	= next_point[rank_index].snap;
					id_to_link 		= next_point[rank_index].id;
					merit_to_link	= next_point[rank_index].merit;

					thisjob	= get_job(vh, tree, key, data, dkey, ind, snap_to_link, id_to_link, merit_to_link, snap_curr);

				}
			}

			// Remove a que for overlapped ones
			check_overlap(NEXT_T, thisjob, thisjob.ind, ind, cut, next_point, rank_index);
			check_overlap(NEXT_T, thisjob, thisjob.dind, ind, cut, next_point, rank_index);
			check_overlap(NEXT_T, thisjob, thisjob.tind, ind, cut, next_point, rank_index);
			check_overlap(NEXT_T, thisjob, thisjob.tind2, ind, cut, next_point, rank_index);
			ncut 	= cut.size(); // update ncut

			// Gather que
			int owner;
			owner = thisjob.ind % size;
			job_ind 	= std::move(commque(LINKJOB_T, thisjob, owner));

			owner = thisjob.dind % size;
			job_dind 	= std::move(commque(LINKJOB_T, thisjob, owner));

			owner = thisjob.tind % size;
			job_tind 	= std::move(commque(LINKJOB_T, thisjob, owner));

			owner = thisjob.tind2 % size;
			job_tind2 	= std::move(commque(LINKJOB_T, thisjob, owner));

			// Do Job1
			// 1a for data[ind]
			if(job_ind.size()>0){
				for(auto j:job_ind){
					if(j.jobnum != 1) continue;
					DoJob1a(vh, data, j, snap_curr);
				}
			}

			// 1b for tree[tind]
			if(job_tind.size()>0){
				for(auto j:job_tind){
					if(j.jobnum != 1) continue;
					DoJob1b(vh, tree, key, data, j);
				}
			}

			// Do Job2
			
			// 2a for data[ind]
			if(job_ind.size()>0){
				for(auto j:job_ind){
					if(j.jobnum != 2) continue;
					DoJob2a(vh, tree, data, j, snap_curr);

				}
			}


			// 2b for data[dind]
			if(job_dind.size()>0){
				for(auto j:job_dind){
					if(j.jobnum != 2) continue;
					DoJob2b(vh, data, j.dind, snap_curr);
				}
			}


			// 2c for tree[tind]
			if(job_tind.size()>0){
				for(auto j:job_tind){
					if(j.jobnum != 2) continue;
					DoJob2c(tree, key, data, j);
				}
			}

			// 2d for tree[tind2]
			if(job_tind2.size()>0){
				for(auto j:job_tind2){
					if(j.jobnum != 2) continue;
					DoJob2d(tree, j);
				}
			}
			

			// Do Job 3

			// 3a for data[ind]
			if(job_ind.size()>0){
				for(auto j:job_ind){
					if(j.jobnum != 3) continue;
					DoJob3a(vh, data, j, snap_curr);
				}
			}

			// 3b for tree[tind]
			if(job_tind.size()>0){
				for(auto j:job_tind){
					if(j.jobnum != 3) continue;
					DoJob3b(tree, j.tind, j.tind2);
				}
			}

			// Do Job 4

			// 4a for data[ind]
			if(job_ind.size()>0){
				for(auto j:job_ind){
					if(j.jobnum != 4) continue;
					DoJob4a(vh, tree, data, j, snap_curr);
				}
			}

			// 4b for tree[tind]
			if(job_tind.size()>0){
				for(auto j:job_tind){
					if(j.jobnum != 4) continue;
					DoJob4b(tree, key, data, j);
				}
			}

			// 4c for tree[tind2]
			if(job_tind2.size()>0){
				for(auto j:job_tind2){
					if(j.jobnum != 4) continue;
					DoJob4c(tree, key, j);
				}
			}


			// 4d for data[dind];
			// 4e for data[dind];
			// TODO-----
			// modtree should be changed in tind2 que
			if(job_dind.size()>0){
				for(auto j:job_dind){
					if(j.jobnum != 4) continue;
					if(j.dind<0) continue;

					Tree::TreeSt modtree = tree[j.tind2];
					CT_I32 lsnap = -1;
					CT_I32 lid = -1;
					for(CT_I32 k=modtree.endind; k>=0; k--){
						if(modtree.snap[k]<=j.snap) break;
						lsnap 	= modtree.snap[k];
					}

					if(lsnap<0){
						DoJob4d(vh, data, j.dind, snap_curr);
					}else{
						CT_I32 job4type[3];

						job4type[1]	= lsnap;
						job4type[2]	= lid;
						if(lsnap <= snap_curr){
							job4type[0] = 1;
						}else{
							CT_I32 index0, index1;

							index0	= wheresnap(sinfo, snap_curr);
							index1 	= wheresnap(sinfo, lsnap);
							if(std::abs(index0-index1) > vh.ctree_n_search){
								job4type[0] = 2;

								// this to be updated in the future
								//modtree.stat 	= -2;
								//modtree.frag_bid= j.tind;
							}else{
								job4type[0] = 3;
							}
						}

						DoJob4e(vh, data, j.dind, snap_curr, job4type);

					}
				}
			}

			// synchronize data & reset dkey
			syn_data(thisjob, data, dkey);
if(myrank == 0){
for(CT_I32 i=0; i<data[0].last_ind+1; i++){
	if(data[i].snap>100 || data[i].snap0>100) LOG()<<" snap !! "<<i<<" / "<<data[i].snap<<" / "<<data[i].snap0<<" / "<<rank_index;
	if(data[i].id < 0 && data[i].snap > 0) LOG()<<" id !! "<<i<<" / "<<data[i].id<<" / "<<data[i].snap0<<" / "<<rank_index;

}
}

			// synchronize tree & reset key
			syn_tree(thisjob, tree, key);

			rank_index 	+= size;
			
			if(rank_index >= ncut){
				jobdone = 1;
			}
	
        	int done_count = 0;
        	MPI_Allreduce(&jobdone, &done_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        	
        	if (done_count == size) break;
		}

		
///---------------Old ver...


//		JobArray linkjob(ncut);
//		
//		// Collect Job Type
//		for(CT_I32 i=0; i<ncut; i++){
//			
//			
//			int owner = i % size;
//			if(owner != rank) continue;
//
//			ind 	= cut[i];
//			
//			if(data[ind].stat == -1 || islink[i] < 0){
//				continue;		// already linked to another
//			}
//
//
//			snap_to_link 	= next_point[i].snap;
//			id_to_link 		= next_point[i].id;
//			merit_to_link	= next_point[i].merit;
//
//
//			linkjob[i]	= get_job(vh, tree, key, data, ind, snap_to_link, id_to_link, merit_to_link, snap_curr);
//
//
//		}
//
//		MPI_Barrier(MPI_COMM_WORLD);
//
//		// Merge Job list
//		for(CT_I32 i=0; i<ncut; i++){
//			if(islink[i] < 0) continue;
//
//			int owner = i % size;
//			
//			std::vector<std::uint8_t> blob_t;
//			if (rank == owner) blob_t = serialize(linkjob[i]);
//			bcast_blob_from_owner(owner, blob_t);
//	        if (rank != owner) deserialize(blob_t, linkjob[i]);
//		}
//
//		MPI_Barrier(MPI_COMM_WORLD);
//
//		// Do Job 1
//		for(CT_I32 i=0; i<ncut; i++){
//			if(islink[i] < 0) continue;
//			if( linkjob[i].jobnum != 1 ) continue;
//
//			int owner = linkjob[i].ind % size;
//			if(owner != rank) continue;
//
//			DoJob1(vh, tree, key, data, linkjob[i], snap_curr);
//		}
//
//		// Do Job 2
//		for(CT_I32 i=0; i<ncut; i++){
//			if(islink[i] < 0) continue;
//			if( linkjob[i].jobnum != 2 ) continue;
//
//			int owner = linkjob[i].ind % size;
//			CT_I32 dind = -1;
//			Tree::Tree_BID tind = -1;
//			Tree::Tree_BID tind2 = -1;
//
//			CT_I32 barray[3];
//			if(owner == rank){
//				DoJob2a(vh, data, linkjob[i], snap_curr);
//
//				dind = dkey[ linkjob[i].c_tree.snap[ linkjob[i].c_tree.endind ] + dkey[0] * linkjob[i].c_tree.id[ linkjob[i].c_tree.endind ] ];
//				tind 	= Tree::getkey(key, data[linkjob[i].ind].snap, data[linkjob[i].ind].id);
//				tind2 	= Tree::getkey(key, linkjob[i].snap, linkjob[i].id);
//				//barray[0]	= dkey[ linkjob[i].c_tree.snap[ linkjob[i].c_tree.endind ] + dkey[0] * linkjob[i].c_tree.id[ linkjob[i].c_tree.endind ] ];
//				//barray[1]	= Tree::getkey(key, data[linkjob[i].ind].snap, data[linkjob[i].ind].id);
//				//barray[2] 	= Tree::getkey(key, linkjob[i].snap, linkjob[i].id);
//			}
//
//			MPI_Bcast(&dind, 1, MPI_INT, owner, MPI_COMM_WORLD);
//			MPI_Bcast(&tind, 1, Tree::Tree_BID, owner, MPI_COMM_WORLD);
//			MPI_Bcast(&tind2, 1, Tree::Tree_BID, owner, MPI_COMM_WORLD);
//			//MPI_Bcast(&barray[0], 3, MPI_INT, owner, MPI_COMM_WORLD);
//			//dind 	= barray[0];
//			//tind 	= barray[1];
//			//tind2 	= barray[2];
//
//			if(dind % size == rank){
//				DoJob2b(vh, data, dind, snap_curr);
//			}
//		
//
//			if(tind % size == rank){
//				DoJob2c(tree, key, data, linkjob[i]);
//			}
//			
//			
//			if(tind2 % size == rank){
//				DoJob2d(tree, key, linkjob[i]);
//			}
//		}
//
//		// Do Job 3
//		for(CT_I32 i=0; i<ncut; i++){
//			if(islink[i] < 0) continue;
//			if( linkjob[i].jobnum != 3 ) continue;
//
//			int owner = linkjob[i].ind % size;
//			Tree::Tree_BID tind = -1;
//			Tree::Tree_BID tind2 = -1;
//
//			if(owner == rank){
//				DoJob3a(vh, data, linkjob[i], snap_curr);
//
//				tind 	= Tree::getkey(key, data[linkjob[i].ind].snap, data[linkjob[i].ind].id);
//				tind2 	= Tree::getkey(key, linkjob[i].snap, linkjob[i].id);
//			}
//
//			MPI_Bcast(&tind, 1, Tree::Tree_BID, owner, MPI_COMM_WORLD);
//			MPI_Bcast(&tind2, 1, Tree::Tree_BID, owner, MPI_COMM_WORLD);
//
//			if(tind%size == rank){
//				DoJob3b(tree, tind, tind2);
//			}
//		}
//
//		// Do Job 4
//		for(CT_I32 i=0; i<ncut; i++){
//			if(islink[i] < 0) continue;
//			if( linkjob[i].jobnum != 4 ) continue;
//
//			int owner = linkjob[i].ind % size;
//			Tree::Tree_BID tind = -1;
//			Tree::Tree_BID tind2 = -1;
//			CT_I32 dind = -1;
//
//			if(owner == rank){
//				DoJob4a(vh, data, linkjob[i], snap_curr);
//
//				dind = dkey[ linkjob[i].c_tree.snap[ linkjob[i].c_tree.endind ] + dkey[0] * linkjob[i].c_tree.id[ linkjob[i].c_tree.endind ] ];
//				tind 	= Tree::getkey(key, data[linkjob[i].ind].snap, data[linkjob[i].ind].id);
//				tind2 	= Tree::getkey(key, linkjob[i].snap, linkjob[i].id);
//			}
//
//			MPI_Bcast(&dind, 1, MPI_INT, owner, MPI_COMM_WORLD);
//			MPI_Bcast(&tind, 1, Tree::Tree_BID, owner, MPI_COMM_WORLD);
//			MPI_Bcast(&tind2, 1, Tree::Tree_BID, owner, MPI_COMM_WORLD);
//
//			if(tind % size == rank){
//				DoJob4b(tree, key, data, linkjob[i]);
//			}
//
//			bool gostop;
//			if(tind2 % size == rank){
//				DoJob4c(tree, key, data, dind, linkjob[i], gostop);
//			}
//			
//
//			if(dind<0) continue;
//
//
//			MPI_Bcast(&gostop, 1, MPI_C_BOOL, owner, MPI_COMM_WORLD);
//
//			if(!gostop){
//				if(dind % size == rank){
//					DoJob4d(vh, data, dind, snap_curr);
//				}
//			}else{
//				CT_I32 job4type[3];
//				if(tind2 % size == rank){
//					Tree::TreeSt& modtree = tree[tind2];					
//
//					job4type[1]	= modtree.snap[0];
//					job4type[2]	= modtree.id[0];
//					if(modtree.snap[0] <= snap_curr){
//						job4type[0] = 1;
//					}else{
//						CT_I32 index0, index1;
//
//						index0	= wheresnap(sinfo, snap_curr);
//						index1 	= wheresnap(sinfo, modtree.snap[0]);
//						if(std::abs(index0-index1) > vh.ctree_n_search){
//							job4type[0] = 2;
//
//							modtree.stat 	= -2;
//							modtree.frag_bid= tind;
//						}else{
//							job4type[0] = 3;
//						}
//
//					}
//				}
//
//
//				MPI_Bcast(&job4type[0], 3, MPI_INT, owner, MPI_COMM_WORLD);
//
//				if(dind % size == rank){
//					DoJob4e(vh, data, dind, snap_curr, job4type);
//				}
//			}
//
//		}
//
//		// merge data, tree
//
//
//		for(CT_I32 i=0; i<data[0].last_ind+1; i++){
//			int owner = i % size;
//
//			std::vector<std::uint8_t> blob_t;
//			if (rank == owner) blob_t = serialize(data[i]);
//			bcast_blob_from_owner(owner, blob_t);
//	        if (rank != owner) deserialize(blob_t, data[i]);
//		}
//
//		MPI_Barrier(MPI_COMM_WORLD);
//
//for(CT_I32 i=0; i<24; i++){
//	if(myrank == i){
//		LOG()<<myrank<<" / "<<tree[1].endind;
//	}
//	MPI_Barrier(MPI_COMM_WORLD);
//}
//
//
//		for(CT_I32 i=0; i<tree[0].lind+1; i++){
//			int owner = i % size;
//
//			std::vector<std::uint8_t> blob_t;
//			if (rank == owner) blob_t = serialize(tree[i]);
//			bcast_blob_from_owner(owner, blob_t);
//	        if (rank != owner) deserialize(blob_t, tree[i]);
//		}
//
//		MPI_Barrier(MPI_COMM_WORLD);
//
//for(CT_I32 i=0; i<24; i++){
//	if(myrank == i){
//		LOG()<<myrank<<" / "<<tree[1].endind;
//	}
//	MPI_Barrier(MPI_COMM_WORLD);
//}
//
//		// reset key
////ifdef CTREE_USE_OMP
////		pragma omp parallel for default(none) 
////			shared(tree, key)
////endif
//		for(CT_I32 i=0; i<tree[0].lind+1; i++){
//
//			for(CT_I32 j=0; j<tree[i].endind+1; j++){
//				key[ tree[i].snap[j] + key[0].key * tree[i].id[j] ].ind = i;
//			}
//		}
//


#else

		// parallelizatin here..
		for(CT_I32 i=0; i<ncut; i++){
			if(islink[i] < 0) continue;

			ind 	= cut[i];

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
#endif



	if(myrank == 0){
		auto t1 = std::chrono::steady_clock::now();
		howlong5 	= std::chrono::duration<double>(t1 - t0).count();
		LOG()<<" Gather Next point      : "<<howlong1;
		LOG()<<" Gather Check Arr       : "<<howlong2;
		LOG()<<" Check Link             : "<<howlong3;
		LOG()<<" Close Link             : "<<howlong4;
		LOG()<<" Do Link                : "<<howlong5;
		
		//savedata(vh, data, snap_curr);
		//savetree_ctree(vh, tree, key, snap_curr);
		
		LOG()<<oldncut<<" / "<<cut.size();
		LOG()<<oldnn<<" / "<<next_point.size();
		ControlArray data2 = loaddata(vh, snap_curr);
		Tree::TreeArray tree2;
		Tree::TreeKeyArray key2;
		loadtree_ctree(vh, tree2, key2, snap_curr);

		validate_data(data, data2);
		validate_tree(tree, tree2);
		//validate_treekey(key, key2);

			

		if(snap_curr == 90) u_stop();
	}

	}

	// Link MPI
//	void DoJob1(const vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& key, ControlArray& data, LinkJob& job, CT_I32 snap_curr){
//
//		// tree[ tind ]
//		expandbr(vh, data, job.ind, tree, key, job.id, job.snap, job.merit);
//		
//		// data [ ind ]
//		ctfree(vh, data, job.ind, job.snap, job.id, snap_curr);
//	}
#ifdef CTREE_USE_MPI
	void init_job(LinkJob& thisjob){
		thisjob.jobnum 	= -1;
		thisjob.ind 	= -1;
		thisjob.id 		= -1;
		thisjob.snap 	= -1;
		thisjob.merit 	= -1.;
		thisjob.dind 	= -1;
		thisjob.tind 	= -1;
		thisjob.tind2 	= -1;
	}
	void filter_sort(std::vector<std::pair<CT_I32, CT_I32>>& pairs){
    	// remove < 0
    	pairs.erase(std::remove_if(pairs.begin(), pairs.end(), [](auto& p){ return p.first < 0; }), pairs.end());

    	// sort
    	std::stable_sort(pairs.begin(), pairs.end(), [](auto& a, auto& b){return a.first < b.first;});

	}

	void syn_data(LinkJob& thisjob, ControlArray& data, ControlKey& dkey){

		int rank = 0, size = 1;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	    MPI_Comm_size(MPI_COMM_WORLD, &size);


	    std::vector<CT_I32> all_ind(size), all_dind(size), all_num(size);

	    MPI_Allgather(&thisjob.ind,  1, MPI_INT, all_ind.data(),  1, MPI_INT, MPI_COMM_WORLD);
    	MPI_Allgather(&thisjob.dind, 1, MPI_INT, all_dind.data(), 1, MPI_INT, MPI_COMM_WORLD);
    	MPI_Allgather(&thisjob.jobnum,  1, MPI_INT, all_num.data(),  1, MPI_INT, MPI_COMM_WORLD);

    	std::vector<std::pair<CT_I32, CT_I32>> pairs_ind;  pairs_ind.reserve(size);
    	std::vector<std::pair<CT_I32, CT_I32>> pairs_dind; pairs_dind.reserve(size);
    	for (int r=0; r<size; ++r) {
        	pairs_ind.emplace_back(all_ind[r],  all_num[r]);
        	pairs_dind.emplace_back(all_dind[r], all_num[r]);
    	}

    	filter_sort(pairs_ind);
    	filter_sort(pairs_dind);

    	// merge for ind
    	for(CT_I32 i=0; i<(CT_I32)pairs_ind.size(); i++){
    		if(pairs_ind[i].second >0){

    			int owner = pairs_ind[i].first % size;

    			std::vector<std::uint8_t> blob_t;
    			if(rank == owner) blob_t = serialize(data[pairs_ind[i].first]);
    			bcast_blob_from_owner(owner, blob_t);
    			if(rank != owner) deserialize(blob_t, data[pairs_ind[i].first]);

    			data[pairs_ind[i].first].Free_plist();
				data[pairs_ind[i].first].n_ptcl 		= -1;

				if(data[pairs_ind[i].first].stat >= 0){
					dkey[data[pairs_ind[i].first].snap0 + dkey[0]*data[pairs_ind[i].first].id0] = pairs_ind[i].first;
				}else{
					dkey[data[pairs_ind[i].first].snap0 + dkey[0]*data[pairs_ind[i].first].id0] = -1;
				}
    		}
    	}

    	// merge for dind
    	for(CT_I32 i=0; i<(CT_I32)pairs_dind.size(); i++){
    		if(pairs_dind[i].second>0){

    			int owner = pairs_dind[i].first % size;

    			std::vector<std::uint8_t> blob_t;
    			if(rank == owner) blob_t = serialize(data[pairs_dind[i].first]);
    			bcast_blob_from_owner(owner, blob_t);
    			if(rank != owner) deserialize(blob_t, data[pairs_dind[i].first]);

    			if(pairs_dind[i].second == 2 || pairs_dind[i].second ==4){
    				data[pairs_dind[i].first].Free_plist();
					data[pairs_dind[i].first].n_ptcl 		= -1;
    			}

    			if(data[pairs_dind[i].first].stat >= 0){
					dkey[data[pairs_dind[i].first].snap0 + dkey[0]*data[pairs_dind[i].first].id0] = pairs_dind[i].first;
				}else{
					dkey[data[pairs_dind[i].first].snap0 + dkey[0]*data[pairs_dind[i].first].id0] = -1;
				}
    		}
    	}
    }

    void check_overlap(MPI_Datatype& NEXT_T, LinkJob& thisjob, CT_I32 jobind, CT_I32 ind, std::vector<CT_I32>& cut, NextArray& next_point, CT_I32 rank_index){
    	int rank, size;
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    	MPI_Comm_size(MPI_COMM_WORLD, &size);

    	

    	// collect jobind
    	std::vector<CT_I32> all_ind(size);
    	MPI_Allgather(&jobind, 1, MPI_INT, all_ind.data(), 1, MPI_INT, MPI_COMM_WORLD);

    	// sort
    	std::vector<std::pair<CT_I32,int>> pairs;
    	pairs.reserve(size);
    	for (int r = 0; r < size; ++r) {
    		//if(all_ind[r]<0) continue;
        	pairs.emplace_back(all_ind[r], r);
    	}

    	pairs.erase(std::remove_if(pairs.begin(), pairs.end(),
                           [](const auto& p){ return p.first < 0; }),
            pairs.end());


    	std::sort(pairs.begin(), pairs.end(),
              [](const auto& a, const auto& b){
                  if (a.first != b.first) return a.first < b.first;
                  return a.second < b.second;
              });

    	// winner or loser
    	std::vector<char> is_loser(size, 0);
    	//bool any_dup = false;

    	for (int i = 0; i < (int) pairs.size(); ) {
        	int j = i + 1;
        	while (j < (int) pairs.size() && pairs[j].first == pairs[i].first) ++j;

        	// [i, j) : same ind group
        	if (j - i > 1) {
            	//any_dup = true;
            	// i is winner, i+1..j-1 loser
            	for (int k = i + 1; k < j; ++k) {
                	int loser_rank = pairs[k].second;
                	is_loser[loser_rank] = 1;
            	}
        	}
        	i = j;
    	}

    	// send ind of loser
    	int sendcnt = is_loser[rank] ? 1 : 0;
    	std::vector<int> recvcnts(size);
    	MPI_Allgather(&sendcnt, 1, MPI_INT, recvcnts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    	std::vector<int> displs(size, 0);
    	int total = 0;
    	for (int i = 0; i < size; ++i) {
        	displs[i] = total;
        	total    += recvcnts[i];
    	}

    	std::vector<CT_I32> gathered_loser_inds(total);
    	MPI_Allgatherv(is_loser[rank] ? &ind : nullptr, sendcnt, MPI_INT,
                   gathered_loser_inds.data(), recvcnts.data(), displs.data(), MPI_INT,
                   MPI_COMM_WORLD);

    	// push back to cut
    	cut.insert(cut.end(), gathered_loser_inds.begin(), gathered_loser_inds.end());

    	// append next_array
    	NextArray recv_next(total);
    	NextSt local_next_elem	= next_point[rank_index];
    	if (total > 0) {
            MPI_Allgatherv(is_loser[rank] ? &local_next_elem : nullptr, sendcnt, NEXT_T,
                recv_next.data(), recvcnts.data(), displs.data(), NEXT_T, MPI_COMM_WORLD);
        	next_point.insert(next_point.end(), recv_next.begin(), recv_next.end());
    	}

    	// change job for loser
    	if(is_loser[rank] == 1){
    		init_job(thisjob);
    	}

    }

    

//		CT_I32 local_ind[2] = {thisjob.ind, thisjob.dind};
//		
//			
//		std::vector<CT_I32> all_inds(size*2);
//		MPI_Allgather(&local_ind, 2, MPI_INT, all_inds.data(), 2, MPI_INT, MPI_COMM_WORLD);
//
//		std::vector<CT_I32> uniq_sorted;
//		uniq_sorted.reserve(all_inds.size());
//
//    		
//   		for (CT_I32 x : all_inds) {
//       		if (x >= 0) uniq_sorted.push_back(x);
//   		}
//
//    		
//   		std::sort(uniq_sorted.begin(), uniq_sorted.end());
//   		uniq_sorted.erase(std::unique(uniq_sorted.begin(), uniq_sorted.end()), uniq_sorted.end());

//   		for(CT_I32 i=0; i<(CT_I32) uniq_sorted.size(); i++){
//
//			int owner = uniq_sorted[i] % size;
//
//			std::vector<std::uint8_t> blob_t;
//			if(rank == owner) blob_t = serialize(data[uniq_sorted[i]]);
//			bcast_blob_from_owner(owner, blob_t);
//        	if(rank != owner) deserialize(blob_t, data[uniq_sorted[i]]);
//		}

		// reset dkey
//		for(CT_I32 i=0; i<(CT_I32) uniq_sorted.size(); i++){
//			if(data[uniq_sorted[i]].stat >= 0){
//				dkey[data[uniq_sorted[i]].snap0 + dkey[0]*data[uniq_sorted[i]].id0] = uniq_sorted[i];
//			}else{
//				dkey[data[uniq_sorted[i]].snap0 + dkey[0]*data[uniq_sorted[i]].id0] = -1;
//			}
//		}
	

	void syn_tree(LinkJob& thisjob, Tree::TreeArray& tree, Tree::TreeKeyArray& key){

		int rank = 0, size = 1;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	    MPI_Comm_size(MPI_COMM_WORLD, &size);

		Tree::Tree_BID local_ind[2] = {thisjob.tind, thisjob.tind2};
			
		std::vector<Tree::Tree_BID> all_inds(size*2);

		MPI_Datatype CT_BID_T;
    	if(sizeof(Tree::Tree_BID)==8){
    		CT_BID_T = MPI_INT64_T;
    	}else if(sizeof(Tree::Tree_BID)==4){
    		CT_BID_T = MPI_INT32_T;
    	}
		
		MPI_Allgather(&local_ind, 2, CT_BID_T, all_inds.data(), 2, CT_BID_T, MPI_COMM_WORLD);

		std::vector<CT_I32> uniq_sorted;
   		uniq_sorted.reserve(all_inds.size());

    		
   		for (CT_I32 x : all_inds) {
       		if (x >= 0) uniq_sorted.push_back(x);
   		}

    		
   		std::sort(uniq_sorted.begin(), uniq_sorted.end());
   		uniq_sorted.erase(std::unique(uniq_sorted.begin(), uniq_sorted.end()), uniq_sorted.end());


   		for(CT_I32 i=0; i<(CT_I32) uniq_sorted.size(); i++){

			int owner = uniq_sorted[i] % size;

			std::vector<std::uint8_t> blob_t;
			if (rank == owner) blob_t = serialize(tree[uniq_sorted[i]]);
			bcast_blob_from_owner(owner, blob_t);
        	if (rank != owner) deserialize(blob_t, tree[uniq_sorted[i]]);

        	if(tree[uniq_sorted[i]].isfree == 1){
				tree[uniq_sorted[i]].id.resize(0);
				tree[uniq_sorted[i]].snap.resize(0);
				tree[uniq_sorted[i]].p_merit.resize(0);
				tree[uniq_sorted[i]].m_id.resize(0);
				tree[uniq_sorted[i]].m_snap.resize(0);
				tree[uniq_sorted[i]].m_merit.resize(0);
				tree[uniq_sorted[i]].m_bid.resize(0);
				tree[uniq_sorted[i]].father_bid = -1;
				tree[uniq_sorted[i]].numprog = 0;
				tree[uniq_sorted[i]].endind = -1;
				tree[uniq_sorted[i]].isfree 	= -1;
			}
        	
		}

		// reset key
		Tree::TreeSt dumtree;
		for(CT_I32 i=0; i<(CT_I32) uniq_sorted.size(); i++){
			dumtree 	= tree[uniq_sorted[i]];

#ifdef CTREE_USE_OMP
			#pragma omp parallel for default(none) shared(dumtree, key, uniq_sorted, i)
#endif
			for(CT_I32 j=0; j<dumtree.endind+1; j++){
				//key[ dumtree.snap[j] + key[0].key * dumtree.id[j] ].ind = uniq_sorted[i];
				key[ dumtree.snap[j] + key[0] * dumtree.id[j] ] = uniq_sorted[i];
			}
		}


	}

	void DoJob1a(const vctree_set::Settings& vh, ControlArray& data, LinkJob& job, CT_I32 snap_curr){
		// data [ ind ]
		ctfree(vh, data, job.ind, job.snap, job.id, snap_curr);
	}

	void DoJob1b(const vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& key, ControlArray& data, LinkJob& job){		// tree[ tind ]
		expandbr(vh, data, job.ind, tree, key, job.id, job.snap, job.merit);
		// tree[ tind ]
	}

	void DoJob2a(const vctree_set::Settings& vh, Tree::TreeArray& tree, ControlArray& data, LinkJob& job, CT_I32 snap_curr){
		Tree::TreeSt c_tree = tree[job.tind2];
		// data[ ind ]
		ctfree(vh, data, job.ind, c_tree.snap[0], c_tree.id[0], snap_curr); // data
	}

	void DoJob2b(const vctree_set::Settings& vh, ControlArray& data, CT_I32 dind, CT_I32 snap_curr){
		// data[ dind ]
		if(dind>=0){
			data[dind].stat 	= -1;
			ctfree(vh, data, dind, -1, -1, snap_curr);
		}
	}

	void DoJob2c(Tree::TreeArray& tree, Tree::TreeKeyArray& key, ControlArray& data, LinkJob& job){
		Tree::TreeSt c_tree = tree[job.tind2];
		//c_tree.p_merit[0]	= job.merit;
		tree[job.tind].p_merit[0] 	= job.merit;

		// tree [ tind ]
		for(CT_I32 i=c_tree.endind; i>=0; i--){
			Tree::treeinput(tree, key, data[job.ind].snap0, data[job.ind].id0, c_tree.snap[i], c_tree.id[i], c_tree.p_merit[i]);
		}
	}

	void DoJob2d(Tree::TreeArray& tree, LinkJob& job){
		// tree [ tind2 ]
		//Tree::treefree(tree, key, job.snap, job.id);
		//Tree::Tree_BID comp_key 	= Tree::getkey(key, job.snap, job.id);
		//Tree::TreeSt& comp_tree 	= tree[comp_key];
		Tree::TreeSt& comp_tree 	= tree[job.tind2];
		comp_tree.isfree 	= 1;

	}

	void DoJob3a(const vctree_set::Settings& vh, ControlArray& data, LinkJob& job, CT_I32 snap_curr){
		// data [ ind ]
		data[job.ind].stat = -1;
		ctfree(vh, data, job.ind, -1, -1, snap_curr);
	}

	void DoJob3b(Tree::TreeArray& tree, Tree::Tree_BID org_bid, Tree::Tree_BID com_bid){
		Tree::TreeSt& tree0 = tree[org_bid];
		// tree [ tind ]
		tree0.stat = -2;
		tree0.frag_bid = com_bid;
	}

	void DoJob4a(const vctree_set::Settings& vh, Tree::TreeArray& tree, ControlArray& data, LinkJob& job, CT_I32 snap_curr){
		Tree::TreeSt c_tree = tree[job.tind2];
		// data[ind]
		ctfree(vh, data, job.ind, c_tree.snap[0], c_tree.id[0], snap_curr);	
	}

	void DoJob4b(Tree::TreeArray& tree, Tree::TreeKeyArray& key, ControlArray& data, LinkJob& job){

		Tree::TreeSt c_tree = tree[job.tind2];
		// tree[ind]
		for(CT_I32 i=c_tree.endind; i>=0; i--){
			if(c_tree.snap[i] > job.snap) continue;
			Tree::treeinput(tree, key, data[job.ind].snap0, data[job.ind].id0, c_tree.snap[i], c_tree.id[i], c_tree.p_merit[i]);
		}
	}

	void DoJob4c(Tree::TreeArray& tree, Tree::TreeKeyArray& key, LinkJob& job){

		//tree[tind2]
		Tree::modifytree(tree, key, job.snap, job.id, job.snap);

		// if dind branch had a single checkpoint, so it is removed after modifytree, gostop = false

		//if(dind>=0){
		//	gostop 	= Tree::istree(key, data[dind].snap0, data[dind].id0);
		//}else{
		//	gostop = false;
		//}
		
	}

	void DoJob4d(const vctree_set::Settings& vh, ControlArray& data, CT_I32 dind, CT_I32 snap_curr){
		data[dind].stat = -1;
		ctfree(vh, data, dind, -1, -1, snap_curr);
	}

	void DoJob4e(const vctree_set::Settings& vh, ControlArray& data, CT_I32 dind, CT_I32 snap_curr, CT_I32* jobtype){
		if(jobtype[0] == 1){
			data[dind].stat = 1;
			ctfree(vh, data, dind, jobtype[1], jobtype[2], snap_curr);
		}else if(jobtype[0] == 2){
			data[dind].stat = -1;
			ctfree(vh, data, dind, -1, -1, snap_curr);
		}else if(jobtype[0] == 3){
			data[dind].stat = 0;
			ctfree(vh, data, dind, jobtype[1], jobtype[2], snap_curr);
		}
	}


	CT_Merit link_commerit(const vctree_set::Settings& vh, Tree::TreeSt tree0, CT_I32 snap_to_link, CT_I32 id_to_link, CT_I32 snap_curr){

		std::vector<CT_I32> brorg_id(tree0.endind+1), brorg_snap(tree0.endind+1);
		CT_I32 org_n = 0;
		

		//// Gather snap & ID
		for(CT_I32 i=0; i<tree0.endind+1; i++){
			if(tree0.snap[i] > snap_curr + vh.ctree_n_step_dn){
				brorg_id[org_n] 	= tree0.id[i];
				brorg_snap[org_n]  	= tree0.snap[i];
				org_n ++;
			}
		}
		if(org_n == 0){
			for(CT_I32 i=0; i<tree0.endind+1; i++){
				if(tree0.snap[i] > snap_curr){
					brorg_id[org_n] 	= tree0.id[i];
					brorg_snap[org_n]  	= tree0.snap[i];
					org_n ++;
				}
			}
		}

		brorg_id.resize(org_n);
		brorg_snap.resize(org_n);

		return brcompare(vh, snap_to_link, id_to_link, brorg_snap, brorg_id);
	}

	LinkJob get_job(const vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& key, ControlArray& data, ControlKey& dkey, CT_I32 ind, CT_I32 snap_to_link, CT_I32 id_to_link, CT_Merit merit_to_link, CT_I32 snap_curr){
		
		LinkJob thisjob;

		thisjob.ind 	= ind;
		thisjob.snap 	= snap_to_link;
		thisjob.id 		= id_to_link;
		thisjob.merit 	= merit_to_link;
		
		thisjob.tind 	= Tree::getkey(key, data[ind].snap, data[ind].id);
		thisjob.tind2 	= Tree::getkey(key, snap_to_link, id_to_link);

		//Job 1: Does connection point has a tree?
		bool istree = Tree::istree(key, snap_to_link, id_to_link);

		if(istree){
			thisjob.dind 	= dkey[tree[thisjob.tind2].snap[tree[thisjob.tind2].endind] + dkey[0] * tree[thisjob.tind2].id[tree[thisjob.tind2].endind]];
		}

		if(!istree){
			thisjob.jobnum 	= 1;
			return thisjob;
		}

		//Job 2: Connection point has a tree. But it can be fully merged to this branch
		Tree::TreeSt this_tree 		= Tree::gettree(tree, key, data[ind].snap, data[ind].id);
		Tree::TreeSt comp_tree		= Tree::gettree(tree, key, snap_to_link, id_to_link);

		// Debug check
		if(thisjob.dind>=0 && (data[thisjob.dind].id0 != tree[thisjob.tind2].id[tree[thisjob.tind2].endind] || data[thisjob.dind].snap0 != tree[thisjob.tind2].snap[tree[thisjob.tind2].endind] )){
			int myrank	= mpi_rank();
			LOG()<<"?? : "<<myrank<<" / "<<thisjob.dind<<" / "<<thisjob.tind<<" / "<<thisjob.tind2;
			LOG()<<"   : "<<data[thisjob.dind].snap0<<" / "<<data[thisjob.dind].id0;
			LOG()<<"   : "<<data[thisjob.dind].snap<<" / "<<data[thisjob.dind].id;
			LOG()<<"   : "<<tree[thisjob.tind2].snap[tree[thisjob.tind2].endind]<<" / "<<tree[thisjob.tind2].id[tree[thisjob.tind2].endind];
			LOG()<<" to: "<<snap_to_link<<" / "<<id_to_link;

			u_stop();
		}

		if(comp_tree.snap[comp_tree.endind] < this_tree.snap[0]){
			thisjob.jobnum = 2;
			return thisjob;
		}

		//Job: Merit comparison
		CT_Merit merit_com = link_commerit(vh, comp_tree, snap_to_link, id_to_link, snap_curr);
		CT_Merit merit_org = link_commerit(vh, this_tree, snap_to_link, id_to_link, snap_curr);
	
		//Job 3: Tree of the connection point is better
		if(merit_com > merit_org){
			thisjob.jobnum = 3;
			return thisjob;
		}


		//Job 4: Save both tree: This tree occupy a part comp_tree
		thisjob.jobnum = 4;
		return thisjob;
	}

	JobArray commque(MPI_Datatype& LINKJOB_T, LinkJob& thisjob, int owner){


		int rank = 0, size = 1;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	    MPI_Comm_size(MPI_COMM_WORLD, &size);

	    const int TAG_DATA = 100;
    	const int TAG_DONE = 101;
		
		JobArray job;
		job.resize(0);
		
		std::vector<int> bucket(size, 0);
		if (owner >= 0 && owner < size) bucket[owner] = 1;

		int my_bucket = 0;
		MPI_Reduce_scatter_block(bucket.data(), &my_bucket, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		int expected_done = my_bucket - ((owner == rank) ? 1 : 0);


		if(owner >= 0 && owner < size && owner != rank){ // send	
			MPI_Send(&thisjob, 1, LINKJOB_T, owner, TAG_DATA, MPI_COMM_WORLD);
			MPI_Send(nullptr, 0, MPI_BYTE, owner, TAG_DONE, MPI_COMM_WORLD);

		}


		if (expected_done > 0) {
			int done = 0;
			while(done < expected_done){
				MPI_Status st;
        		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &st);

            	if (st.MPI_TAG == TAG_DATA) {
		           	LinkJob getjob;
                	MPI_Recv(&getjob, 1, LINKJOB_T, st.MPI_SOURCE, TAG_DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                	job.push_back(getjob);
            	} else if (st.MPI_TAG == TAG_DONE) {
                	MPI_Recv(nullptr, 0, MPI_BYTE, st.MPI_SOURCE, TAG_DONE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                	++done;
            	}

			}

		}

		if(owner == rank){
			job.push_back(thisjob);
		}

		return job;
	}
#endif	
					

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

	//
	void delgal(const vctree_set::Settings& vh, ControlArray& data, ControlKey& dkey){

		std::vector<CT_I32> cut = get_control(data, -1);
		CT_I32 ncut = cut.size();
		if(ncut == data[0].last_ind+1) return;

		CT_I32 nn 	= ncut + vh.ctree_nstep;
		ControlArray data2 = allocate(vh, nn);

#ifdef CTREE_USE_OMP
		#pragma omp parallel for default(none) \
			shared(ncut, dkey, data, data2, cut)
#endif
		for(CT_I32 i=0; i<ncut; i++){
			data2[i]	= data[cut[i]];
			dkey[  data[cut[i]].snap0 + dkey[0]*data[cut[i]].id0  ] = i;
		}
		data2[0].last_ind = ncut-1;

		data 	= std::move(data2);

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
		double dt_classify, dt_readsnap, dt_commerit, dt_link, dt_addgal, dt_delgal;//, dt_collectpid, 

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

		//-----
		// For Null tree input
		//-----
		if(tree.size()==0 || key.size()==0){
			tree.resize(1);
			//key.resize(1);


			tree[0].lind 	= 0;

			// find maximum snapshot
			CT_I32 max_snap = -1;
			for(CT_I32 i=0; i<(CT_I32) sinfo.size(); i++){
				if(sinfo[i].snum > max_snap) max_snap = sinfo[i].snum;
			}

			// find maximum ID
			CT_I32 max_id = -1;
			for(CT_I32 i=0; i<(CT_I32) gal.size(); i++){
				if(gal[i].id > max_id) max_id = gal[i].id;
			}
			max_id 	+= 10000; //buffer

			CT_I32 max_key = max_snap+1;

			key.resize(max_snap + max_key*max_id, {-1});
			//key[0].key 	= max_key;
			key[0]			= max_key;
			//vh.treekey 	= max_snap+1;
		}

		//-----
		// Data Key
		//-----
		ControlKey dkey(key.size(), -1);
		//dkey[0] 	= key[0].key;
		dkey[0] 	= key[0];
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
				//LOG()<<"		"<<tree[0].lind<<" / "<<tree.size();

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
			//t0 = std::chrono::steady_clock::now();
			//PIDArray pid 	= collectpid(vh, data, tree, key);
			//if(myrank == 0){
			//	auto t1 = std::chrono::steady_clock::now();
			//	dt_collectpid = std::chrono::duration<double>(t1 - t0).count();
			//}

			//----- Read particles at this snapshot
			t0 = std::chrono::steady_clock::now();
			SnapPT pid0 	= readsnap(vh, data, sinfo[i].snum);
			if(myrank == 0){
				auto t1 = std::chrono::steady_clock::now();
				dt_readsnap = std::chrono::duration<double>(t1 - t0).count();
			}

			//----- Compute Merit
			t0 = std::chrono::steady_clock::now();
			commerit2(vh, data, tree, key, pid0, sinfo[i].snum);
			if(myrank == 0){
				auto t1 = std::chrono::steady_clock::now();
				dt_commerit = std::chrono::duration<double>(t1 - t0).count();
			}

			//----- Compute Merit
			//t0 = std::chrono::steady_clock::now();
			//commerit(vh, data, pid, pid0, sinfo[i].snum);
			//if(myrank == 0){
			//	auto t1 = std::chrono::steady_clock::now();
			//	dt_commerit = std::chrono::duration<double>(t1 - t0).count();
			//}

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

			//----- Delete Galaxies with a broken tree
			t0 = std::chrono::steady_clock::now();
			delgal(vh, data, dkey);
			if(myrank == 0){
				auto t1 = std::chrono::steady_clock::now();
				dt_delgal = std::chrono::duration<double>(t1 - t0).count();
			}


			//----- Time Log
			if(myrank == 0){
				LOG()<<"        Time report";
				LOG()<<"          Classify Galaxies in   "<<dt_classify;
				//LOG()<<"          Collect PID in         "<<dt_collectpid;
				LOG()<<"          Read Snapshot ptcls in "<<dt_readsnap;
				LOG()<<"          Compute merit in       "<<dt_commerit;
				LOG()<<"          Link branch in         "<<dt_link;
				LOG()<<"          Add new galaxies in    "<<dt_addgal;
				LOG()<<"          Delete broken gals in  "<<dt_delgal;
				LOG()<<" 	";
			}

#ifdef CTREE_USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);	
#endif
		}

		classify(vh, data, sinfo, vh.snapi, number);
		finalize(vh, data, dkey, tree, key, sinfo, vh.snapi, number);
		LOG()<<" finialize";
		LOG()<<" merging ?";

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
	std::vector<std::uint8_t> serialize_d1(const ControlSt& x) {
	    std::vector<std::uint8_t> buf;

	    //append_pod<CT_ID>(buf, x.id0);
	    //append_pod<CT_snap>(buf, x.snap0);
	    //append_pod<CT_ID>(buf, x.id);
	    //append_pod<CT_snap>(buf, x.snap);
	    
	    //append_vec<CT_double>(buf, x.pos);
	    //append_vec<CT_double>(buf, x.vel);
	    
	    //append_pod<CT_I32>(buf, x.stat);
	    //append_pod<CT_I32>(buf, x.detstat);

	    append_vec<PIDSt>(buf, x.p_list);
		append_pod<CT_I32>(buf, x.n_ptcl);

		//append_vec<ListSt>(buf, x.list);
		//append_pod<CT_I32>(buf, x.list_n);

		//append_pod<CT_I32>(buf, x.last_ind);
	    return buf;


	}

	void deserialize_d1(const std::vector<std::uint8_t>& buf, ControlSt& out) {
	    const std::uint8_t* p   = buf.data();
	    const std::uint8_t* end = p + buf.size();

	    //out.id0 		= read_pod<CT_ID>(p, end);
	    //out.snap0 		= read_pod<CT_snap>(p, end);
	    //out.id 			= read_pod<CT_ID>(p, end);
	    //out.snap 		= read_pod<CT_snap>(p, end);

	    //out.pos 		= read_vec<CT_double>(p, end);
	    //out.vel 		= read_vec<CT_double>(p, end);

	    //out.stat 		= read_pod<CT_I32>(p, end);
	    //out.detstat		= read_pod<CT_I32>(p, end);

	    out.p_list 		= read_vec<PIDSt>(p, end);
	    out.n_ptcl 		= read_pod<CT_I32>(p, end);

	    //out.list 		= read_vec<ListSt>(p, end);
	    //out.list_n 		= read_pod<CT_I32>(p, end);

	    //out.last_ind 	= read_pod<CT_I32>(p, end);

	    if (p != end) throw std::runtime_error("TFSt deserialize: trailing bytes");
	}

	// Serialize & Deserialize for Control Array2
	std::vector<std::uint8_t> serialize_d2(const ControlSt& x) {
	    std::vector<std::uint8_t> buf;

	    //append_pod<CT_ID>(buf, x.id0);
	    //append_pod<CT_snap>(buf, x.snap0);
	    //append_pod<CT_ID>(buf, x.id);
	    //append_pod<CT_snap>(buf, x.snap);
	    
	    //append_vec<CT_double>(buf, x.pos);
	    //append_vec<CT_double>(buf, x.vel);
	    
	    //append_pod<CT_I32>(buf, x.stat);
	    //append_pod<CT_I32>(buf, x.detstat);

	    //append_vec<PIDSt>(buf, x.p_list);
		//append_pod<CT_I32>(buf, x.n_ptcl);

		append_vec<ListSt>(buf, x.list);
		append_pod<CT_I32>(buf, x.list_n);

		//append_pod<CT_I32>(buf, x.last_ind);
	    return buf;


	}

	void deserialize_d2(const std::vector<std::uint8_t>& buf, ControlSt& out) {
	    const std::uint8_t* p   = buf.data();
	    const std::uint8_t* end = p + buf.size();

	    //out.id0 		= read_pod<CT_ID>(p, end);
	    //out.snap0 		= read_pod<CT_snap>(p, end);
	    //out.id 			= read_pod<CT_ID>(p, end);
	    //out.snap 		= read_pod<CT_snap>(p, end);

	    //out.pos 		= read_vec<CT_double>(p, end);
	    //out.vel 		= read_vec<CT_double>(p, end);

	    //out.stat 		= read_pod<CT_I32>(p, end);
	    //out.detstat		= read_pod<CT_I32>(p, end);

	    //out.p_list 		= read_vec<PIDSt>(p, end);
	    //out.n_ptcl 		= read_pod<CT_I32>(p, end);

	    out.list 		= read_vec<ListSt>(p, end);
	    out.list_n 		= read_pod<CT_I32>(p, end);

	    //out.last_ind 	= read_pod<CT_I32>(p, end);

	    if (p != end) throw std::runtime_error("TFSt deserialize: trailing bytes");
	}

	// Serialize & Deserialize for Next Array
	std::vector<std::uint8_t> serialize(const NextSt& x) {
	    std::vector<std::uint8_t> buf;

	    append_pod<CT_ID>(buf, x.id);
	    append_pod<CT_snap>(buf, x.snap);
	    append_pod<CT_Merit>(buf, x.merit);
	    return buf;
	}

	void deserialize(const std::vector<std::uint8_t>& buf, NextSt& out) {
	    const std::uint8_t* p   = buf.data();
	    const std::uint8_t* end = p + buf.size();

	    out.id 			= read_pod<CT_ID>(p, end);
	    out.snap 		= read_pod<CT_snap>(p, end);
	    out.merit		= read_pod<CT_Merit>(p, end);

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

	// Serialize & Deserialize for LinkJob
//	std::vector<std::uint8_t> serialize(const LinkJob& x) {
//		std::vector<std::uint8_t> buf;
//
//		append_pod<CT_I32>(buf, x.jobnum);
//		append_pod<CT_I32>(buf, x.ind);
//		append_pod<CT_I32>(buf, x.snap);
//		append_pod<CT_I32>(buf, x.id);
//		append_pod<CT_Merit>(buf, x.merit);
//		
//		append_vec<Tree::Tree_GID>(buf, x.c_tree.id);
//		append_vec<Tree::Tree_Snap>(buf, x.c_tree.snap);
//		append_vec<Tree::Tree_merit>(buf, x.c_tree.p_merit);
//		append_pod<Tree::Tree_I32>(buf, x.c_tree.endind);
//
//		return buf;
//    }
//
//    void deserialize(const std::vector<std::uint8_t>& buf, LinkJob& out) {
//	    const std::uint8_t* p   = buf.data();
//	    const std::uint8_t* end = p + buf.size();
//
//	    out.jobnum 	= read_pod<CT_I32>(p, end);
//	    out.ind  	= read_pod<CT_I32>(p, end);
//	    out.snap  	= read_pod<CT_I32>(p, end);
//	    out.id  	= read_pod<CT_I32>(p, end);
//	    out.merit  	= read_pod<CT_Merit>(p, end);
//
//	    Tree::TreeSt dumtree;
//
//	    dumtree.id  	= read_vec<Tree::Tree_GID>(p, end);
//	    dumtree.snap  	= read_vec<Tree::Tree_Snap>(p, end);
//	    dumtree.p_merit	= read_vec<Tree::Tree_merit>(p, end);
//	    dumtree.endind 	= read_pod<Tree::Tree_I32>(p, end);
//
//	    out.c_tree 	= std::move(dumtree);
//
//	    if (p != end) throw std::runtime_error("TFSt deserialize: trailing bytes");
//	}

	// Serialize & Deserialize for Data
	std::vector<std::uint8_t> serialize(const ControlSt& x) {
		std::vector<std::uint8_t> buf;

		append_pod<CT_I32>(buf, x.detstat);
		append_pod<CT_I32>(buf, x.stat);
		append_pod<CT_I32>(buf, x.n_ptcl);
		append_pod<CT_I32>(buf, x.id);
		append_pod<CT_I32>(buf, x.id0);
		append_pod<CT_I32>(buf, x.snap);
		append_pod<CT_I32>(buf, x.snap0);

		append_vec<ListSt>(buf, x.list);
		append_pod<CT_I32>(buf, x.list_n);
	
		return buf;	
	}
	void deserialize(const std::vector<std::uint8_t>& buf, ControlSt& out) {
	    const std::uint8_t* p   = buf.data();
	    const std::uint8_t* end = p + buf.size();

	    out.detstat 	= read_pod<CT_I32>(p, end);
	    out.stat 		= read_pod<CT_I32>(p, end);
	    out.n_ptcl 		= read_pod<CT_I32>(p, end);

	    out.id 			= read_pod<CT_I32>(p, end);
	    out.id0 		= read_pod<CT_I32>(p, end);
	    out.snap 		= read_pod<CT_I32>(p, end);
	    out.snap0 		= read_pod<CT_I32>(p, end);

	    out.list 		= read_vec<ListSt>(p, end);
	    out.list_n 		= read_pod<CT_I32>(p, end);
		
		//out.Free_plist();
		//out.n_ptcl 		= -1;
		if (p != end) throw std::runtime_error("TFSt deserialize: trailing bytes");
	}

	// Serialize & Deserialize for Tree
	std::vector<std::uint8_t> serialize(const Tree::TreeSt& x) {
		std::vector<std::uint8_t> buf;

		append_vec<Tree::Tree_GID>(buf, x.id);
		append_vec<Tree::Tree_Snap>(buf, x.snap);
		append_vec<Tree::Tree_merit>(buf, x.p_merit);

		append_pod<Tree::Tree_I32>(buf, x.father_bid);
		append_pod<Tree::Tree_I32>(buf, x.frag_bid);
		append_pod<Tree::Tree_I32>(buf, x.numprog);
		append_pod<Tree::Tree_I32>(buf, x.endind);

		append_pod<Tree::Tree_I32>(buf, x.stat);
		append_pod<Tree::Tree_I32>(buf, x.isfree);

		return buf;	
	}
	void deserialize(const std::vector<std::uint8_t>& buf, Tree::TreeSt& out) {
	    const std::uint8_t* p   = buf.data();
	    const std::uint8_t* end = p + buf.size();

	    out.id 		= read_vec<Tree::Tree_GID>(p, end);
	    out.snap 	= read_vec<Tree::Tree_Snap>(p, end);
	    out.p_merit = read_vec<Tree::Tree_merit>(p, end);

	    out.father_bid 	= read_pod<Tree::Tree_I32>(p, end);
	    out.frag_bid 	= read_pod<Tree::Tree_I32>(p, end);
	    out.numprog 	= read_pod<Tree::Tree_I32>(p, end);
	    out.endind 		= read_pod<Tree::Tree_I32>(p, end);

	    out.stat 		= read_pod<Tree::Tree_I32>(p, end);
	    out.isfree 		= read_pod<Tree::Tree_I32>(p, end);
		
//		if(out.isfree == 1){
//			out.id.resize(0);
//			out.snap.resize(0);
//			out.p_merit.resize(0);
//
//			out.m_id.resize(0);
//			out.m_snap.resize(0);
//			out.m_merit.resize(0);
//			out.m_bid.resize(0);
//
//			out.father_bid = -1;
//			out.numprog = 0;
//			out.endind = -1;
//
//			out.isfree 	= -1;
//		}

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

// For Debugging
	void savetree_ctree(const vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& treekey, CT_I32 snap_curr){
	    namespace fs = std::filesystem;

	    constexpr std::int32_t TAG_GID = savetree_gettype<Tree::Tree_GID>();
	    constexpr std::int32_t TAG_BID = savetree_gettype<Tree::Tree_BID>();
	    constexpr std::int32_t TAG_SNAP = savetree_gettype<Tree::Tree_Snap>();
	    constexpr std::int32_t TAG_MERIT = savetree_gettype<Tree::Tree_merit>();



	    //-----
	    // Save TreeKey
	    //-----
	    const std::string path = vh.out_dir + "/ctree_key_" + i4(snap_curr) + "_s2.dat";
	    LOG()<<"    Writing TreeKey in "<<path;

	    std::ofstream out(path, std::ios::binary);
	    if (!out) throw std::runtime_error("save_treekey_bin: cannot open " + path);

	    // First element as type specifer
	    //      I32 type
	    //      1 I32
	    //      2 I64
	    //      3 float
	    //      4 double
	    out.write(reinterpret_cast<const char*>(&TAG_BID), sizeof(TAG_BID));

	    // Second elements Key size in I64
	    std::int64_t n = static_cast<std::int64_t>(treekey.size());
	    out.write(reinterpret_cast<const char*>(&n), sizeof(n));

	    // Third elements as Tree Key
	    //std::int64_t tk = treekey[0].key;
	    std::int64_t tk = treekey[0];
	    out.write(reinterpret_cast<const char*>(&tk), sizeof(tk));

	    // Left for elements
	    std::vector<Tree::Tree_BID> inds;
	    inds.reserve(treekey.size());
	    //for (const auto& k : treekey) inds.push_back(k.ind);
	    for (const auto& k : treekey) inds.push_back(k);

	    if (!inds.empty()) {
	        out.write(reinterpret_cast<const char*>(inds.data()),
	                  inds.size() * sizeof(Tree::Tree_BID));
	    }


	    //-----
	    // Save Tree
	    //-----
	    const std::string path2 = vh.out_dir + "/ctree_tree_" + i4(snap_curr) + "_s2.dat";
	    LOG()<<"    Writing Tree in "<<path;

	    std::ofstream out2(path2, std::ios::binary);
	    if (!out2) throw std::runtime_error("save_tree_bin: cannot open " + path);

	    // 1 GID specifer
	    out2.write(reinterpret_cast<const char*>(&TAG_GID), sizeof(TAG_GID));
	    // 2 Snap specifer
	    out2.write(reinterpret_cast<const char*>(&TAG_SNAP), sizeof(TAG_SNAP));
	    // 3 BID specifer
	    out2.write(reinterpret_cast<const char*>(&TAG_BID), sizeof(TAG_BID));
	    // 4 merit specifer
	    out2.write(reinterpret_cast<const char*>(&TAG_MERIT), sizeof(TAG_MERIT));


	    // 5 Tree Size in I64
	    std::int64_t n2 = static_cast<std::int64_t>(tree.size());
	    out2.write(reinterpret_cast<const char*>(&n2), sizeof(n2));

	    // 6 End ind in I64
	    std::int64_t n3 = static_cast<std::int64_t>(tree[0].lind);
	    out2.write(reinterpret_cast<const char*>(&n3), sizeof(n3));

	    // loop for tree
	    
	    for(auto t : tree){
	        
	        Tree::Tree_I32 n_branch = t.endind + 1;
	        Tree::Tree_I32 n_numprg = t.numprog;
	        Tree::Tree_I32 father_bid = t.father_bid;
	        Tree::Tree_I32 tstat 	= t.stat;

	      
	        
	        // Size for the main branch
	        out2.write(reinterpret_cast<const char*>(&n_branch), sizeof(n_branch));

	        // skip empty tree
	        if(n_branch <= 0) continue;

	        
	        // Size for the merged branch
	        out2.write(reinterpret_cast<const char*>(&n_numprg), sizeof(n_numprg));

	        // Branch ID of father
	        out2.write(reinterpret_cast<const char*>(&father_bid), sizeof(father_bid));

	        // Branch Stat
	        out2.write(reinterpret_cast<const char*>(&tstat), sizeof(tstat));

	        // Branch_ID
	        for(std::int32_t i=0; i<n_branch; i++){
	            out2.write(reinterpret_cast<const char*>(&t.id[i]), sizeof(t.id[i]));
	        }

	        // Branch_snap
	        for(std::int32_t i=0; i<n_branch; i++){
	            out2.write(reinterpret_cast<const char*>(&t.snap[i]), sizeof(t.snap[i]));
	        }

	        // Branch_ID (progenitor)
	        for(std::int32_t i=0; i<n_branch; i++){
	            out2.write(reinterpret_cast<const char*>(&t.p_id[i]), sizeof(t.p_id[i]));
	        }

	        // Branch_Snap (progenitor)
	        for(std::int32_t i=0; i<n_branch; i++){
	            out2.write(reinterpret_cast<const char*>(&t.p_snap[i]), sizeof(t.p_snap[i]));
	        }

	        // Branch_Merit (progenitor)
	        for(std::int32_t i=0; i<n_branch; i++){
	            out2.write(reinterpret_cast<const char*>(&t.p_merit[i]), sizeof(t.p_merit[i]));
	        }

	        if(n_numprg>0){

	            // Merged Branch ID
	            for(std::int32_t i=0; i<n_numprg; i++){
	                out2.write(reinterpret_cast<const char*>(&t.m_id[i]), sizeof(t.m_id[i]));
	            }

	            // Merged Branch Snap
	            for(std::int32_t i=0; i<n_numprg; i++){
	                out2.write(reinterpret_cast<const char*>(&t.m_snap[i]), sizeof(t.m_snap[i]));
	            }

	            // Merged Branch Merit
	            for(std::int32_t i=0; i<n_numprg; i++){
	                out2.write(reinterpret_cast<const char*>(&t.m_merit[i]), sizeof(t.m_merit[i]));
	            }

	            // Merged Branch BID
	            for(std::int32_t i=0; i<n_numprg; i++){
	                out2.write(reinterpret_cast<const char*>(&t.m_bid[i]), sizeof(t.m_bid[i]));
	            }
	        }
	    }  
	}


	void loadtree_ctree(const vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& treekey, CT_I32 snap_curr){

	    int myrank  = mpi_rank();
	    // Simple version used
	    {
	        // file check
	        

	        const std::string path = vh.out_dir + "/ctree_key_" + i4(snap_curr) + "_s2.dat";
	        if(myrank==0){
	            LOG() << "    Reading TreeKey from " << path;
	        }

	        std::ifstream in(path, std::ios::binary);
	        if (!in) {
	            LOG()<<"load_treekey_bin: cannot open " + path;
	            u_stop();
	        }

	        // Read Data Type
	        std::int32_t bidtag;
	        loadtree_read(in, bidtag);
	    
	        // Read Nkey
	        std::int64_t nbid;
	        loadtree_read(in, nbid);

	        // TreeKey
	        std::int64_t keyval;
	        loadtree_read(in, keyval);

	        // Allocate
	        //auto treekey_load = loadtree_mkvector(bidtag, nbid);
	        std::vector<std::int32_t> treekey_load;
	        treekey_load.resize(nbid);
	        loadtree_vecread<std::int32_t>(in, treekey_load, nbid);

	        // Copy
	        //Tree::TreeKeyArray treekey;
	        treekey.resize(nbid);

	        for(std::int32_t i=0; i<nbid; i++){
	 
	            //treekey[i].ind  = treekey_load[i];
	            treekey[i]  = treekey_load[i];
	        }

	        //treekey[0].key = keyval;
	        treekey[0] = keyval;


	    } 


	    {
	        
	        const std::string path = vh.out_dir + "/ctree_tree_" + i4(snap_curr) + "_s2.dat";
	        if(myrank==0){
	            LOG() << "    Reading TreeKey from " << path;
	        }

	        std::ifstream in(path, std::ios::binary);
	        if (!in) {
	            LOG()<<"load_treekey_bin: cannot open " + path;
	            u_stop();
	        }

	        // Read Data Type
	        std::int32_t gidtag, snaptag, bidtag, mertag;
	        std::int32_t nbranch, nmerge, fid, stat;



	        loadtree_read(in, gidtag);
	        loadtree_read(in, snaptag);
	        loadtree_read(in, bidtag);
	        loadtree_read(in, mertag);
	    
	        // Read Ntree
	        std::int64_t ntree, lastind;
	        loadtree_read(in, ntree);
	        loadtree_read(in, lastind);

	        //tree.resize(ntree);
	        
	        tree.resize(lastind+vh.ctree_nstep);
	        tree[0].lind = lastind;

	        //for(std::int64_t i=0; i<ntree; i++){
	        for(std::int64_t i=0; i<tree[0].lind+1; i++){

	            loadtree_read(in, nbranch);

	            tree[i].endind = nbranch-1;

	            if(nbranch <= 0) continue;

	            loadtree_read(in, nmerge);
	            loadtree_read(in, fid);
	            loadtree_read(in, stat);

	            tree[i].father_bid = fid;
	            tree[i].numprog = nmerge;
	            tree[i].stat 	= stat;

	            tree[i].id.resize(nbranch);
	            tree[i].snap.resize(nbranch);
	            tree[i].p_id.resize(nbranch);
	            tree[i].p_snap.resize(nbranch);
	            tree[i].p_merit.resize(nbranch);

	            loadtree_vecread<std::int32_t>(in, tree[i].id, nbranch);
	            loadtree_vecread<std::int32_t>(in, tree[i].snap, nbranch);
	            loadtree_vecread<std::int32_t>(in, tree[i].p_id, nbranch);
	            loadtree_vecread<std::int32_t>(in, tree[i].p_snap, nbranch);
	            loadtree_vecread<double>(in, tree[i].p_merit, nbranch);

	            if(nmerge >= 1){
	                tree[i].m_id.resize(nmerge);
	                tree[i].m_snap.resize(nmerge);
	                tree[i].m_merit.resize(nmerge);
	                tree[i].m_bid.resize(nmerge);

	                loadtree_vecread<std::int32_t>(in, tree[i].m_id, nmerge);
	                loadtree_vecread<std::int32_t>(in, tree[i].m_snap, nmerge);
	                loadtree_vecread<double>(in, tree[i].m_merit, nmerge);
	                loadtree_vecread<std::int32_t>(in, tree[i].m_bid, nmerge);
	            }


	        }




	    }
	}
	void savedata(const vctree_set::Settings& vh, ControlArray& data, CT_I32 snap_curr){


	    namespace fs = std::filesystem;

	    //-----
	    // Save TreeKey
	    //-----
	    const std::string path = vh.out_dir + "/data_" + i4(snap_curr) + "_s2.dat";
	    LOG()<<"    Writing Data in "<<path;

	    std::ofstream out(path, std::ios::binary);
	    if (!out) throw std::runtime_error("save_treekey_bin: cannot open " + path);

	    // Data Size
	    std::int64_t n = static_cast<std::int64_t>(data.size());
	    out.write(reinterpret_cast<const char*>(&n), sizeof(n));

	    for(auto d : data){
	        out.write(reinterpret_cast<const char*>(&d.id), sizeof(d.id));
	        out.write(reinterpret_cast<const char*>(&d.id0), sizeof(d.id0));

	        out.write(reinterpret_cast<const char*>(&d.snap), sizeof(d.snap));
	        out.write(reinterpret_cast<const char*>(&d.snap0), sizeof(d.snap0));

	        out.write(reinterpret_cast<const char*>(&d.stat), sizeof(d.stat));

	        out.write(reinterpret_cast<const char*>(&d.n_ptcl), sizeof(d.n_ptcl));
	        out.write(reinterpret_cast<const char*>(&d.list_n), sizeof(d.list_n));

	        out.write(reinterpret_cast<const char*>(&d.last_ind), sizeof(d.last_ind));
	    }
	      
	}

	ControlArray loaddata(const vctree_set::Settings& vh, CT_I32 snap_curr){


	    const std::string path = vh.out_dir + "/data_" + i4(snap_curr) + "_s2.dat";
	    std::ifstream in(path, std::ios::binary);
	    if (!in) {
	        LOG()<<"load_treekey_bin: cannot open " + path;
	        u_stop();
	    }

	    std::int64_t ndata;

	    loadtree_read(in, ndata);

	    
	    ControlArray data = allocate(vh, ndata);

	    CT_I32 dumint;
	    for(auto &d: data){

	        loadtree_read(in, dumint);
	        d.id    = dumint;

	        loadtree_read(in, dumint);
	        d.id0    = dumint;

	        loadtree_read(in, dumint);
	        d.snap    = dumint;

	        loadtree_read(in, dumint);
	        d.snap0    = dumint;

	        loadtree_read(in, dumint);
	        d.stat    = dumint;

	        loadtree_read(in, dumint);
	        d.n_ptcl    = dumint;

	        loadtree_read(in, dumint);
	        d.list_n    = dumint;

	        loadtree_read(in, dumint);
	        d.last_ind    = dumint;
	    }

	    return data;
	}

	void validate_data(ControlArray& data, ControlArray& data2){

		LOG()<<"----- Size Check : "<<data.size()<<" / "<<data2.size();
		if(data.size() != data2.size()){
			LOG()<<"		size corrupted";
			u_stop();
		}

		LOG()<<"----- Lastind Check : "<<data[0].last_ind<<" / "<<data2[0].last_ind;
		if(data[0].last_ind != data2[0].last_ind){
			LOG()<<"		Last index corrupted";
			u_stop();
		}

		LOG()<<"----- Get into the element inspection";
		for(CT_I32 i=0; i<data[0].last_ind+1; i++){

			CT_I32 v1, v2;

			v1 	= data[i].id;
			v2 	= data2[i].id;

			if(v1 != v2){
				LOG()<<"		ID corrupted : "<<v1<<" / "<<v2;
				LOG()<<"		@ : "<<i;
				//u_stop();
			}

			v1 	= data[i].id0;
			v2 	= data2[i].id0;

			if(v1 != v2){
				LOG()<<"		ID0 corrupted : "<<v1<<" / "<<v2;
				LOG()<<"		@ : "<<i;
				//u_stop();
			}

			v1 	= data[i].snap;
			v2 	= data2[i].snap;

			if(v1 != v2){
				LOG()<<"		Snap corrupted : "<<v1<<" / "<<v2;
				LOG()<<"		@ : "<<i;
				//u_stop();
			}

			v1 	= data[i].snap0;
			v2 	= data2[i].snap0;

			if(v1 != v2){
				LOG()<<"		Snap0 corrupted : "<<v1<<" / "<<v2;
				LOG()<<"		@ : "<<i;
				//u_stop();
			}

			v1 	= data[i].stat;
			v2 	= data2[i].stat;

			if(v1 != v2){
				LOG()<<"		stat corrupted : "<<v1<<" / "<<v2;
				LOG()<<"		@ : "<<i;
				//u_stop();
			}

			v1 	= data[i].n_ptcl;
			v2 	= data2[i].n_ptcl;

			if(v1 != v2){
				LOG()<<"		n_ptcl corrupted : "<<v1<<" / "<<v2;
				LOG()<<"		@ : "<<i;
				//u_stop();
			}

			v1 	= data[i].list_n;
			v2 	= data2[i].list_n;

			if(v1 != v2){
				LOG()<<"		list_n corrupted : "<<v1<<" / "<<v2;
				LOG()<<"		@ : "<<i;
				//u_stop();
			}
		}

		LOG()<<"Data are purity ^-^";
		
	}

	void validate_tree(Tree::TreeArray& tree, Tree::TreeArray& tree2){

		LOG()<<"----- Size Check : "<<tree.size()<<" / "<<tree2.size();
		if(tree.size() != tree2.size()){
			LOG()<<"		size corrupted";
			//u_stop();
		}

		LOG()<<"----- Lastind Check : "<<tree[0].lind<<" / "<<tree2[0].lind;
		if(tree[0].lind != tree2[0].lind){
			LOG()<<"		Last index corrupted";
			u_stop();
		}

		LOG()<<"----- Get into the element inspection";
		for(CT_I32 i=0; i<tree2[0].lind+1; i++){

			CT_I32 v1, v2;

			v1 	= tree[i].endind;
			v2 	= tree2[i].endind;

			if(v1 != v2){
				LOG()<<"		endind corrupted : "<<v1<<" / "<<v2;
				LOG()<<"		@ : "<<i;
				//u_stop();
			}

			if(v1 < 0) continue;

			CT_I32 w1, w2;
			
			for(CT_I32 j=0; j<v1+1; j++){
				w1 	= tree[i].id[j];
				w2 	= tree2[i].id[j];

				if(w1 != w2){
					LOG()<<"		ID list corrupted : "<<w1<<" / "<<w2;
					LOG()<<"		@ : "<<i<<" / "<<j;
					//u_stop();
				}

				w1 	= tree[i].snap[j];
				w2 	= tree2[i].snap[j];

				if(w1 != w2){
					LOG()<<"		Snap list corrupted : "<<w1<<" / "<<w2;
					LOG()<<"		@ : "<<i<<" / "<<j;
					//u_stop();
				}
			}

			w1 	= tree[i].father_bid;
			w2 	= tree2[i].father_bid;

			if(w1 != w2){
				LOG()<<"		father_bid corrupted : "<<w1<<" / "<<w2;
				LOG()<<"		@ : "<<i;
				//u_stop();
			}

			w1 	= tree[i].frag_bid;
			w2 	= tree2[i].frag_bid;

			if(w1 != w2){
				//LOG()<<"		frag corrupted : "<<w1<<" / "<<w2;
				//LOG()<<"		@ : "<<i;
				//u_stop();
			}

			w1 	= tree[i].numprog;
			w2 	= tree2[i].numprog;

			if(w1 != w2){
				LOG()<<"		numprog corrupted : "<<w1<<" / "<<w2;
				LOG()<<"		@ : "<<i;
				//u_stop();
			}

			
			w1 	= tree[i].stat;
			w2 	= tree2[i].stat;

			if(w1 != w2){
				//LOG()<<"		stat_bid corrupted : "<<w1<<" / "<<w2;
				//LOG()<<"		@ : "<<i;
				//u_stop();
			}

			w1 	= tree[i].isfree;
			w2 	= tree2[i].isfree;

			if(w1 != w2){
				LOG()<<"		isfree corrupted : "<<w1<<" / "<<w2;
				LOG()<<"		@ : "<<i;
				//u_stop();
			}
		}

		LOG()<<"Tree are purity ^-^";
		
	}

	void validate_treekey(Tree::TreeKeyArray& key, Tree::TreeKeyArray& key2){

		LOG()<<"----- Size Check : "<<key.size()<<" / "<<key2.size();
		if(key.size() != key2.size()){
			LOG()<<"		size corrupted";
			u_stop();
		}

		for(CT_I32 i=0; i<(CT_I32) key.size(); i++){

			CT_I32 v1, v2;

			v1 	= key[i];
			v2 	= key2[i];

			if(v1 != v2){
				LOG()<<"		key corrupted : "<<v1<<" / "<<v2;
				LOG()<<"		@ : "<<i;
				//u_stop();
			}
		}

		LOG()<<"Key are purity ^-^";
		
	}
}