#include "global/io.h"
#include "global/allvar.h"
#include "makebr/makebr.h"
#include "utilities/utilities.h"
#include <numeric>
#include <limits>





namespace Makebr{
	// Some Helpers
	//template <typename T>
	//auto g_maxid(const std::vector<T>& ids){
	//	if (ids.empty()) throw std::runtime_error("g_maxid: empty container");
    //	return *std::max_element(ids.begin(), ids.end());
	//}

	void mainloop(vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& key){

		
		MBRHeader mbh;
		Treelog log;
		int myrank  = mpi_rank();

		mbh 	= setheader(vh);

		//----- Change the TF output name to synchronize the exact snapshot list
		if(myrank == 0) tfcname(vh);

		//----- Gather snapshot list
		std::vector<IO_VR::VRT_Snap> snaplist;
		IO_VR::VRT_Snap maxsnap = -1;

		for(IO_VR::VRT_Snap i=vh.snapi; i<vh.snapf+1; i++){
			if(is_snap(vh, i)){
				snaplist.push_back(i);
				if(i > maxsnap) maxsnap = i;
			}
		}

		
		//// Read Galaxy First
		std::vector<IO_VR::GalArray> g_all;
		std::vector<TFSt> t_all;

		
		g_all.resize(maxsnap+1);
		t_all.resize(maxsnap+1);



		auto t0 = std::chrono::steady_clock::now();
		if(myrank == 0){
			LOG() <<"    Makebr) Read All TreeFrog and Galaxies first";
		}

#ifdef CTREE_USE_MPI
	    int rank = 0, size = 1;
	    IO_VR::VRT_GID local_max = -1;
	    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	    MPI_Comm_size(MPI_COMM_WORLD, &size);

	    for(IO_VR::VRT_Snap s : snaplist){
	    	if(s % size == rank){
	    		//LOG()<<rank<<" / "<<s;
				g_all[s] 	= IO_VR::r_gal(vh, s, -1, false);
				if(s<maxsnap){ t_all[s] 	= readtreefrog(vh, mbh, s); }
				else{ t_all[s] = readtreefrog(vh, mbh, snaplist[0]);} // dummy read for the last snapshot

				local_max 	= std::max(local_max, getmaxid(g_all[s]));
			}
			//log.max_id 	= std::max(log.max_id, getmaxid(g_curr));
			
		}
    
    	MPI_Barrier(MPI_COMM_WORLD);

    	IO_VR::VRT_GID global_max = -1;
    	MPI_Allreduce(&local_max, &global_max, 1,
              mpi_type::type<IO_VR::VRT_GID>(), MPI_MAX, MPI_COMM_WORLD);

    	log.max_id = global_max; 

    	//Synchronize
    	for (IO_VR::VRT_Snap s : snaplist) {
        	int owner = s % size;

	        // ---- TF ----
	        std::vector<std::uint8_t> blob_t;
	        if (rank == owner) blob_t = serialize(t_all[s]);
	        bcast_blob_from_owner(owner, blob_t);
	        if (rank != owner) deserialize(blob_t, t_all[s]);

	        // ---- Gal ----
	        std::vector<std::uint8_t> blob_g;
	        if (rank == owner) blob_g = serialize(g_all[s]);
	        bcast_blob_from_owner(owner, blob_g);
	        if (rank != owner) deserialize(blob_g, g_all[s]);
    	}
    	MPI_Barrier(MPI_COMM_WORLD);
#else
		for(IO_VR::VRT_Snap s : snaplist){
			//LOG()<<s;
			g_all[s] 	= IO_VR::r_gal(vh, s, -1, false);
			if(s<maxsnap) t_all[s] 	= readtreefrog(vh, mbh, s);
			log.max_id 	= std::max(log.max_id, getmaxid(g_all[s]));
		}
#endif

		if(myrank == 0){
			auto t1 = std::chrono::steady_clock::now();
			double dt_sec = std::chrono::duration<double>(t1 - t0).count();
			LOG() << "      Done in : " << dt_sec << " [sec]";
		}
		

		//----- Update Tree Key
		//		key = minimum of 10^n greater than maxsnap
		t0 = std::chrono::steady_clock::now();
		if(myrank == 0){
			LOG() <<"    Makebr) Make Tree Key";
		}

		//Tree::Tree_I64 p = 1;
    	//while (p <= (Tree::Tree_I64) maxsnap) {
        //	p *= 10;
    	//}
    	//Tree::Tree_I64 p = maxsnap+1;
		vh.treekey = maxsnap+1;


		
		//----- Allocate TreeArray
		Tree::Tree_BID treesize;
		if(treesize > 2147483647 && sizeof(Tree::Tree_BID)==4){
			if(myrank == 0){
				LOG()<<" key will be overflow: change Tree_BID to I64";
				u_stop();
			}
		}

		treesize 	= maxsnap + vh.treekey * log.max_id + 1;

		tree.resize(treesize);
		key.resize(treesize, {-1});
		//key[0].key 	= vh.treekey;		// key also stored in the first element
		key[0]		= vh.treekey;
		tree[0].lind = 0;		// last index stored in the first element

		if(myrank == 0){
			auto t1 = std::chrono::steady_clock::now();
			double dt_sec = std::chrono::duration<double>(t1 - t0).count();
			LOG() << "      Done in : " << dt_sec << " [sec]";
		}

		//----- Initialize Tree
		t0 = std::chrono::steady_clock::now();
		if(myrank == 0){
			LOG() <<"    Makebr) Tree Initialize";
		}

		for(IO_VR::GalSt& g : g_all[vh.snapi]){
			Tree::treeinit(tree, key, (Tree::Tree_Snap) vh.snapi, (Tree::Tree_GID) g.id);
		}

		if(myrank == 0){
			auto t1 = std::chrono::steady_clock::now();
			double dt_sec = std::chrono::duration<double>(t1 - t0).count();
			LOG() << "      Done in : " << dt_sec << " [sec]";
		}

#ifdef CTREE_USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif

		//// Adjust snapshot order
		if(vh.treedir == "DES"){
			// nothing done
		} else if (vh.treedir == "PRG"){
			// reverse
			std::reverse(snaplist.begin(), snaplist.end());
		} else {
			LOG()<<"Wrong Tree direction: treedir";
			u_stop();
		}
		
		//// Main Loop
		t0 = std::chrono::steady_clock::now();
		if(myrank == 0){
			LOG() <<"    Makebr) Main Loop starts";
		}

		for (std::size_t i = 0; i < snaplist.size()-1; ++i) {
    		IO_VR::VRT_Snap s = snaplist[i];


    		// End loop if this is the last snapshot
			//if(s == vh.snapf){
			//	// 123123 TODO HERE
			//	LOG()<<"HERE SHOULD BE DONE";
			//	u_stop();
			//}

			//// Make New Tree
			TFSt& t_curr = t_all[s];
			IO_VR::GalArray& g_next = g_all[snaplist[i+1]];

			for(TF_id l=0; l<(TF_id) t_curr.id.size();l++){
				if(t_curr.num[l] == 0){
					if(!Tree::istree(key, (Tree::Tree_Snap) s, (Tree::Tree_GID) t_curr.id[l])){
						Tree::treeinit(tree, key, (Tree::Tree_Snap) s, (Tree::Tree_GID) t_curr.id[l]);
					}

					log.n_new ++;
				}
				
			}

#ifdef CTREE_USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif		
	
			
			// Find Connection point 2
			{
			EvolArray Evoldum;
			Evoldum.resize(t_curr.id.size());

#ifdef CTREE_USE_MPI
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	    	MPI_Comm_size(MPI_COMM_WORLD, &size);

	    	std::vector<std::pair<TF_id, EvolSt>> local_updates;

	    	int local_nall  = 0;


			for(TF_id l=0; l<(TF_id) t_curr.id.size();l++){
				if(t_curr.num[l] == 0) continue;
				
				int owner = l % size;
				if(owner != rank) continue;

				local_nall ++;

				TF_off ind1 = t_curr.off[l];
				TF_off ind2 = t_curr.off[l] + t_curr.num[l] - 1;

				std::vector<TF_res> dum_id = slice<TF_res>(t_curr.res, ind1, ind2);
				std::vector<TF_npart> dum_part = slice<TF_npart>(t_curr.npart, ind1, ind2);
				std::vector<TF_merit> dum_merit = slice<TF_merit>(t_curr.merit, ind1, ind2);

				//----- Connection candidate with the maximum merit
				auto dum_maxmerit 	= std::max_element(dum_merit.begin(), dum_merit.end());
    			std::size_t idx = static_cast<std::size_t>(std::distance(dum_merit.begin(), dum_maxmerit));

    			IO_VR::VRT_GID best_id   = dum_id[idx];
    			IO_VR::VRT_I32 best_part = dum_part[idx];
    			TF_merit best_merit= *dum_maxmerit;

    			//----- End connection due to a too low merit
    			if(best_merit < vh.meritlimit) continue;


    			std::vector<IO_VR::VRT_I32> where_idx;
				
    			for (std::size_t m = 0; m < g_next.size(); ++m) {
				    const auto& g = g_next[m];
				    if (g.id == best_id && g.npart == best_part) {
				        where_idx.push_back(m);
				    }
				}

				if(where_idx.size()==0){
					// Not matched
					continue;
				}else if(where_idx.size()>2){
					// n_link > 1 or error?
					LOG()<<"Not implemented yet or error occurs";
					u_stop();
				}

				Evoldum[l].idc 	= t_curr.id[l];
				Evoldum[l].idn 	= g_next[where_idx[0]].id;
				Evoldum[l].snapc= s;
				Evoldum[l].snapn= snaplist[i+1];
				Evoldum[l].merit= best_merit;

    			local_updates.emplace_back(l, Evoldum[l]);

			}

			MPI_Barrier(MPI_COMM_WORLD);

			int global_nall = 0;
    		MPI_Allreduce(&local_nall, &global_nall, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    		log.n_all += global_nall;

			//----- Synchronize Evoldum
			for (int owner = 0; owner < size; ++owner) {
    			std::vector<std::uint8_t> blob;
    			if (rank == owner) {
		       		blob = serialize(local_updates);
			    }
			    bcast_blob_from_owner(owner, blob);
			    if (rank != owner) {
			        deserialize(blob, Evoldum);
			    }
			}
			MPI_Barrier(MPI_COMM_WORLD);
	
//for(TF_id l=0; l<(TF_id) t_curr.id.size();l++){
//	if(l<100 && myrank==0 && s==49){
//	LOG()<<Evoldum[l].snapc<<" & "<<Evoldum[l].idc<<" --> "<<Evoldum[l].snapn<<" / "<<Evoldum[l].idn<<" with "<<Evoldum[l].merit;
//}
//}
#else	

			for(TF_id l=0; l<(TF_id) t_curr.id.size();l++){
				if(t_curr.num[l] == 0) continue;
				
				log.n_all ++;

				TF_off ind1 = t_curr.off[l];
				TF_off ind2 = t_curr.off[l] + t_curr.num[l] - 1;

				std::vector<TF_res> dum_id = slice<TF_res>(t_curr.res, ind1, ind2);
				std::vector<TF_npart> dum_part = slice<TF_npart>(t_curr.npart, ind1, ind2);
				std::vector<TF_merit> dum_merit = slice<TF_merit>(t_curr.merit, ind1, ind2);

				//----- Connection candidate with the maximum merit
				auto dum_maxmerit 	= std::max_element(dum_merit.begin(), dum_merit.end());
    			std::size_t idx = static_cast<std::size_t>(std::distance(dum_merit.begin(), dum_maxmerit));

    			IO_VR::VRT_GID best_id   = dum_id[idx];
    			IO_VR::VRT_I32 best_part = dum_part[idx];
    			TF_merit best_merit= *dum_maxmerit;

    			//----- End connection due to a too low merit
    			if(best_merit < vh.meritlimit) continue;


    			std::vector<IO_VR::VRT_I32> where_idx;
				
    			for (std::size_t m = 0; m < g_next.size(); ++m) {
				    const auto& g = g_next[m];
				    if (g.id == best_id && g.npart == best_part) {
				        where_idx.push_back(m);
				    }
				}

				if(where_idx.size()==0){
					// Not matched
					continue;
				}else if(where_idx.size()>2){
					// n_link > 1 or error?
					LOG()<<"Not implemented yet or error occurs";
					u_stop();
				}

				Evoldum[l].idc 	= t_curr.id[l];
				Evoldum[l].idn 	= g_next[where_idx[0]].id;
				Evoldum[l].snapc= s;
				Evoldum[l].snapn= snaplist[i+1];
				Evoldum[l].merit= best_merit;

				//if(l<100 && myrank==0 && s==49){
				//	LOG()<<Evoldum[l].snapc<<" & "<<Evoldum[l].idc<<" --> "<<Evoldum[l].snapn<<" / "<<Evoldum[l].idn<<" with "<<Evoldum[l].merit;
				//}
			}
			//LOG()<<log.n_all;
#endif
			// here from using evoldum, determine the branch

			//----- Sort Evoldum
			std::sort(Evoldum.begin(), Evoldum.end(), [](const EvolSt& a, const EvolSt& b) {
              	if (a.idn != b.idn) return a.idn < b.idn;
              	return a.idc < b.idc;
          	});

			std::vector<IO_VR::VRT_I32> uni_idx;
          	for(IO_VR::VRT_GID l=1; l<(IO_VR::VRT_GID) Evoldum.size(); l++){
          		if(Evoldum[l].idn == 0) continue;
          		if(Evoldum[l].idn > Evoldum[l-1].idn) uni_idx.push_back(l);
          	}
          	uni_idx.push_back(Evoldum.size());

          	//----- Determine Connection point with the maximum merit

          	// MPI implemented here?
          	for(IO_VR::VRT_GID l=0; l<(IO_VR::VRT_GID) uni_idx.size()-1; l++){
          		IO_VR::VRT_GID ind1 = uni_idx[l];
          		IO_VR::VRT_GID ind2 = uni_idx[l+1]-1;
          		log.n_link ++;


          		//----- Single connection
          		if(ind1==ind2){
          			Tree::treeinput(tree, key, (Tree::Tree_Snap) Evoldum[ind1].snapc, (Tree::Tree_GID) Evoldum[ind1].idc, 
							(Tree::Tree_Snap) Evoldum[ind1].snapn, (Tree::Tree_Snap) Evoldum[ind1].idn, Evoldum[ind1].merit);
          			continue;
          		}


          		//----- Multiple connection
          		EvolArray Evoltmp(Evoldum.begin() + ind1, Evoldum.begin() + (ind2+1));


				//----- Sort by merit
				std::sort(Evoltmp.begin(), Evoltmp.end(), [](const EvolSt& a, const EvolSt& b) {
              		if (a.merit != b.merit) return a.merit > b.merit;
              		return a.idc < b.idc;
          		});

				//----- Primary connection
				Tree::treeinput(tree, key, (Tree::Tree_Snap) Evoltmp[0].snapc, (Tree::Tree_GID) Evoltmp[0].idc, 
							(Tree::Tree_Snap) Evoltmp[0].snapn, (Tree::Tree_Snap) Evoltmp[0].idn, Evoltmp[0].merit);


				//----- Merge other branch to this
				Tree::Tree_BID keyind0, keyind;

				//keyval = Evoltmp[0].snapc + key[0].key * Evoltmp[0].idc;
				//keyind0 = key[keyval].ind;
				keyind0 	= Tree::getkey(key, Evoltmp[0].snapc, Evoltmp[0].idc);

				for(IO_VR::VRT_GID m=1; m<(IO_VR::VRT_GID) Evoltmp.size(); m++){

					if(!Tree::istree(key, Evoltmp[m].snapc, Evoltmp[m].idc)){
						Tree::treeinit(tree, key, Evoltmp[m].snapc, Evoltmp[m].idc);
					}

					keyind = Tree::getkey(key, Evoltmp[m].snapc, Evoltmp[m].idc);
					//keyind = key[keyval].ind;

					if((std::int32_t) tree[keyind].id.size() < vh.minbranchlength){
						continue;
					}

					tree[keyind0].m_id.push_back(Evoltmp[m].idc);
					tree[keyind0].m_snap.push_back(Evoltmp[m].snapc);
					tree[keyind0].m_merit.push_back(Evoltmp[m].merit);
					tree[keyind0].m_bid.push_back(keyind);
					tree[keyind0].numprog ++;
					tree[keyind].father_bid = keyind0;

					log.n_broken ++;
				}
	
          	}


			} // For localizing EvolDum

			if(myrank == 0) LOG() <<"      Connection done for Snapshot = "<<s;
			if(myrank == 0) LOG() <<"        (n_all = "<<log.n_all<<" / n_link = "<<log.n_link<<" / n_broken = "<<log.n_broken<<" )";
			if(myrank == 0){
				LOG()<<"        Memory report";
				LOG()<<"			"<<how_big<Tree::TreeSt>(tree)<<" GB for tree";
				LOG()<<"			"<<how_big<Tree::TreeKey>(key)<<" GB for key";
				LOG()<<" ";
			}
		}

		if(myrank == 0){
			auto t1 = std::chrono::steady_clock::now();
			double dt_sec = std::chrono::duration<double>(t1 - t0).count();
			LOG() << "      Done in : " << dt_sec << " [sec]";
		}
		//return tree;
	}

	// Some utilities
	IO_VR::VRT_GID getmaxid(const IO_VR::GalArray& gal) {
    	if (gal.empty()) throw std::runtime_error("gal is empty");
    	auto it = std::max_element(gal.begin(), gal.end(),
                               [](const IO_VR::GalSt& a, const IO_VR::GalSt& b){
                                   return a.id < b.id;
                               });
    	return it->id;
	}


	MBR_Snap findnextsnap(const vctree_set::Settings& vh, const MBR_Snap& snap_curr){
		MBR_Snap snap_next 	= snap_curr;
		while (true){
			if(vh.treedir == "DES"){
				snap_next ++;
			}else if(vh.treedir == "PRG"){
				snap_next --;
			}else{
				LOG()<<		"Wrong Tree Direction";
				u_stop();
			}

			// Filename
        	char buf[256];
        	std::snprintf(buf, sizeof(buf), "tree.snapshot_%04dVELOCIraptor.tree", snap_next);
        	std::string dumfname = vh.vr_dir_tree + "/" + buf;
        		
        	// Is File?
        	if (fs::exists(dumfname)) {
            	return snap_next;
        	}

        	// Range Check
        	if (snap_next < vh.snapi || snap_next > vh.snapf) {
            	return -1;
        	}
		}
	}

	bool is_snap(const vctree_set::Settings& vh, const IO_VR::VRT_Snap snap_curr){
		char buf[256];
    	std::snprintf(buf, sizeof(buf), "tree.snapshot_%04dVELOCIraptor.tree", snap_curr);
    	std::string dumfname = vh.vr_dir_tree + "/" + buf;
    		
    	// Is File?
    	if (fs::exists(dumfname)) {
        	return true;
    	} else{
    		return false;
    	}
	}

	MBRHeader setheader(const vctree_set::Settings& vh){
		MBRHeader MB_header;

		int32_t snapmin = std::min(vh.snapi, vh.snapf);
  		int32_t snapmax = std::max(vh.snapi, vh.snapf);

		if(vh.treedir == "DES"){
      		MB_header.snapi  = snapmin;
      		MB_header.snapf  = snapmax;
  
      		MB_header.tag_num    = "NumDesc";
      		MB_header.tag_off    = "DescOffsets";
      		MB_header.tag_result = "Descendants";
      		MB_header.tag_npart  = "DescNpart";
      		MB_header.tag_merit  = "Merits";
      		MB_header.tag_nlink  = "Nsteps_search_new_links";
    	}
    	else if(vh.treedir == "PRG"){
      		MB_header.snapi  = snapmax;
      		MB_header.snapf  = snapmin;
  
      		MB_header.tag_num    = "NumProgen";
      		MB_header.tag_off    = "ProgenOffsets";
      		MB_header.tag_result = "Progenitors";
      		MB_header.tag_npart  = "ProgenNpart";
      		MB_header.tag_merit  = "Merits";
      		MB_header.tag_nlink  = "Nsteps_search_new_links";
  
    	}
    	else{
    		LOG()<<"	Wrong Tree direction 'treedir' ";
    		u_stop();
      		
    	}
    	return MB_header;
	}

	bool tfcname(const vctree_set::Settings& vh){

		try {
    		const fs::path tfout = fs::path(vh.vr_dir_tree);

    		// Read snaplist
    		std::vector<int> slist;
    		{
      			std::ifstream in(tfout / "tree.snaplist.txt");
      			if (!in) {
        			std::cerr << "[tfcname] cannot open: " << (tfout / "tree.snaplist.txt") << "\n";
        			return false;
      			}
      			std::string line;
      			while (std::getline(in, line)) {
        			int v;
        			if (extract_snap_from_line(line, v)) slist.push_back(v);
      			}
    		}

    		// list tfile
    		struct Item { int num; fs::path path; };
    		std::vector<Item> items;
    		for (const auto& de : fs::directory_iterator(tfout)) {
      			if (!de.is_regular_file()) continue;
      			const auto fn = de.path().filename().string();
      			if (fn.rfind("tree.snapshot", 0) == 0 && de.path().extension() == ".tree") {
        			int num;
        			if (extract_number_from_filename(fn, num)) {
          				items.push_back({num, de.path()});
        			}
      			}
    		}

    		if (items.size() != slist.size()) {
      			LOG() << "[tfcname] count mismatch: "
                	<< "files=" << items.size()
                	<< ", snaplist=" << slist.size() << "\n";
      			return false;
    		}

    		// Sort
    		std::sort(items.begin(), items.end(),
              [](const Item& a, const Item& b){ return a.num < b.num; });

    		// rename
    		for (std::size_t i = 0; i < items.size(); ++i) {
      			const std::string newbase = "tree.snapshot_" + i4(slist[i]) + "VELOCIraptor";
      			const fs::path target = tfout / newbase;

      			if (items[i].path.filename() == target.filename()) continue; // 이미 같으면 스킵

      			// 충돌 방지: 동일 이름 파일이 있다면 먼저 제거/백업 등 필요 시 처리
      			if (fs::exists(target)) {
        			LOG() << "[tfcname] target already exists, skipping: " << target << "\n";
        			continue;
      			}
      			fs::rename(items[i].path, target);  // mv org -> no-ext name
      			items[i].path = target;             // 경로 갱신
    		}

    		// add .tree
    		for (auto& it : items) {
      			const fs::path target = fs::path(it.path.string() + ".tree");
      			if (fs::exists(target)) {
        			std::cerr << "[tfcname] target already exists, skipping: " << target << "\n";
        			continue;
      			}
      			fs::rename(it.path, target);        // mv name -> name.tree
      			it.path = target;
    		}

    		return true;
  		} catch (const std::exception& e) {
    		std::cerr << "[tfcname] exception: " << e.what() << "\n";
    		return false;
  		}
	}
	
	//-----
	// Read TreeFrog
	//-----
	TFSt readtreefrog(const vctree_set::Settings& vh, const MBRHeader& mbh, const MBR_Snap& snap_curr){


		std::string tfname 	= vh.vr_dir_tree + "/tree.snapshot_" + i4(snap_curr) + "VELOCIraptor.tree";
	
	    hid_t file_id = H5Fopen(tfname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	    if (file_id < 0) throw std::runtime_error("Failed to open file: " + tfname);
	
	    try {
	        //std::vector<TF_num> 	v_num 	= IO_VR::VR_HDF5_rdbyname<TF_num>(file_id, mbh.tag_num);
	        //std::vector<TF_off> 	v_off 	= IO_VR::VR_HDF5_rdbyname<TF_off>(file_id, mbh.tag_off);
	        //std::vector<TF_res> 	v_res 	= IO_VR::VR_HDF5_rdbyname<TF_res>(file_id, mbh.tag_result);
	        //std::vector<TF_merit> 	v_merit	= IO_VR::VR_HDF5_rdbyname<TF_merit>(file_id, mbh.tag_merit);
	        //std::vector<TF_npart> 	v_npart	= IO_VR::VR_HDF5_rdbyname<TF_npart>(file_id, mbh.tag_npart);
	        //std::vector<TF_nlink> 	v_nlink	= IO_VR::VR_HDF5_rdbyname<TF_nlink>(file_id, mbh.tag_nlink);
	        //std::vector<TF_id> 		v_id 	= IO_VR::VR_HDF5_rdbyname<TF_id>(file_id, mbh.tag_id);
	

	        TFSt TFData;

	       	TFData.num 			= IO_VR::VR_HDF5_rdbyname<TF_num>(file_id, mbh.tag_num);
	       	TFData.off 			= IO_VR::VR_HDF5_rdbyname<TF_off>(file_id, mbh.tag_off);
	       	TFData.res 			= IO_VR::VR_HDF5_rdbyname<TF_res>(file_id, mbh.tag_result);
	       	TFData.merit 		= IO_VR::VR_HDF5_rdbyname<TF_merit>(file_id, mbh.tag_merit);
	       	TFData.npart 		= IO_VR::VR_HDF5_rdbyname<TF_npart>(file_id, mbh.tag_npart);
	       	TFData.nlink 		= IO_VR::VR_HDF5_rdbyattr<TF_nlink>(file_id, mbh.tag_nlink);
	       	TFData.id 			= IO_VR::VR_HDF5_rdbyname<TF_id>(file_id, mbh.tag_id);



	        H5Fclose(file_id);
	        return TFData;
	    } catch (...) {
	        H5Fclose(file_id);
	        throw;
	    }
	}

	template <typename T> std::vector<T> slice(const std::vector<T>& v, std::size_t ind1, std::size_t ind2_inclusive){
    	if (ind1 >= v.size() || ind1 > ind2_inclusive) return {};
    	std::size_t last_exclusive = std::min(ind2_inclusive + 1, v.size());
    	return std::vector<T>(v.begin() + ind1, v.begin() + last_exclusive);
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

    // Serialize & Deserialize for TFSt
	std::vector<std::uint8_t> serialize(const Makebr::TFSt& x) {
	    std::vector<std::uint8_t> buf;
	    append_vec(buf, x.num);
	    append_vec(buf, x.off);
	    append_vec(buf, x.res);
	    append_vec(buf, x.merit);
	    append_vec(buf, x.npart);
	    append_vec(buf, x.nlink);
	    append_vec(buf, x.id);
	    return buf;
	}

	void deserialize(const std::vector<std::uint8_t>& buf, Makebr::TFSt& out) {
	    const std::uint8_t* p   = buf.data();
	    const std::uint8_t* end = p + buf.size();
	    out.num   = read_vec<Makebr::TF_num>(p, end);
	    out.off   = read_vec<Makebr::TF_off>(p, end);
	    out.res   = read_vec<Makebr::TF_res>(p, end);
	    out.merit = read_vec<Makebr::TF_merit>(p, end);
	    out.npart = read_vec<Makebr::TF_npart>(p, end);
	    out.nlink = read_vec<Makebr::TF_nlink>(p, end);
	    out.id    = read_vec<Makebr::TF_id>(p, end);
	    if (p != end) throw std::runtime_error("TFSt deserialize: trailing bytes");
	}
	


	// Serialize & Deserialize for GalArray
	std::vector<std::uint8_t> serialize(const IO_VR::GalArray& A) {
	    std::vector<std::uint8_t> buf;
	    std::uint64_t n = static_cast<std::uint64_t>(A.size());
	    append_pod(buf, n);
	    for (const auto& g : A) {
	        append_pod(buf, g.id);
	        append_pod(buf, g.snap);
	        append_pod(buf, g.npart);
	        append_vec(buf, g.pid); // pid is vector<VRT_PID>
	    }
	    return buf;
	}

	void deserialize(const std::vector<std::uint8_t>& buf, IO_VR::GalArray& A) {
	    const std::uint8_t* p   = buf.data();
	    const std::uint8_t* end = p + buf.size();
	    std::uint64_t n = read_pod<std::uint64_t>(p, end);
	    A.resize(static_cast<size_t>(n));
	    for (size_t i = 0; i < A.size(); ++i) {
	        A[i].id    = read_pod<IO_VR::VRT_GID>(p, end);
	        A[i].snap  = read_pod<IO_VR::VRT_Snap>(p, end);
	        A[i].npart = read_pod<IO_VR::VRT_I32>(p, end);
	        A[i].pid   = read_vec<IO_VR::VRT_PID>(p, end);
	    }
	    if (p != end) throw std::runtime_error("GalArray deserialize: trailing bytes");
	}
	
	// Serialize & Deserialize for Evoldum
	std::vector<std::uint8_t> serialize(const std::vector<std::pair<TF_id, EvolSt>>& upd) {
	    std::vector<std::uint8_t> buf;
	    std::uint64_t n = upd.size();
	    append_pod(buf, n);
	    for (const auto& pr : upd) {
	        append_pod(buf, pr.first);   // index
	        append_pod(buf, pr.second);  // EvolSt
	    }
	    return buf;
	}

	void deserialize(const std::vector<std::uint8_t>& blob, EvolArray& evoldum) {
    	const std::uint8_t* p = blob.data();
   		const std::uint8_t* e = p + blob.size();
    	std::uint64_t n = read_pod<std::uint64_t>(p, e);
	    for (std::uint64_t k = 0; k < n; ++k) {
	        TF_id idx = read_pod<TF_id>(p, e);
	        EvolSt it = read_pod<EvolSt>(p, e);
	        if (idx >= (TF_id) evoldum.size()) evoldum.resize(static_cast<std::size_t>(idx) + 1);
	        evoldum[idx] = it;
	    }
    	if (p != e) throw std::runtime_error("deserialize: trailing bytes");
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


















