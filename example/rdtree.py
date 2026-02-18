import numpy as np
import pickle as pickle
import argparse
from typing import Optional
from pathlib import Path

class rdctree:
	def __init__(self):
		self.keyname 	= "./ctree_key.dat"
		self.treename 	= "./ctree_tree.dat"
		self.varname 	= "./ctree.pkl"

	def get_dtype(self, tag):
		"""
		Tag meaning (as in your writer):
		  1 -> int32
		  2 -> int64
		  3 -> float32
		  4 -> float64
		"""
		
		if(tag == 1):
			return np.int32
		if(tag == 2):
			return np.int64
		if(tag == 3):
			return np.dtype('<f4')
		if(tag == 4):
			return np.dtype('<f8')

	def get_array(self, type, nn):

		return np.zeros(nn, dtype=self.get_dtype(type))


	def rdkey(self):
		with open(self.keyname, "rb") as f:
			bidtag 	= np.fromfile(f, dtype=np.int32, count=1)
			nkey 	= np.fromfile(f, dtype=np.int64, count=1)
			keyv 	= np.fromfile(f, dtype=np.int64, count=1)

			self.key 	= np.fromfile(f, dtype=np.int32, count=nkey[0])

		self.key[0]	= keyv[0]

	def rdtree(self):
		with open(self.treename, "rb") as f:

			gidtag 	= np.fromfile(f, dtype=np.int32, count=1)
			snaptag = np.fromfile(f, dtype=np.int32, count=1)
			bidtag 	= np.fromfile(f, dtype=np.int32, count=1)
			mertag 	= np.fromfile(f, dtype=np.int32, count=1)

			#ntree 	= np.fromfile(f, dtype=np.int64, count=1)
			lind 	= np.fromfile(f, dtype=np.int64, count=1)


			tree: List[Dict[str, Any]] = []

			for _i in range(lind[0]+1):
				
				nbranch 	= (np.fromfile(f, dtype=np.int32, count=1))[0]

				if nbranch <= 0:
					tree.append({"stat": -1})
					continue

				nmerge 	= (np.fromfile(f, dtype=np.int32, count=1))[0]
				fid 	= (np.fromfile(f, dtype=self.get_dtype(bidtag), count=1))[0]
				tstat 	= np.fromfile(f, dtype=np.int32, count=1)

				idtmp	= np.fromfile(f, dtype=self.get_dtype(gidtag), count=nbranch)
				snaptmp = np.fromfile(f, dtype=self.get_dtype(snaptag), count=nbranch)
				#pidtmp 	= np.fromfile(f, dtype=self.get_dtype(gidtag), count=nbranch)
				#psnaptmp= np.fromfile(f, dtype=self.get_dtype(snaptag), count=nbranch)
				mertmp 	= np.fromfile(f, dtype=self.get_dtype(mertag), count=nbranch)

				if(nmerge>=1):
					midtmp 	= np.fromfile(f, dtype=self.get_dtype(gidtag), count=nmerge)
					msnaptmp= np.fromfile(f, dtype=self.get_dtype(snaptag), count=nmerge)
					mmertmp = np.fromfile(f, dtype=self.get_dtype(mertag), count=nmerge)
					mbidtmp = np.fromfile(f, dtype=self.get_dtype(bidtag), count=nmerge)
				else:
					midtmp = self.get_array(gidtag, 1)
					msnaptmp = self.get_array(snaptag, 1)
					mmertmp = self.get_array(mertag, 1)
					mbidtmp = self.get_array(bidtag, 1)

				tree.append(
				{
					"br_len": np.int32(nbranch),
					"father_ID": fid,
					"n_mergebr": np.int32(nmerge),
					"id": idtmp,
					"snap": snaptmp,
					#"p_id": pidtmp,
					#"p_snap": psnaptmp,
					"merit": mertmp,
					"m_id": midtmp,
					"m_snap": msnaptmp,
					"m_merit": mmertmp,
					"m_bid": mbidtmp,
					"stat": np.int32(tstat),
				}
			)

		self.tree 	= tree

	def run(self):
		#print(f"[Ctree Python reader] a path for the keyfile (default: './ctree_key.dat': {self.keyname}")
		#print(f"[Ctree Python reader] a path for the treefile (default: './ctree_tree.dat': {self.treename}")
		#print(f"[Ctree Python reader] a path for the output pickle (default: './ctree.pkl': {self.varname}")
		
		self.rdkey()
		self.rdtree()

		data = {"tree":self.tree, "key":self.key}

		fname 	= self.varname

		with open(fname, 'wb') as f:
			pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)

		return True



def prompt_path(prompt, default):
	"""Ask user for a path. If empty input -> default."""
	s = input(f"{prompt} [default: {default}] > ").strip()
	return s if s else default

def main(argv=None):

	reader 	= rdctree()

	DEFAULT_KEY 		= reader.keyname
	DEFAULT_TREE 		= reader.treename
	DEFAULT_OUT 		= reader.varname

	parser = argparse.ArgumentParser(
		description="Read CTree binary outputs (ctree_key.dat, ctree_tree.dat)."
	)
	parser.add_argument("--key", default=None, help=f"Path to key file (default: {DEFAULT_KEY})")
	parser.add_argument("--tree", default=None, help=f"Path to tree file (default: {DEFAULT_TREE})")
	parser.add_argument("--pickle", default=None, help=f"Path to output pickle file (default: {DEFAULT_OUT})")

	args = parser.parse_args(argv)

	
	reader.keyname = args.key if args.key is not None else prompt_path("Enter key filename", DEFAULT_KEY)
	reader.treename = args.tree if args.tree is not None else prompt_path("Enter tree filename", DEFAULT_TREE)
	reader.varname = args.pikcle if args.tree is not None else prompt_path("Enter output varname", DEFAULT_OUT)

	ok 	= reader.run()

	return 0 if ok else 1

if __name__ == "__main__":
	raise SystemExit(main())

