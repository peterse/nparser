#scratch pytho
import sys, os
#print sys.path
import ROOT as R
import dctROOTv9 as dR
from collections import OrderedDict
import math
import numpy as np
from decimal import Decimal



class N_event:
	"""Characterizes a given set of events for a specified root file"""
	#Rotating frame -3deg => rotating all vectors +3deg
	correction_angle = np.radians(3) #degrees
	traits_list = ("Name", "ID", "4vec", "reconMass", "KE")
	


	def __init__(self, filename, export=False, printing=True):
		current_tree = None
		self.sFT = dR.fileTools(filename)
		#Exporting
		self.export = export
		self.export_case = "angle"
		#Printing
		self.printing = printing

		#Useful parameters and attributes
		self.rot_matrix = self.rot2D_matrix(self.correction_angle)

		#Get the first tree in the dct
		for tree_name, handl in self.sFT.tree_handle_dct.iteritems():
			self.treeinfo = (tree_name, handl)
			self.tot_evt = handl.GetEntries()
			break
		
		#! ! ! ! Initialize Framework (hists) here ! ! ! ! 


		#no, I don't want a generator go away Python 3
		self.rnge = range(0, 20)
		#dictionary of values to be exported
		#Tuple of (list, running index)
		self.export_dct = {"proton": [[0]*4*len(self.rnge), 0], "neutron": [[0]*4*len(self.rnge), 0], "nu_px": [[0]*len(self.rnge), 0], "nu_py": [[0]*len(self.rnge), 0]}

		
		for i in self.rnge:

			(my_parts, summ) = self.event_splitter(i)
			if self.printing:
				print "\n\n"
				dR.dctTools(my_parts).printer()
				print "\n"
				dR.dctTools(summ).printer()


		if self.export:		
		#Cleanup lists in export_dct
			for k, lst in self.export_dct.iteritems():
				lst[0] = filter(lambda x: x!=0, lst[0])

			if self.export_case == "angle":
				out_lst = [(self.export_dct["neutron"][0], "n_trans_angle"), (self.export_dct["proton"][0], "p_trans_angle")]
				self.write_hist_angle("nparser_angles", out_lst)
			elif self.export_case == "momentum":
				out_lst = [(self.export_dct["nu_px"][0], "nu_x_momentum"), (self.export_dct["nu_py"][0], "nu_y_momentum")]
				self.write_hist_momentum("nparser_nu_momenta", out_lst)


	def template_odct(self):
		dct = OrderedDict({})
		for trait in self.traits_list:
			dct[trait] = None
		return dct

	def event_splitter(self, event):
	#Splits an event into a dct of particles and their 4vecs IN DETECTOR COORDS
		
				

		self.treeinfo[1].GetEvent(event)
		#Base number of particles on nFSPart
		leaf_obj = self.treeinfo[1].GetLeaf("mc_nFSPart")
		N_parts = int(leaf_obj.GetValue())
		print "%i particles found in event %i" %(N_parts, event)

		#Initialize dictionary of particles in interaction
		part_dct = OrderedDict({})
		for n_part in range(N_parts):
			part_dct[n_part] = self.template_odct()
				
		#Get 4vectors and populate
		vecs = self.__get_4vecs(event, N_parts)
		for part_num, vec in enumerate(vecs):
			part_dct[part_num]["4vec"] = vec
	
		#Get names, populate
		for part_num, nth_dct in part_dct.iteritems():
			(ID, name) = self.__get_name(event, part_num)
			nth_dct["ID"] = ID
			nth_dct["Name"] = name

			#Remember where we put the lepton
			#! ! ! ! GENERALIZE !! ! ! !
			if name == "antimuon":
				self.i_big_lep = part_num

			#Get mass, correct 4vec
			nth_dct["reconMass"] = self.vec_mass(nth_dct["4vec"])
			#Correct the yz coordinates of 4-vectors from detector -> beamline coords
			nth_dct["4vec"] = self.yz_rotation(nth_dct["4vec"])
			#mass = dR.PDGTools.fetch_mass(nth_dct["ID"])
			nth_dct["KE"] = nth_dct["4vec"][0] - nth_dct["reconMass"]


		my_lepton = self.lepton_check(part_dct)

		#Get particle group info by summing up total momentum before adding neutrino
		P_out = self.get_P(part_dct)
		nu_dct = self.get_neutrino(part_dct)

		#Get the summary of the reaction
		summary = self.event_summary(event, part_dct, nu_dct)
		#Add the neutrino now
		part_dct[N_parts] = nu_dct

		#Parse over each neutron and find angles
		for part_num, nth_dct in part_dct.iteritems():
			if nth_dct["Name"] in ["neutron", "proton"]:
				name = nth_dct["Name"]
				v_b = nth_dct["4vec"]
				v_mu = part_dct.get(self.i_big_lep).get("4vec")
				v_nu = part_dct.get(N_parts).get("4vec")
				#angle between batryon and the mu-nu plane in degrees
				nth_dct["mu-nu_angle"] = self.trans_plane_angle(v_b, v_mu, v_nu)
				#mu-baryon angle in xy-plane
				azi = self.azi_angle(v_b, v_mu)
				nth_dct["azi_angle"] = azi
				
				#Export angles to hist; these are NOT one-to-one with the event
				if self.export and self.export_case == "angle":
					if name=="neutron" and summary["out"]["neutron"] == 1 and summary["out"]["proton"] == 0:
						self.export2lst( azi, self.export_dct[name])



		return part_dct, summary


	def get_neutrino(self, part_dct):
	#Get nu info by combining mc results and info on outgoing particles
		nu_dct = self.template_odct()

		#nu identity
		ID = self.treeinfo[1].GetLeaf("mc_incoming").GetValue()
		nu_dct["ID"] = ID
		nu_dct["Name"] = dR.PDGTools.decode_ID(ID)
		#ASSUME NU CARRIES ALL THE P_OUT TO THE RXN
		tot_P =  self.get_P(part_dct)
		E_nu = self.treeinfo[1].GetLeaf("mc_incomingE").GetValue()
		m_p = dR.PDGTools.fetch_mass("proton")
		#N will be the number of at-rest protons expected to have reacted with nu_in
		N_guess = E_nu - tot_P[0]
		#4-vector constructed from E_nu, p_out
		nu_dct["4vec"] = (E_nu, tot_P[1], tot_P[2], tot_P[3])
		nu_dct["reconMass"] = self.vec_mass(nu_dct["4vec"])

		if self.export and self.export_case == "momentum":
			self.export2lst(tot_P[1], self.export_dct["nu_px"])
			self.export2lst(tot_P[2], self.export_dct["nu_py"])

		return nu_dct


	def __get_name(self, event, index, incoming=False):
	#Get the name of a particle based on FS PDG or incoming info
		if incoming:
			leaf_obj = self.treeinfo[1].GetLeaf("mc_incoming")
			evt_id = 0
		else:
			leaf_obj = self.treeinfo[1].GetLeaf("mc_FSPartPDG")
			evt_id = leaf_obj.GetValue(index)
		return int(evt_id), dR.PDGTools.decode_ID(evt_id) #(ID, name)
	

	def __get_4vecs(self, event, N_parts):
	#Return a list of tuples whose position corresponds to part_num

		#NOTE: [[0,0,0,0]]*N_parts is a list of N-parts copies of [0,0,0,0], the correct format is below:
		my_4vecs = [ [0,0,0,0] for i in range(N_parts)]
		vals = [0]*N_parts
		#populate Momenta
		for dim_i, dim in enumerate(["x", "y", "z"]):
			dimleaf = "mc_FSPartP" + dim
			P_obj = self.treeinfo[1].GetLeaf(dimleaf)
			#print my_4vecs
			for part_num in range(len(my_4vecs)):
				my_4vecs[part_num][dim_i+1] = P_obj.GetValue(part_num)
	
	
		#populate Energies
		E_obj = self.treeinfo[1].GetLeaf("mc_FSPartE")
		for part_num in range(N_parts):
			my_4vecs[part_num][0] = E_obj.GetValue(part_num)
		#print my_4vecs
		my_4vecs = [tuple(lst) for lst in my_4vecs]
		return my_4vecs #(E, px, py, pz)


	def lepton_check(self, dct):
	# ! ! ! IN PROGRESS ! ! ! 
	#Find the identity of the lepton in an event, and throw errors if many present
	#Takes the part dict (pre-nu) as a parameter
		lepton_check = set([(num, nth_dct["Name"]) for num, nth_dct in dct.iteritems()])
		#! ! ! ! Setting up uniqueness stuff for deciding incoming neutrino ! ! ! 
		return None
		
	def get_P(self, part_dct):
		#GET ENERGY OF INC NEUTRINO
 		#INCLUDE PROTON INTERXN IN THE INITIAL NEUTRINO SCATTER
		tot_P = [0,0,0,0]

		for part_num, nth_dct in part_dct.iteritems():
			tot_P = [sum(x) for x in zip(nth_dct["4vec"], tot_P)]

		return tuple(tot_P)

	def event_summary(self, event, dct, nu_dct):
	#Summarizes the full event 
		summ_dct = OrderedDict({"event": event})

		all_parts = [nth_dct["Name"] for __, nth_dct in dct.iteritems()]
		
		#Outgoing, ingoing baryons
		summ_dct["out"] = {}
		summ_dct["in"] = {}
		for part in ["proton", "neutron"]:
			summ_dct["out"][part] = all_parts.count(part)


		#Particles involved in initial reaction are deterministic
		summ_dct["in"]["proton"] = summ_dct.get("out").get("proton") + 1
		summ_dct["in"]["neutron"] = summ_dct.get("out").get("neutron") - 1

		#Total outgoing momentum, mass
		summ_dct["out"]["P_out"] = self.get_P(dct)
		summ_dct["out"]["mass_out"] = sum([nth_dct["reconMass"] for __, nth_dct in dct.iteritems()]).real
		summ_dct["out"]["KE_out"] = sum([nth_dct["KE"] for __, nth_dct in dct.iteritems()]).real
		summ_dct["Z"] = self.treeinfo[1].GetLeaf("mc_targetZ").GetValue()
		#! ! ! ! ! ! ! 
		#Q = (KE_out + m_mu) - (E_in)
		summ_dct["reconQ"] = summ_dct["out"]["KE_out"] + dct.get(self.i_big_lep).get("reconMass") - nu_dct["4vec"][0] 
		summ_dct["mc_Q2"] = self.treeinfo[1].GetLeaf("mc_Q2").GetValue()
		return summ_dct
	
	def trans_plane_angle(self, v_b, v_mu, v_nu):
	#Given the 4-vectors of nu and mu (v_nu, v_mu), find the angle (deg) of neutron out of this plane

		
		#print "input vectors: ", v_b (baryon), v_mu, v_nu
		v_nu = v_nu[1:]
		v_mu = v_mu[1:]
		
		#set up vector normal to mu - nu plane
		v_plane = np.cross(v_nu, v_mu)
		mag = sum(map(lambda x: x*x, v_plane)) ** .5
		v_plane = v_plane / mag

		#Find angle between the vectors (which is the complement of the angle we're looking for)
		v_b = v_b[1:]
		nmag = sum(map(lambda x: x*x, v_b)) ** .5
		v_b = np.array(v_b) / nmag
		not_phi = np.dot(v_plane, v_b)

		#print mag, nmag
		#print v_plane,"\n", v_b

		return 90 - np.arccos(not_phi)

	def azi_angle(self, v_b, v_mu):
	#Passed a baryon and muon 4vec, return the angle (deg) between them in the transverse plane
		azi_b = np.degrees(np.arctan2(v_b[2], v_b[1]))			
		azi_mu = np.degrees(np.arctan2(v_mu[2], v_mu[1]))	
		diff = azi_b - azi_mu
		#Range: 0-360
		if diff < 0:
			diff = diff + 360
		return diff


	#Exporting methods
	def export2lst(self, source, target_dct):
	#send source to target[i_export] - inst var
		try:
			target_dct[0][target_dct[1]] = source
		except IndexError:
			#Didn't initialize large enough
			target_dct[0].append(source)
		target_dct[1] += 1
	
	def write_hist_angle(self, title, tupl_lst):
	#Pass this fxn a list of tuples [(lst, "name"), (lst2, "name")...]
		hout = R.TFile(title+".root", "RECREATE")
		for (lst, name) in tupl_lst:
			h1 = R.TH1F(name, name, 100, 100, 260)
			for val in lst:
				h1.Fill(val)
			h1.Write()

	def write_hist_momentum(self, title, tupl_lst):
	#see write_hist_angle
		hout = R.TFile(title+".root", "RECREATE")
		for (lst, name) in tupl_lst:
			h1 = R.TH1F(name, name, 100, -500, 500)
			for val in lst:
				h1.Fill(val)
			h1.Write()

	@staticmethod
	def rot2D_matrix(theta):
	#Initialize a rotation matrix for a given theta, in radians
		return ( (math.cos(theta), -math.sin(theta)), (math.sin(theta), math.cos(theta)) )


	@staticmethod
	def vec_mass(tupl):
	#Return the mass of a 4vec (E, px, py, pz)
		sqrs = map(lambda x: x*x, tupl)
		try:
			return (sqrs[0] - sum(sqrs[1:]))**.5
		#Return imaginary number if norm(E) < norm(p)
		except ValueError:
			return complex(0,abs(sqrs[0] - sum(sqrs[1:]))**.5)

	def yz_rotation(self, p_vec):

		zpyp = (p_vec[3], p_vec[2]) #(z, y)
		#code golf!
		#cor_zy2 = tuple([sum([a*b for a,b in zip(zpyp, i)]) for i in self.rot_matrix])
		#Explicit switch
		cor_zy = [0, 0]
		for i, rot_tup in enumerate(self.rot_matrix):
			#Matrix multiplication / linear alg
			cor_zy[i] = sum([a*b for a,b in zip(zpyp, rot_tup)])
		cor_zy = tuple(cor_zy) #(z', y')
		

		return (p_vec[0], p_vec[1], cor_zy[1], cor_zy[0])


if __name__ == "__main__":
	N_event("smallmc.root", export=False, printing=True)



