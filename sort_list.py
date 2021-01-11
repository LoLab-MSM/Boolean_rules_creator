import numpy as np
import random
import itertools

def check_and_swap_bit(rule_list, pos):
	available_states = rule_list
	result_list = []
	for b in rule_list:
		newtup = []
		for bit in range(8):
			if bit == pos and b[bit] == 0:
				newtup.append(1)
			elif bit == pos and b[bit] == 1:
				newtup.append(0)
			else:
				newtup.append(b[bit])
		tup_newtup = tuple(newtup)
		if tup_newtup in rule_list:
			result_list.append(b)
	return result_list

def get_pattern(element, not_elim):
	pattern = []
	for i in not_elim:
		pattern.append(element[i])
	return  pattern

def sublist(a, idx):
    return [a[i] for i in idx]

def check_and_swap_bit_same(rule_list, pos, not_elim):
	result_list = []
	patterns = []
	tmp_lst = []
	tmp2_lst = []
	result_lists = []
	for b in rule_list:
		newtup = []
		
		for bit in range(8):
			if bit == pos and b[bit] == 0:
				newtup.append(1)
			elif bit == pos and b[bit] == 1:
				newtup.append(0)
			else:
				newtup.append(b[bit])
		tup_newtup = tuple(newtup)
		if tup_newtup not in rule_list:   
			result_list.append(b)
		else:
			p = get_pattern(b, not_elim)
			if p not in patterns:
				patterns.append(p)
				tmp_lst.append([a for a in rule_list if sublist(a, not_elim) == p])
	

	possible_combinations = len(tmp_lst)	
	all_combinations = [a for a in itertools.product([0, 1], repeat=possible_combinations)]	

	#======= asymmetric: eliminate random
	tmp=[]
	coin = random.randint(0, 1)
	for i in range(len(tmp_lst)):
		tmp.append([a for a in tmp_lst[i] if a[pos] == coin]) 
	tmp.append(result_list)
	tmp2_lst.append([val for sublist in tmp for val in sublist])


	return tmp2_lst
	
#the elimination of transitions starts with the full backward model transition lists	
for rule_nr in [1,2,3,4,5,6]: 
	if rule_nr == 1:
		workfile = 'EMT_incbw_rule1.txt'
	elif rule_nr ==2:
		workfile = 'EMT_incbw_rule2.txt'
	elif rule_nr ==3:
		workfile = 'EMT_incbw_rule3.txt'
	elif rule_nr ==4:
		workfile = 'EMT_incbw_rule4.txt'
	elif rule_nr ==5:
		workfile = 'EMT_incbw_rule5.txt'
	elif rule_nr ==6:
		workfile = 'EMT_incbw_rule6.txt'

	with open(workfile, 'r') as f:
		read_data = f.read()

	rule = eval(read_data)

	list_NICD = []
	list_Notch = []
	list_TP53 = []
	list_TP63 = []
	list_miRNA = []
	list_EMTreg = []
	list_ECM = []
	list_DNAdam = []

			
	if rule_nr == 1:
		# NICD depends on Notch (1), TP53 (2), TP63 (3), 
		# no NICD (0), miRNA (4), EMTreg (5), ECM (6), DNAdam(7)
		list_miRNA = check_and_swap_bit(rule, 4)
		list_EMTreg = check_and_swap_bit(list_miRNA, 5)
		list_ECM = check_and_swap_bit(list_EMTreg, 6)
		list_DNAdam = check_and_swap_bit(list_ECM, 7)

		new_set = list_DNAdam
		
	elif rule_nr == 2:
		# Notch depends on ECM (6), miRNA (4), 
		# no NICD (0), Notch(1), TP53 (2), TP63 (3), EMTreg(5), DNAdam(7)
		list_NICD = check_and_swap_bit(rule, 0)
		list_TP53 = check_and_swap_bit(list_NICD, 2)
		list_TP63 = check_and_swap_bit(list_TP53, 3)
		list_EMTreg = check_and_swap_bit(list_TP63, 5)
		list_DNAdam = check_and_swap_bit(list_EMTreg, 7)

		new_set = list_DNAdam

	elif rule_nr == 3:
		# TP53 depends on DNAdam (7), NICD (0), miRNA (4), EMTreg (5), TP63_TP73 (3)
		# no Notch (1), TP53 (2), ECM (6)
		list_Notch = check_and_swap_bit(rule, 1)
		list_ECM = check_and_swap_bit(list_Notch, 6)

		new_set = list_ECM

	elif rule_nr == 4:
		# TP63 depends on DNAdam (7), miRNA (4), NICD (0), TP53 (2)
		# no Notch (1), TP63 (3), EMTreg(5), ECM (6)
		list_Notch = check_and_swap_bit(rule, 1)
		list_EMTreg = check_and_swap_bit(list_Notch, 5)
		list_ECM = check_and_swap_bit(list_EMTreg, 6)

		new_set = list_ECM

	elif rule_nr == 5:
		# miRNA depends on  TP53 (2), TP63_TP73 (3), EMTreg (5)
		# no NICD (0), Notch (1), miRNA (4), ECM (6), DNAdam (7)
		list_NICD = check_and_swap_bit(rule, 0)
		list_Notch = check_and_swap_bit(list_NICD, 1)
		list_ECM = check_and_swap_bit(list_Notch, 6)
		list_DNAdam = check_and_swap_bit(list_ECM, 7)

		new_set = list_DNAdam		

	elif rule_nr == 6:
		# EMTreg depends on  NICD (0), miRNA (4)
		# no Notch (1), TP53 (2), TP63 (3), EMTreg (5), ECM (6), DNAdam (7)
		list_Notch = check_and_swap_bit(rule, 1)
		list_TP53 = check_and_swap_bit(list_Notch, 2)
		list_TP63 = check_and_swap_bit(list_TP53, 3)
		list_ECM = check_and_swap_bit(list_TP63, 6)
		list_DNAdam = check_and_swap_bit(list_ECM, 7)

		new_set = list_DNAdam

	with open('EMT_userguided_rule'+str(rule_nr)+'.txt', 'w') as f:
    		for item in new_set:
    	    		f.write(str(item)+", ")
