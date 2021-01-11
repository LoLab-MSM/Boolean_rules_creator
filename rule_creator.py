import itertools as it
import tempfile

import numpy as np
import json
import ast
import copy
import re
import random
import math
from sympy.logic.boolalg import to_dnf
from sympy.logic.boolalg import Xor
from sympy.parsing.sympy_parser import parse_expr
from time import time
import os

# ==================rule creator==================

def is_in_rule_list(elem, rulelist):
    for rule in rulelist:
        if elem in rule:
            return True
    return False


def rule_creator(worklist, active_SS, backward_paths, exclude=[]):
    n = len(active_SS)
    distance = [None]*len(worklist)
    # create a list from which the rules will be taken:
    # first entry contains the states, where x1 needs to be changed
    # second entry contains the states, where x2 needs to be changed...
    rulelist = [[] for _ in range(n)]

    # fill list of distances from active steady state
    d = 0
    for j in worklist:
        distance[d] = sum(abs(np.subtract(active_SS, j)))
        d = d+1

    # create a mask to sort all elements according to their distances to the main steady state
    mask = []
    for i in range(n):
        mask.append(np.where(np.array(distance) == i+1)[0])

    # fill the rulelist with the elements that are one flip away from the main steady state
    # here: mask[0] chooses the states that are one flip away from the main steady state
    old_element = active_SS

    for i in range(len(np.array(worklist)[mask[0]])):
        active_element = worklist[mask[0][i]]
        # compute the difference from the active element to the main steady state
        diff = tuple(abs(np.subtract(active_element, old_element)))
        # get index of difference - i.e., which node flip its state
        # this is done by determining the index of where the difference !=0, i.e. 1
        node_to_flip = diff.index(1)
        if (active_element, node_to_flip) not in exclude:
            rulelist[node_to_flip].append(active_element)

    # choose the active element from worklist
    # use elements that are one flip closer to the steady state
    # and one flip further away to include possible pathways that might only be available via a detour
    for k in range(n-1):
        wl_one_step_closer = [worklist[mask[k][i]] for i in range(len(np.array(worklist)[mask[k]]))]
        # determine if there exists a state one flip further away from the SS
        if len(mask) > k+2:
            wl_one_step_farther = [worklist[mask[k+2][i]] for i in range(len(np.array(worklist)[mask[k+2]]))]
        else:
            wl_one_step_farther = []
        for target_element in wl_one_step_closer+wl_one_step_farther:
            for j in range(len(np.array(worklist)[mask[k+1]])):
                second_element = np.array(worklist)[mask[k+1][j]]
                diff_flipstate = abs(np.subtract(second_element, target_element))
                # if the distance between the active one-flip-state and the two-flip-state is 1,
                # associate the two-flip-state with the one-flip state and eliminate from list
                if sum(diff_flipstate) == 1:
                    # state is assigned to the active state -> get index of flip
                    node_to_flip = tuple(diff_flipstate).index(1)
                    if (tuple(second_element), node_to_flip) not in exclude:
                        rulelist[node_to_flip].append(tuple(second_element))
                        if backward_paths == 1 and (tuple(target_element), node_to_flip) not in exclude:
                            rulelist[node_to_flip].append(tuple(target_element))

    return rulelist


# ==================file check and prep==================


# -------read the input file and put it in a sorted order---------------------
class IvInfo:
    def __init__(self, iv, frequency):
        self.iv = iv
        self.frequency = frequency


def read_steady_states(fn):
    with open(fn) as fs:
        data = json.load(fs)
    fs.close()

    print('\t')
    print("data: ", data)

    SS = {}
    for striv, sslist in data.items():
        iv = ast.literal_eval(striv)
        for ssfreq in sslist:
            ss    = repr(ssfreq[0])
            freq  = ssfreq[1]
            if ss not in SS:
                SS[ss] = []
            SS[ss].append(IvInfo(tuple(iv), freq))
    return data, SS


def print_SS(SS):
    for ss, ivilist in SS.items():
        print('\t')
        print('steady state: {}'.format(ss))
        for ivi in ivilist:
            print('\t{} <- {}'.format(ivi.iv, ivi.frequency))
        print('\t')


# check, if a state leads to a non-steady-state
def check_ss(SS):
    flag = 0
    ss_list = []
    for ss, ivilist in SS.items():
        ss_list.append(tuple(ast.literal_eval(ss)))

    for ss, ivilist in SS.items():
        for ivi in ivilist:
            if ivi.iv in ss_list and ivi.iv != tuple(ast.literal_eval(ss)):
                print("steady state ", ivi.iv, "is not a real steady state")
                flag = 1
    if flag == 1:
        print("the steady states in your input file do not make sense - please change your input!")
        print('\t')
        exit()
    elif flag == 0:
        print("your steady states are fine, please continue to reachability check")
        print('\t')


# -------check if every target can be reached by the full pathway----------------

def is_transition_valid(x, y, blacklist):
    node_to_flip = list(abs(np.subtract(x, y))).index(1)
    return (x, node_to_flip) not in blacklist


def get_accessible_states(x, lst, blacklist):
    # here we have to use blacklist so that we do not generate transitions that are forbidden
    return [y for y in lst if sum(abs(np.subtract(x, y))) == 1 and is_transition_valid(x, y, blacklist)]


def call_check_reachability(iv, ss, available_states, blacklist):
    active_list = [iv]
    print("check reachability for steady state ", ss, '\n')
    maxdepth = 2 ** (len(ss))
    for i in range(maxdepth):
        if min([sum(abs(np.subtract(x, ss))) for x in active_list]) <= 1:
            return True
        active_list = [get_accessible_states(x, available_states, blacklist) for x in active_list]
        active_list = [item for sublist in active_list for item in sublist] #flatten the list
        if len(active_list) == 0:
            return False
    return False


def apply_rule(state, state_list):
    pathway = []
    rule_idx = []

    for i in range(len(state_list)):
        if state in state_list[i]:
            rule_idx.append(i)

    # print(rule)
    for r in rule_idx:
        ns = list(state)
        if state[r] == 0:
            ns[r] = 1
            pathway.append(tuple(ns))
        elif state[r] == 1:
            ns[r] = 0
            pathway.append(tuple(ns))
    #print(state, pathway)
    return pathway


def check_validity_of_removal(initial_state, states_list, ss):
    maxdepth = 2 ** (len(ss))
    #state_list_flat = [item for sublist in states_list for item in sublist]  # flatten the list
    #state_list_flat = list(dict.fromkeys(state_list_flat))  # remove duplicates

    paths = [initial_state]
    for i in range(1, maxdepth):
        paths = [apply_rule(state, states_list) for state in paths]
        paths = [item for sublist in paths for item in sublist]  # flatten the list
        paths = list(dict.fromkeys(paths))  # remove duplicates
        #print(initial_state, paths, ss)
        if ss in paths:
            #print("ss in paths")
            return True
    return False


def get_freq(iv, SS, data):
    for str_iv, ssfreq in data.items():
        iv_list = tuple(ast.literal_eval(str_iv))
        if iv == iv_list:
            for freqlist in ssfreq:
                if tuple(freqlist[0]) == SS:
                    freq = freqlist[1]
    return freq


# -------create a hierarchy of objects to deal with targets (all states
# that should be able to reach), CombSS (objects that are capable of holding
# SS and targets with their corresponding initial values), and SSInfo (the class
# that holds the full list of all objects)

# Store specific ss and all additional states that eventually end up
# in that steady states
class Targets:
    def __init__(self, _ss, _targets):
        self.ss = _ss
        self.targets = _targets

    def check_reachability(self, iv, available_states, blacklist): #LE
        reachable = False
        for tgt in self.targets:
            reachable = reachable or call_check_reachability(iv, tgt, available_states, blacklist) #LE
        return reachable

    def eliminate_transitions(self, rules, ivs, freq, blacklist):
        print('eliminate_transitions calling')
        ct = 0
        maxrules = 6
        perc = math.floor(100/maxrules)
        keep_rules = 0
        freq_pos = 0
        new_rules = []
        check_rules = []

        for i in range(len(rules)):
            new_rules.append(list(set(rules[i])))

        dbg_ct=0
        for state in ivs:
            dbg_ct+=1
            #print('looking at state',state,dbg_ct,len(ivs))
            ct = sum([x.count(state) for x in new_rules])
            keep_rules = math.ceil(freq[freq_pos]/perc)#*0.4)
            #print(keep_rules, freq[freq_pos])
            #exit()
            if keep_rules > maxrules:
                keep_rules = maxrules
            if ct > keep_rules:
                i = 0
                rule_choice = list(range(len(symbols)))
                rule_choice = [m for m in range(len(symbols)) if state in new_rules[m]]
                #print(rule_choice)
                #exit()
                #for i in range(ct-keep_rules):
                while i < ct-keep_rules and len(rule_choice)>0:
                    check_rules = copy.deepcopy(new_rules)
                    choose_rule = random.choice(rule_choice)
                    #print(choose_rule)
                    #while state not in new_rules[choose_rule]:
                    #    rule_choice.remove(choose_rule)
                    #    choose_rule = random.choice(rule_choice)
                    check_rules[choose_rule].remove(state)
                    for tgt in self.targets:
                        #print("chose: ", state, choose_rule)
                        #print("chose: ", state, choose_rule, check_validity_of_removal(state, check_rules, tgt))
                        if check_validity_of_removal(state, check_rules, tgt) == True:
                            new_rules[choose_rule].remove(state)
                            i = i+1
                            #rule_choice.remove(choose_rule)
                            break
                    rule_choice.remove(choose_rule)

            freq_pos = freq_pos+1
        return new_rules


    def create_rules(self, available_states, frequencies, backwardpaths, blacklist):
        rule_list = [[] for i in symbols]
        #print(len(available_states), frequencies)
        #exit()
        a = time()
        for ss in self.targets:
            print('tgts: {}', len(self.targets), len(available_states), len(ss), backwardpaths, blacklist)
            rules = rule_creator(worklist=available_states, active_SS=ss, backward_paths=backwardpaths, exclude=blacklist)
            for k in range(len(symbols)):
                rule_list[k] += rules[k]
        b = time()
        rule_list = self.eliminate_transitions(rule_list, available_states, frequencies, blacklist)
        print('cri {} {}'.format(b-a,time()-b))
        return rule_list


class CombSS:
    def __init__(self, _sss, _ivs):
        self.ivs = _ivs
        self.targets = [Targets(ss, [ss]) for ss in _sss]

    def sss(self):
        return set([x.ss for x in self.targets])

    def is_strict_subset(self, css):
        return self.sss().issubset(css.sss()) and self.sss() != css.sss()

    def append_target(self, ss, ivs):
        for tgt in self.targets:
            if tgt.ss == ss:
                tgt.targets.extend(ivs)
                return
        print('ERRROR: ss={} could not be found!'.format(ss))

    def check_reachability(self, blacklist): #LE
        for iv in self.ivs:
            reachability = True
            for tgt in self.targets:
                reachability = reachability and tgt.check_reachability(iv, self.ivs, blacklist)
            if reachability == False:
                print('ERROR: iv={} could not reach ss={}'.format(iv, self.sss()))
                print('Please consider changing your input!')
                print('\t')
                exit()

    def create_rules(self, backwardpaths, blacklist, data):
        rule_list = [[] for i in symbols]
        for tgt in self.targets:
            frequencies = [get_freq(iv, tgt.ss, data) for iv in self.ivs]
            # print(frequencies)
            rules = tgt.create_rules(self.ivs, frequencies, backwardpaths, blacklist)
            for k in range(len(symbols)):
                rule_list[k] += rules[k]
        #exit()
        return rule_list


class SSInfo:
    def __init__(self, data):
        self.comb_steadystates = []
        self.data = data
        # create list of single and multiple steady states
        print('\t\t')
        print('============================')
        print('Start creating SSinfo object')
        print('\t')
        for str_iv, ssfreq in data.items():
            iv = tuple(ast.literal_eval(str_iv))
            sss = set([tuple(x[0]) for x in ssfreq])
            if iv in sss:
                continue
            self.add_or_create_entry(sss, iv)  # append entries to comb_steadystates

        # add targets for double and higher SS
        print('Populating the target entries.')
        for css in self.comb_steadystates:
            if len(css.sss()) >= 2:
                for css2 in self.comb_steadystates:
                    if css2.is_strict_subset(css):
                        print('{} is a subset of {}'.format(css2.sss(), css.sss()))
                        for tgts2 in css2.targets:
                            print('\tto target ss={} appending {}'.format(tgts2.ss, css2.ivs))
                            css.append_target(tgts2.ss, css2.ivs)

        print('Finished creating SSInfo object.')
        print('================================')

    #def print(self):  # print function that gives the info of the object
    #    print('\t\t')
    #    print('your SSInfo object created by your data:')
    #    for css in self.comb_steadystates:
    #        print(' ss={}   ivs={}'.format(css.sss(), css.ivs))
    #        for target in css.targets:
    #            print('\t ss={}    targets={}'.format(target.ss, target.targets))
    #    print('\t')

    def check_reachability(self, blacklist): #LE
        for css in self.comb_steadystates:
            css.check_reachability(blacklist) #LE

    def create_rules(self, backwardpaths, blacklist): #LE
        print(len(self.data),len(symbols))
        a = time()
        #random.seed(128)
        rule_list = [[] for i in symbols]
        for css in self.comb_steadystates:
            rules = css.create_rules(backwardpaths, blacklist, self.data)
            for k in range(len(symbols)):
                rule_list[k] += rules[k]
        print('rk took: {}'.format(time()-a))
        return rule_list

    def add_or_create_entry(self, sss, iv):
        idx = self.find(sss)
        if idx == -1:
            self.comb_steadystates.append(CombSS(sss, [iv]))
        else:
            self.comb_steadystates[idx].ivs.append(iv)

    def find(self, sss):
        idx = 0
        for x in self.comb_steadystates:
            if x.sss() == sss:
                return idx
            idx += 1
        return -1


# =============================rule manipulation: =====================
# fct to check if all rules have reached the desired length
def lst_or(a,b):
    result = False
    for i in range(len(a)):
        result = result or (a[i]>b[i])
        #print(a[i], b[i], a[i]>b[i], result)

    #print('result: ', result)
    return result

# start with the full backward-pathway system & eliminate states from
# given rule list
def rule_manipulation(bw_rules, percentages, simulations):
    n = len(bw_rules)
    rulefilelist = ['']*simulations

    for sim in range(simulations):

        bw_rulelist = copy.deepcopy(bw_rules)
        transition_states = []
        states_count = []
        rule_length = [0]*n
        #print(rule_length)

        for i in range(n):
            rule_length[i] = len(bw_rulelist[i])

        new_rule_length = [math.ceil(a*b) for a,b in zip(rule_length, percentages)]

        # keep track how many times the transition state is still left
        # it must not be totally eliminated from the system, otherwise
        # it would create an additional steady state
        for k in range(n):
            for i in bw_rulelist[k]:
                #print(i)
                if i not in transition_states:
                    transition_states.append(i)
                    states_count.append(1)

                else:
                    states_count[transition_states.index(i)] = states_count[transition_states.index(i)]+1

        # as long as not every rule is shortened to the according length, continue eliminating states
        while lst_or(rule_length, new_rule_length):
            choose_rule = random.choice(range(n))
            if len(bw_rulelist[choose_rule]) == 0:
                continue

            choose_transition = random.choice(bw_rulelist[choose_rule])
            idx = transition_states.index(choose_transition)
            if states_count[idx] > 1:
                bw_rulelist[choose_rule].remove(choose_transition)
                states_count[idx] = states_count[idx]-1
                rule_length[choose_rule] = rule_length[choose_rule]-1

        # translate rules to string
        rules = ['' for i in range(n)]
        for k in range(n):
            ruletext = []
            #rulelist[k] = list(set(rulelist[k]))
            for j in range(len(bw_rulelist[k])):
                ruletext.append(
                    ' & '.join(['{}{}'.format('' if bw_rulelist[k][j][i] == 1 else ' ~', symbols[i]) for i in range(n)]))
                ruletext[j] = "(" + ruletext[j] + ")"
            if ruletext != []:
                # rules[k] = 'Xor(' + ' | '.join(ruletext) + r',{})'.format(symbols[k])
                rules[k] = ' | '.join(ruletext)
                # print(len(rulelist[k]))
                rules[k] = "Xor((" + rules[k] + "), " + symbols[k] + ")"
                # print("rules[k] ", rules[k])
                rules[k] = parse_expr(rules[k])
                rules[k] = to_dnf(rules[k], True)
                rules[k] = str(rules[k]).replace('&', 'and').replace('|', 'or').replace('~', 'not ')
            else:
                rules[k] = symbols[k]
            rulefilelist[sim] = rulefilelist[sim]+'1: {}* = {}'.format(symbols[k], rules[k])+'\n'

    with open("rule_sets.json", "w") as fs:
        fs.write(json.dumps(rulefilelist, indent=2))
    fs.close()


symbols = []

# =================create rules with human guided input======================

def get_pattern(element, not_elim):
    pattern = []
    for i in not_elim:
        pattern.append(element[i])
    return  pattern


def slist(a, idx):
    return [a[i] for i in idx]


def self_elimination(rule_list, pos, not_elim, perc):
    result_list = []
    patterns = []
    tmp_lst = []
    tmp2_lst = []
    #print(rule_list)
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
        #print(tup_newtup)
        if tup_newtup not in rule_list:
            result_list.append(b)
        else:
            p = get_pattern(b, not_elim)
            if p not in patterns:
                patterns.append(p)
                #print(patterns)
                tmp_lst.append([a for a in rule_list if slist(a, not_elim) == p])


        #print(newtup)
    #print(tmp_lst)
    #exit()
    possible_combinations = len(tmp_lst)
    all_combinations = [a for a in it.product([0, 1], repeat=possible_combinations)]

    print("poss comb:", possible_combinations)
    print("perc: ", perc)
    acceptance_state = int(np.floor(len(all_combinations)*perc))
    print("acceptance state, out of: ", acceptance_state, len(all_combinations))
    print("acceptance state, state: ", acceptance_state, all_combinations[acceptance_state])

    tmp = []
    for j in range(len(tmp_lst)):
        tmp.append([a for a in tmp_lst[j] if a[pos] == all_combinations[acceptance_state][j]])
        #print(tmp)
    #print(tmp)
    #exit()
    tmp.append(result_list)
    #print(tmp)
    #exit()
    tmp2_lst.append([val for sublist in tmp for val in sublist])
    #print(tmp2_lst)
    #exit()
    return tmp2_lst[0]


def select_rules_human_guided(perc, _symbols):
    print("perc full: ", perc)
    for i in range(len(perc)):
        if perc[i]==1:
            perc[i]=0.99
    global symbols
    symbols = _symbols
    n = len(symbols)
    a = time()

    workfile1 = 'EMT_userguided_rule1.txt'
    workfile2 = 'EMT_userguided_rule2.txt'
    workfile3 = 'EMT_userguided_rule3.txt'
    workfile4 = 'EMT_userguided_rule4.txt'
    workfile5 = 'EMT_userguided_rule5.txt'
    workfile6 = 'EMT_userguided_rule6.txt'

    with open(workfile1, 'r') as f:
        read_data = f.read()

    rule1 = list(eval(read_data))

    with open(workfile2, 'r') as f:
        read_data = f.read()

    rule2 = list(eval(read_data))

    with open(workfile3, 'r') as f:
        read_data = f.read()

    rule3 = list(eval(read_data))

    with open(workfile4, 'r') as f:
        read_data = f.read()

    rule4 = list(eval(read_data))

    with open(workfile5, 'r') as f:
        read_data = f.read()

    rule5 = list(eval(read_data))

    with open(workfile6, 'r') as f:
        read_data = f.read()

    rule6 = list(eval(read_data))

    #rule = [rule1, rule2, rule3, rule4, rule5, rule6, [], []]
    #print(rule)
    #print("perc: ", perc)
    #exit()

    # NICD (0) depends on Notch (1), TP53 (2), TP63 (3)
    #print(rule1)
    #exit()
    rule_NICD = self_elimination(rule1, 0, [1, 2, 3], perc[0])
    #Sprint(rule_NICD)
    #exit()

    # Notch (1) depends on ECM (6), miRNA (4)
    rule_Notch = self_elimination(rule2, 1, [4, 6], perc[1])

    # TP53 (2) depends on DNAdam (7), NICD (0), miRNA (4), EMTreg (5), TP63_TP73 (3)
    rule_TP53 = self_elimination(rule3, 2, [0, 3, 4, 5, 7], perc[2])

    # TP63 (3) depends on DNAdam (7), miRNA (4), NICD (0), TP53 (2)
    rule_TP63 = self_elimination(rule4, 3, [0, 2, 4, 7], perc[3])

    # miRNA (4) depends on  TP53 (2), TP63_TP73 (3), EMTreg (5)
    rule_miRNA = self_elimination(rule5, 4, [2, 3, 5], perc[4])

    # EMTreg (5) depends on  NICD (0), miRNA (4)
    rule_EMTreg = self_elimination(rule6, 5, [0, 4], perc[5])

    simple_rulelist = [rule_NICD, rule_Notch, rule_TP53, rule_TP63, rule_miRNA, rule_EMTreg, [], []]

    rules = ['' for i in range(len(symbols))]
    #print(simple_rulelist)
    #exit()

    #f = open("rules.txt", "w+")

    #for k in range(n):
    #    f.write(symbols[k]+ ' = False\n')


    # Sympy
    str_rules = ""
    str_rules_cpp = ""

    for k in range(n):
        ruletext = []
        #print(simple_rulelist)
        #exit()
        for j in range(len(simple_rulelist[k])):
            ruletext.append(' & '.join(['{}{}'.format('' if simple_rulelist[k][j][i] == 1 else ' ~', symbols[i]) for i in range(n)]))
            ruletext[j] = "(" + ruletext[j] + ")"
        #print(ruletext)
        #exit()
        if ruletext!=[]:
            #rules[k] = 'Xor(' + ' | '.join(ruletext) + r',{})'.format(symbols[k])
            rules[k] = ' | '.join(ruletext)
            #print(len(rulelist[k]))
            rules[k] = "Xor((" + rules[k] + "), " + symbols[k] + ")"
            #print(rules[k])
            #exit()

            # C++ output
            str_rules_cpp += '    if(k=={})\n'.format(k)
            str_rule_cpp = str(rules[k]).replace('&','&&').replace('|','||').replace('~','!')
            for sym, sym_idx in zip(symbols,range(len(symbols))):
                str_rule_cpp = str_rule_cpp.replace(sym,"x[{}]".format(sym_idx))
            str_rules_cpp += '        return {};\n'.format(str_rule_cpp)

            # print("rules[k] ", rules[k])
            #rules[k] = parse_expr(rules[k])
            #rules[k] = to_dnf(rules[k], True)
            rules[k] = str(rules[k]).replace('&', 'and').replace('|', 'or').replace('~', 'not ')
        else:
            rules[k] = symbols[k]
            str_rules_cpp += '    if(k=={0})\n        return x[{0}];\n'.format(k)
        #print('1: {}* = {}'.format(symbols[k], rules[k]))
        #exit()
        str_rule = '1: {}* = {}\n'.format(symbols[k], rules[k])
        #f.write(str_rule)
        str_rules += str_rule
    #f.close()

    #print(str_rules_cpp)
    with open("c_simulator.cpp_template") as fs_cpp:
        str_cpp = fs_cpp.read()

    fs_cpp = tempfile.NamedTemporaryFile(mode='w', delete=False)
    #with open("/tmp/c_simulator-{}.cpp", "w") as fs_cpp:
    fs_cpp.write(str_cpp.replace('{0}',str(len(symbols))).replace('{1}',str_rules_cpp))
    fs_cpp.flush()
    os.fsync(fs_cpp.fileno())
    fs_cpp.close()
    print("generating rules took ", time()-a)

    return str_rules, simple_rulelist, fs_cpp.name



# =================create rules with read-in data======================

def creating_rules(fn, _symbols):
    global symbols
    symbols = _symbols
    a = time()

    [data, SS] = read_steady_states(fn)
    #[data, SS] = read_steady_states("SCLC_10_nodes.json")
    #[data, SS] = read_steady_states("EMT_incbw.json")
    #[data, SS] = read_steady_states("EMT_nobw.json")

    print_SS(SS)
    #check_ss(SS)

    ss_info = SSInfo(data)
    #ss_info.print()
    #ss_info.check_reachability(blacklist=[((1, 1, 1), 0)])
    #ss_info.check_reachability(blacklist=[])

    #print(ss_info.create_rules())
    #print(ss_info.create_rules())
    #rulelist = ss_info.create_rules(backwardpaths=1, blacklist=[((1, 1, 1), 0)])
    print('create rules \n')
    rulelist = ss_info.create_rules(backwardpaths=1, blacklist=[])

    n = len(symbols)
    blacklisted  = [(x,6) for x in it.product(range(2), repeat=n)]
    blacklisted += [(x,7) for x in it.product(range(2), repeat=n)]
    #print(blacklisted)

    simple_rulelist = [[] for _ in range(n)]
    for i in range(len(simple_rulelist)):
        simple_rulelist[i] = list(dict.fromkeys(rulelist[i]))

    print(simple_rulelist)

    for i in range(len(simple_rulelist)):
        print("rule ", i, len(simple_rulelist[i]))

    # add expert knowledge pathways
    #rulelist[0].append((1, 1, 1, 0))
    #rulelist[2].append((1, 1, 0, 0))
    #print(rulelist)
    #print(rulelist[0][0])
    #print(rulelist[0][0][3])


    # ============== create multiple rulesets =============

    #print(simple_rulelist)
    #rule_manipulation(simple_rulelist, [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], 10)


    #exit()

    # ============== reorder/simplify rules ===============

    #n = len(active_SS)
    rules = ['' for i in range(n)]

    #f = open("rules.txt", "w+")

    #for k in range(n):
    #    f.write(symbols[k]+ ' = False\n')

    # Sympy
    str_rules = ""
    str_rules_cpp = ""
    for k in range(n):
        ruletext = []
        rulelist[k] = simple_rulelist[k]
        for j in range(len(rulelist[k])):
            ruletext.append(' & '.join(['{}{}'.format('' if rulelist[k][j][i] == 1 else ' ~', symbols[i]) for i in range(n)]))
            ruletext[j] = "(" + ruletext[j] + ")"
        if ruletext!=[]:
            #rules[k] = 'Xor(' + ' | '.join(ruletext) + r',{})'.format(symbols[k])
            rules[k] = ' | '.join(ruletext)
            #print(len(rulelist[k]))
            rules[k] = "Xor((" + rules[k] + "), " + symbols[k] + ")"

            # C++ output
            str_rules_cpp += '    if(k=={})\n'.format(k)
            str_rule_cpp = str(rules[k]).replace('&','&&').replace('|','||').replace('~','!')
            for sym, sym_idx in zip(symbols,range(len(symbols))):
                str_rule_cpp = str_rule_cpp.replace(sym,"x[{}]".format(sym_idx))
            str_rules_cpp += '        return {};\n'.format(str_rule_cpp)

            # print("rules[k] ", rules[k])
            #rules[k] = parse_expr(rules[k])
            #rules[k] = to_dnf(rules[k], True)
            rules[k] = str(rules[k]).replace('&', 'and').replace('|', 'or').replace('~', 'not ')
        else:
            rules[k] = symbols[k]
            str_rules_cpp += '    if(k=={0})\n        return x[{0}];\n'.format(k)
        print('1: {}* = {}'.format(symbols[k], rules[k]))
        str_rule = '1: {}* = {}\n'.format(symbols[k], rules[k])
        #f.write(str_rule)
        str_rules += str_rule
    #f.close()

    print(str_rules_cpp)
    with open("c_simulator.cpp_template") as fs_cpp:
        str_cpp = fs_cpp.read()

    fs_cpp = tempfile.NamedTemporaryFile(mode='w', delete=False)
    #with open("/tmp/c_simulator-{}.cpp", "w") as fs_cpp:
    fs_cpp.write(str_cpp.replace('{0}',str(len(symbols))).replace('{1}',str_rules_cpp))
    fs_cpp.flush()
    os.fsync(fs_cpp.fileno())
    fs_cpp.close()
    print("generating rules took ", time()-a)

    return str_rules, simple_rulelist, fs_cpp.name


if __name__ == '__main__':
    creating_rules('modified_freq.json', ["x1","x2","x3","x4"])

    #print(ruletext)

