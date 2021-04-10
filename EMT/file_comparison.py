import json
import ast
import math
import itertools as it


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
            #print('ss: ',ss)
            freq  = ssfreq[1]
            #print('freq: ', freq)
            if ss not in SS:
                SS[ss] = []
            SS[ss].append(IvInfo(tuple(iv), freq))
    return data, SS


def comparator(freq_multiple_ss, fn_comparison):
    with open(fn_comparison) as f2:
        orig_data = json.load(f2)
    f2.close()

    modified_data = {str(x[0]) : x[1] for x in freq_multiple_ss}

    N=0
    diff = 0

    for key in modified_data.keys():
        if len(modified_data[key]) > 1:
            #print('new IC')
            for e_orig in orig_data[key]:
                ss_orig   = e_orig[0]
                freq_orig = e_orig[1]
                # find in list
                for e_mod in modified_data[key]:
                    ss_mod   = e_mod[0]
                    freq_mod = float(e_mod[1])#*0.1
                    if ss_mod == ss_orig:
                        break
                N = N + 1
                diff = diff + abs(freq_orig - freq_mod)**2

    RMS = math.sqrt(diff/float(N))
    #print('RMS: ', RMS, diff, N)
    return RMS #diff #RMS


def map_freqs(fn):

    with open(fn) as f2:
        data = json.load(f2)
    f2.close()

    N=0
    idx=[]
    all_idx = 0

    n = len(ast.literal_eval(next(iter(data))))
    #print(n)

    key_list = list(it.product([0, 1], repeat=n))
    key_list = [list(x) for x in key_list]

    for ordered_idx in key_list:
        idx.append((all_idx, len(data[str(ordered_idx)])))
        if len(data[str(ordered_idx)]) > 1:
            N=N+len(data[str(ordered_idx)])
        all_idx = all_idx+1

    return N, idx


def comparator_single_test(sim_data, fn_comparison):
    with open(fn_comparison) as f2:
        orig_data = json.load(f2)
    f2.close()

    with open(sim_data) as fs:
        freq_multiple_ss = json.load(fs)
    fs.close()

    modified_data = {str(x[0]) : x[1] for x in freq_multiple_ss}

    N=0
    diff = 0

    for key in modified_data.keys():
        if len(modified_data[key]) > 1:
            for e_orig in orig_data[key]:
                ss_orig   = e_orig[0]
                freq_orig = e_orig[1]
                for e_mod in modified_data[key]:
                    ss_mod   = e_mod[0]
                    freq_mod = float(e_mod[1])
                    if sum([x1 - x2 for (x1, x2) in zip(ss_orig, ss_mod)]) == 0:
                        if freq_orig - freq_mod > 5:
                            if freq_orig > 50 and freq_mod < 50 or freq_orig < 50 and freq_mod > 50:
                                print("IC: ", key, "ss", ss_orig, "orig freq", freq_orig, "mod freq", freq_mod, "!! reverse states")
                            else:
                                print("IC: ", key, "ss", ss_orig, "orig freq", freq_orig, "mod freq", freq_mod)

                N = N + 1
                diff = diff + abs(freq_orig - freq_mod)**4
            print('new IC \n')

    RMS = math.pow(diff/float(N),0.25)
    return RMS


if __name__ == '__main__':
    comparator_single_test('simulation_cpp_testing.json', 'EMT_paper.json')

