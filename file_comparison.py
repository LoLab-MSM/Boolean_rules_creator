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

#d = read_steady_states('testing.json')
"""
def comparator(fn):
    with open(fn) as fs:
        modified_data = json.load(fs)
    fs.close()

    with open('original_data.json') as f2:
        orig_data = json.load(f2)
    f2.close()

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
                #print('original file, IC: ', key, 'reaches steady state ', ss_orig, 'freq: ', freq_orig)
                #print('original file, IC: ', key, 'reaches steady state ', ss_mod, 'freq: ', freq_mod)
                #print('diff: ', abs(freq_orig - freq_mod))


    RMS = math.sqrt(diff/float(N))
    #print('RMS: ', RMS, diff, N)
    return RMS
"""

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
                #diff = max(abs(freq_orig - freq_mod), diff)
                #print('original file, IC: ', key, 'reaches steady state ', ss_orig, 'freq: ', freq_orig)
                #print('original file, IC: ', key, 'reaches steady state ', ss_mod, 'freq: ', freq_mod)
                #print('diff: ', abs(freq_orig - freq_mod))


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

    #n = int(math.sqrt(len(data)))
    n = len(ast.literal_eval(next(iter(data))))
    #print(n)

    key_list = list(it.product([0, 1], repeat=n))
    key_list = [list(x) for x in key_list]
    #print(key_list[0])

    #print(orig_data[str(list(key_list[0]))])

    for ordered_idx in key_list:
        #print(ordered_idx)
        #print(data[str(ordered_idx)])
        idx.append((all_idx, len(data[str(ordered_idx)])))
        if len(data[str(ordered_idx)]) > 1:
            N=N+len(data[str(ordered_idx)])
        all_idx = all_idx+1
        #    print(orig_data[str(list(key_list[ordered_idx]))])

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
            #print('new IC')
            for e_orig in orig_data[key]:
                #print('IC', key)
                ss_orig   = e_orig[0]
                freq_orig = e_orig[1]
                #print('orig data: ', ss_orig, freq_orig)
                # find in list
                for e_mod in modified_data[key]:
                    ss_mod   = e_mod[0]
                    freq_mod = float(e_mod[1])#*0.1
                    #print('sim data: ', ss_mod, freq_mod)
                    if sum([x1 - x2 for (x1, x2) in zip(ss_orig, ss_mod)]) == 0:
                        if freq_orig - freq_mod > 5:
                            if freq_orig > 50 and freq_mod < 50 or freq_orig < 50 and freq_mod > 50:
                                print("IC: ", key, "ss", ss_orig, "orig freq", freq_orig, "mod freq", freq_mod, "!! reverse states")
                            else:
                                print("IC: ", key, "ss", ss_orig, "orig freq", freq_orig, "mod freq", freq_mod)
                            #print("IC: ", key, "orig ss", ss_orig, "freq ss", ss_mod, "orig freq", freq_orig, "mod freq", freq_mod)
                    #if ss_mod == ss_orig:
                    #    break

                N = N + 1
                #diff = diff + abs(freq_orig - freq_mod)**2
                diff = diff + abs(freq_orig - freq_mod)**4
                #print('original file, IC: ', key, 'reaches steady state ', ss_orig, 'freq: ', freq_orig)
                #print('original file, IC: ', key, 'reaches steady state ', ss_mod, 'freq: ', freq_mod)
                #print('diff: ', abs(freq_orig - freq_mod))
            print('new IC \n')


    #RMS = math.sqrt(diff/float(N))
    RMS = math.pow(diff/float(N),0.25)
    #print('RMS: ', RMS, diff, N)
    return RMS

if __name__ == '__main__':
    comparator_single_test('simulation_cpp_testing.json', 'EMT_paper.json')


#print(comparator("modified_run.json"))


#x = comparator()
#print(x)

        # for i in range(len(modified_data[key])):
        #     for j in range(len(orig_data[key])):
        #         if
        #
        #         print(modified_data[key][i][1])
        #         print('original file, IC: ', key, 'reaches steady state ', orig_data[key][i][0], 'freq: ', orig_data[key][i][1])
        #         print('modified file, IC: ', key, 'reaches steady state ', modified_data[key][i][0], 'freq: ', float(modified_data[key][i][1])*0.1)
        #         print('diff: ', abs(float(orig_data[key][i][1])-float(modified_data[key][i][1])*0.1))
