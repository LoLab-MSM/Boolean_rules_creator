

import rule_creator

fn = 'MM_steady_states.json'
symbols = ['EN', 'S', 'ES', 'P']

str_rules, simple_rulelist, fs_cpp_name = creating_rules(fn,symbols,0)



print('\n\nString rules:')
print(str_rules)

print('simple rulelist')
print(simple_rulelist)
