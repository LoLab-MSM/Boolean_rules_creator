# Boolean_rules_creator

The optimization.py file is the most important run file that includes
  * the execution of the rule creation from the file rule_creator.py using the function creationg_rules()
  * the translation of the newly created rules to C++ and the compilation of the resulting file based on the template c_simluator.cpp_template
  * the execution of the DEAP based genetic algorithm to find an optimal model
  
To call on the expert knowledge optimization, run the file opt_human.py, which bases its model selection on the elimination of dependencies
