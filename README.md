# Boolean_rules_creator
![](https://img.shields.io/badge/Python-v3.7-informational?style=flat&logo=python&logoColor=white&color=2bbc8a?link=http://https://www.python.org/left&link=http://right) ![](https://img.shields.io/badge/JsonCPP-v11-informational?style=flat&logo=C++&logoColor=white&color=2bbc8a)


## Overview
Welcome to the Boolean_rules_creator, a generalizable unsupervised approach to generate parameter-free, logic-based mechanistic hypotheses of cellular processes, described by multiple discrete states. For reference, see the preprint linked below. <div style="page-break-after: always"></div>
Unsupervised logic-based mechanism inference for network-driven biological processes 
(M. Prugger et. al 2020), now available through bioRxiv:
 https://www.biorxiv.org/content/10.1101/2020.12.15.422874v1.full
 <div style="page-break-after: always"></div>
 This tool provides a method of mechanism inference in biological processes where the input states and attractors (steady states) are known.
 The repository contains two different examples, as described in the paper, looking at Enzyme-Substrate Kinetics (ES) kinetics as a simpler model with 4 nodes, and an Epithelial-to-mesenchymal transition (EMT) model with 7 nodes. Descriptions of these models and their Python scripts are outlined below.

----------------------------------------------------------------------------------------
##  &#128295; Installation

Boolean_rules_creator is dependent upon Python libraries DEAP, Numpy, and Sympy, which can be installed usuing pip or conda. 
```shell
pip install numpy sympy deap 
```
Once dependencies are installed, clone the repository using Git:

```shell
> git clone https://github.com/LoLab-VU/Boolean_rules_creator.git
> cd Boolean_rules_creator
```


In addition, you will need a current version of g++ on your computer to run the C++ compiler. We recommend using JsonCpp, available at https://github.com/open-source-parsers/jsoncpp. Follow the installation instructions there and make sure to also follow the link for the amalgamated source in order add the correct path to optimization.py.

----------------------------------------------------------------------------------------

## &#128196; Examples from M. Prugger et. al. 2020
### Enzyme-Substrate Kinetics

For the ES kinetics model, there is only one possibly steady state outcome, so the Boolean ruleset can be generated without optimization. In this example, only the creating_rules function from rule_creator.py is used to recreate the forward, backwards, and expert-guided rulesets.
<div style="page-break-after: always"></div>
The ES folder contains the following files:

* Jupyter Notebook ES_rules.ipynb which contains a full walkthrough of the example (references the file fig1b.jpg)
* rule_creator.py which contains the creating_rules function
* ES_rules.py is a short Python script that calls the creating function for the ES data
* ES_steady_states.json is the input file with the possible initial and steady states in the following format:

```
"[initial state1]:
    [[[steady state2],100]],
[initial state2]:
    [[[steady state3],100]],
  ...
```
The number accompanying each pair of states is the frequency with which the intial state will lead to the steady state. Since there is only one steady state for this system, this frequency is 100 for all intial states.


#
### EMT Transition
For models with more than one steady state, as seen with the EMT model in the paper, unsupervised model optimization is used to generate the correct rule list.

In order to run the optimization, use the 'optimization.py' file, which includes:
* the execution of the rule creation from the file rule_creator.py using the function creationg_rules()
* the translation of the newly created rules to C++ and the compilation of the resulting file based on the template c_simluator.cpp_template
* the execution of the DEAP based genetic algorithm to find an optimal model

In order for this file to run, installation of a C++ compiler is necessary. See the installation instructions above for more details. The path to this compiler must be included in the optimization.py script, as mentioned in the information about the amalgamated source from Jsoncpp. Modify this to your system within the script in line 43, as shown below:
```python
os.system('time g++ -DOUTPUT_FILE=\\"{2}\\" /path/to/jsoncpp-master/dist/jsoncpp.cpp -I/path/to/jsoncpp-master/dist/json  -O3 -fopenmp -x c++ {0} -o {1} && {1}'.format(file_name,exe_file,json_file))
```

To call on the expert knowledge optimization, run the file opt_human.py, which bases its model selection on the elimination of dependencies.

The EMT file in the repository contains the following files:

* EMT_incbw_ruleX.txt: A list of each transition in the full backward model -> these files are are a starting point to setup the expert knowledge optimization.

* EMT_paper.json: This is the data file that describes the end point of the simulation results from the original model. This is the file to which we compare our automatically  generated models.

* EMT_userguided_ruleX.txt: These files are obtained by running sort_list.py with the EMT_incbw_ruleX.txt files. Since the elimination of dependencies that are not self-dependencies are unique, these files are already sorted out in the according way. The optimizer therefore only reads the transitions where it needs to still decide, whether transitions to remove or not, according to the self-dependency selection.

* c_simulator.cpp_template: A template file on which the corresponding C++ file gets created: Each model creates a different rule set. To quickly execute the asynchronous updating simulation, the new ruleset gets compiled into its own temporary C++ file and runs outside of the python environment. The results are used for the optimizer in the optimization.py file or the opt_human.py

* file_comparison.py: This file performs the comparison between two json files: the EMT_paper.json file and the files generated by the temp files created from c_simulator.cpp_template. It returns the root mean square error between the dynamics given in those files

* opt_human.py: Runfile that performs the optimization based on the expert knowledge. This code needs the input files EMT_userguided_ruleX.txt, since the removal of non-self dependencies are unique, and the optimizer only has to decide on the removal of transitions regarding self-dependencies.

* optimization.py: Runfile that performs the unsupervised optimization over the whole possible search space

* rule_creator.py: Executes the creation of the Boolean rules and also includes the algorithm for the expert knowledge inclusion, where the states are filtered by their appearance in the rule formulation

* sort_list.py: Takes the files EMT_incbw_ruleX.txt as input and removes the desired dependencies (excluding self-dependencies) which are then saved in the files EMT_userguided_ruleX.txt. If we translate these files into rules, the resulting string would no longer include the species that are eliminated in this file. Self-dependencies, however, are still possible in these lists and are then eliminated in the opt_human.py file.

#

The output is a folder "optimization_output" that stores the results for each created model: 
 * the files "fitness-XXX-XXX.txt" contains the information for each created model. The first three numbers determine the generation of the genetic algorithm, the second three numbers determine the individual of the population. The file itself contains the various random executions of the genetic algorithm (equivalent to the variable parallel_sims).
Each entry in these files starts with the computed fitness value, then it displays the determined rules that result in this fitness value in a way such that it can be used to be read by Boolean2Rules. Next, the simulation result for every initial condition is given (this is the same information as is stored and used by the json-files). Finally, we also include the transition list for every rule.
 * the files "run-XXX-XXX.json" are the corrsponding json files.
outside of the "optimization_output" folder, there will also be created an additional file "rms_fit.txt", that stores the rms values for each simulation for quick access.
