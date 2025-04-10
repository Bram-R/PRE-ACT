Probabilistic PRE-ACT State-Transition Model for Cost-Effectiveness Analysis
Project Overview
This repository contains an implementation of a Probabilistic State-Transition Model (STM) for evaluating the cost-effectiveness of healthcare interventions. The model simulates the transitions of a cohort of patients across different health states over time and estimates costs, quality-adjusted life years (QALYs), and life-years (LYs). The analyses use probabilistic sensitivity analysis (PSA) to account for parameter uncertainty and generate cost-effectiveness outputs, including incremental cost-effectiveness ratios (ICERs), cost-effectiveness acceptability curves (CEACs), expected value of perfect information (EVPI), and cost-effectiveness planes.
Key Features
•	State-Transition modeling: implements a discrete-time Markov model with a defined set of health states and transitions.
•	Comprehensive incorporation of toxicity of time
•	Probabilistic Sensitivity Analysis (PSA): Monte Carlo simulations (n = 5000) sample from probability distributions to account for parameter uncertainty.
•	Cost-effectiveness metrics: computes ICERs, cost-effectiveness planes, CEACs, EVPI, and Expected Loss Curves (ELC).
•	Graphical Representation: Generates key plots to visualize results, including density distributions, cost-effectiveness planes, and CEACs.
Installation & dependencies
The model requires the following R packages:
•	docstring, 
•	DiagrammeR, 
•	dampack, 
•	BCEA
Project structure
/project_directory
│── f_input.R                       	# Function to generate input parameters
│── f_model.R                     	# State-transition model function
│── f_gen_pop_utility.R          	# Population utility function
│── f_gen_pop_mortality.R    	# Population mortality function
│── f_interpolate_toxicity.R 	# Function to interpolate toxicity over time
│── main_script.R            	# Main script for execution and analysis
│── README.docx                	# This documentation
Methods

Model structure
A state-transition model (STM) was developed to assess the cost-effectiveness of integrating AI-assisted decision-making in breast cancer treatment. The model consists of four health states:
•	Event-free survival
•	Locoregional recurrence (LRR)
•	Distant metastasis (DM)
•	Death (absorbing state)
Patients transition between states in monthly cycles based on transition probabilities derived from clinical trials and epidemiological data.
The model compares:
•	Current Practice (CP) – Standard of care
•	CP + AI-Assisted Decision-Making (AI) – AI-based screening patients at high-risk for toxicity
The primary outcomes include:
•	Total costs (£)
•	Quality-adjusted life-years (QALYs)
•	Life-years (LYs)
•	Incremental cost-effectiveness ratios (ICERs)

 
