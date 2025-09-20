from src.defineReactionNetwork import generate_reaction_network

rule generate_reaction_network:
    input: 
        data/case{case_number}/enzymatic_reactions.csv data/case{case_number}/spontaneous_reactions.csv
    output:
        data/case{case_number}/reaction_network
    run:
        generate_reaction_network(input)