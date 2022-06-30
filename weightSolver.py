# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 10:26:07 2019

@author: NickolasT
"""

import numpy as np
import mapper as mp
import time
import random

# This file takes array of solver objects at individual angles and adjusts which
# leaf position solution and weights to use to optimize full 3D fluence result

def genetic(idealVolume, solverArray, initWeights,solverParams,logFile):
    # solverArray[up/low/avg][angle][row][x]
    # initWeights[angle] (INITIAL GUESS)
    
    # Algorithm is a genetic evolution by reproducing from pairs of individuals
    # of the previous generation, using either cross-over (retain genes)
    # or randomize any gene. Choosing reproductive individuals is based on
    # proportionate selection with a generic fitness bias factor
    # ================PARAMETERS==================
    GENERATIONS = 3                # Number of generations to evolve
    POP_SIZE = 10                   # Number of individuals per generation
                                    # inclusive of two best performing parents
    
    CROSS_FREQ = .8                 # Frequency to generate cross-over child
                                    # else, then generate mutated child
    
    CROSS_WEIGHT_THRESH = 5E-2      # Threshold width to match weights of two parents
                                    # to do cross-over
    MUT_WEIGHT_AMOUNT = .3          # Width of weight mutation allowed +/-
    MUT_GENE_FREQ = .1              # Frequency to mutate a gene in mutated child
    FITNESS_BIAS = 3                # Factor to favor selecting high fitness children
                                    # to reproduce for next gen
    #================INITIALIZE====================
    nChoices = len(solverArray)     # Number of choices to choose from
    nAngles = len(initWeights)      # Number of gantry angles
    bestScore = 1E10                # Tracks the best score
    
    # Get fluenceArray from solvers (indices are [up/low/avg][angle][row][x])
    fluencesArray = []
    for i in range(len(solverArray)):
        fluencesArray.append(mp.getSolverArrayFluences(solverArray[i]))
    fluencesArray = np.array(fluencesArray)
    
    choices = np.ones((POP_SIZE, nAngles))*(nChoices-1)     # 2D array of indices for which leaf pos solution to use
    weights = np.array(initWeights)                         # 1D array of weights for respective leaf solutions
    for i in range(POP_SIZE):
        weights = np.vstack((weights,initWeights))
        
    # =============EVOLVE========================
    for g in range(GENERATIONS):
        t = time.time()
        
        # =================REPRODUCE GEN=====================
        for p in range(1,POP_SIZE):
            # indexA and B denote which 2 individuals will reproduce a child
            indexA = 0
            indexB = 1
            # ***Do proportionate selection for parents after 0th gen***
            if g > 0:
                parentA = random.random()
                parentB = random.random()
                # Find index for parents given CDF of previous gen scores
                for i in range(POP_SIZE):
                    if scoresCDF[i] < parentA < scoresCDF[i+1]:
                        indexA = i
                    if scoresCDF[i] < parentB < scoresCDF[i+1]:
                        indexB = i
            # Reproduce child using either cross-over or mutation
            if random.random() < CROSS_FREQ:
                # Produce cross-over child
                for i in range(nAngles):
                    # Weights (if weights are about equal, then keep it, else mutate it)
                    if (abs(weights[indexB][i] - weights[indexA][i]) < CROSS_WEIGHT_THRESH):
                        weights[p][i] = weights[indexA][i]
                    else:
                        weights[p][i] += random.uniform(-1,1) * MUT_WEIGHT_AMOUNT
                        #need to make sure there are no negative weights 
                        weights[p][i] = max([weights[p][i], 0])         # !!!ADD CODE!!!
                    # Choices
                    if (choices[indexB][i] != choices[indexA][i]):
                        if random.random() < .5:
                            choices[p][i] = (choices[p][i]+1) % nChoices
                        else:
                            choices[p][i] = (choices[p][i]-1) % nChoices
                    else:
                        choices[p][i] = choices[indexA][i]
            else:
                # Produce mutated child
                for i in range(nAngles):
                    # Weights
                    if random.random() < MUT_GENE_FREQ:
                        weights[p][i] += random.uniform(-1,1) * MUT_WEIGHT_AMOUNT
                        weights[p][i] = max([weights[p][i], 0])         # !!!ADD CODE!!!
                    # Choices
                    if random.random() < MUT_GENE_FREQ:
                        if random.random() < .5:
                            choices[p][i] = (choices[p][i]+1) % nChoices
                        else:
                            choices[p][i] = (choices[p][i]-1) % nChoices
        # ====================FITNESS FUNCTION====================
        # Calculate the fitness of all individuals in the population
        bestScoreIndex = 0
        scores = np.ones(POP_SIZE)
        scores[0] = 1/bestScore**FITNESS_BIAS   # Generic FITNESS_BIAS to favor
                                                # high fitness individuals for reproduction
        for p in range(1,POP_SIZE):
            outFluences = []
            for c,choice in enumerate(choices[p]):
                outFluences.append(fluencesArray[int(choice)][c]*weights[p][c])
            volume = mp.mapFluences(outFluences,solverParams)
            score = mp.diffVolume(volume,idealVolume)[1]
            scores[p] = 1/score**FITNESS_BIAS
            # Update best score
            if score <= bestScore:
                bestScoreIndex = p
                bestScore = score
        # ===================PREP NEXT GEN=============================
        # Keep best child and move to index 0
        choices = np.vstack((choices[bestScoreIndex],
                             choices[:bestScoreIndex],
                             choices[bestScoreIndex+1:]))
        weights = np.vstack((weights[bestScoreIndex],
                             weights[:bestScoreIndex],
                             weights[bestScoreIndex+1:]))
        # Determine selection reproduction by calculating likelihood of reproduction
        # based on fitness factor
        # Calculate probability array for proportionate selection
        scoresCDF = [0]
        cummulScore = 0
        for i in scores:
            cummulScore += i/sum(scores)
            scoresCDF.append(cummulScore)
        print("Generation: %d/%d\tError: %f\tGenTime: %f s\tEstTime: %f s" % (g,GENERATIONS,
                                                                                bestScore, time.time() - t,
                                                                                (time.time() - t)*(GENERATIONS - g)))  
        logFile.write("Generation: %d/%d\tError: %f\tGenTime: %f s\tEstTime: %f s\n" % (g,GENERATIONS,
                                                                                bestScore, time.time() - t,
                                                                                (time.time() - t)*(GENERATIONS - g)))  
    # Return the best result
    # Assume best individual is always index 0
    outWeights = weights[0]
    outSolverArray = []
    for i in range(nAngles):
        outSolverArray.append(solverArray[int(choices[0][i])][i])
    return outSolverArray, outWeights
def geneticOLD(idealVolume, solverArray, initWeights,solverParams, logFile):
    # solverArray[up/low/avg][angle][row][x]
    # initWeights[angle] (INITIAL GUESS)
    
    # Algorithm is a genetic evolution by reproducing from the top two individuals
    # of the previous generation, using either cross-over (retain genes)
    # or randomize any gene.
    # ================PARAMETERS==================
    GENERATIONS = 80                # Number of generations to evolve
    POP_SIZE = 20                   # Number of individuals per generation
                                    # inclusive of two best performing parents
    
    CROSS_FREQ = .6                 # Frequency to generate cross-over child
                                    # else, then generate mutated child
    
    CROSS_WEIGHT_THRESH = 1E-2      # Width of weight mutation allowed +/-
    MUT_WEIGHT_AMOUNT = .3          # Frequency to mutate a gene in mutated child
    MUT_GENE_FREQ = .2              # Factor to favor selecting high fitness children
                                    # to reproduce for next gen
    #================INITIALIZE====================
    nChoices = len(solverArray)     # Number of choices to choose from
    nAngles = len(initWeights)      # Number of gantry angles
    bestScore = 1E10                # Tracks the best score
    bestScore2 = 1E10               # Tracks the second best score
    
    # Get fluenceArray from solvers (indices are [up/low/avg][angle][row][x])
    fluencesArray = []
    for i in range(len(solverArray)):
        fluencesArray.append(mp.getSolverArrayFluences(solverArray[i]))
    fluencesArray = np.array(fluencesArray)
    
    choices = np.ones((POP_SIZE, nAngles))*(nChoices-1)
    weights = np.array(initWeights)
    for i in range(POP_SIZE):
        weights = np.vstack((weights,initWeights))
        
    # =============EVOLVE========================
    for g in range(GENERATIONS):
        t = time.time()
        
        # =================REPRODUCE GEN=====================
        for p in range(2,POP_SIZE):
            # Reproduce child using either cross-over or mutation
            if random.random() < CROSS_FREQ:
                # Produce cross-over child
                for i in range(nAngles):
                    # Weights
                    if (abs(weights[1][i] - weights[0][i]) < CROSS_WEIGHT_THRESH):
                        weights[p][i] = weights[0][i]
                    else:
                        weights[p][i] += random.uniform(-1,1) * MUT_WEIGHT_AMOUNT
                    # Choices
                    if (choices[1][i] != choices[0][i]):
                        if random.random() < .5:
                            choices[p][i] = (choices[p][i]+1) % nChoices
                        else:
                            choices[p][i] = (choices[p][i]-1) % nChoices
                    else:
                        choices[p][i] = choices[0][i]
            else:
                # Produce mutated child
                for i in range(nAngles):
                    # Weights
                    if random.random() < MUT_GENE_FREQ:
                        weights[p][i] += random.uniform(-1,1) * MUT_WEIGHT_AMOUNT
                    # Choices
                    if random.random() < MUT_GENE_FREQ:
                        if random.random() < .5:
                            choices[p][i] = (choices[p][i]+1) % nChoices
                        else:
                            choices[p][i] = (choices[p][i]-1) % nChoices
        # ====================FITNESS FUNCTION====================
        bestScoreIndex = 0
        bestScore2Index = 1
        for p in range(2,POP_SIZE):
            outFluences = []
            for c,choice in enumerate(choices[p]):
                outFluences.append(fluencesArray[int(choice)][c]*weights[p][c])
            volume = mp.mapFluences(outFluences,solverParams)
            score = mp.diffVolume(volume,idealVolume)[1]
            # Update best score and its position
            if score <= bestScore:
                bestScoreIndex = p
                bestScore = score
            # Update second best score and its position
            elif bestScore < score < bestScore2:
                bestScore2Index = p
                bestScore2 = score
        # ===================PREP NEXT GEN=============================
        # Move top two to front of population to be parents of next generation
        
        # May need to offset bestScoreIndex if the second best moves it
        if (bestScore2Index > bestScoreIndex):
            bestScoreIndex += 1
        
        choices = np.vstack((choices[bestScore2Index],
                             choices[:bestScore2Index],
                             choices[bestScore2Index+1:]))
        weights = np.vstack((weights[bestScore2Index],
                             weights[:bestScore2Index],
                             weights[bestScore2Index+1:])) 
        choices = np.vstack((choices[bestScoreIndex],
                             choices[:bestScoreIndex],
                             choices[bestScoreIndex+1:]))
        weights = np.vstack((weights[bestScoreIndex],
                             weights[:bestScoreIndex],
                             weights[bestScoreIndex+1:]))
        
        
        print("Generation: %d/%d\tError: %f\tGenTime: %f s\tEstTime: %f s" % (g,GENERATIONS,
                                                                                bestScore, time.time() - t,
                                                                                (time.time() - t)*(GENERATIONS - g)))  
        logFile.write("Generation: %d/%d\tError: %f\tGenTime: %f s\tEstTime: %f s\n" % (g,GENERATIONS,
                                                                                bestScore, time.time() - t,
                                                                                (time.time() - t)*(GENERATIONS - g)))  
    # Return the best result
    outWeights = weights[0]
    outSolverArray = []
    for i in range(nAngles):
        outSolverArray.append(solverArray[int(choices[0][i])][i])
    return outSolverArray, outWeights


def randSearch(idealVolume, solverArray, initWeights,solverParams, logFile):
    # solverArray[up/low/avg][angle][row][x]
    # initWeights[angle] (INITIAL GUESS)
    
    # Algorithm is a randomizing search by reproducing from the best individual
    # of the previous generation, using either cross-over (retain genes)
    # or randomize any gene.
    # ================PARAMETERS==================
    GENERATIONS = 20                # Number of generations to evolve
    POP_SIZE = 10                   # Number of individuals per generation
                                    # inclusive of one best performing parent
    
    CROSS_FREQ_WEIGHT = .7          # Frequency for cross-over of weight genes
    CROSS_FREQ_CHOICE = .7          # Frequency for cross-over of choice genes
    MUT_FREQ_WEIGHT = .8            # Frequence for mutation of weight genes
    MUT_FREQ_CHOICE = .8            # Frequency for mutation of choice genes
    
    CROSS_WEIGHT_THRESH = 1E-2      # Difference threshold for weights equaling
    MUT_WEIGHT_AMOUNT = .1          # Width of allowed weight mutation +/-
    MUT_GENE_FREQ = .8              # Frequency to mutate a gene
    
    #================INITIALIZE====================
    nChoices = len(solverArray)     # Number of choices to choose from
    nAngles = len(initWeights)      # Number of gantry angles
    bestScore = 1E10                # Tracks the best score
    
    # Get fluenceArray from solvers (indices are [up/low/avg][angle][row][x])
    fluencesArray = []
    for i in range(len(solverArray)):
        fluencesArray.append(mp.getSolverArrayFluences(solverArray[i]))
    fluencesArray = np.array(fluencesArray)
    
    # Choice array (indices are [individual][angle])
    choices = np.ones((POP_SIZE, nAngles))*(nChoices-2)
    # weight array (indices are [individual][angle])
    weights = np.array(initWeights)
    for i in range(POP_SIZE):
        weights = np.vstack((weights,initWeights))
        
    # =============EVOLVE========================
    for g in range(GENERATIONS):
        t = time.time()
        
        # =================REPRODUCE GEN=====================
        for p in range(2,POP_SIZE):
            # ----------Weight----------
            if random.random() < CROSS_FREQ_WEIGHT:
                for i in range(nAngles):
                    if (abs(weights[p][i] - weights[0][i]) > CROSS_WEIGHT_THRESH) and random.random() < MUT_FREQ_WEIGHT:
                        weights[p][i] += random.uniform(-1,1) * MUT_WEIGHT_AMOUNT 
            elif random.random() < MUT_FREQ_WEIGHT:
                for i in range(nAngles):
                    if random.random() < MUT_GENE_FREQ:
                        weights[p][i] += random.uniform(-1,1) * MUT_WEIGHT_AMOUNT 
                    
            # ----------Choice---------    
            if random.random() < CROSS_FREQ_CHOICE:
                for i in range(nAngles):
                    if (choices[p][i] != choices[0][i]):
                        if random.random() < .5:
                            choices[p][i] = (choices[p][i]+1) % nChoices
                        else:
                            choices[p][i] = (choices[p][i]-1) % nChoices
            
            elif random.random() < MUT_FREQ_CHOICE:
                for i in range(nAngles):
                    if random.random() < MUT_GENE_FREQ:
                        if random.random() < .5:
                            choices[p][i] = (choices[p][i]+1) % nChoices
                        else:
                            choices[p][i] = (choices[p][i]-1) % nChoices
            
        # ====================FITNESS FUNCTION====================
        bestScoreIndex = 0
        for p in range(2,POP_SIZE):
            outFluences = []
            for c,choice in enumerate(choices[p]):
                outFluences.append(fluencesArray[int(choice)][c]*weights[p][c])
            volume = mp.mapFluences(outFluences,solverParams)
            score = mp.diffVolume(volume,idealVolume)[1]
            # Update best score and its position
            if score <= bestScore:
                bestScoreIndex = p
                bestScore = score
        
        # Prepare for next Generation
        # Move top two to front of population to be parents of next generation
        choices = np.vstack((choices[bestScoreIndex],
                             choices[:bestScoreIndex],
                             choices[bestScoreIndex+1:]))
        weights = np.vstack((weights[bestScoreIndex],
                             weights[:bestScoreIndex],
                             weights[bestScoreIndex+1:]))  
        print("Generation: %d/%d\tError: %f\tGenTime: %f s\tEstTime: %f s" % (g,GENERATIONS,
                                                                                bestScore, time.time() - t,
                                                                                (time.time() - t)*(GENERATIONS - g)))  
        logFile.write("Generation: %d/%d\tError: %f\tGenTime: %f s\tEstTime: %f s\n" % (g,GENERATIONS,
                                                                                bestScore, time.time() - t,
                                                                                (time.time() - t)*(GENERATIONS - g)))  
    # Return the best result
    outWeights = weights[0]
    outSolverArray = []
    for i in range(nAngles):
        outSolverArray.append(solverArray[int(choices[0][i])][i])
    return outSolverArray, outWeights
