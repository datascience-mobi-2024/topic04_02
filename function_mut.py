def diff_weighted(feature_pos, feature_neg, aa:str, ideal_pos:dict, ideal_neg:dict, sort = True):
    """
    Calculates the weighted sum of deviations for a given amino acid ('aa') based on positive and negative features.

    Args:
        feature_pos (list): List representing the positive features.
        feature_neg (List): List representing the negative features.
        aa (str): String of the amino acid sequence.
        ideal_pos (dict): Dictionary containing ideal values for positive features.
        ideal_neg (dict): Dictionary containing ideal values for negative features.
        sort (bool, optional): Flag indicating whether to sort the features by deviation (default: True).

    Returns:
        tuple: A tuple containing two elements:
            - sum_dev (float): The total weighted sum of deviations for all features.
            - sorted_keys (list, optional): If sort=True, a list of features sorted by their weighted deviation (highest first).
                - WT_weight (dict, optional): If sort=False, a dictionary containing the weighted deviation for each feature.
    """
    from function import rel_aa_comp
    import operator
    
    WT_weight = {}
    sum_dev = 0
    
    #positive features
    for key in feature_pos:
        if len(key) <=2:
            
            weighted_diff = abs(rel_aa_comp(aa, key) - (ideal_pos[key] * feature_pos[key]))
            WT_weight[key] = weighted_diff
            sum_dev += weighted_diff
            
        elif 'motif' in key:
            weighted_diff = abs(rel_aa_comp(aa, key[0:2]) - (ideal_pos[key] * feature_pos[key]))
            WT_weight[key] = weighted_diff
            sum_dev += weighted_diff
            
    
    #negative features
    for key in feature_neg:
        if len(key) <=2:
            weighted_diff = abs(rel_aa_comp(aa, key) - ideal_neg[key] * feature_neg[key])
            WT_weight[key] = weighted_diff
            sum_dev += weighted_diff
            
        elif 'motif' in key:
            weighted_diff = abs(rel_aa_comp(aa, key[0:2]) - (ideal_pos[key] * feature_pos[key]))
            WT_weight[key] = weighted_diff
            sum_dev += weighted_diff    

    if sort:
        sorted_keys = sorted(WT_weight.items(), key=operator.itemgetter(1), reverse=True)
        return sum_dev, sorted_keys
    else:
        return sum_dev, WT_weight #sum_dev is a positive value of all deviations, higher values indicate a worse fit
    
    

def mut_apply(AA_list, Mut_list):
    """
    Applies a list of mutations to a list of amino acids.

    Args:
        AA_list (list): A list containing the original amino acid sequence.
        Mut_list (list): A list of mutation strings defining the substitutions to be applied in the form ('WT-POS-MUT')

    Returns:
        list: A new list containing the amino acid sequence after applying the mutations.
    """
    if len(Mut_list) > 0:
        if 'M' in Mut_list[0]:
            Mut_list = Mut_list[1:]
        for n in Mut_list:
            AA_pos = int(n.split('-')[1]) - 1
            AA_mut = n.split('-')[2]
            AA_list[AA_pos] = AA_mut
        return AA_list
    else:
        return(print('Mut list is empty'))

def mut_live_test (AA_list, Mut_list, pos_corr, neg_corr, ideal_pos_value, ideal_neg_value):
    from function_mut import mut_apply
    from function_mut import diff_weighted
    AAs = ''.join(AA_list)
    WT_sum, WT_diff = diff_weighted(pos_corr, neg_corr, AAs, ideal_pos_value, ideal_neg_value)
    AA_mut = mut_apply(AA_list, Mut_list)
    AA_muts = ''.join(AA_mut)
    MUT_sum, WT_diff = diff_weighted(pos_corr, neg_corr, AA_muts, ideal_pos_value, ideal_neg_value)
    Diff = abs(WT_sum) - abs(MUT_sum) # calculates the difference between the weighted sum of deviations before and after the mutation
    
    # if the difference is positive the mutation is beneficial
    # if Diff is negative the mutation is not beneficial
    return Diff #The higher the difference the better the mutation
    
    
def mutator2(AA_list:list, free_AA, deviation,pos_corr:dict, sorted_freq_pos, neg_corr:dict, conserv_substitution, ideal_pos_value, ideal_neg_value):
    """
    Generates a list of potential mutations based on deviations and correlations.

    Args:
        AAs (str): The amino acid sequence.
        free_AA (np.ndarray): Array containing information about free amino acids Col1: Prot name, Col2: aminoacid position, Col3: Aminoacid (one letter code).
        deviation (list): List containing deviations from ideal values for features.
        pos_corr_list (list): List of amino acids that positively correlate with desired features.
        sorted_freq_pos (list): Possibly sorted list of frequencies for amino acids contributing to positive features (usage unclear).
        neg_corr_list (list): List of amino acids that negatively correlate with desired features.
        conserv_substitution (dict): Dictionary containing a list of possible conservative substitutions for each amino acid.

    Returns:
        list: A list of mutation strings defining potential substitutions.
    """
    
    import itertools
    from function import rel_aa_comp
    from function_mut import mut_apply
    from function_mut import mut_live_test
    
    pos_corr_list = list(pos_corr.keys())
    neg_corr_list = list(neg_corr.keys())
    AAs = ''.join(AA_list)
    free_AA_dict = {a: b for a, b in zip(free_AA[:,2], free_AA[:,1] )} # create dictionary from array the key is the absolute aminoacid position and the value the aminoacid
    first_entry = deviation[0][0]
    mut_list = []
    # determine possible substitutions if the first entry is a single amino acid
    if len(first_entry) == 1:
        
        if first_entry in pos_corr_list: #checks if amount of aminoacid should be increaed
            for key in free_AA_dict:
                aminoacid = free_AA_dict[key]
                AA_subst = conserv_substitution[aminoacid]
                
                
                for k in pos_corr_list: # Check if the current amino acid is in the positive correlation list
                    if len(k) == 1:     #selects the first feature with one aminoacid
                        if rel_aa_comp(AAs, k) < ideal_pos_value[k]:    #checks if the relative composition is suboptimal
                            if k in AA_subst:
                                mut_aa = (aminoacid + '-' + key + '-' + k)
                                Diff = mut_live_test(AA_list, [mut_aa], pos_corr, neg_corr, ideal_pos_value, ideal_neg_value) # live test if mutation is benefitical
                                if Diff > 0:
                                    mut_list.append(mut_aa)
                                    AA_list = mut_apply(AA_list, [mut_aa]) #applies mutation to the protein sequence
                                    break


        elif first_entry in neg_corr_list: #checks if amount of aminoacid should be decreased
            for key in free_AA_dict:
                aminoacid = free_AA_dict[key]
                if aminoacid == first_entry:
                    AA_subst = conserv_substitution[aminoacid]
                    
                    #tries to substitute the aminoacid to the aminoacids that comes first in the sorted_freq_pos list
                    for k in pos_corr_list:
                        if aminoacid == k:
                            if rel_aa_comp(AAs, k) < ideal_pos_value[k]:
                                mut_aa = aminoacid + '-' + key + '-' + k
                                Diff = mut_live_test(AA_list, [mut_aa], pos_corr, neg_corr, ideal_pos_value, ideal_neg_value) # live test if mutation is benefitical
                                if Diff > 0: # live test if mutation is benefitical
                                    AA_list = mut_apply(AA_list, [mut_aa]) #applies mutation to the protein sequence
                                    mut_list.append(mut_aa)
                                    break
                            
                            
 
    # determine possible substitutions if the first entry is a pair of amino acids
    elif len(first_entry) == 2:
        
        if first_entry in pos_corr_list: # checks if the amount of aminoacids should be increased
            for key in free_AA_dict:
                aminoacid = free_AA_dict[key]
                #if aminoacid not in sorted_freq_pos: #checks if the aminoacid overall positively contributes to one of the pos_corr features
                AA_subst = conserv_substitution[aminoacid] # list of possible substitutions for the current amino acid
                
                #mutation
                for n in range(len(deviation)): # iterates through all deviations and takes the first deviation (from pos corr) which can be increaed
                    entry = deviation[n][0]
                    if entry in pos_corr_list:
                        subst = [] #creates a list of possible substitutions that increase one of the amino acids in the highest entry that positively correlates
                        for k in AA_subst:
                            if k in first_entry:
                                subst = k
                                
                        if len(subst) == 1: # if only one substitution increases one of the aminoacids this substitution will be used
                            mut_aa = aminoacid + '-' + key + '-' + subst
                            Diff = mut_live_test(AA_list, [mut_aa], pos_corr, neg_corr, ideal_pos_value, ideal_neg_value) # live test if mutation is benefitical
                            if Diff > 0:
                                AA_list = mut_apply(AA_list, [mut_aa]) #applies mutation to the protein sequence
                                mut_list.append(mut_aa)
                                break
                            
                        elif len(subst) == 2: # if the aminoacid can be substituted to both aminoacids in first_entry choose the one that has the least frequency
                            comp = [(aa, rel_aa_comp(AAs, aa)) for aa in subst] # calculate the relative amino acid composition of the possible substitutions
                            lowest_comp = min(comp, key=lambda pair: pair[1]) # find the amino acid with the lowest relative composition
                            mut_aa = aminoacid + '-' + key + '-' + lowest_comp[0]
                            Diff = mut_live_test(AA_list, [mut_aa], pos_corr, neg_corr, ideal_pos_value, ideal_neg_value) # live test if mutation is benefitical
                            if Diff > 0:
                                AA_list = mut_apply(AA_list, [mut_aa])
                                mut_list.append(mut_aa) # append the mutation to the mutation list
                                break

         
                    
        elif first_entry in neg_corr_list: # checks if the amount of aminoacids should be decreased
            for key in free_AA_dict:
                aminoacid = free_AA_dict[key]
                if aminoacid in first_entry: #checks if the aminoacid is present in the first_entry
                    AA_subst = conserv_substitution[aminoacid] # list of possible substitutions for the current amino acid
                    
                    #mutation               
                    for n in range(len(deviation)): # iterates through all deviations and takes the first deviation (from pos corr) which can be increaed
                        entry = deviation[n][0]
                        if entry in pos_corr_list: 
                            subst = [] #creates a list of possible substitutions that increase one of the amino acids in the first_entry
                            for k in AA_subst:
                                if k in entry:
                                    subst.append(k)
                            if len(subst) == 1: # if only one substitution increases one of the aminoacids this substitution will be used
                                mut_aa = aminoacid + '-' + key + '-' + subst[0]
                                Diff = mut_live_test(AA_list, [mut_aa], pos_corr, neg_corr, ideal_pos_value, ideal_neg_value) # live test if mutation is benefitical
                                if Diff > 0:
                                    AA_list = mut_apply(AA_list, [mut_aa])
                                    mut_list.append(mut_aa) 
                                    break
                                
                            elif len(subst) == 2: # if the aminoacid can be substituted to both aminoacids in first_entry choose the one that has the least frequency
                                comp = [(aa, rel_aa_comp(AAs, aa)) for aa in subst] # calculate the relative amino acid composition of the possible substitutions
                                lowest_comp = min(comp, key=lambda pair: pair[1]) # find the amino acid with the lowest relative composition
                                mut_aa = aminoacid + '-' + key + '-' + lowest_comp[0]
                                Diff = mut_live_test(AA_list, [mut_aa], pos_corr, neg_corr, ideal_pos_value, ideal_neg_value) # live test if mutation is benefitical
                                if Diff > 0:
                                    AA_list = mut_apply(AA_list, [mut_aa])
                                    mut_list.append(mut_aa) # append the mutation to the mutation list
                                    break
    
    
    elif 'motif' in first_entry:
        if first_entry in pos_corr_list: # checks if the amount of aminoacids should be increased
            for key in free_AA_dict:
                    aa_back = AAs[key-1]    # aminoacid before the current aminoacid
                    aa_current = AAs[key]   # current aminoacid
                    aa_for = AAs[key+1]     # aminoacid after the current aminoacid
                    current_list = [aa_back, aa_current, aa_for] # list of the three aminoacids
                    
                    if (aa_current in first_entry[0]) and (aa_for not in first_entry[1]):    #checks if current aminoacid is part of the motif, and next aminoacid not
                        for_subst = conserv_substitution[aa_back] # possible substitutions for the aminoacid before the motif
                        for s in for_subst: #iterates through the possible substitutions
                            if s in first_entry[1]:
                                mut_aa = f'{aa_back}-{key+1}-{s}'
                                Diff = mut_live_test(AA_list, [mut_aa], pos_corr, neg_corr, ideal_pos_value, ideal_neg_value) # live test if mutation is benefitical
                                if Diff > 0:
                                    AA_list = mut_apply(AA_list, [mut_aa])
                                    mut_list.append(mut_aa)
                                    break  
                                          
                    if (aa_current in first_entry[1]) and (aa_back not in first_entry[0]):   #checks if current aminoacid is part of the motif, and previous aminoacid not
                        back_subst = conserv_substitution[aa_for] # possible substitutions for the aminoacid after the motif
                        for b in back_subst: #iterates through the possible substitutions
                            if b in first_entry[0]:
                                mut_aa = f'{aa_for}-{key-1}-{b}'
                                Diff = mut_live_test(AA_list, [mut_aa], pos_corr, neg_corr, ideal_pos_value, ideal_neg_value) # live test if mutation is benefitical
                                if Diff > 0:
                                    AA_list = mut_apply(AA_list, [mut_aa])
                                    mut_list.append(mut_aa)
                                    break
    
    return AA_list, mut_list