def diff_weighted(feature_pos, feature_neg, aa:str, ideal_pos:float, ideal_neg:float, sort = True):
    from function import rel_aa_comp
    import operator
    WT_weight = {}
    sum_dev = 0
    for key in feature_pos:
        if len(key) <=2:
            weighted_diff = abs(rel_aa_comp(aa, key) - ideal_pos * feature_pos[key])
            WT_weight[key] = weighted_diff
            sum_dev += weighted_diff
    
    for key in feature_neg:
        if len(key) <=2:
            weighted_diff = abs(rel_aa_comp(aa, key) - ideal_neg * feature_neg[key])
            WT_weight[key] = weighted_diff
            sum_dev += weighted_diff

    if sort:
        sorted_keys = sorted(WT_weight.items(), key=operator.itemgetter(1), reverse=True)
        return sum_dev, sorted_keys
    else:
        return sum_dev, WT_weight