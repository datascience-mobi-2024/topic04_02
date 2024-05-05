#function that calculates the relative amino acid composition based on sequence entry and definded amino acid property
def rel_aa_comp(Sequence:str, AA_property:str) -> list: 
    count = 0
    for element in Sequence:  
        if element in AA_property:
            count += 1
    return count

def rel_aa(Sequence:str, AA_property:str) -> list: 
    count = 0
    for element in Sequence:  
        if element in AA_property:
            count += 1
    return count