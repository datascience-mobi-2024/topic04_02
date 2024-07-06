class Protein:
    '''
    Protein class that takes a sequence and returns a SPARC prediction, optionally with a PDB file for mutational improvement of melting point. Computes features for the sequence and stores them as attributes of the object. Full list of features can be accessed by calling the object's .features attribute.
    
    Initial attributes:
    seq: str, protein sequence
    PDB: str, name of PDB file in ./data/pdbs directory, optional for SPARC but required for .mutate() method
    dpath: str, path to directory where temp files will be stored, defaults to ./data/fastas
    S4path: str, path to directory where S4pred is stored, required for SPARC, defaults to ./data/s4pred
    
    Methods:
    mutate( ): performs mutational improvement of melting point, requires PDB file. Returns mutated melting point, sequence, and difference in melting point and adds them to the object as .mutatedmp, .mutatedseq, and .mutatedmpdiff
    mutreduce( ): reduces amount of mutations while still maintaining improved melting point. Requires mutate() method to be used first. Returns reduced sequence and melting point and adds them to the object as .mut_seq_reduced and .mut_mp_reduced
    '''
    import os
    import numpy as np
    def __init__(self, seq:str, PDB = None, dpath = os.path.abspath('./data/fastas'), S4path = os.path.abspath('./data/s4pred')):
        self.PDB = PDB
        from SPARC import SPARC
        self.sequence = seq.upper()
        result = SPARC(self.sequence, 'test', dpath, S4path)
        self.SPARC = result[0][0]
        self.features = result[1]
        self.features_scaled = result[2]
        self.features_pca = result[3]
        for i in self.features.index:
            setattr(self, i, self.features[i])
    def __str__(self):
        return f"Protein with predicted melting point of {str(self.SPARC).strip('[]')} °C"
    def mutate(self,pqr_output_path='./data/pqrs',iterations = 100, threshhold = 1000):
        from function_mut import prot_mut
        output_mut = prot_mut(pdb_path = './data/pdbs', pdb_file=self.PDB,pqr_output_path=pqr_output_path,iterations = iterations, threshhold = threshhold)
        self.mutatedmp = output_mut[0][1][0][0]
        self.mutatedseq = output_mut[1][0]
        self.mutatedmpdiff = self.mutatedmp - self.SPARC
        mutationlist = []
        for n in range(len(self.mutatedseq)):
            if self.mutatedseq[n] != self.sequence[n]:
                mutationlist.append(f'{self.sequence[n]}{n+1}{self.mutatedseq[n]}')
        self.mutationlist = mutationlist
        return [self.mutatedmp, self.mutatedseq, self.mutatedmpdiff, self.mutationlist]
    def mutreduce(self,name = 'proteinxyz'):
        from function_mut import mutation_decreaser
        outputdecreaser = mutation_decreaser(mut_temp = self.mutatedmp, wt_temp = self.SPARC, wt_protein = self.sequence, mut_protein = self.mutatedseq, name = name)
        self.mut_mp_reduced = outputdecreaser[2]
        self.mut_seq_reduced = outputdecreaser[0]
        mutationlist_reduced = []
        for n in range(len(self.mut_seq_reduced)):
            if self.mut_seq_reduced[n] != self.sequence[n]:
                mutationlist_reduced.append(f'{self.sequence[n]}{n+1}{self.mut_seq_reduced[n]}')
        self.mutationlist_reduced = mutationlist_reduced
        removed_mutations = []
        for n in range(len(self.mutatedseq)):
            if self.mutatedseq[n] != self.mut_seq_reduced[n]:
                removed_mutations.append(f'{self.mutatedseq[n]}{n+1}')
        self.removed_mutations = removed_mutations
        return [self.mut_mp_reduced, self.mut_seq_reduced, self.mutationlist_reduced]