#fasparse: function to parse the output of S4pred
def fasparse(faspath):
    import numpy as np
    helixl = []
    sheetl = []
    data = np.loadtxt(faspath, dtype={'names': ('index', 'col1', 'col2', 'val1', 'val2', 'val3'),'formats': ('i4', 'S1', 'S1', 'f4', 'f4', 'f4')})
    indiceshelix = data['index'][data['col2'] == b'H']
    indicessheet = data['index'][data['col2'] == b'E']
    helix = np.split(indiceshelix, np.where(np.diff(indiceshelix) != 1)[0]+1)
    sheet = np.split(indicessheet, np.where(np.diff(indicessheet) != 1)[0]+1)
    return [helix, sheet]

# SPARC: "Sequence-based Protein Attribute-deRived Melt Point Calculator"
def SPARC(data:str, name, datapath, directoryS4, removefiles = True):
    """
    data: string, protein sequence
    name: arbitrary name for the protein
    datapath: string, path to data directory, only "/" allowed
    directoryS4: string, compltete path to directory of S4pred, only "/" allowed
    """
    import joblib
    import numpy as np
    import pandas as pd
    import os
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
    from sklearn.metrics import mean_squared_error, r2_score
    model1 = joblib.load('./data/gbr_model1.joblib')
    scaler1 = joblib.load('./data/scaler1.joblib')
    pca1 = joblib.load('./data/pca1.joblib')
    prokaryotes1 = joblib.load('./data/prokaryotes1.joblib')
    features = pd.Series(np.zeros(prokaryotes1.shape[1]-1))
    features.index = prokaryotes1.drop(columns='meltPoint').columns
    with open(os.path.join(datapath,f'{name}.fasta'), "w") as fasta_file:
        fasta_file.write(f">{name}\n{data}\n")
    
    fastapath = os.path.join(datapath,f'{name}.fasta')
    faspath = os.path.join(datapath,f'{name}.fas')
    if datapath.startswith('.'):                #added to work for dynamic paths
        fastapath = f'{os.path.abspath(fastapath)}'
        faspath = f'{os.path.abspath(faspath)}'
    if os.path.isfile(faspath):
        pass
    else:
        os.chdir(directoryS4)#replace with AA2s4pred funciton
        os.system(f'python3 run_model.py "{fastapath}" > "{faspath}"')
        os.chdir('../../')
    Protlen = len(data)
    features.iloc[0] = Protlen
    HelixSparc = fasparse(faspath)[0]
    SheetSparc = fasparse(faspath)[1]
    features.iloc[1] = len(HelixSparc)  #counts
    features.iloc[2] = len(SheetSparc)
    helixind = np.concatenate(HelixSparc)-1
    helixseq = np.array(list(data))[list(helixind)]
    sheetind = np.concatenate(SheetSparc)-1
    sheetseq = np.array(list(data))[list(sheetind)]
    features.iloc[3] = len(helixind) / Protlen #percs
    features.iloc[4] = len(sheetind) / Protlen
    features.iloc[5] = np.mean(np.array([len(arr) for arr in HelixSparc])) #avg lengths
    features.iloc[6] = np.mean(np.array([len(arr) for arr in SheetSparc]))
    features.iloc[7] = features.iloc[3] + features.iloc[4]
    aapairs = features.index[8:87]
    for aa in aapairs:
        features[aa] = (data.count(aa[0]) + data.count(aa[1])) / Protlen
    for acid in ['A', 'V', 'I', 'L', 'M', 'F', 'W','N', 'Q', 'S', 'T', 'Y','D', 'E','R', 'H', 'K', 'C', 'P', 'G' ]:
        if len(helixseq) != 0:
            features[f'{acid}helix'] = helixseq.tolist().count(acid) / len(helixseq)
        if len(sheetseq) != 0:
            features[f'{acid}sheet'] = sheetseq.tolist().count(acid) / len(sheetseq)
        features[acid] = data.count(acid) / Protlen
    motifs = [motif[:2] for motif in np.array(features.iloc[107:243].index)]
    for motif in motifs:
        features[f'{motif}motif'] = data.count(motif) / Protlen
    features['EALRmotif'] = data.count('EALR') / Protlen
    features['LEALmotif'] = data.count('LEAL') / Protlen
    features['HydrophobicAA'] = features['A'] + features['V'] + features['I'] + features['L'] + features['M'] + features['F'] + features['W']
    features['ChargedAA'] = features['R'] + features['H'] + features['K'] + features['D'] + features['E']
    features['PolarAA'] = features['N'] + features['Q'] + features['S'] + features['T'] + features['Y']
    featuresdf = pd.DataFrame(features).T
    features_scaled = scaler1.transform(featuresdf)
    features_pca = pca1.transform(features_scaled)
    prediction = model1.predict(features_pca)
    if removefiles == True:
        os.remove(fastapath)
        os.remove(faspath)
    
    return [prediction,features,features_scaled,features_pca]