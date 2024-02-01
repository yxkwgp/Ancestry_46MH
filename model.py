######################Part1: Build model######################
## Import packages
import numpy as np
import pandas as pd
from sklearn import preprocessing
from sklearn.naive_bayes import CategoricalNB
from scipy.special import logsumexp
from sklearn import metrics

## Load training set
dataStr = pd.read_csv('training_set.txt', sep = ' ')
MHsNames = list(dataStr)[2:]

## LabelEncoder coding
# Y: superpop
lePop = preprocessing.LabelEncoder()
lePop.fit(dataStr['superPop'].drop_duplicates())
superPop = lePop.transform(dataStr['superPop'])
# X: genotypes of MHs
LabelEncoders = []
dataNum = [[1] * (dataStr.shape[1] - 2)] * dataStr.shape[0]
dataNum = np.array(dataNum)
for i in range(2, dataStr.shape[1]):
    exec(f'le{i - 2} = preprocessing.LabelEncoder()')
    exec(f'le{i - 2}.fit(dataStr.iloc[:, {i}].drop_duplicates())')
    exec(f'dataNum[:, {i} - 2] = le{i - 2}.transform(dataStr.iloc[:, {i}])')
    exec(f'LabelEncoders.append(le{i - 2})')
    # The coding mode of each column of the training set is stored in le0, le1, le2, ..., le45 for subsequent coding of the genotypes from samples to be predicted into the serial numbers of the corresponding categories. 

## Build model and update
CNBmodel = CategoricalNB()
CNBmodel.fit(dataNum, superPop)
for i in range(dataNum.shape[1]):
    CNBmodel.category_count_[i] = np.append(CNBmodel.category_count_[i], [[0],[0],[0],[0],[0]], axis=1)
CNBmodel._update_feature_log_prob(alpha=1)
# In rare cases, some of the genotypes of the samples to be predicted are not in the training set, and in order to avoid errors, a new category with a count of 0 in any class (superpopulation) is added to each MH and the probabilities are recalculated. And the resulting zero-probability will be eliminated by default Laplace smoothing step. The genotypes out of vocabulary will be seen as the new catogory when being predicted. 

#####################Part2: predict input######################
## Function for calculating AUC
def calcAUC(keptIndex):
    # Extract kept colums of training set
    dataPart = dataNum[:, keptIndex]
    # Leave-one-out cross validation to get probability table from training set
    probTable = np.zeros(shape=(2504,5))
    for i in range(2504):
        trainX = np.delete(dataPart, [2 * i, 2 * i + 1], axis=0)
        testX = dataPart[2 * i: 2 * i + 2,]
        trainY = np.append(superPop[:2 * i], superPop[2 * i + 2:])
        CNB = CategoricalNB()
        CNB.fit(trainX, trainY)
        for j in range(dataPart.shape[1]):
            CNB.category_count_[j] = np.append(CNB.category_count_[j], [[0],[0],[0],[0],[0]], axis=1)
        CNB._update_feature_log_prob(alpha=1)
        jll = np.zeros(CNB.class_count_.shape[0])
        for k in range(CNB.n_features_in_):
            indice1 = testX[0, k]
            indice2 = testX[1, k]
            jll += CNB.feature_log_prob_[k][:, indice1]
            jll += CNB.feature_log_prob_[k][:, indice2]
        total_ll = jll + CNB.class_log_prior_
        log_prob_all = logsumexp(total_ll)
        normalizedLogProb = total_ll - log_prob_all
        probTable[i,] = np.exp(normalizedLogProb)
    # Calculate AUC based on target y and probabilities
    yPred = CNB.classes_[np.argmax(probTable, axis=1)]
    yTrue = superPop[::2]
    lb = preprocessing.LabelBinarizer()
    yTrueforROC = lb.fit_transform(yTrue)
    AUC = []
    for l in range(5):
        AUC.append(metrics.roc_auc_score(yTrueforROC[:,l], probTable[:,l]))
    return AUC

## Check MHs integrity of input file
dataInput = pd.read_csv('input.txt', sep = ' ', keep_default_na = False, dtype = 'string')
IDs = list(dataInput['IDs'])[::2]
dataInput = dataInput.reindex(columns = MHsNames)
dataInputNa = dataInput.isna()
for m in range(dataInput.shape[0]):
    for n in range(dataInput.shape[1]):
        if not dataInputNa.iloc[m, n]:
            SNPs = dataInput.iloc[m, n].split('-')
            if len(SNPs) == 1:
                dataInput.iloc[m, n] = np.nan
            else:
                for SNP in SNPs:
                    for character in SNP:
                        if character not in ['A', 'T', 'C', 'G']:
                            dataInput.iloc[m, n] = np.nan
                            break
                    else:
                        continue
                    break
dataAvail = pd.DataFrame(columns = MHsNames, dtype = bool)
dataInputNotNa = dataInput.notna()
for i in range(int(dataInputNotNa.shape[0] / 2)):
    newLine = dataInputNotNa.iloc[2 * i, :] & dataInputNotNa.iloc[2 * i + 1, :]
    dataAvail.loc[len(dataAvail)] = newLine.astype('int')

## Predict superpopulation based on genotypes of input
probInput = np.zeros(shape=(dataAvail.shape[0],5)) # For saving robabilities in superpopulations
for i in range(dataAvail.shape[0]):
    jll = np.zeros(CNBmodel.class_count_.shape[0])
    for j in range(dataAvail.shape[1]):
        if dataAvail.iloc[i, j]:
            A1 = dataInput.iloc[2 * i, j]
            A2 = dataInput.iloc[2 * i + 1, j]
            # LabelEncoder coding
            exec(f"""
if A1 in le{j}.classes_:
    A1Num = int(le{j}.transform([A1]))
else:
    A1Num = len(le{j}.classes_)
if A2 in le{j}.classes_:
    A2Num = int(le{j}.transform([A2]))
else:
    A2Num = len(le{j}.classes_)
            """)
            jll += CNBmodel.feature_log_prob_[j][:, A1Num]
            jll += CNBmodel.feature_log_prob_[j][:, A2Num]
    jll += jll + CNBmodel.class_log_prior_
    log_prob_all = logsumexp(jll)
    normalizedLogProb = jll - log_prob_all
    probInput[i,] = np.exp(normalizedLogProb)
# Predict based on probabilities in probInput
yPredNum = list(CNBmodel.classes_[np.argmax(probInput, axis=1)])
yPred = lePop.inverse_transform(yPredNum)

## Calculate AUC
dataMode = dataAvail.groupby(MHsNames)
AUCdict = dict()
for mode in dataMode.groups.keys():
    keptIndex = []
    for index in range(len(mode)):
        if mode[index]:
            keptIndex.append(index)
    AUC = calcAUC(keptIndex)
    IDindices = list(dataMode.groups[mode])
    for IDindex in IDindices:
        AUCdict[IDindex] = AUC

## Calculate likelihood ratio
LRs = []
for i in range(len(yPredNum)):
    otherSum = 0
    for j in range(5):
        if j != yPredNum[i]:
            otherSum += probInput[i, j]
    LR = np.max(probInput[i]) / otherSum
    LRs.append(LR)

#################### Part3: Data outputing ###################
# IDs yPred Probabilities... LR AUC... N_46 N_12  LR(likelihood ratio)=P(y=yPred)/P(y=others)
outfile = open('output.txt', 'w')
outfile.write('IDs\tyPred\tProb_AFR\tProb_AMR\tProb_EAS\tProb_EUR\tProb_SAS\tLR\tAUC_AFR\tAUC_AMR\tAUC_EAS\tAUC_EUR\tAUC_SAS\tN_46\tN_12\n')
for i in range(len(IDs)):
    # ID
    outfile.write(IDs[i] + '\t')
    # yPred
    outfile.write(str(yPred[i]) + '\t')
    # Probabilities
    probRound = ['%.2e'%j for j in probInput[i]]
    outfile.write('\t'.join(probRound) + '\t')
    # LR
    outfile.write('%.2e'%LRs[i] + '\t')
    # AUC...
    AUCround = ['%.4f'%j for j in AUCdict[i]]
    outfile.write('\t'.join(AUCround) + '\t')
    # Numbers of available markers in 46-MH panel and 12-MH panel
    N_46 = sum(dataAvail.iloc[i])
    N_12 = sum(dataAvail.iloc[i, [3, 7, 8, 10, 14, 26, 27, 29, 32, 34, 38, 41]])
    outfile.write(str(N_46) + '\t' + str(N_12) + '\n')
outfile.close()
