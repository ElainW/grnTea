import os
import sys
import pandas as pd
from deepforest import CascadeForestClassifier
from sklearn.model_selection import train_test_split
from scipy.io import arff
import pickle
from sklearn.metrics import balanced_accuracy_score, precision_score, recall_score, confusion_matrix

class GntpClassifier_DF21():
    def __init__(self):
        self.n_feature = 15
        return
####
    ####gnerate the arff file for training data (with label)
    def gnrt_training_arff_from_xTEA_output(self, sf_01_list, sf_11_list, sf_arff, b_balance=False):
        with open(sf_arff, "w") as fout_arff:
            fout_arff.write("@RELATION\tinsgntp\n\n")
            fout_arff.write("@ATTRIBUTE\tAluclipcns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tL1clipcns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tSVAclipcns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tAluldisccns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tAlurdisccns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tL1ldisccns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tL1rdisccns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tSVAldisccns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tSVArdisccns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tpolyA\tNUMERIC\n")  # total polyA

            # fout_arff.write("@ATTRIBUTE\trpolyA\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tlcov\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\trcov\tNUMERIC\n")

            fout_arff.write("@ATTRIBUTE\tclip\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tfullmap\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tclipratio\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tdiscratio\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\trawlclip\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\trawrclip\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tdiscordant\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tconcordant\tNUMERIC\n")
            #fout_arff.write("@ATTRIBUTE\tindel\tNUMERIC\n")
            #fout_arff.write("@ATTRIBUTE\tclass\t{0,1,2}\n\n")
            fout_arff.write("@ATTRIBUTE\tclass\t{1,2}\n\n")
            fout_arff.write("@DATA\n")
            # with open(sf_00_list) as fin_00:
            #     for line in fin_00:
            #         sf_rslt = line.rstrip()
            #         l_features = self.load_in_feature_from_xTEA_output(sf_rslt)
            #         for rcd in l_features:
            #             fout_arff.write(",".join(rcd) + "\n") ##

            n_11=0
            with open(sf_11_list) as fin_11:
                for line in fin_11:
                    sf_rslt, sf_gt_ft = line.rstrip().split()
                    l_features = self.load_in_feature_from_xTEA_output(sf_rslt, sf_gt_ft, b_train=True)
                    for rcd in l_features:
                        fout_arff.write(",".join(rcd) + ",2,\n")
                    n_11+=1

            n_01=0
            with open(sf_01_list) as fin_01:
                for line in fin_01:
                    if n_01>n_11 and b_balance==True:
                        break
                    sf_rslt, sf_gt_ft = line.rstrip().split()
                    l_features = self.load_in_feature_from_xTEA_output(sf_rslt, sf_gt_ft, b_train=True)
                    for rcd in l_features:
                        fout_arff.write(",".join(rcd) + ",1,\n")
                    n_01+=1

####
    ####train the model
    def train_model(self, sf_arff, sf_model, f_ratio=0.2):
        data = arff.loadarff(sf_arff)
        df = pd.DataFrame(data[0])
        xVar = df.iloc[:, :self.n_feature]
        yVar = df.iloc[:, self.n_feature]
        yVar = yVar.astype('int')
    
        X_train, X_test, y_train, y_test = train_test_split(xVar, yVar, test_size=f_ratio, random_state=42, shuffle=True)
        clf = CascadeForestClassifier(n_jobs=-1, random_state=0, n_estimators=20)
        # clf = svm.SVC(kernel='linear')
        # clf=svm.SVC(kernel='rbf') #Gaussian Kernel
        clf.fit(X_train, y_train)
        with open(sf_model + "n20", 'wb') as file:
            pickle.dump(clf, file)
        preds = clf.predict(X_test)
        bal_accuracy = round(balanced_accuracy_score(y_test, preds), 3)
        precision = round(precision_score(y_test, preds), 3)
        recall = round(recall_score(y_test, preds), 3)
        mat = confusion_matrix(y_test, preds)
        print('Balanced Accuracy: {} / Precision: {} / Recall: {}'.format(bal_accuracy, precision, recall))
        print("Confusion matrix:")
        print(mat)
        print("-------------------------------------------------------------------------------------------")
        # tab = pd.crosstab(y_test, preds, rownames=['Actual Result'], colnames=['Predicted Result'])
        # print(tab)
        
        clf = CascadeForestClassifier(n_jobs=-1, random_state=0, n_estimators=10)
        # clf = svm.SVC(kernel='linear')
        # clf=svm.SVC(kernel='rbf') #Gaussian Kernel
        clf.fit(X_train, y_train)
        with open(sf_model + "n10", 'wb') as file:
            pickle.dump(clf, file)
        preds = clf.predict(X_test)
        bal_accuracy = round(balanced_accuracy_score(y_test, preds), 3)
        precision = round(precision_score(y_test, preds), 3)
        recall = round(recall_score(y_test, preds), 3)
        mat = confusion_matrix(y_test, preds)
        print('Balanced Accuracy: {} / Precision: {} / Recall: {}'.format(bal_accuracy, precision, recall))
        print("Confusion matrix:")
        print(mat)
        print("-------------------------------------------------------------------------------------------")
        # tab = pd.crosstab(y_test, preds, rownames=['Actual Result'], colnames=['Predicted Result'])
        # print(tab)
        
        clf = CascadeForestClassifier(n_jobs=-1, random_state=0, n_estimators=15)
        # clf = svm.SVC(kernel='linear')
        # clf=svm.SVC(kernel='rbf') #Gaussian Kernel
        clf.fit(X_train, y_train)
        with open(sf_model + "n15", 'wb') as file:
            pickle.dump(clf, file)
        preds = clf.predict(X_test)
        bal_accuracy = round(balanced_accuracy_score(y_test, preds), 3)
        precision = round(precision_score(y_test, preds), 3)
        recall = round(recall_score(y_test, preds), 3)
        mat = confusion_matrix(y_test, preds)
        print('Balanced Accuracy: {} / Precision: {} / Recall: {}'.format(bal_accuracy, precision, recall))
        print("Confusion matrix:")
        print(mat)
        print("-------------------------------------------------------------------------------------------")
        # tab = pd.crosstab(y_test, preds, rownames=['Actual Result'], colnames=['Predicted Result'])
        # print(tab)

    ####
    ####
    ####Given xTEA output, generate the prepared arff file
    def prepare_arff_from_xTEA_output(self, sf_xTEA, sf_genotype_ft, sf_arff, b_train):
        with open(sf_arff, "w") as fout_arff:
            fout_arff.write("@RELATION\tinsgntp\n\n")
            fout_arff.write("@ATTRIBUTE\tAluclipcns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tL1clipcns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tSVAclipcns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tAluldisccns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tAlurdisccns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tL1ldisccns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tL1rdisccns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tSVAldisccns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tSVArdisccns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tpolyA\tNUMERIC\n")  # total polyA

            # fout_arff.write("@ATTRIBUTE\trpolyA\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tlcov\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\trcov\tNUMERIC\n")

            fout_arff.write("@ATTRIBUTE\tclip\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tfullmap\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tclipratio\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tdiscratio\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\trawlclip\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\trawrclip\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tdiscordant\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tconcordant\tNUMERIC\n")
            #fout_arff.write("@ATTRIBUTE\tindel\tNUMERIC\n")
            #fout_arff.write("@ATTRIBUTE\tclass\t{0,1,2}\n\n")
            fout_arff.write("@ATTRIBUTE\tclass\t{1,2}\n\n")
            fout_arff.write("@DATA\n")

            l_features = self.load_in_feature_from_xTEA_output(sf_xTEA, sf_genotype_ft, b_train)
            for rcd in l_features:
                fout_arff.write(",".join(rcd) + "\n")

        data = arff.loadarff(sf_arff)
        df = pd.DataFrame(data[0])
        xVar = df.iloc[:, :self.n_feature]
        #yVar = df.iloc[:, self.n_feature]
        #yVar = yVar.astype('int')

        X=xVar.to_numpy()
        # X_train, X_test, y_train, y_test = train_test_split(xVar, yVar, test_size=0.999999)
        # return X_test
        return X
####

    ####clf is the trained model
    def predict_for_site(self, sf_model, sf_xTEA, sf_genotype_ft, sf_new):
        rf_model_df21 = CascadeForestClassifier()
        rf_model_df21.load(sf_model)

        sf_arff = sf_genotype_ft + ".arff"
        # site_features=self.prepare_arff_from_xTEA_output_two_category(sf_xTEA, sf_arff)

        site_features = self.prepare_arff_from_xTEA_output(sf_xTEA, sf_genotype_ft, sf_arff, b_train=False)
        preds=None
        if len(site_features)>0:
            preds = rf_model_df21.predict(site_features)
        with open(sf_xTEA) as fin_xTEA, open(sf_new, "w") as fout_new: # may need to update this!!!
            if None is preds:
                return
            i_idx = 0
            for line in fin_xTEA: # need to deal with this!!!
                sinfo = line.rstrip()
                s_gntp = "0/0"
                if preds[i_idx] == 1:
                    s_gntp = "0/1"
                elif preds[i_idx] == 2:
                    s_gntp = "1/1"
                sinfo += ("\t" + s_gntp + "\n")
                fout_new.write(sinfo)
                i_idx += 1

    ####
    def load_in_feature_from_xTEA_output(self, sf_xtea, sf_genotype_ft, b_train=True):
        '''
        Significant changes in grnTea as the output file is different
        here we merge the original features with features from x_genotype_feature.py
        '''
        l_all_features = []
        if os.path.isfile(sf_xtea)==False or os.path.isfile(sf_genotype_ft)==False:
            return l_all_features
        gt_ft_dict = {} # store features extracted from x_genotype_feature.py to merge with those in sf_xtea
        with open(sf_genotype_ft) as fin_ft:
            for line in fin_ft:
                fields = line.split()
                chrm = fields[0]
                ins_pos = int(fields[1])
                n_af_clip = int(fields[2])
                n_full_map = int(fields[3])
                n_l_raw_clip = int(fields[4])
                n_r_raw_clip = int(fields[5])
                n_disc_pairs = int(fields[6])
                n_concd_pairs = int(fields[7])
                n_disc_large_indel = int(fields[8])
                n_polyA = int(fields[10])
                n_disc_chrms = int(fields[11])
                if chrm not in gt_ft_dict:
                    gt_ft_dict[chrm] = {}
                if ins_pos in gt_ft_dict[chrm]:
                    sys.exit(f"{chrm}:{pos} duplicate entry")
                gt_ft_dict[chrm][ins_pos] = [n_af_clip, n_full_map, n_l_raw_clip, n_r_raw_clip, n_disc_pairs, n_concd_pairs, n_disc_large_indel, n_polyA, n_disc_chrms]
        with open(sf_xtea) as fin_xtea:
            for line in fin_xtea:
                fields = line.split()
                chrm = fields[0]
                ins_pos = int(fields[1])
                if chrm not in gt_ft_dict:
                    continue
                if ins_pos not in gt_ft_dict[chrm]:
                    continue
                l_features = self._parser_features(fields, gt_ft_dict[chrm][ins_pos])
                # if b_train == True: # currently don't have this yet
                #     s_gntp = fields[-1]
                #     if s_gntp == "FP" or s_gntp == "0/0":
                #         s_gntp = "0"
                #     elif s_gntp == "0/1" or s_gntp == "1/0":
                #         s_gntp = "1"
                #     elif s_gntp == "1/1":
                #         s_gntp = "2"
                #     l_features.append(s_gntp)
                # else:
                #     l_features.append("1")
                l_feature2 = [str(x) for x in l_features]
                l_all_features.append(l_feature2)
        return l_all_features
####
    ####
    # parse out the features
    def _parser_features(self, l_fields, gt_ft):  #
        l_features = []
        f_lcov = float(l_fields[25])
        if f_lcov == 0:
            f_lcov = 0.0000000001
        f_rcov = float(l_fields[16])
        if f_rcov == 0:
            f_rcov = 0.0000000001
        n_polyA = gt_ft[8]
        
        l_features.append(float(l_fields[7])/(f_lcov+f_rcov))#Alu-clip-algn-on-consensus
        l_features.append(float(l_fields[8])/(f_lcov+f_rcov))#L1-clip-algn-on-consensus
        l_features.append(float(l_fields[9])/(f_lcov+f_rcov))#SVA-lclip-algn-on-consensus
        l_features.append(float(l_fields[10])/f_lcov)#Alu-ldisc-algn-on-consensus
        l_features.append(float(l_fields[11])/f_rcov)#Alu-rdisc-algn-on-consensus
        l_features.append(float(l_fields[12])/f_lcov)#L1-ldisc-algn-on-consensus
        l_features.append(float(l_fields[13])/f_rcov)#L1-rdisc-algn-on-consensus
        l_features.append(float(l_fields[14])/f_lcov)#SVA-ldisc-algn-on-consensus
        l_features.append(float(l_fields[15])/f_lcov)#SVA-rdisc-algn-on-consensus
        l_features.append(n_polyA/(f_lcov+f_rcov))#polyA        

        # l_features.append()#right-polyA
        l_features.append(f_lcov)#left-local-coverage
        l_features.append(f_rcov)#right-local-coverage
        
        n_af_clip = gt_ft[0]
        n_full_map = gt_ft[1]
        l_features.append(n_af_clip/(f_lcov+f_rcov))#effective clip
        l_features.append(n_full_map/(f_lcov+f_rcov))#effective fully mapped

        f_clip_ratio=0
        if n_af_clip + n_full_map > 0:
            f_clip_ratio = n_af_clip / (n_af_clip + n_full_map)  # clip ration
        l_features.append(f_clip_ratio)
        
        n_disc_pairs, n_concd_pairs = gt_ft[4], gt_ft[5]
        f_disc_ratio=0
        if n_disc_pairs + n_concd_pairs>0:
            f_disc_ratio = n_disc_pairs / (n_disc_pairs + n_concd_pairs)  # disc ratio
        l_features.append(f_disc_ratio)
        
        n_l_raw_clip, n_r_raw_clip = gt_ft[2], gt_ft[3]
        l_features.append(n_l_raw_clip/f_lcov)#raw left-clip
        l_features.append(n_r_raw_clip/f_rcov)# raw right-clip
        l_features.append(n_disc_pairs/(f_lcov+f_rcov))#discordant pairs
        l_features.append(n_concd_pairs/(f_lcov+f_rcov))# conordant pairs
        #l_features.append(float(l_fields[41]) / (f_lcov + f_rcov))  # indels reads

        self.n_feature = len(l_features)
        return l_features
#pkl_filename="/n/data1/hms/dbmi/park
