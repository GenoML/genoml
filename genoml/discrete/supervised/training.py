# Importing the necessary packages 
import pandas as pd
from sklearn.metrics import accuracy_score, balanced_accuracy_score, log_loss, roc_auc_score, confusion_matrix, roc_curve, auc
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier, GradientBoostingClassifier, BaggingClassifier
from sklearn.svm import SVC
from sklearn.naive_bayes import ComplementNB
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis
from xgboost import XGBClassifier
import sklearn
import time
from joblib import dump, load
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.feature_selection import RFE
from sklearn.model_selection import train_test_split

# Define the train class 
class train:
    def __init__(self, df):
        #code that will prepare the data
       
        y = df.PHENO
        X = df.drop(columns=['PHENO'])
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42) # 70:30
        IDs_train = X_train.ID
        IDs_test = X_test.ID
        X_train = X_train.drop(columns=['ID'])
        X_test = X_test.drop(columns=['ID'])

        #saving the preped data the other classes will need
        self.df = df
        self.X_train = X_train
        self.X_test = X_test
        self.y_train = y_train
        self.y_test = y_test
        self.IDs_train = IDs_train
        self.IDs_test = IDs_test

        #Where the results will be stored 
        self.log_table = None
        self.best_algo = None
        self.algo = None
        self.rfe_df = None

        #The methods we will use
        self.algorithms = [
        LogisticRegression(),
        RandomForestClassifier(),
        AdaBoostClassifier(),
        GradientBoostingClassifier(),
        SGDClassifier(loss='modified_huber'),
        SVC(probability=True),
        MLPClassifier(),
        KNeighborsClassifier(),
        LinearDiscriminantAnalysis(),
        QuadraticDiscriminantAnalysis(),
        BaggingClassifier(),
        XGBClassifier()
        ]
    
    #report and data summary you want 
    def summary(self):
        print("Your data looks like this (showing the first few lines of the left-most and right-most columns)...")
        print("#"*70)
        print(self.df.describe())
        print("#"*70)

    def compete(self, verbose= False):
        log_cols=["Algorithm", "AUC_Percent", "Accuracy_Percent", "Balanced_Accuracy_Percent", "Log_Loss", "Sensitivity", "Specificity", "PPV", "NPV", "Runtime_Seconds"]
        log_table = pd.DataFrame(columns=log_cols)

        for algo in self.algorithms:
            
            start_time = time.time()
            
            algo.fit(self.X_train, self.y_train)
            name = algo.__class__.__name__

            print("")
            print("#"*70)
            print("")
            print(name)

            test_predictions = algo.predict_proba(self.X_test)
            test_predictions = test_predictions[:, 1]
            rocauc = roc_auc_score(self.y_test, test_predictions)
            print("AUC: {:.4%}".format(rocauc))

            test_predictions = algo.predict(self.X_test)
            acc = accuracy_score(self.y_test, test_predictions)
            print("Accuracy: {:.4%}".format(acc))

            test_predictions = algo.predict(self.X_test)
            balacc = balanced_accuracy_score(self.y_test, test_predictions)
            print("Balanced Accuracy: {:.4%}".format(balacc))
            
            CM = confusion_matrix(self.y_test, test_predictions)
            TN = CM[0][0]
            FN = CM[1][0]
            TP = CM[1][1]
            FP = CM[0][1]
            sensitivity = TP/(TP+FN)
            specificity = TN/(TN+FP)
            PPV = TP/(TP+FP)
            NPV = TN/(TN+FN)
            

            test_predictions = algo.predict_proba(self.X_test)
            ll = log_loss(self.y_test, test_predictions)
            print("Log Loss: {:.4}".format(ll))
            
            end_time = time.time()
            elapsed_time = (end_time - start_time)
            print("Runtime in seconds: {:.4}".format(elapsed_time))

            log_entry = pd.DataFrame([[name, rocauc*100, acc*100, balacc*100, ll, sensitivity, specificity, PPV, NPV, elapsed_time]], columns=log_cols)
            log_table = log_table.append(log_entry)

        print("#"*70)

        print("")

        """
        log_outfile = self.run_prefix + '.training_withheldSamples_performanceMetrics.csv'

        print(f"This table below is also logged as {log_outfile} and is in your current working directory...")
        print("#"*70)
        print(log_table)
        print("#"*70)

        log_table.to_csv(log_outfile, index=False)
        """

        self.log_table = log_table

        return log_table

    def results(self):
        best_performing_summary = self.log_table[self.log_table.AUC_Percent == self.log_table.AUC_Percent.max()]
        best_algo = best_performing_summary.at[0,'Algorithm']
    
        """
        best_algo_name_out = self.run_prefix + ".best_algorithm.txt"
        file = open(best_algo_name_out,'w')
        file.write(best_algo)
        file.close() 
        """

        self.best_algo = best_algo

        return best_algo

    def feature_ranking(self):
        best_algo = self.best_algo
        X_train = self.X_train
        y_train = self.y_train
        if (best_algo == 'SVC') or (best_algo == 'ComplementNB') or (best_algo == 'KNeighborsClassifier') or (best_algo == 'QuadraticDiscriminantAnalysis') or (best_algo == 'BaggingClassifier'):
        
            print("Even if you selected to run feature ranking, you can't generate feature ranks using SVC, ComplementNB, KNeighborsClassifier, QuadraticDiscriminantAnalysis, or BaggingClassifier... it just isn't possible.")
        
        else:

            top_ten_percent = (len(X_train)//10)
            # core_count = args.n_cores
            names = list(X_train.columns)
            rfe = RFE(estimator=self.algo)
            rfe.fit(X_train, y_train)
            rfe_out = zip(rfe.ranking_, names)
            rfe_df = pd.DataFrame(rfe_out, columns = ["RANK","FEATURE"])
            #table_outfile = self.run_prefix + '.trainedModel_trainingSample_featureImportance.csv'
            #rfe_df.to_csv(table_outfile, index=False)

            self.rfe_df = rfe_df

            return rfe_df

    def AUC(self, save = False):
        plot_out = self.run_prefix + '.trainedModel_withheldSample_ROC.png'

        test_predictions = self.algo.predict_proba(self.X_test)
        test_predictions = test_predictions[:, 1]

        fpr, tpr, thresholds = roc_curve(self.y_test, test_predictions)
        roc_auc = auc(fpr, tpr)

        plt.figure()
        plt.plot(fpr, tpr, color='purple', label='ROC curve (area = %0.2f)' % roc_auc)
        plt.plot([0, 1], [0, 1], color='cyan', linestyle='--', label='Chance (area = %0.2f)' % 0.5)
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False positive rate')
        plt.ylabel('True positive rate')
        plt.title('Receiver operating characteristic (ROC) - ' + self.best_algo)
        plt.legend(loc="lower right")
        if (save):
            plt.savefig(plot_out, dpi = 600)

        #print()
        #print(f"We are also exporting a ROC curve for you here {plot_out} this is a graphical representation of AUC in the withheld test data for the best performing algorithm.")
    
    def export_prob_hist(self):
        # Exporting withheld test data
        test_predicteds_probs = self.algo.predict_proba(self.X_test)
        test_case_probs = test_predicteds_probs[:, 1]
        test_predicted_cases = self.algo.predict(self.X_test)

        test_case_probs_df = pd.DataFrame(test_case_probs)
        test_predicted_cases_df = pd.DataFrame(test_predicted_cases)
        y_test_df = pd.DataFrame(self.y_test)
        IDs_test_df = pd.DataFrame(self.IDs_test)

        test_out = pd.concat([IDs_test_df.reset_index(), y_test_df.reset_index(drop=True), test_case_probs_df.reset_index(drop=True), test_predicted_cases_df.reset_index(drop=True)], axis = 1, ignore_index=True)
        test_out.columns=['INDEX','ID',"CASE_REPORTED","CASE_PROBABILITY","CASE_PREDICTED"]
        test_out = test_out.drop(columns=['INDEX'])

        test_outfile = self.run_prefix + '.trainedModel_withheldSample_Predictions.csv'
        test_out.to_csv(test_outfile, index=False)

        print("")
        print("Preview of the exported predictions for the withheld test data that has been exported as", test_outfile, "these are pretty straight forward.")
        print("They generally include the sample ID, the previously reported case status (1 = case), the case probability from the best performing algorithm and the predicted label from that algorithm")
        print("")
        print("#"*70)
        print(test_out.head())
        print("#"*70)


        # Exporting training data, which is by nature overfit.
        train_predicteds_probs = self.algo.predict_proba(self.X_train)
        train_case_probs = train_predicteds_probs[:, 1]
        train_predicted_cases = self.algo.predict(self.X_train)

        train_case_probs_df = pd.DataFrame(train_case_probs)
        train_predicted_cases_df = pd.DataFrame(train_predicted_cases)
        y_train_df = pd.DataFrame(self.y_train)
        IDs_train_df = pd.DataFrame(self.IDs_train)

        train_out = pd.concat([IDs_train_df.reset_index(), y_train_df.reset_index(drop=True), train_case_probs_df.reset_index(drop=True), train_predicted_cases_df.reset_index(drop=True)], axis = 1, ignore_index=True)
        train_out.columns=['INDEX','ID',"CASE_REPORTED","CASE_PROBABILITY","CASE_PREDICTED"]
        train_out = train_out.drop(columns=['INDEX'])

        train_outfile = self.run_prefix + '.trainedModel_trainingSample_Predictions.csv'
        train_out.to_csv(train_outfile, index=False)

        print("")
        print("Preview of the exported predictions for the training samples which is naturally overfit and exported as", train_outfile, "in the similar format as in the withheld test dataset that was just exported.")
        print("#"*70)
        print(train_out.head())
        print("#"*70)

        # Export historgrams of probabilities
        genoML_colors = ["cyan","purple"]

        g = sns.FacetGrid(train_out, hue="CASE_REPORTED", palette=genoML_colors, legend_out=True,)
        g = (g.map(sns.distplot, "CASE_PROBABILITY", hist=False, rug=True))
        g.add_legend()

        plot_out = self.run_prefix + '.trainedModel_withheldSample_probabilities.png'
        g.savefig(plot_out, dpi=600)

        print("")
        print("We are also exporting probability density plots to the file", plot_out, "this is a plot of the probability distributions of being a case, stratified by case and control status in the withheld test samples.")

    def export_model(self):
        best_algo = self.best_algo

        if best_algo == 'LogisticRegression':
            algo = getattr(sklearn.linear_model, best_algo)()

        if  best_algo == 'SGDClassifier':
            algo = getattr(sklearn.linear_model, best_algo)(loss='modified_huber')

        if (best_algo == 'RandomForestClassifier') or (best_algo == 'AdaBoostClassifier') or (best_algo == 'GradientBoostingClassifier') or  (best_algo == 'BaggingClassifier'):
            algo = getattr(sklearn.ensemble, best_algo)()

        if best_algo == 'SVC':
            algo = getattr(sklearn.svm, best_algo)(probability=True)

        if best_algo == 'ComplementNB':
            algo = getattr(sklearn.naive_bayes, best_algo)()

        if best_algo == 'MLPClassifier':
            algo = getattr(sklearn.neural_network, best_algo)()

        if best_algo == 'XGBClassifier':
            algo = getattr(xgboost, best_algo)()

        if best_algo == 'KNeighborsClassifier':
            algo = getattr(sklearn.neighbors, best_algo)()

        if (best_algo == 'LinearDiscriminantAnalysis') or (best_algo == 'QuadraticDiscriminantAnalysis'):
            algo = getattr(sklearn.discriminant_analysis, best_algo)()

        algo.fit(self.X_train, self.y_train)
        name = algo.__class__.__name__

        print("...remember, there are occasionally slight fluctuations in model performance on the same withheld samples...")

        print("#"*70)

        print(name)

        test_predictions = algo.predict_proba(self.X_test)
        test_predictions = test_predictions[:, 1]
        rocauc = roc_auc_score(self.y_test, test_predictions)
        print("AUC: {:.4%}".format(rocauc))

        test_predictions = algo.predict(self.X_test)
        acc = accuracy_score(self.y_test, test_predictions)
        print("Accuracy: {:.4%}".format(acc))

        test_predictions = algo.predict(self.X_test)
        balacc = balanced_accuracy_score(self.y_test, test_predictions)
        print("Balanced Accuracy: {:.4%}".format(balacc))

        test_predictions = algo.predict_proba(self.X_test)
        ll = log_loss(self.y_test, test_predictions)
        print("Log Loss: {:.4}".format(ll))

        ### Save it using joblib
       
        algo_out = self.run_prefix + '.trainedModel.joblib'
        dump(algo, algo_out)

        print("#"*70)

        print(f"... this model has been saved as {algo_out} for later use and can be found in your working directory.")

        self.algo = algo

        return algo

    def save_results(self, path, algorithmResults = False, bestAlgorithm = False, featureRankings = False):
        if(algorithmResults):
            log_table = self.log_table
            log_outfile = path + '.training_withheldSamples_performanceMetrics.csv'

            print(f"This table below is also logged as {log_outfile} and is in your current working directory...")
            print("#"*70)
            print(log_table)
            print("#"*70)

            log_table.to_csv(log_outfile, index=False)

        if(bestAlgorithm):
            best_algo_name_out = path + ".best_algorithm.txt"
            file = open(best_algo_name_out,'w')
            file.write(self.best_algo)
            file.close() 

        if(featureRankings):
            table_outfile = path + '.trainedModel_trainingSample_featureImportance.csv'
            self.rfe_df.to_csv(table_outfile, index=False)   
        