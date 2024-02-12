import numpy as np
import pandas as pd
import scipy
from sklearn.decomposition import PCA
from .perturbationobject import Perturbome
import seaborn as sns
import matplotlib.pyplot as plt

flatten = lambda *n: (e for a in n
                      for e in (flatten(*a) if isinstance(a, (tuple, list)) else (a,)))


class CalculateInteractions:
    """
    Class for calculating the interactions between perturbations in a pertubation object.

    Attributes:
    perturbation_object : Perturbome or None
        The perturbation object.

    precision : float or None
        Precision calculated based on the number of treatment samples.

    interaction_values : pd.DataFrame or None
        Interaction values between perturbations.
    """

    def __init__(self, perturbation_object=None, no_treatment_samples=None):
        self.precision = None
        self.perturbation_object = None
        self.interaction_values = None
        self.interaction_categories = None
        self.no_treatment_samples = no_treatment_samples
        if perturbation_object is not None and no_treatment_samples is not None:
            self.add_perturbome(perturbation_object, no_treatment_samples)

    def add_perturbome(self, perturbation_object, no_treatment_samples):
        """
        Add a perturbome object and calculate precision.

        :param perturbation_object: The perturbome object to be added.
        :type perturbation_object: Perturbome

        :param no_treatment_samples: dataframe containing all the no treatment samples row by row
        :type no_treatment_samples: pandas.Dataframe

        :raises ValueError: If the perturbation_object is not a Perturbome object.

        """
        if not isinstance(perturbation_object, Perturbome):
            raise ValueError("Perturbation object is not a perturbome object defined in perturbationobject.py")
        self.perturbation_object = perturbation_object
        self.precision = self.calculate_precision(no_treatment_samples)

    def multiple_testing_corr(self, df):
        """
        Perform multiple testing correction (Bonferroni) on all columns in the dataframe which begin with 'p_val'.

        :param df: Dataframe containing the interaction values.
        :type df: pandas.Dataframe

        :return: Dataframe containing the corrected interaction values.
        :rtype: pandas.Dataframe
        """

        pvals_corr = df.copy()
        # Bonferroni correction
        pvals_corr = pvals_corr.filter(like='p_val')
        pvals_corr = pvals_corr.apply(lambda x: x * len(pvals_corr.columns))
        # replace the p_val columns with the corrected p_val columns
        for col in pvals_corr.columns:
            df[col] = pvals_corr[col]
        return df

    def pca(self, n_comp=5):
        """
        Perform Principal Component Analysis (PCA) on the perturbations, interactions and no_treatment samples. Used to reduce # of features.

        :param n_comp: Number of components to keep.
        :type n_comp: int

        :return: Components resulting from PCA.
        :rtype: numpy.ndarray
        """
        pca = PCA(n_components=n_comp)
        self.no_treatment_samples.columns = range(len(self.no_treatment_samples.columns))
        matrix = pd.concat([pd.DataFrame.from_dict(self.perturbation_object.perturbations, orient='index'),
                            pd.DataFrame.from_dict(self.perturbation_object.interactions, orient='index'),
                            self.no_treatment_samples], ignore_index=True)
        pca.fit(matrix.T)
        # seperate perturbations, interactions and no_treatment samples
        perts = pca.components_[:, :len(self.perturbation_object.perturbations)]
        self.perturbation_object.perturbations = {key: perts[:, i] for i, key in
                                                  enumerate(self.perturbation_object.perturbations.keys())}
        inters = pca.components_[:,
                 len(self.perturbation_object.perturbations):len(self.perturbation_object.perturbations) +
                                                             len(self.perturbation_object.interactions)]
        self.perturbation_object.interactions = {key: inters[:, i] for i, key in
                                                 enumerate(self.perturbation_object.interactions.keys())}
        no_treat = pca.components_[:, len(self.perturbation_object.perturbations) +
                                      len(self.perturbation_object.interactions):]
        self.no_treatment_samples = pd.DataFrame(no_treat, columns=self.no_treatment_samples.index).T

    def calculate_precision(self, no_treatment, feature_reduction=None):
        """
        Calculate the precision matrix of the perturbations.

        First, a Yeojohnson transformation is applied to the data.
        Then the covariance matrix is calculated, and the precision matrix is calculated as the inverse of the covariance matrix.
        If the covariance matrix is not invertible (more features than samples), a PCA is performed and the precision matrix is calculated from the PCA components.

        :param no_treatment: Matrix of the no treatment samples. In shape of (n_samples, n_features).
        :type no_treatment: pandas.Dataframe

        :param feature_reduction: Method for feature reduction. Default is PCA. options: PCA, PINV
        :type feature_reduction: str

        :return: Precision matrix.
        :rtype: numpy.ndarray
        """
        if feature_reduction is None:
            feature_reduction = "PINV"

        if len(no_treatment) < len(no_treatment.columns) and feature_reduction == "PCA":
            # if there are more features than samples, perform PCA
            self.pca(n_comp=len(no_treatment))
        # create empty dataframe for the transformed data
        transformed_data = pd.DataFrame(columns=self.no_treatment_samples.columns)
        for index, row in self.no_treatment_samples.iterrows():
            transformed_data.loc[len(transformed_data)] = pd.Series(scipy.stats.yeojohnson(row)[0],
                                                                    index=self.no_treatment_samples.columns)
        # covariance matrix
        cov = transformed_data.cov()
        # precision matrix
        if feature_reduction == "PINV":
            precision = np.linalg.pinv(cov)
        else:
            precision = np.linalg.inv(cov)
        return precision

    def interaction_value(self, pert1, pert2, combination):
        """
        Calculate the interaction value between two perturbations using Mahalanobis distance.

        :param pert1: 1-D array of perturbation 1.
        :type pert1: numpy.ndarray

        :param pert2: 1-D array of perturbation 2.
        :type pert2: numpy.ndarray

        :param combination: 1-D array of the combination of perturbation 1 and 2.
        :type combination: numpy.ndarray

        :return: Tuple containing the values of alpha, beta, gamma, deviation_magnitude, emergent_effect, deviation.
        :rtype: tuple
        """
        dim = len(pert1)
        assert len(pert2) == dim
        assert len(combination) == dim

        # get the linear estimate of the perturbation
        estimate_mean = pert1 + pert2

        # deviation of the real combination form the estimate
        deviation = combination - estimate_mean

        # decompose into alpha, beta and gamma
        alpha_beta, _, _, _ = np.linalg.lstsq(np.vstack([pert1, pert2]).T, deviation.T, rcond=None)
        alpha, beta = alpha_beta

        # closest vector to deviation within the plane spanned by pert_mean1 and pert__mean2
        in_span = alpha * pert1 + beta * pert2
        emergent_effect = deviation - in_span
        gamma = scipy.spatial.distance.mahalanobis(emergent_effect, np.zeros(dim), self.precision)
        deviation_magnitude = scipy.spatial.distance.mahalanobis(deviation, np.zeros(dim), self.precision)
        dev_pert1 = scipy.spatial.distance.mahalanobis((alpha - 1) * pert1, np.zeros(dim), self.precision)
        dev_pert2 = scipy.spatial.distance.mahalanobis((beta - 1) * pert2, np.zeros(dim), self.precision)
        p_val_emergent = scipy.stats.chi2.sf(gamma, dim)  # apparently using chi2 distribution to calculate p-value
        p_val_deviation = scipy.stats.chi2.sf(deviation_magnitude, dim)
        p_val_pert1 = scipy.stats.chi2.sf(dev_pert1, dim)
        p_val_pert2 = scipy.stats.chi2.sf(dev_pert2, dim)
        return alpha, beta, gamma, deviation_magnitude, emergent_effect, deviation, p_val_deviation, p_val_emergent, p_val_pert1, p_val_pert2

    def get_interaction_values(self, perturbome=None):
        """
        Calculate the interaction values for all perturbations in the perturbation object.

        :param perturbome: A perturbation object from perturbationobject.py.
                           Defaults to None if not provided.
        :type perturbome: Perturbome or None

        :return: DataFrame with the interaction values.
        :rtype: pandas.Dataframe
        """
        if perturbome is None:
            perturbome = self.perturbation_object
        df = pd.DataFrame(
            columns=["alpha", "beta", "gamma", "deviation_magnitude", "emergent_effect", "deviation", "p_val_deviation",
                     "p_val_emergent", "p_val_pert1", "p_val_pert2"])
        for key in perturbome.perturbations:
            for key1 in perturbome.perturbations:
                if key != key1:
                    pert1 = perturbome.perturbations[key]
                    pert2 = perturbome.perturbations[key1]
                    combination = perturbome.interactions[(key, key1)]
                    alpha, beta, gamma, deviation_magnitude, emergent_effect, deviation, p_val_deviation, p_val_emergent, p_val_pert1, p_val_pert2 = self.interaction_value(
                        pert1,
                        pert2,
                        combination)
                    df.loc[key + "_" + key1] = [alpha, beta, gamma, deviation_magnitude, emergent_effect, deviation,
                                                p_val_deviation, p_val_emergent, p_val_pert1, p_val_pert2]

        self.interaction_values = self.multiple_testing_corr(df)
        return df

    def categorize_interactions(self, interaction_values=None, alpha_threshold=0.05):
        """
        Categorize interactions based on pvalues and interaction values.
        Interaction values are categorized as:
        - 000: Non-interacting
        - 001: Emergent effect
        - 100: uni-directional negative without emergent effect
        - 200: uni-directional positive without emergent effect
        - 010: uni-directional negative without emergent effect
        - 020: uni-directional positive without emergent effect
        - 101: uni-directional negative with emergent effect
        - 201: uni-directional positive with emergent effect
        - 011: uni-directional negative with emergent effect
        - 021: uni-directional positive with emergent effect
        - 110: bi-directional negative without emergent effect
        - 120: bi-directional negative and positive without emergent effect
        - 210: bi-directional positive and negative without emergent effect
        - 220: bi-directional positive without emergent effect
        - 111: bi-directional negative with emergent effect
        - 121: bi-directional negative and positive with emergent effect
        - 211: bi-directional positive and negative with emergent effect
        - 221: bi-directional positive with emergent effect

        First digit: perturbation 1 (e.g. 1 = pertubation 2 has a negative effect on perturbation 1)
        Second digit: perturbation 2 (e.g. 2 = pertubation 1 has a positive effect on perturbation 2)
        Third digit: emergent effect (e.g. 1 = emergent effect)

        :param alpha_threshold: Threshold for the p-value.
        :type alpha_threshold: float

        :param interaction_values: DataFrame containing the interaction values.
        :type interaction_values: pandas.Dataframe

        """
        if interaction_values is None:
            interaction_values = self.interaction_values

        interaction_cat = {}
        # categorize interactions
        for i, row in enumerate(interaction_values.iterrows()):
            row = row[1]
            if row[
                "p_val_deviation"] <= alpha_threshold:
                if row["p_val_emergent"] <= alpha_threshold and row["p_val_pert1"] <= alpha_threshold and row[
                    "p_val_pert2"] <= alpha_threshold:
                    if row["alpha"] < 1 and row["beta"] < 1:
                        interaction_cat[
                            interaction_values.index[i]] = "111"  # bi-directional negative with emergent effect
                    elif row["alpha"] > 1 and row["beta"] > 1:
                        interaction_cat[
                            interaction_values.index[i]] = "221"  # bi-directional positive with emergent effect
                    elif row["alpha"] < 1 and row["beta"] > 1:
                        interaction_cat[interaction_values.index[
                            i]] = "211"  # bi-directional negative and positive with emergent effect
                    elif row["alpha"] > 1 and row["beta"] < 1:
                        interaction_cat[interaction_values.index[
                            i]] = "121"  # bi-directional positive and negative with emergent effect
                elif row["p_val_emergent"] >= alpha_threshold and row["p_val_pert1"] <= alpha_threshold and row[
                    "p_val_pert2"] <= alpha_threshold:
                    if row["alpha"] < 1 and row["beta"] < 1:
                        interaction_cat[interaction_values.index[i]] = "110"
                    elif row["alpha"] > 1 and row["beta"] > 1:
                        interaction_cat[interaction_values.index[i]] = "220"
                    elif row["alpha"] < 1 and row["beta"] > 1:
                        interaction_cat[interaction_values.index[i]] = "120"
                    elif row["alpha"] > 1 and row["beta"] < 1:
                        interaction_cat[interaction_values.index[i]] = "210"
                elif row["p_val_emergent"] <= alpha_threshold and row["p_val_pert1"] >= alpha_threshold and row[
                    "p_val_pert2"] <= alpha_threshold:
                    if row["beta"] < 1:
                        interaction_cat[interaction_values.index[i]] = "011"
                    elif row["beta"] > 1:
                        interaction_cat[interaction_values.index[i]] = "021"
                elif row["p_val_emergent"] <= alpha_threshold and row["p_val_pert1"] <= alpha_threshold and row[
                    "p_val_pert2"] >= alpha_threshold:
                    if row["alpha"] < 1:
                        interaction_cat[interaction_values.index[i]] = "101"
                    elif row["alpha"] > 1:
                        interaction_cat[interaction_values.index[i]] = "201"
                elif row["p_val_emergent"] >= alpha_threshold and row["p_val_pert1"] >= alpha_threshold and row[
                    "p_val_pert2"] <= alpha_threshold:
                    if row["beta"] < 1:
                        interaction_cat[interaction_values.index[i]] = "010"
                    elif row["beta"] > 1:
                        interaction_cat[interaction_values.index[i]] = "020"
                elif row["p_val_emergent"] >= alpha_threshold and row["p_val_pert1"] <= alpha_threshold and row[
                    "p_val_pert2"] >= alpha_threshold:
                    if row["alpha"] < 1:
                        interaction_cat[interaction_values.index[i]] = "100"
                    elif row["alpha"] > 1:
                        interaction_cat[interaction_values.index[i]] = "200"
            else:
                interaction_cat[interaction_values.index[i]] = "000"  # non-interacting
        self.interaction_categories = interaction_cat
        return interaction_cat

    def plot_interactions_histogram(self, interaction_cat=None):
        """
        Plot heatmap of interaction categories.
        """
        if interaction_cat is None:
            interaction_cat = self.interaction_categories
        interaction_cat = pd.DataFrame(interaction_cat, index=["interaction_category"]).T
        interaction_cat = interaction_cat.groupby("interaction_category").size()
        interaction_cat = interaction_cat / interaction_cat.sum()
        interaction_cat.plot(kind="bar")
        plt.xlabel("Interaction category")
        plt.ylabel("Frequency")
        plt.show()
