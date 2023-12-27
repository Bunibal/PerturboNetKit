import numpy as np
import pandas as pd
import scipy
from sklearn.decomposition import PCA
from pertubationobject import Perturbome


class CalculateInteractions:
    """
    Class for calculating the interactions between perturbations in a pertubation object.
    """

    def __init__(self, perturbation_object, no_treatment_samples):
        if not isinstance(perturbation_object, Perturbome):
            raise ValueError("Perturbation object is not a perturbome object defined in pertubationobject.py")
        self.perturbation_object = perturbation_object
        self.precision = self.calculate_precision(no_treatment_samples)

    def pca(self, n_comp=5):
        """
        Perform a PCA on the perturbations
        :param n_comp: number of components to keep
        :return: components
        """
        pca = PCA(n_components=n_comp)
        return pca.fit(self.perturbation_object.perturbations)

    def calculate_precision(self, no_treatment):
        """
        Calculate the precision matrix of the perturbations. First a Yeojohnson transformation is applied to the data.
        Then the covariance matrix is calculated and the precision matrix is calculated as the inverse of the covariance matrix.
        :param no_treatment: matrix of the no treatment samples
        :return: precision matrix
        """
        # create empty dataframe for the transformed data
        transformed_data = pd.DataFrame(columns=no_treatment.columns)
        for row in no_treatment:
            transformed_data.loc[len(transformed_data)] = scipy.stats.yeojohnson(row)
        # covariance matrix
        cov = transformed_data.cov()
        # precision matrix
        precision = np.linalg.inv(cov)
        return precision

    def interaction_value(self, pert1, pert2, combination):
        """
        Calculate the interaction value between two perturbations using the mahanalobis distance
        :param pert1: 1-d array of perturbation 1
        :param pert2: 1-d array of perturbation 2
        :param combination: 1-d array of the combination of perturbation 1 and 2
        :return: the values of alpha, beta, gamma, deviation_magnitude, emergent_effect, deviation
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
        return alpha, beta, gamma, deviation_magnitude, emergent_effect, deviation

    def get_interaction_values(self, perturbome):
        """
        Calculate the interaction values for all perturbations in the pertubation object
        :param perturbome: A pertubation object from pertubationobject.py
        :return: dataframe with the interaction values
        """
        df = pd.DataFrame(columns=["alpha", "beta", "gamma", "deviation_magnitude", "emergent_effect", "deviation"])
        for key in perturbome.perturbations:
            for key1 in perturbome.perturbations:
                if key != key1:
                    pert1 = perturbome.perturbations[key]
                    pert2 = perturbome.perturbations[key1]
                    combination = perturbome.interactions[(key, key1)]
                    alpha, beta, gamma, deviation_magnitude, emergent_effect, deviation = self.interaction_value(pert1,
                                                                                                                 pert2,
                                                                                                                 combination)
                    df.loc[key + "_" + key1] = [alpha, beta, gamma, deviation_magnitude, emergent_effect, deviation]
        return df
