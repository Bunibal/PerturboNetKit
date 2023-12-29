import numpy as np
import pandas as pd
import scipy
from sklearn.decomposition import PCA
from .perturbationobject import Perturbome

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
    """

    def __init__(self):
        self.precision = None
        self.perturbation_object = None

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

    def pca(self, n_comp=5):
        """
        Perform Principal Component Analysis (PCA) on the perturbations.

        :param n_comp: Number of components to keep.
        :type n_comp: int

        :return: Components resulting from PCA.
        :rtype: numpy.ndarray
        """
        pca = PCA(n_components=n_comp)
        return pca.fit(self.perturbation_object.perturbations)

    def calculate_precision(self, no_treatment):
        """
        Calculate the precision matrix of the perturbations.

        First, a Yeojohnson transformation is applied to the data.
        Then the covariance matrix is calculated, and the precision matrix is calculated as the inverse of the covariance matrix.

        :param no_treatment: Matrix of the no treatment samples.
        :type no_treatment: pandas.Dataframe

        :return: Precision matrix.
        :rtype: numpy.ndarray
        """
        # create empty dataframe for the transformed data
        transformed_data = pd.DataFrame(columns=no_treatment.columns)
        for row in no_treatment:
            transformed_data.loc[len(transformed_data)] = scipy.stats.yeojohnson(row)  # TODO: like this or sth else?

        # covariance matrix
        cov = transformed_data.cov()
        # precision matrix
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
        return alpha, beta, gamma, deviation_magnitude, emergent_effect, deviation

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
