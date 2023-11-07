def interaction_values(pert1, pert2, combination, precision, dim):

    # get the linear estimate of the perturbation
    estimate_mean = pert1 + pert2

    # deviation of the real combination form the estimate
    deviation = combination - estimate_mean

    # decompose into alpha, beta and gamma
    alpha_beta,_,_,_ = np.linalg.lstsq(np.vstack([pert1, pert2]).T, deviation.T,rcond=None)
    alpha, beta = alpha_beta

    # closest vector to deviation within the plane spanned by pert_mean1 and pert__mean2
    in_span = alpha*pert1 + beta*pert2
    emergent_effect = deviation - in_span
    gamma = scipy.spatial.distance.mahalanobis(emergent_effect,np.zeros(dim), precision)
    deviation_magnitude = scipy.spatial.distance.mahalanobis(deviation,      np.zeros(dim), precision)
    return alpha, beta, gamma, deviation_magnitude, emergent_effect, deviation  
