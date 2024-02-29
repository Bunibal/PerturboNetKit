class Perturbome:
    """
    Class for storing pertubation vectors and their interactions
    example:
    perturbations = {"A": [1, 2, 3], "B": [4, 5, 6]}
    interactions = {("A", "B"): [1, 2, 3]}
    """

    def __init__(self, perturbations, interactions):
        self.perturbations = perturbations
        self.interactions = interactions
        self.check_object()

    def check_object(self):
        """
        Check if the pertubation object is valid
        """
        if self.perturbations is None:
            raise ValueError("Perturbations are not defined")
        if self.interactions is None:
            raise ValueError("Interactions are not defined")
        if type(self.perturbations) is not dict:
            raise ValueError("Perturbations are not a dictionary")
        if type(self.interactions) is not dict:
            raise ValueError("Interactions are not a dictionary")
        if len(self.perturbations) == 0:
            raise ValueError("Perturbations are empty")
        if len(self.interactions) == 0:
            raise ValueError("Interactions are empty")
        for key in self.perturbations:
            for key1 in self.perturbations:
                if key != key1:
                    if (key, key1) not in self.interactions or (key1, key) not in self.interactions:
                        raise ValueError("Interaction between " + str(key) + " and " + str(key1) + " is not defined")
