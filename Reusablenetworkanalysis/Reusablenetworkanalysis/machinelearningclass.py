import sklearn
from .Analysisclass import Analysis
from .Vectorarithmetic import CalculateInteractions


class Machinelearning(CalculateInteractions):
    def __init__(self, analysis_object):
        if not isinstance(analysis_object, Analysis):
            raise ValueError("Analysis object is not a Analysis object defined in Analysisclass.py")
        super().__init__()
