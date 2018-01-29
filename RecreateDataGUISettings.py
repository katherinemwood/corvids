__author__ = 'sean'

class RecreateDataGUISettings(object):

    def __init__(self, debug, min_score, max_score, poss_vals, mean, mean_precision, variance, variance_precision, num_samples, check_vals, use_SD):
        self.debug = debug
        self.min_score = min_score
        self.max_score = max_score
        self.poss_vals = poss_vals
        self.mean = mean
        self.mean_precision = mean_precision
        self.variance = variance
        self.variance_precision = variance_precision
        self.num_samples = num_samples
        self.check_vals = check_vals
        self.use_SD = use_SD

